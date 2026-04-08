#define _CRT_SECURE_NO_WARNINGS 1
#include <stdio.h>
#include <vector>
#include <cmath>
#include <algorithm>

#include "raylib/src/raylib.h"
#include "dsp/optimizer.h"
#include <complex>

struct PoleModal
{
	const float Ts = 1.0f / 48000.0f;

	float a1 = 0, a2 = 0;
	float z1 = 0, z2 = 0;

	float pre = 0, pim = 0, rre = 0, rim = 0;

	std::complex<float> pole = { 0, 0 };
	std::complex<float> residue = { 0, 0 };
	std::complex<float> step1 = { 1, 0 };

	inline float ProcessSample()
	{
		float y = z1;
		z1 = -a1 * y + z2;
		z2 = -a2 * y;
		return y;
	}

	void CalcPole(float pre, float pim, float rre, float rim)
	{
		this->pre = pre;
		this->pim = pim;
		this->rre = rre;
		this->rim = rim;

		float R = expf(pre * Ts);
		float O = pim * Ts;
		a1 = -2.0f * R * cosf(O);
		a2 = R * R;

		pole = std::complex<float>(pre, pim);
		residue = std::complex<float>(rre, rim);
		step1 = std::exp(pole * Ts);
	}

	void InjectImpulse(float tau, float v)
	{
		if (tau < 0.0f) tau = 0.0f;
		if (tau >= 1.0f) tau = 1.0f;
		float dt1 = (1.0f - tau) * Ts;
		std::complex<float> shift = std::exp(pole * dt1);
		std::complex<float> A1 = v * residue * shift;
		std::complex<float> A2 = A1 * step1;
		float g1 = 2.0f * A1.real();
		float g2 = 2.0f * A2.real();
		z1 += g1;
		z2 += g2 + a1 * g1;
	}

	void Reset()
	{
		z1 = 0;
		z2 = 0;
	}
};

void fft(float re[], float im[], int N, int inv)
{
	int i, j, k, m;
	int len = N;
	j = 0;
	for (i = 1; i < len; i++) {
		int bit = len >> 1;
		while (j & bit) {
			j ^= bit;
			bit >>= 1;
		}
		j ^= bit;
		if (i < j) {
			float tmp;
			tmp = re[i]; re[i] = re[j]; re[j] = tmp;
			tmp = im[i]; im[i] = im[j]; im[j] = tmp;
		}
	}
	for (m = 2; m <= len; m <<= 1) {
		float angle = -inv * 2.0 * 3.1415926535897932384626f / m;
		float wm_re = cosf(angle);
		float wm_im = sinf(angle);
		for (k = 0; k < len; k += m) {
			float w_re = 1.0;
			float w_im = 0.0;

			for (j = 0; j < m / 2; j++) {
				int t = k + j;
				int u = t + m / 2;
				float tr = w_re * re[u] - w_im * im[u];
				float ti = w_re * im[u] + w_im * re[u];
				float ur = re[t];
				float ui = im[t];
				re[t] = ur + tr;
				im[t] = ui + ti;
				re[u] = ur - tr;
				im[u] = ui - ti;
				float next_w_re = w_re * wm_re - w_im * wm_im;
				float next_w_im = w_re * wm_im + w_im * wm_re;
				w_re = next_w_re;
				w_im = next_w_im;
			}
		}
	}
	if (inv == -1) {
		for (i = 0; i < len; i++) {
			re[i] /= len;
			im[i] /= len;
		}
	}
}


class SystemModal
{
public:
	constexpr static int NumOrders = 6;
private:
	PoleModal poles[NumOrders];
	std::vector<float> params;
public:
	void CalcPoles(std::vector<float>& params)
	{
		if (params.size() != NumOrders * 4) return;
		this->params = params;
		for (int i = 0; i < NumOrders; i++) {
			float pre = params[i * 4 + 0];
			float pim = params[i * 4 + 1];
			float rre = params[i * 4 + 2];
			float rim = params[i * 4 + 3];
			poles[i].CalcPole(pre, pim, rre, rim);
		}
	}
	void GetParams(std::vector<float>& params) const
	{
		params = this->params;
	}
	void InjectImpulse(float tau, float v)
	{
		for (int i = 0; i < NumOrders; i++) {
			poles[i].InjectImpulse(tau, v);
		}
	}
	float ProcessSample()
	{
		float y = 0;
		for (int i = 0; i < NumOrders; i++) {
			y += poles[i].ProcessSample();
		}
		return y;
	}
	void Reset()
	{
		for (int i = 0; i < NumOrders; i++) {
			poles[i].Reset();
		}
	}
};

class ModalLPFOptimizer
{
private:
	constexpr static int NumOrders = SystemModal::NumOrders;
	SystemModal modal1, modal2;
	AdamOptimizer optimizer;

	float randf()
	{
		return ((float)rand() / (float)RAND_MAX) * ((float)rand() / (float)RAND_MAX) * (rand() % 2 ? 1 : -1);
	}

	constexpr static int TestLen = 32;
	constexpr static int TestTauN = 8;
	float testre[TestLen];
	float testim[TestLen];
	float lastloss = -1;
	std::vector<float> params2;
	float Error(std::vector<float>& params)
	{
		float errv2 = 0;
		for (int i = 0; i < TestTauN; ++i)
		{
			modal1.Reset();
			params2 = params;
			Regularization(params2);
			modal1.CalcPoles(params2);
			float tau = (float)i / (float)TestTauN;
			modal1.InjectImpulse(tau, 1.0);
			for (int i = 0; i < TestLen; ++i)
			{
				float x = (float)i / (float)TestLen;
				float w = (1.0 - x * x);
				float window = w * w;
				testre[i] = modal1.ProcessSample() * window;
				testim[i] = 0;
			}
			fft(testre, testim, TestLen, 1);
			float errv = 0;
			for (int i = 0; i < TestLen / 2; ++i)
			{
				float x = (float)i / (float)(TestLen / 2 - 1);
				float target = x > 0.95 ? 0.05 : 1.0;
				float mag1 = sqrtf(testre[i] * testre[i] + testim[i] * testim[i]);
				float e1 = fabsf(logf(mag1) - logf(target));
				errv += e1 * e1;
			}
			errv /= TestLen;
			errv2 += errv;
		}
		lastloss = errv2;
		return errv2;
	}
public:
	ModalLPFOptimizer()
	{
		std::vector<float> initParams(NumOrders * 4);
		for (int i = 0; i < NumOrders; i++) {
			initParams[i * 4 + 0] = -fabsf(randf()) - 0.1; // pre
			initParams[i * 4 + 1] = fabsf(randf()) * 24000 + 0.1; // pim
			initParams[i * 4 + 2] = randf(); // rre
			initParams[i * 4 + 3] = randf(); // rim
		}
		modal1.CalcPoles(initParams);
		modal2.CalcPoles(initParams);
		optimizer.SetupOptimizer(NumOrders * 4, initParams,0.05);
		optimizer.SetErrorFunc([this](std::vector<float>& params) { return Error(params); });
	}
	float RunOptimizer()
	{
		optimizer.RunOptimizer(100);
		return lastloss;
	}
	void GetNowParams(std::vector<float>& params)
	{
		optimizer.GetNowVec(params);
	}


	float sigmoid(float x)
	{
		if (x >= 0.0f)
		{
			float e = expf(-x);
			return 1.0f / (1.0f + e);
		}
		else
		{
			float e = expf(x);
			return e / (1.0f + e);
		}
	}
	float softplus(float x)
	{
		if (x > 20.0f) return x;
		if (x < -20.0f) return expf(x);
		return logf(1.0f + expf(x));
	}
	void Regularization(std::vector<float>& input)
	{
		const float Fs = 48000.0f;
		const float M_PI = 3.14159265358979323846f;
		const float NyquistRad = M_PI * Fs; // 因为 O = pim * Ts，所以 pim 必须是 rad/s
		for (int i = 0; i < NumOrders; i++)
		{
			float rawSigma = input[i * 4 + 0];
			float rawOmega = input[i * 4 + 1];
			float rawRre = input[i * 4 + 2];
			float rawRim = input[i * 4 + 3];
			float sigma = -(10.0f + 12000.0f * softplus(rawSigma));
			float omega = 0.45f * NyquistRad * sigmoid(rawOmega);
			float gainScale = 2.5f / (float)NumOrders;
			float rre = gainScale * tanhf(rawRre);
			float rim = gainScale * tanhf(rawRim);
			input[i * 4 + 0] = sigma;
			input[i * 4 + 1] = omega;
			input[i * 4 + 2] = rre;
			input[i * 4 + 3] = rim;
		}
	}
};


Color HSVtoRGB(float h, float s, float v)
{
	float r, g, b;

	int i = int(h * 6.0f);
	float f = h * 6.0f - i;
	float p = v * (1.0f - s);
	float q = v * (1.0f - f * s);
	float t = v * (1.0f - (1.0f - f) * s);

	switch (i % 6)
	{
	case 0: r = v; g = t; b = p; break;
	case 1: r = q; g = v; b = p; break;
	case 2: r = p; g = v; b = t; break;
	case 3: r = p; g = q; b = v; break;
	case 4: r = t; g = p; b = v; break;
	case 5: r = v; g = p; b = q; break;
	}

	return Color{
		(unsigned char)(r * 255),
		(unsigned char)(g * 255),
		(unsigned char)(b * 255),
		255
	};
}

int main2()
{
	constexpr int ScreenW = 1800;
	constexpr int ScreenH = 960;
	constexpr int FFTLenResp = 8192;
	constexpr float SampleRate = 48000.0f;

	// ---- Sweep / STFT 参数 ----
	constexpr float SweepDurSec = 3.0f;
	constexpr int SweepSamples = (int)(SweepDurSec * SampleRate);
	constexpr float SweepF0 = 2400.0f;
	constexpr float SweepF1 = 24000.0f;

	constexpr int STFTSize = 1024;
	int hopSize = 512;                 // 可调
	constexpr float DbMin = -16.0f;
	constexpr float DbMax = 16.0f;

	InitWindow(ScreenW, ScreenH, "IIR Modal BLIT Optimizer + Sweep Spectrogram");
	SetTargetFPS(60);

	ModalLPFOptimizer opt;

	auto Clamp = [](float x, float a, float b) -> float
		{
			if (x < a) return a;
			if (x > b) return b;
			return x;
		};

	auto Lerp = [](float a, float b, float t) -> float
		{
			return a + (b - a) * t;
		};

	auto FreqToX = [&](float f, const Rectangle& rc) -> float
		{
			float fmin = 20.0f;
			float fmax = 24000.0f;
			f = Clamp(f, fmin, fmax);
			float t = (log10f(f) - log10f(fmin)) / (log10f(fmax) - log10f(fmin));
			return rc.x + t * rc.width;
		};

	auto DbToY = [&](float db, const Rectangle& rc) -> float
		{
			float dbMin = -30.0f;
			float dbMax = 30.0f;
			db = Clamp(db, dbMin, dbMax);
			float t = (db - dbMin) / (dbMax - dbMin);
			return rc.y + rc.height * (1.0f - t);
		};

	auto LogFreqToY = [&](float f, const Rectangle& rc) -> float
		{
			float fmin = 20.0f;
			float fmax = 24000.0f;
			f = Clamp(f, fmin, fmax);
			float t = (log10f(f) - log10f(fmin)) / (log10f(fmax) - log10f(fmin));
			return rc.y + rc.height * (1.0f - t);
		};

	auto ColorMapDb = [&](float db) -> Color
		{
			float t = (db - DbMin) / (DbMax - DbMin);
			t = Clamp(t, 0.0f, 1.0f);

			// 深色 -> 蓝 -> 青 -> 黄 -> 白
			float h, s, v;
			if (t < 0.25f)
			{
				float u = t / 0.25f;
				h = Lerp(0.72f, 0.62f, u);
				s = 1.0f;
				v = Lerp(0.08f, 0.45f, u);
			}
			else if (t < 0.55f)
			{
				float u = (t - 0.25f) / 0.30f;
				h = Lerp(0.62f, 0.50f, u);
				s = 1.0f;
				v = Lerp(0.45f, 0.95f, u);
			}
			else if (t < 0.82f)
			{
				float u = (t - 0.55f) / 0.27f;
				h = Lerp(0.50f, 0.15f, u);
				s = Lerp(1.0f, 0.85f, u);
				v = 1.0f;
			}
			else
			{
				float u = (t - 0.82f) / 0.18f;
				h = Lerp(0.15f, 0.0f, u);
				s = Lerp(0.85f, 0.0f, u);
				v = 1.0f;
			}
			return HSVtoRGB(h, s, v);
		};

	auto DrawGrid = [&](const Rectangle& rc)
		{
			DrawRectangleLinesEx(rc, 1.0f, GRAY);

			const float freqMarks[] = {
				20, 50, 100, 200, 500,
				1000, 2000, 5000, 10000, 20000,
			};

			for (float f : freqMarks)
			{
				float x = FreqToX(f, rc);
				DrawLine((int)x, (int)rc.y, (int)x, (int)(rc.y + rc.height), Fade(GRAY, 0.35f));

				char txt[32];
				if (f >= 1000.0f) sprintf(txt, "%.0fk", f / 1000.0f);
				else sprintf(txt, "%.0f", f);

				DrawText(txt, (int)x - 12, (int)(rc.y + rc.height + 8), 18, LIGHTGRAY);
			}

			for (int db = -30; db <= 30; db += 10)
			{
				float y = DbToY((float)db, rc);
				DrawLine((int)rc.x, (int)y, (int)(rc.x + rc.width), (int)y,
					db == 0 ? Fade(WHITE, 0.55f) : Fade(GRAY, 0.30f));

				char txt[32];
				sprintf(txt, "%d dB", db);
				DrawText(txt, (int)rc.x - 70, (int)y - 10, 18, LIGHTGRAY);
			}
		};

	auto DrawSpecGrid = [&](const Rectangle& rc, float durationSec)
		{
			DrawRectangleLinesEx(rc, 1.0f, GRAY);

			// 时间刻度：始终铺满整个框
			for (int i = 0; i <= 6; ++i)
			{
				float t = durationSec * (float)i / 6.0f;
				float x = rc.x + rc.width * ((float)i / 6.0f);
				DrawLine((int)x, (int)rc.y, (int)x, (int)(rc.y + rc.height), Fade(GRAY, 0.25f));

				char txt[32];
				sprintf(txt, "%.2fs", t);
				DrawText(txt, (int)x - 18, (int)(rc.y + rc.height + 8), 18, LIGHTGRAY);
			}

			// 线性频率刻度
			const float freqMarks[] = { 0, 4000, 8000, 12000, 16000, 20000, 24000 };
			for (float f : freqMarks)
			{
				float y = rc.y + rc.height * (1.0f - f / 24000.0f);
				DrawLine((int)rc.x, (int)y, (int)(rc.x + rc.width), (int)y, Fade(GRAY, 0.25f));

				char txt[32];
				if (f >= 1000.0f) sprintf(txt, "%.0fk", f / 1000.0f);
				else sprintf(txt, "%.0f", f);
				DrawText(txt, (int)rc.x - 48, (int)y - 10, 18, LIGHTGRAY);
			}
		};

	auto BuildSpectrumDb = [&](SystemModal& modal, std::vector<float>& outDb, float tau)
		{
			static float re[FFTLenResp];
			static float im[FFTLenResp];

			modal.Reset();
			modal.InjectImpulse(tau, 1.0f);

			for (int i = 0; i < FFTLenResp; ++i)
			{
				re[i] = modal.ProcessSample();
				im[i] = 0.0f;
			}

			fft(re, im, FFTLenResp, 1);

			outDb.resize(FFTLenResp / 2 + 1);
			for (int i = 0; i <= FFTLenResp / 2; ++i)
			{
				float mag = sqrtf(re[i] * re[i] + im[i] * im[i]);
				if (mag < 1.0e-20f) mag = 1.0e-20f;
				outDb[i] = 20.0f * log10f(mag);
			}
		};

	auto DrawSpectrum = [&](const std::vector<float>& db, const Rectangle& rc, Color color)
		{
			if (db.size() < 2) return;

			const int halfN = (int)db.size() - 1;
			for (int px = 0; px < (int)rc.width - 1; ++px)
			{
				float t0 = (float)px / rc.width;
				float t1 = (float)(px + 1) / rc.width;

				float f0 = powf(10.0f, log10f(20.0f) + t0 * (log10f(24000.0f) - log10f(20.0f)));
				float f1 = powf(10.0f, log10f(20.0f) + t1 * (log10f(24000.0f) - log10f(20.0f)));

				float bin0 = f0 * (float)FFTLenResp / SampleRate;
				float bin1 = f1 * (float)FFTLenResp / SampleRate;

				int i0 = (int)bin0;
				int i1 = (int)bin1;

				i0 = std::max(0, std::min(i0, halfN));
				i1 = std::max(0, std::min(i1, halfN));

				float y0 = DbToY(db[i0], rc);
				float y1 = DbToY(db[i1], rc);
				DrawLineV({ rc.x + (float)px, y0 }, { rc.x + (float)px + 1.0f, y1 }, color);
			}
		};
	auto BuildSweepSignal = [&](SystemModal& modal, std::vector<float>& out)
		{
			out.assign(SweepSamples, 0.0f);
			modal.Reset();

			float phase = 0.0f;

			for (int n = 0; n < SweepSamples; ++n)
			{
				float u = (float)n / (float)(SweepSamples - 1);

				// 指数扫频：100 Hz -> 24000 Hz
				float freq = SweepF0 * powf(SweepF1 / SweepF0, u);
				float incr = freq / SampleRate;

				float prevPhase = phase;
				phase += incr;

				// 一个 sample 内可能跨越多次（高频端接近 Nyquist 时仍然安全）
				while (phase >= 1.0f)
				{
					float frac = (1.0f - prevPhase) / incr;   // crossing position in this sample
					float tau = Clamp(frac, 0.0f, 0.999999f);
					modal.InjectImpulse(tau, 1.0f);

					phase -= 1.0f;
					prevPhase = 0.0f;
				}

				out[n] = modal.ProcessSample();
			}
		};

	auto BuildSpectrogramDb = [&](const std::vector<float>& signal, std::vector<float>& outDb,
		int& outFrames, int& outBins)
		{
			const int bins = STFTSize / 2 + 1;
			outBins = bins;

			if ((int)signal.size() < STFTSize)
			{
				outFrames = 0;
				outDb.clear();
				return;
			}

			// hopSize 越大，这里的 frames 越少，STFT 计算量越低
			const int frames = 1 + ((int)signal.size() - STFTSize) / hopSize;
			outFrames = frames;
			outDb.assign(frames * bins, DbMin);

			// Hann window
			static std::vector<float> win;
			if ((int)win.size() != STFTSize)
			{
				win.resize(STFTSize);
				for (int i = 0; i < STFTSize; ++i)
				{
					win[i] = 0.5f - 0.5f * cosf(2.0f * 3.14159265358979323846f * (float)i / (float)(STFTSize - 1));
				}
			}

			std::vector<float> re(STFTSize, 0.0f);
			std::vector<float> im(STFTSize, 0.0f);

			for (int frame = 0; frame < frames; ++frame)
			{
				int base = frame * hopSize;

				for (int i = 0; i < STFTSize; ++i)
				{
					re[i] = signal[base + i] * win[i];
					im[i] = 0.0f;
				}

				fft(re.data(), im.data(), STFTSize, 1);

				for (int k = 0; k < bins; ++k)
				{
					float mag = sqrtf(re[k] * re[k] + im[k] * im[k]);
					if (mag < 1.0e-20f) mag = 1.0e-20f;
					outDb[frame * bins + k] = 20.0f * log10f(mag);
				}
			}
		};

	auto DrawSpectrogram = [&](const std::vector<float>& specDb, int frames, int bins, const Rectangle& rc)
		{
			if (frames <= 0 || bins <= 1) return;

			// 不管 frames 有多少，都把它们拉伸铺满整个绘图区
			for (int px = 0; px < (int)rc.width; ++px)
			{
				float tx0 = (float)px / (float)rc.width;
				float tx1 = (float)(px + 1) / (float)rc.width;

				int f0 = (int)floorf(tx0 * frames);
				int f1 = (int)floorf(tx1 * frames);

				if (f0 < 0) f0 = 0;
				if (f0 >= frames) f0 = frames - 1;
				if (f1 <= f0) f1 = f0 + 1;
				if (f1 > frames) f1 = frames;

				for (int py = 0; py < (int)rc.height; ++py)
				{
					float ty0 = (float)py / (float)rc.height;
					float ty1 = (float)(py + 1) / (float)rc.height;

					// 线性频率轴：顶部 Nyquist，底部 0 Hz
					int k1 = (int)floorf((1.0f - ty0) * (bins - 1));
					int k0 = (int)floorf((1.0f - ty1) * (bins - 1));

					if (k0 < 0) k0 = 0;
					if (k1 < 0) k1 = 0;
					if (k0 >= bins) k0 = bins - 1;
					if (k1 >= bins) k1 = bins - 1;
					if (k1 < k0) std::swap(k0, k1);

					// 一个屏幕像素对应若干 STFT 帧时，取峰值，避免细节被平均抹掉
					float peakDb = DbMin;
					for (int fr = f0; fr < f1; ++fr)
					{
						for (int k = k0; k <= k1; ++k)
						{
							float v = specDb[fr * bins + k];
							if (v > peakDb) peakDb = v;
						}
					}

					Color c = ColorMapDb(peakDb);
					DrawRectangle((int)rc.x + px, (int)rc.y + py, 1, 1, c);
				}
			}
		};


	std::vector<float> nowParams;
	std::vector<float> nowDb;
	std::vector<float> sweepSignal;
	std::vector<float> specDb;
	int specFrames = 0;
	int specBins = 0;
	float lastLoss = 0.0f;
	int runCount = 0;

	while (!WindowShouldClose())
	{
		// 调 hopSize
		if (IsKeyPressed(KEY_ONE))   hopSize = 32;
		if (IsKeyPressed(KEY_TWO))   hopSize = 64;
		if (IsKeyPressed(KEY_THREE)) hopSize = 128;
		if (IsKeyPressed(KEY_FOUR))  hopSize = 256;

		lastLoss = opt.RunOptimizer();
		++runCount;

		opt.GetNowParams(nowParams);
		std::vector<float> drawParams = nowParams;
		opt.Regularization(drawParams);

		SystemModal nowModal;
		nowModal.CalcPoles(drawParams);

		// 左侧：多 tau 频响
		std::vector<std::vector<float>> multiTauDb;
		for (float tau = 0.0f; tau < 1.0f; tau += 0.1f)
		{
			multiTauDb.emplace_back();
			BuildSpectrumDb(nowModal, multiTauDb.back(), tau);
		}

		// 右侧：扫频冲激串 + STFT 瀑布
		BuildSweepSignal(nowModal, sweepSignal);
		BuildSpectrogramDb(sweepSignal, specDb, specFrames, specBins);

		BeginDrawing();
		ClearBackground(Color{ 18, 18, 22, 255 });

		// 左侧：多 tau 频响
		Rectangle leftPlot = { 90, 60, 760, 760 };
		Rectangle rightPlot = { 980, 60, 720, 760 };

		DrawGrid(leftPlot);
		DrawSpecGrid(rightPlot, SweepDurSec);

		for (int i = 0; i < (int)multiTauDb.size(); ++i)
		{
			float tau = 0.1f * (float)i;
			Color col = HSVtoRGB(tau, 1.0f, 1.0f);
			DrawSpectrum(multiTauDb[i], leftPlot, col);
		}

		// 右侧：扫频冲激串 + STFT 瀑布
		DrawSpectrogram(specDb, specFrames, specBins, rightPlot);

		DrawText("Static Response (multi-tau)", 90, 20, 28, RAYWHITE);
		DrawText("Sweep Waterfall / STFT", 980, 20, 28, RAYWHITE);

		char info[512];
		sprintf(info,
			"Run: %d   Loss: %.6f   Sweep: %.0f Hz -> %.0f Hz / %.1f s   STFT: %d   Hop: %d   Keys: [1]=32 [2]=64 [3]=128 [4]=256",
			runCount, lastLoss, SweepF0, SweepF1, SweepDurSec, STFTSize, hopSize);
		DrawText(info, 90, ScreenH - 40, 20, LIGHTGRAY);

		// 右下角画个 dB 色条
		{
			Rectangle cb = { rightPlot.x + rightPlot.width - 22, rightPlot.y + 20, 16, 220 };
			for (int i = 0; i < (int)cb.height; ++i)
			{
				float t = 1.0f - (float)i / (float)(cb.height - 1);
				float db = Lerp(DbMin, DbMax, t);
				DrawRectangle((int)cb.x, (int)(cb.y + i), (int)cb.width, 1, ColorMapDb(db));
			}
			DrawRectangleLinesEx(cb, 1.0f, LIGHTGRAY);

			char txt[32];
			sprintf(txt, "%.0f dB", DbMax);
			DrawText(txt, (int)cb.x - 70, (int)cb.y - 6, 18, LIGHTGRAY);
			sprintf(txt, "%.0f dB", DbMin);
			DrawText(txt, (int)cb.x - 70, (int)(cb.y + cb.height - 12), 18, LIGHTGRAY);
		}

		EndDrawing();
	}

	CloseWindow();
	return 0;
}
int main1()
{
	constexpr int ScreenW = 1800;
	constexpr int ScreenH = 960;
	constexpr float SampleRate = 48000.0f;

	// ---- Sweep / STFT 参数 ----
	constexpr float SweepDurSec = 6.0f;
	constexpr int SweepSamples = (int)(SweepDurSec * SampleRate);
	constexpr float SweepF0 = 2400.0f;
	constexpr float SweepF1 = 24000.0f;

	constexpr int STFTSize = 512;
	int hopSize = 64;
	constexpr float DbMin = -60.0f;
	constexpr float DbMax = 30.0f;

	InitWindow(ScreenW, ScreenH, "SystemModal Sweep Spectrogram");
	SetTargetFPS(60);

	auto Clamp = [](float x, float a, float b) -> float
		{
			if (x < a) return a;
			if (x > b) return b;
			return x;
		};

	auto Lerp = [](float a, float b, float t) -> float
		{
			return a + (b - a) * t;
		};

	auto DrawSpecGrid = [&](const Rectangle& rc, float durationSec)
		{
			DrawRectangleLinesEx(rc, 1.0f, GRAY);

			for (int i = 0; i <= 6; ++i)
			{
				float t = durationSec * (float)i / 6.0f;
				float x = rc.x + rc.width * ((float)i / 6.0f);
				DrawLine((int)x, (int)rc.y, (int)x, (int)(rc.y + rc.height), Fade(GRAY, 0.25f));

				char txt[32];
				sprintf(txt, "%.2fs", t);
				DrawText(txt, (int)x - 18, (int)(rc.y + rc.height + 8), 18, LIGHTGRAY);
			}

			const float freqMarks[] = { 0, 4000, 8000, 12000, 16000, 20000, 24000 };
			for (float f : freqMarks)
			{
				float y = rc.y + rc.height * (1.0f - f / 24000.0f);
				DrawLine((int)rc.x, (int)y, (int)(rc.x + rc.width), (int)y, Fade(GRAY, 0.25f));

				char txt[32];
				if (f >= 1000.0f) sprintf(txt, "%.0fk", f / 1000.0f);
				else sprintf(txt, "%.0f", f);
				DrawText(txt, (int)rc.x - 48, (int)y - 10, 18, LIGHTGRAY);
			}
		};

	auto ColorMapDb = [&](float db) -> Color
		{
			float t = (db - DbMin) / (DbMax - DbMin);
			t = Clamp(t, 0.0f, 1.0f);

			float h, s, v;
			if (t < 0.25f)
			{
				float u = t / 0.25f;
				h = Lerp(0.72f, 0.62f, u);
				s = 1.0f;
				v = Lerp(0.08f, 0.45f, u);
			}
			else if (t < 0.55f)
			{
				float u = (t - 0.25f) / 0.30f;
				h = Lerp(0.62f, 0.50f, u);
				s = 1.0f;
				v = Lerp(0.45f, 0.95f, u);
			}
			else if (t < 0.82f)
			{
				float u = (t - 0.55f) / 0.27f;
				h = Lerp(0.50f, 0.15f, u);
				s = Lerp(1.0f, 0.85f, u);
				v = 1.0f;
			}
			else
			{
				float u = (t - 0.82f) / 0.18f;
				h = Lerp(0.15f, 0.0f, u);
				s = Lerp(0.85f, 0.0f, u);
				v = 1.0f;
			}
			return HSVtoRGB(h, s, v);
		};

	auto BuildSweepSignal = [&](SystemModal& modal, std::vector<float>& out)
		{
			out.assign(SweepSamples, 0.0f);
			modal.Reset();

			float phase = 0.0f;

			for (int n = 0; n < SweepSamples; ++n)
			{
				float u = (float)n / (float)(SweepSamples - 1);

				// 指数扫频
				float freq = SweepF0 * powf(SweepF1 / SweepF0, u);
				float incr = freq / SampleRate;

				float prevPhase = phase;
				phase += incr;

				while (phase >= 1.0f)
				{
					float frac = (1.0f - prevPhase) / incr;
					float tau = Clamp(frac, 0.0f, 0.999999f);
					modal.InjectImpulse(tau, 1.0f);

					phase -= 1.0f;
					prevPhase = 0.0f;
				}

				out[n] = modal.ProcessSample();
			}
		};

	auto BuildSpectrogramDb = [&](const std::vector<float>& signal, std::vector<float>& outDb,
		int& outFrames, int& outBins)
		{
			const int bins = STFTSize / 2 + 1;
			outBins = bins;

			if ((int)signal.size() < STFTSize)
			{
				outFrames = 0;
				outDb.clear();
				return;
			}

			const int frames = 1 + ((int)signal.size() - STFTSize) / hopSize;
			outFrames = frames;
			outDb.assign(frames * bins, DbMin);

			static std::vector<float> win;
			if ((int)win.size() != STFTSize)
			{
				win.resize(STFTSize);
				for (int i = 0; i < STFTSize; ++i)
				{
					win[i] = 0.5f - 0.5f * cosf(2.0f * 3.14159265358979323846f * (float)i / (float)(STFTSize - 1));
				}
			}

			std::vector<float> re(STFTSize, 0.0f);
			std::vector<float> im(STFTSize, 0.0f);

			for (int frame = 0; frame < frames; ++frame)
			{
				int base = frame * hopSize;

				for (int i = 0; i < STFTSize; ++i)
				{
					re[i] = signal[base + i] * win[i];
					im[i] = 0.0f;
				}

				fft(re.data(), im.data(), STFTSize, 1);

				for (int k = 0; k < bins; ++k)
				{
					float mag = sqrtf(re[k] * re[k] + im[k] * im[k]);
					if (mag < 1.0e-20f) mag = 1.0e-20f;
					outDb[frame * bins + k] = 20.0f * log10f(mag);
				}
			}
		};

	auto DrawSpectrogram = [&](const std::vector<float>& specDb, int frames, int bins, const Rectangle& rc)
		{
			if (frames <= 0 || bins <= 1) return;

			for (int px = 0; px < (int)rc.width; ++px)
			{
				float tx0 = (float)px / (float)rc.width;
				float tx1 = (float)(px + 1) / (float)rc.width;

				int f0 = (int)floorf(tx0 * frames);
				int f1 = (int)floorf(tx1 * frames);

				if (f0 < 0) f0 = 0;
				if (f0 >= frames) f0 = frames - 1;
				if (f1 <= f0) f1 = f0 + 1;
				if (f1 > frames) f1 = frames;

				for (int py = 0; py < (int)rc.height; ++py)
				{
					float ty0 = (float)py / (float)rc.height;
					float ty1 = (float)(py + 1) / (float)rc.height;

					int k1 = (int)floorf((1.0f - ty0) * (bins - 1));
					int k0 = (int)floorf((1.0f - ty1) * (bins - 1));

					if (k0 < 0) k0 = 0;
					if (k1 < 0) k1 = 0;
					if (k0 >= bins) k0 = bins - 1;
					if (k1 >= bins) k1 = bins - 1;
					if (k1 < k0) std::swap(k0, k1);

					float peakDb = DbMin;
					for (int fr = f0; fr < f1; ++fr)
					{
						for (int k = k0; k <= k1; ++k)
						{
							float v = specDb[fr * bins + k];
							if (v > peakDb) peakDb = v;
						}
					}

					Color c = ColorMapDb(peakDb);
					DrawRectangle((int)rc.x + px, (int)rc.y + py, 1, 1, c);
				}
			}
		};

	// ------------------------------------------------------------
	// 这里直接写入你求得的“4组上半平面极点 + 对应留数”
	// 注意：不要把共轭负虚部那4组再放进来
	// ------------------------------------------------------------
	std::vector<float> filterParams =
	{
		-21063.386066f, 23870.9014559f, 0.0583185920643f, -0.528378346331f,
		-17374.1743263f, 67089.1514001f, -0.146510371181f, 0.410903309151f,
		-12151.9956308f, 99676.7704908f, 0.172563663704f, -0.246988256165f,
		-7411.46851302f, 120847.503268f, -0.143387503768f, 0.107064051987f,
		-3835.25322485f, 132839.027666f, 0.0847367893874f, -0.0220047341762f,
		-1170.4463834f, 138130.842319f, -0.02572248283f, -0.00271585828356f
	};


	SystemModal modal;
	modal.CalcPoles(filterParams);

	std::vector<float> sweepSignal;
	std::vector<float> specDb;
	int specFrames = 0;
	int specBins = 0;

	Rectangle plot = { 120, 60, 1500, 760 };

	// 先生成一次
	BuildSweepSignal(modal, sweepSignal);
	BuildSpectrogramDb(sweepSignal, specDb, specFrames, specBins);

	while (!WindowShouldClose())
	{
		if (IsKeyPressed(KEY_ONE))   hopSize = 32;
		if (IsKeyPressed(KEY_TWO))   hopSize = 64;
		if (IsKeyPressed(KEY_THREE)) hopSize = 128;
		if (IsKeyPressed(KEY_FOUR))  hopSize = 256;

		// 如果改了 hopSize，就重算一次 STFT
		if (IsKeyPressed(KEY_ONE) || IsKeyPressed(KEY_TWO) || IsKeyPressed(KEY_THREE) || IsKeyPressed(KEY_FOUR))
		{
			BuildSpectrogramDb(sweepSignal, specDb, specFrames, specBins);
		}

		BeginDrawing();
		ClearBackground(Color{ 18, 18, 22, 255 });

		DrawSpecGrid(plot, SweepDurSec);
		DrawSpectrogram(specDb, specFrames, specBins, plot);

		DrawText("Sweep Waterfall / STFT (Fixed Elliptic Poles + Residues)", 120, 20, 28, RAYWHITE);

		char info[512];
		sprintf(info,
			"Sweep: %.0f Hz -> %.0f Hz / %.1f s   STFT: %d   Hop: %d   Keys: [1]=32 [2]=64 [3]=128 [4]=256",
			SweepF0, SweepF1, SweepDurSec, STFTSize, hopSize);
		DrawText(info, 120, ScreenH - 40, 20, LIGHTGRAY);

		// 色条
		{
			Rectangle cb = { plot.x + plot.width - 22, plot.y + 20, 16, 220 };
			for (int i = 0; i < (int)cb.height; ++i)
			{
				float t = 1.0f - (float)i / (float)(cb.height - 1);
				float db = Lerp(DbMin, DbMax, t);
				DrawRectangle((int)cb.x, (int)(cb.y + i), (int)cb.width, 1, ColorMapDb(db));
			}
			DrawRectangleLinesEx(cb, 1.0f, LIGHTGRAY);

			char txt[32];
			sprintf(txt, "%.0f dB", DbMax);
			DrawText(txt, (int)cb.x - 70, (int)cb.y - 6, 18, LIGHTGRAY);
			sprintf(txt, "%.0f dB", DbMin);
			DrawText(txt, (int)cb.x - 70, (int)(cb.y + cb.height - 12), 18, LIGHTGRAY);
		}

		EndDrawing();
	}

	CloseWindow();
	return 0;
}

int main()
{
	main1();
}