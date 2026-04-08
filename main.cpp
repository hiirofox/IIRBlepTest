#define _CRT_SECURE_NO_WARNINGS 1
#include <stdio.h>
#include <vector>
#include <cmath>
#include <algorithm>

#include "raylib/src/raylib.h"
#include "dsp/optimizer.h"

struct PoleModal
{
	const float Ts = 1.0f / 48000.0f;
	float a1 = 0, a2 = 0;
	float z1 = 0, z2 = 0;
	float pre = 0, pim = 0, rre = 0, rim = 0;
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
	}
	void InjectImpulse(float tau, float v)
	{
		if (tau < 0.0f) tau = 0.0f;
		if (tau >= 1.0f) tau = 1.0f;
		float t1 = (1.0f - tau) * Ts;
		float t2 = (2.0f - tau) * Ts;
		float e1 = expf(pre * t1);
		float e2 = expf(pre * t2);
		float c1 = cosf(pim * t1);
		float s1 = sinf(pim * t1);
		float c2 = cosf(pim * t2);
		float s2 = sinf(pim * t2);
		float g1 = 2.0f * v * e1 * (rre * c1 - rim * s1);
		float g2 = 2.0f * v * e2 * (rre * c2 - rim * s2);
		float dz1 = g1;
		float dz2 = g2 + a1 * g1;
		z1 += dz1;
		z2 += dz2;
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
	constexpr static int NumOrders = 4;
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

	constexpr static int TestLen = 128;
	float testre1[TestLen];
	float testim1[TestLen];
	float testre2[TestLen];
	float testim2[TestLen];
	float lastloss = -1;
	std::vector<float> params2;
	float Error(std::vector<float>& params)
	{
		// ---------- 可调参数 ----------
		constexpr float Fs = 48000.0f;
		constexpr float M_PI = 3.14159265358979323846f;

		// 训练时频谱长度要远大于原来的 128，否则你看到的大量东西只是截断/窗造成的
		constexpr int FFTLen = 32;

		// 多个 tau 一起训练，覆盖 [0, 1)
		constexpr int TauCount = 17;

		// 目标低通边界（按 Nyquist 归一化）
		// pass <= 0.45*Nyquist, stop >= 0.49*Nyquist
		constexpr float PassEdge = 0.45f;
		constexpr float StopEdge = 0.49f;

		// 各项 loss 权重
		constexpr float WPassFlat = 1.0f;   // 平均频响通带平坦
		constexpr float WStop = 2.0f;   // 阻带抑制
		constexpr float WConsistency = 40.0f;   // 不同 tau 去相位后的复频响一致
		constexpr float WGain = 2.0f;   // 通带平均增益锚定在 1
		constexpr float WTail = 0.02f;  // 尾部能量，减轻“靠超慢极点拖尾骗频响”

		auto sqr = [](float x) -> float { return x * x; };
		auto clamp01 = [](float x) -> float
			{
				if (x < 0.0f) return 0.0f;
				if (x > 1.0f) return 1.0f;
				return x;
			};
		auto smoothstep = [&](float a, float b, float x) -> float
			{
				float t = clamp01((x - a) / (b - a));
				return t * t * (3.0f - 2.0f * t);
			};

		// 先做参数规整，再算极点
		params2 = params;
		Regularization(params2);
		modal1.CalcPoles(params2);

		// 保存每个 tau 的“去掉分数延时后的复频响”
		std::vector<float> meanRe(FFTLen / 2 + 1, 0.0f);
		std::vector<float> meanIm(FFTLen / 2 + 1, 0.0f);

		std::vector<float> allRe(TauCount * (FFTLen / 2 + 1), 0.0f);
		std::vector<float> allIm(TauCount * (FFTLen / 2 + 1), 0.0f);

		float tailLoss = 0.0f;

		for (int kt = 0; kt < TauCount; ++kt)
		{
			float tau = (float)kt / (float)TauCount; // [0,1)

			std::vector<float> re(FFTLen, 0.0f);
			std::vector<float> im(FFTLen, 0.0f);

			modal1.Reset();
			modal1.InjectImpulse(tau, 1.0f);

			for (int n = 0; n < FFTLen; ++n)
			{
				float y = modal1.ProcessSample();
				re[n] = y;
				im[n] = 0.0f;

				// 惩罚尾部能量，避免用极慢衰减模态 + 截断 FFT 伪造频响
				if (n >= FFTLen / 4)
					tailLoss += y * y;
			}

			fft(re.data(), im.data(), FFTLen, 1);

			// 你的离散输出从“下一个采样点”开始，
			// 所以不同 tau 之间应只差 delay d = 1 - tau 的线性相位
			float d = 1.0f - tau;

			for (int k = 0; k <= FFTLen / 2; ++k)
			{
				float w = 2.0f * PI * (float)k / (float)FFTLen;
				float c = cosf(w * d);
				float s = sinf(w * d);

				// 去掉 e^{-jwd}，即乘 e^{+jwd}
				float rr = re[k] * c - im[k] * s;
				float ii = re[k] * s + im[k] * c;

				allRe[kt * (FFTLen / 2 + 1) + k] = rr;
				allIm[kt * (FFTLen / 2 + 1) + k] = ii;

				meanRe[k] += rr;
				meanIm[k] += ii;
			}
		}

		for (int k = 0; k <= FFTLen / 2; ++k)
		{
			meanRe[k] /= (float)TauCount;
			meanIm[k] /= (float)TauCount;
		}

		float passFlatLoss = 0.0f;
		float stopLoss = 0.0f;
		float consistencyLoss = 0.0f;
		float gainLoss = 0.0f;

		int passCount = 0;
		int stopCount = 0;

		for (int k = 0; k <= FFTLen / 2; ++k)
		{
			float fn = (float)k / (float)(FFTLen / 2); // 0..1 对应 0..Nyquist

			// 平滑低通模板，而不是硬台阶
			float stopMix = smoothstep(PassEdge, StopEdge, fn); // 0=pass, 1=stop
			float targetMag = 1.0f - stopMix;

			float mr = meanRe[k];
			float mi = meanIm[k];
			float mm = sqrtf(mr * mr + mi * mi + 1.0e-30f);

			if (fn <= PassEdge)
			{
				// 通带平坦且接近 1
				passFlatLoss += sqr(mm - 1.0f);
				gainLoss += sqr(mm - 1.0f);
				++passCount;
			}
			else if (fn >= StopEdge)
			{
				// 阻带能量要小
				stopLoss += mm * mm;
				++stopCount;
			}
			else
			{
				// 过渡带用软目标
				passFlatLoss += 0.25f * sqr(mm - targetMag);
			}

			// 不同 tau 去掉线性相位后，复频响应一致
			for (int kt = 0; kt < TauCount; ++kt)
			{
				float rr = allRe[kt * (FFTLen / 2 + 1) + k];
				float ii = allIm[kt * (FFTLen / 2 + 1) + k];

				float dr = rr - mr;
				float di = ii - mi;

				// 通带里更强调一致性；阻带里稍微弱一点
				float w = (fn <= PassEdge) ? 1.0f : 0.2f;
				consistencyLoss += w * (dr * dr + di * di);
			}
		}

		if (passCount > 0)
		{
			passFlatLoss /= (float)passCount;
			gainLoss /= (float)passCount;
		}
		if (stopCount > 0)
			stopLoss /= (float)stopCount;

		consistencyLoss /= (float)(TauCount * (FFTLen / 2 + 1));
		tailLoss /= (float)(TauCount * FFTLen);

		float errv =
			WPassFlat * passFlatLoss +
			WStop * stopLoss +
			WConsistency * consistencyLoss +
			WGain * gainLoss +
			WTail * tailLoss;

		lastloss = errv;
		return errv;
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
		optimizer.SetupOptimizer(NumOrders * 4, initParams, 1);
		optimizer.SetErrorFunc([this](std::vector<float>& params) { return Error(params); });
	}
	float RunOptimizer()
	{
		optimizer.RunOptimizer(10);
		return lastloss;
	}
	void GetNowParams(std::vector<float>& params)
	{
		optimizer.GetNowVec(params);
	}



	void Regularization(std::vector<float>& input)
	{
		const float Fs = 48000.0f;
		const float M_PI = 3.14159265358979323846f;
		const float NyquistRad = M_PI * Fs; // 因为 O = pim * Ts，所以 pim 必须是 rad/s

		auto sigmoid = [](float x) -> float
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
			};

		auto softplus = [](float x) -> float
			{
				if (x > 20.0f) return x;
				if (x < -20.0f) return expf(x);
				return logf(1.0f + expf(x));
			};

		for (int i = 0; i < NumOrders; i++)
		{
			float rawSigma = input[i * 4 + 0];
			float rawOmega = input[i * 4 + 1];
			float rawRre = input[i * 4 + 2];
			float rawRim = input[i * 4 + 3];

			// 1) 实部必须 < 0，且不要太接近 0，否则靠长拖尾 + 截断频谱作弊
			//    这里下界给一个最小衰减，再留足动态范围
			float sigma = -(80.0f + 12000.0f * softplus(rawSigma));

			// 2) 虚部映射到 (0, 0.98*Nyquist) 的 rad/s
			//    你原来 *24000 的量纲不对；CalcPole 里 O = pim*Ts，所以 pim 应该是 rad/s
			float omega = 0.98f * NyquistRad * sigmoid(rawOmega);

			// 3) 输出模态系数做有界化，避免数值爆炸
			//    总幅度按阶数缩放，让 Adam 更容易收敛
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
int main()
{
	constexpr int ScreenW = 1600;
	constexpr int ScreenH = 900;
	constexpr int FFTLen = 8192;
	constexpr float SampleRate = 48000.0;

	InitWindow(ScreenW, ScreenH, "Modal Optimizer Frequency Response");
	SetTargetFPS(60);

	ModalLPFOptimizer opt;


	auto BuildSpectrumDb = [&](SystemModal& modal, std::vector<float>& outDb, float tau)
		{
			static float re[FFTLen];
			static float im[FFTLen];

			modal.Reset();
			modal.InjectImpulse(tau, 1.0f);/////!

			for (int i = 0; i < FFTLen; ++i)
			{
				float x = (float)i / (float)(FFTLen - 1);
				float w = (1.0f - x * x);
				float window = w * w;
				re[i] = modal.ProcessSample() * window;
				im[i] = 0.0f;
			}

			fft(re, im, FFTLen, 1);

			outDb.resize(FFTLen / 2 + 1);
			for (int i = 0; i <= FFTLen / 2; ++i)
			{
				float mag = sqrtf(re[i] * re[i] + im[i] * im[i]);
				if (mag < 1.0e-20f) mag = 1.0e-20f;
				outDb[i] = 20.0f * log10f(mag);
			}
		};

	auto FreqToX = [&](float f, const Rectangle& rc) -> float
		{
			float fmin = 20.0f;
			float fmax = 24000.0f;
			if (f < fmin) f = fmin;
			if (f > fmax) f = fmax;

			float t = (log10f(f) - log10f(fmin)) / (log10f(fmax) - log10f(fmin));
			return rc.x + t * rc.width;
		};

	auto DbToY = [&](float db, const Rectangle& rc) -> float
		{
			float dbMin = -30.0f;
			float dbMax = 30.0f;
			if (db < dbMin) db = dbMin;
			if (db > dbMax) db = dbMax;

			float t = (db - dbMin) / (dbMax - dbMin);
			return rc.y + rc.height * (1.0f - t);
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

				float bin0 = f0 * (float)FFTLen / SampleRate;
				float bin1 = f1 * (float)FFTLen / SampleRate;

				int i0 = (int)bin0;
				int i1 = (int)bin1;

				if (i0 < 0) i0 = 0;
				if (i0 > halfN) i0 = halfN;
				if (i1 < 0) i1 = 0;
				if (i1 > halfN) i1 = halfN;

				float y0 = DbToY(db[i0], rc);
				float y1 = DbToY(db[i1], rc);
				DrawLineV({ rc.x + (float)px, y0 }, { rc.x + (float)px + 1.0f, y1 }, color);
			}
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
				if (f >= 1000.0f)
					sprintf(txt, "%.0fk", f / 1000.0f);
				else
					sprintf(txt, "%.0f", f);

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

	std::vector<float> nowParams;
	std::vector<float> targetDb;
	std::vector<float> nowDb;
	float lastLoss = 0.0f;
	int runCount = 0;

	while (!WindowShouldClose())
	{
		// 每 run 一次优化，就重算一次频响并绘制
		lastLoss = opt.RunOptimizer();
		++runCount;

		opt.GetNowParams(nowParams);

		// 按优化时的约束方式把 pre 转成稳定极点
		std::vector<float> drawParams = nowParams;
		opt.Regularization(drawParams);

		SystemModal nowModal;
		nowModal.CalcPoles(drawParams);


		BeginDrawing();
		ClearBackground(Color{ 18, 18, 22, 255 });

		Rectangle plot = { 100, 60, (float)ScreenW - 160.0f, (float)ScreenH - 160.0f };
		DrawGrid(plot);

		// 参考/目标曲线
		DrawSpectrum(targetDb, plot, Color{ 80, 170, 255, 255 });

		// 当前训练结果曲线
		for (float tau = 0.0; tau < 1.0; tau += 0.1)
		{
			BuildSpectrumDb(nowModal, nowDb, tau);
			auto col = HSVtoRGB(tau, 1.0f, 1.0f);
			DrawSpectrum(nowDb, plot, col);
		}

		DrawText("Frequency Response", 100, 20, 28, RAYWHITE);

		char info[256];
		sprintf(info, "Run: %d    Loss: %.6f    X: 20 ~ 48000 Hz (log)    Y: -30 ~ 30 dB",
			runCount, lastLoss);
		DrawText(info, 100, ScreenH - 42, 20, LIGHTGRAY);

		EndDrawing();
	}

	CloseWindow();
	return 0;
}