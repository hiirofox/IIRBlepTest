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
	constexpr static int NumOrders = 8;
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
		modal1.Reset();
		modal2.Reset();
		params2 = params;
		Regularization(params2);
		modal1.CalcPoles(params2);
		modal2.CalcPoles(params2);
		modal1.InjectImpulse(0.0, 1.0);
		modal2.InjectImpulse(0.5, 1.0);
		for (int i = 0; i < TestLen; ++i)
		{
			float x = (float)i / (float)TestLen;
			float w = (1.0 - x * x);
			float window = w * w;
			testre1[i] = modal1.ProcessSample() * window;
			testre2[i] = modal2.ProcessSample() * window;
			testim1[i] = 0;
			testim2[i] = 0;
		}
		fft(testre1, testim1, TestLen, 1);
		fft(testre2, testim2, TestLen, 1);
		float errv = 0;
		for (int i = 0; i < TestLen / 2; ++i)
		{
			float x = (float)i / (float)(TestLen / 2 - 1);
			float target = x > 0.75 ? 0.05 : 1.0;
			float mag1 = sqrtf(testre1[i] * testre1[i] + testim1[i] * testim1[i]);
			float mag2 = sqrtf(testre2[i] * testre2[i] + testim2[i] * testim2[i]);
			float e1 = fabsf(logf(mag1) - logf(target));
			float e2 = fabsf(logf(mag2) - logf(target));
			errv += e1 * e1 + e2 * e2;
		}
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
		optimizer.SetupOptimizer(NumOrders * 4, initParams, 0.01);
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
		for (int i = 0; i < NumOrders; i++) {
			float pre = input[i * 4 + 0];
			float pim = input[i * 4 + 1];
			float rre = input[i * 4 + 2];
			float rim = input[i * 4 + 3];
			input[i * 4 + 0] = -pre * pre * 1000.0;
			input[i * 4 + 1] = pim * 24000.0;
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