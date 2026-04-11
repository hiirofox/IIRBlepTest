#define _CRT_SECURE_NO_WARNINGS 1
#include <stdio.h>
#include <vector>
#include <cmath>
#include <algorithm>

#include "raylib/src/raylib.h"

#include "dsp/optimizer.h"
#include "dsp/IIRBlep.h"

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

int main1()
{
	constexpr int ScreenW = 1800;
	constexpr int ScreenH = 960;
	constexpr float SampleRate = 48000.0f;

	// ---- Sweep / STFT 参数 ----
	constexpr float SweepDurSec = 20.0f;
	constexpr int SweepSamples = (int)(SweepDurSec * SampleRate);
	constexpr float SweepF0 = 1.0f;
	constexpr float SweepF1 = 24000.0f;
	constexpr int ResponseSamples = 200;
	constexpr float ResponseTau = 0.0f;

	constexpr int STFTSize = 1024;
	int hopSize = 128;
	constexpr float DbMin = -30.0f;
	constexpr float DbMax = 10.0;

	InitWindow(ScreenW, ScreenH, "IIRBlep2 Sweep Spectrogram");
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

	auto BuildSweepSignal = [&](IIRBlep2::IIRBlep& modal, std::vector<float>& out)
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

				phase += incr;

				while (phase >= 1.0f)
				{
					phase -= (int)phase;
					float frac = phase / incr;
					float tau = Clamp(frac, 0.0f, 0.999999f);
					modal.Add(-1.0f, tau, 1);

				}
				float naive = phase - 0.5;
				modal.Step();
				float residue = modal.Get();
				out[n] = naive + residue;
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

	auto BuildResidualResponse = [&](std::vector<float>& out, int mode)
		{
			out.assign(ResponseSamples, 0.0f);
			IIRBlep2::IIRBlep modal;
			modal.Reset();
			modal.Add(1.0f, ResponseTau, mode);

			for (int n = 0; n < ResponseSamples; ++n)
			{
				modal.Step();
				out[n] = modal.Get();
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

	auto DrawResponsePlot = [&](const std::vector<float>& signal, const Rectangle& rc, const char* title, Color lineColor)
		{
			DrawRectangleLinesEx(rc, 1.0f, GRAY);

			if (signal.empty()) return;

			const float viewAbs = 2.0f;
			const float yScale = rc.height * 0.42f / viewAbs;
			float maxAbs = 1.0e-12f;
			for (float v : signal) {
				maxAbs = std::max(maxAbs, fabsf(v));
			}

			const float centerY = rc.y + rc.height * 0.5f;
			DrawLine((int)rc.x, (int)centerY, (int)(rc.x + rc.width), (int)centerY, Fade(GRAY, 0.45f));

			auto DrawDashedLevel = [&](float level, Color color)
				{
					float y = centerY - level * yScale;
					if (y < rc.y || y > rc.y + rc.height) return;

					const float dashLen = 10.0f;
					const float gapLen = 6.0f;
					for (float x = rc.x; x < rc.x + rc.width; x += dashLen + gapLen)
					{
						float x2 = std::min(x + dashLen, rc.x + rc.width);
						DrawLineEx(Vector2{ x, y }, Vector2{ x2, y }, 1.0f, color);
					}
				};

			DrawDashedLevel(1.0f, Fade(GREEN, 0.7f));
			DrawDashedLevel(-1.0f, Fade(RED, 0.7f));

			for (int i = 0; i <= 5; ++i)
			{
				float t = (float)i / 5.0f;
				float x = rc.x + rc.width * t;
				DrawLine((int)x, (int)rc.y, (int)x, (int)(rc.y + rc.height), Fade(GRAY, 0.18f));

				char txt[32];
				sprintf(txt, "%d", (int)roundf(t * (ResponseSamples - 1)));
				DrawText(txt, (int)x - 10, (int)(rc.y + rc.height + 6), 16, LIGHTGRAY);
			}

			float prevX = rc.x;
			float prevY = centerY - signal[0] * yScale;
			for (size_t i = 1; i < signal.size(); ++i)
			{
				float t = (float)i / (float)(signal.size() - 1);
				float x = rc.x + rc.width * t;
				float y = centerY - signal[i] * yScale;
				DrawLineEx(Vector2{ prevX, prevY }, Vector2{ x, y }, 1.75f, lineColor);
				prevX = x;
				prevY = y;
			}

			char info[128];
			sprintf(info, "%s  |  N=%d  |  range=[-2,2]  |  peak=%.4g", title, ResponseSamples, maxAbs);
			DrawText(info, (int)rc.x + 10, (int)rc.y + 8, 18, RAYWHITE);
		};

	// ------------------------------------------------------------
	// 这里直接写入你求得的“4组上半平面极点 + 对应留数”
	// 注意：不要把共轭负虚部那4组再放进来
	// ------------------------------------------------------------
	IIRBlep2::IIRBlep blep;

	std::vector<float> sweepSignal;
	std::vector<float> specDb;
	std::vector<float> blitResponse;
	std::vector<float> blepResponse;
	std::vector<float> blampResponse;
	int specFrames = 0;
	int specBins = 0;

	Rectangle plot = { 120, 60, 1500, 500 };
	const float responseGap = 20.0f;
	const float responseY = plot.y + plot.height + 40.0f;
	const float responseH = ScreenH - responseY - 76;
	const float responseW = (plot.width - responseGap * 2.0f) / 3.0f;
	Rectangle blitPlot = { plot.x, responseY, responseW, responseH };
	Rectangle blepPlot = { plot.x + responseW + responseGap, responseY, responseW, responseH };
	Rectangle blampPlot = { plot.x + (responseW + responseGap) * 2.0f, responseY, responseW, responseH };

	BuildSweepSignal(blep, sweepSignal);
	BuildSpectrogramDb(sweepSignal, specDb, specFrames, specBins);
	BuildResidualResponse(blitResponse, IIRBlepUtils::BLIT_MODE);
	BuildResidualResponse(blepResponse, IIRBlepUtils::BLEP_MODE);
	BuildResidualResponse(blampResponse, IIRBlepUtils::BLAMP_MODE);

	while (!WindowShouldClose())
	{

		BeginDrawing();
		ClearBackground(Color{ 18, 18, 22, 255 });

		DrawSpecGrid(plot, SweepDurSec);
		DrawSpectrogram(specDb, specFrames, specBins, plot);
		DrawResponsePlot(blitResponse, blitPlot, "BLIT residue: Add(0, 1, 0)", SKYBLUE);
		DrawResponsePlot(blepResponse, blepPlot, "BLEP residue: Add(0, 1, 1)", ORANGE);
		DrawResponsePlot(blampResponse, blampPlot, "BLAMP residue: Add(0, 1, 2)", LIME);

		DrawText("Sweep Waterfall / STFT (IIRBlep2 Shared Poles + Residue Modes)", 120, 20, 28, RAYWHITE);

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
