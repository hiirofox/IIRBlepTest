#define _CRT_SECURE_NO_WARNINGS 1
#include <stdio.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <string>

#include "raylib/src/raylib.h"

#include "dsp/IIRBlep.h"
#include "dsp/TableBlep.h"

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
		float angle = -inv * 2.0f * 3.1415926535897932384626f / m;
		float wm_re = cosf(angle);
		float wm_im = sinf(angle);
		for (k = 0; k < len; k += m) {
			float w_re = 1.0f;
			float w_im = 0.0f;

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

static Color LerpColor(Color a, Color b, float t)
{
	t = std::max(0.0f, std::min(1.0f, t));
	return Color{
		(unsigned char)roundf((float)a.r + ((float)b.r - (float)a.r) * t),
		(unsigned char)roundf((float)a.g + ((float)b.g - (float)a.g) * t),
		(unsigned char)roundf((float)a.b + ((float)b.b - (float)a.b) * t),
		(unsigned char)roundf((float)a.a + ((float)b.a - (float)a.a) * t)
	};
}

struct BlepEvalResult
{
	const char* name = "";
	float freqHz = 0.0f;
	int sampleCount = 0;
	float aliasSuppressionDb = 0.0f;
	float generationMs = 0.0f;
};

static constexpr float EvalSampleRate = 48000.0f;
static constexpr float EvalPi = 3.1415926535897932384626f;

static bool IsPowerOfTwo(int x)
{
	return x > 0 && (x & (x - 1)) == 0;
}

static float IdealIntegratedDiracSaw(float phase, int harmonicCount)
{
	double sum = 0.0;
	for (int h = 1; h <= harmonicCount; ++h)
	{
		sum -= sin(2.0 * (double)EvalPi * (double)h * (double)phase) / ((double)EvalPi * (double)h);
	}
	return (float)sum;
}

static void RemoveMean(std::vector<float>& x)
{
	if (x.empty()) return;

	double sum = 0.0;
	for (float v : x) sum += v;
	const float mean = (float)(sum / (double)x.size());
	for (float& v : x) v -= mean;
}

static void BuildMagnitudeSpectrum(std::vector<float> signal, std::vector<float>& mag)
{
	RemoveMean(signal);

	const int n = (int)signal.size();
	std::vector<float> re(n, 0.0f);
	std::vector<float> im(n, 0.0f);

	for (int i = 0; i < n; ++i)
	{
		const float w = 0.5f - 0.5f * cosf(2.0f * EvalPi * (float)i / (float)(n - 1));
		re[i] = signal[i] * w;
	}

	fft(re.data(), im.data(), n, 1);

	const int bins = n / 2 + 1;
	mag.assign(bins, 0.0f);
	for (int k = 0; k < bins; ++k)
	{
		mag[k] = sqrtf(re[k] * re[k] + im[k] * im[k]);
	}
}

static float MagnitudeErrorSuppressionDb(const std::vector<float>& refMag, const std::vector<float>& testMag)
{
	const int bins = (int)std::min(refMag.size(), testMag.size());
	double refEnergy = 0.0;
	double errEnergy = 0.0;

	for (int k = 1; k + 1 < bins; ++k)
	{
		const double ref = refMag[k];
		const double diff = (double)testMag[k] - ref;
		refEnergy += ref * ref;
		errEnergy += diff * diff;
	}

	if (refEnergy <= 0.0) return 0.0f;
	errEnergy = std::max(errEnergy, refEnergy * 1.0e-16);
	return (float)(10.0 * log10(refEnergy / errEnergy));
}

template<typename Blep>
BlepEvalResult Eval(float freqHz, int sampleCount, const char* name = "Blep")
{
	BlepEvalResult result;
	result.name = name;
	result.freqHz = freqHz;
	result.sampleCount = sampleCount;

	if (!IsPowerOfTwo(sampleCount) || freqHz <= 0.0f || freqHz >= EvalSampleRate * 0.5f)
	{
		result.aliasSuppressionDb = 0.0f;
		result.generationMs = 0.0f;
		return result;
	}

	const float phaseInc = freqHz / EvalSampleRate;
	const int harmonicCount = std::max(1, (int)floorf((EvalSampleRate * 0.5f) / freqHz));
	const int warmupSamples = std::max(4096, (int)ceilf(EvalSampleRate / freqHz) * 16);

	Blep blep;
	blep.Reset();
	float phase = 0.0f;

	auto advanceBlep = [&](float& phaseState, Blep& blepState) -> float
		{
			phaseState += phaseInc;
			while (phaseState >= 1.0f)
			{
				phaseState -= 1.0f;
				const float where = std::max(0.0f, std::min(0.999999f, phaseState / phaseInc));
				blepState.Add(-1.0f, where, BLEP_MODE);
			}

			const float naive = phaseState - 0.5f;
			blepState.Step();
			return naive + blepState.Get();
		};

	for (int i = 0; i < warmupSamples; ++i)
	{
		(void)advanceBlep(phase, blep);
	}

	std::vector<float> blepSignal(sampleCount, 0.0f);
	const auto t0 = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < sampleCount; ++i)
	{
		blepSignal[i] = advanceBlep(phase, blep);
	}
	const auto t1 = std::chrono::high_resolution_clock::now();
	result.generationMs = std::chrono::duration<float, std::milli>(t1 - t0).count();

	std::vector<float> refSignal(sampleCount, 0.0f);
	float refPhase = 0.0f;
	for (int i = 0; i < warmupSamples; ++i)
	{
		refPhase += phaseInc;
		refPhase -= floorf(refPhase);
	}
	for (int i = 0; i < sampleCount; ++i)
	{
		refPhase += phaseInc;
		refPhase -= floorf(refPhase);
		refSignal[i] = IdealIntegratedDiracSaw(refPhase, harmonicCount);
	}

	std::vector<float> refMag;
	std::vector<float> testMag;
	BuildMagnitudeSpectrum(refSignal, refMag);
	BuildMagnitudeSpectrum(blepSignal, testMag);
	result.aliasSuppressionDb = MagnitudeErrorSuppressionDb(refMag, testMag);
	return result;
}

template<typename Blep>
int main1()
{
	constexpr int ScreenW = 1800;
	constexpr int ScreenH = 960;
	constexpr float SampleRate = 48000.0f;

	constexpr float SweepDurSec = 20.0f;
	constexpr int SweepSamples = (int)(SweepDurSec * SampleRate);
	constexpr float SweepF0 = 480.0f;
	constexpr float SweepF1 = 48000.0f;
	constexpr int ResponseSamples = 120;
	constexpr float ResponseTau = 0.0f;

	constexpr int STFTSize = 1024;
	int hopSize = 128;
	constexpr float DbMin = -40.0f;
	constexpr float DbMax = 10.0f;

	InitWindow(ScreenW, ScreenH, "BLEP Residual Diagnostics");
	SetTargetFPS(60);

	Font uiFont = GetFontDefault();
	bool hasCustomFont = false;
	if (FileExists("C:/Windows/Fonts/segoeui.ttf"))
	{
		uiFont = LoadFontEx("C:/Windows/Fonts/segoeui.ttf", 32, nullptr, 0);
		hasCustomFont = uiFont.texture.id != 0;
	}

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

	const Color Bg = Color{ 18, 18, 22, 255 };
	const Color PlotBorder = Color{ 150, 150, 156, 255 };
	const Color GridColor = Color{ 128, 128, 134, 255 };
	const Color TextMain = Color{ 238, 238, 242, 255 };
	const Color TextMuted = Color{ 178, 178, 188, 255 };
	const Color DashedRef = Color{ 172, 172, 180, 255 };

	auto DrawLabel = [&](const char* text, float x, float y, float size, Color color)
		{
			DrawTextEx(uiFont, text, Vector2{ x, y }, size, 1.0f, color);
		};

	auto MeasureLabel = [&](const char* text, float size) -> Vector2
		{
			return MeasureTextEx(uiFont, text, size, 1.0f);
		};

	auto DrawRightLabel = [&](const char* text, float rightX, float y, float size, Color color)
		{
			Vector2 m = MeasureLabel(text, size);
			DrawLabel(text, rightX - m.x, y, size, color);
		};

	auto DrawDashedLine = [&](Vector2 a, Vector2 b, float dashLen, float gapLen, float thick, Color color)
		{
			float dx = b.x - a.x;
			float dy = b.y - a.y;
			float len = sqrtf(dx * dx + dy * dy);
			if (len <= 0.0f) return;

			float ux = dx / len;
			float uy = dy / len;
			for (float d = 0.0f; d < len; d += dashLen + gapLen)
			{
				float d2 = std::min(d + dashLen, len);
				DrawLineEx(
					Vector2{ a.x + ux * d, a.y + uy * d },
					Vector2{ a.x + ux * d2, a.y + uy * d2 },
					thick,
					color);
			}
		};

	auto ColorMapDb = [&](float db) -> Color
		{
			float t = (db - DbMin) / (DbMax - DbMin);
			t = Clamp(t, 0.0f, 1.0f);

			const Color stops[] = {
				Color{ 0, 0, 0, 255 },
				Color{ 5, 12, 45, 255 },
				Color{ 26, 57, 158, 255 },
				Color{ 107, 55, 190, 255 },
				Color{ 205, 46, 105, 255 },
				Color{ 247, 156, 44, 255 },
				Color{ 255, 242, 186, 255 },
				Color{ 255, 255, 255, 255 },
			};
			constexpr int StopCount = sizeof(stops) / sizeof(stops[0]);

			float scaled = t * (float)(StopCount - 1);
			int index = (int)floorf(scaled);
			if (index < 0) index = 0;
			if (index >= StopCount - 1) index = StopCount - 2;
			float frac = scaled - (float)index;
			return LerpColor(stops[index], stops[index + 1], frac);
		};

	auto BuildSweepSignal = [&](Blep& modal, std::vector<float>& out)
		{
			out.assign(SweepSamples, 0.0f);
			modal.Reset();

			float phase = 0.0f;

			for (int n = 0; n < SweepSamples; ++n)
			{
				float u = (float)n / (float)(SweepSamples - 1);
				float freq = SweepF0 * powf(SweepF1 / SweepF0, u);
				float incr = freq / SampleRate;

				phase += incr;

				while (phase >= 1.0f)
				{
					phase -= (int)phase;
					float frac = phase / incr;
					float where = Clamp(frac, 0.0f, 0.999999f);
					modal.Add(-1.0f, where, BLEP_MODE);
				}

				float naive = phase - 0.5f;
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
			Blep modal;
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

			BeginScissorMode((int)rc.x, (int)rc.y, (int)rc.width, (int)rc.height);
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
			EndScissorMode();
		};

	auto DrawSpecGrid = [&](const Rectangle& rc, float durationSec)
		{
			DrawRectangleLinesEx(rc, 1.0f, PlotBorder);

			for (int i = 0; i <= 5; ++i)
			{
				float t = durationSec * (float)i / 5.0f;
				float x = rc.x + rc.width * ((float)i / 5.0f);
				DrawLineEx(Vector2{ x, rc.y }, Vector2{ x, rc.y + rc.height }, 1.0f, Fade(GridColor, 0.24f));

				char txt[32];
				sprintf(txt, "%.0fs", t);
				Vector2 m = MeasureLabel(txt, 14.0f);
				DrawLabel(txt, x - m.x * 0.5f, rc.y + rc.height + 8.0f, 14.0f, TextMuted);
			}

			const float freqMarks[] = { 0, 4000, 8000, 12000, 16000, 20000, 24000 };
			for (float f : freqMarks)
			{
				float y = rc.y + rc.height * (1.0f - f / 24000.0f);
				DrawLineEx(Vector2{ rc.x, y }, Vector2{ rc.x + rc.width, y }, 1.0f, Fade(GridColor, 0.24f));

				char txt[32];
				if (f >= 1000.0f) sprintf(txt, "%.0fk", f / 1000.0f);
				else sprintf(txt, "%.0f", f);
				DrawRightLabel(txt, rc.x - 10.0f, y - 7.0f, 14.0f, TextMuted);
			}
		};

	auto DrawColorBar = [&](const Rectangle& rc)
		{
			for (int i = 0; i < (int)rc.height; ++i)
			{
				float t = 1.0f - (float)i / (float)(rc.height - 1);
				float db = Lerp(DbMin, DbMax, t);
				DrawRectangle((int)rc.x, (int)(rc.y + i), (int)rc.width, 1, ColorMapDb(db));
			}

			DrawRectangleLinesEx(rc, 1.0f, PlotBorder);
			char txt[32];
			DrawLabel("dB", rc.x, rc.y - 22.0f, 14.0f, TextMuted);
			sprintf(txt, "%.0f", DbMax);
			DrawLabel(txt, rc.x + rc.width + 8.0f, rc.y - 4.0f, 14.0f, TextMuted);
			sprintf(txt, "%.0f", DbMin);
			DrawLabel(txt, rc.x + rc.width + 8.0f, rc.y + rc.height - 12.0f, 14.0f, TextMuted);
		};

	auto DrawResponsePlot = [&](const std::vector<float>& signal, const Rectangle& rc,
		const char* title, const char* subtitle, Color lineColor)
		{
			DrawRectangleLinesEx(rc, 1.0f, PlotBorder);
			DrawLabel(title, rc.x, rc.y - 24.0f, 18.0f, TextMain);
			DrawLabel(subtitle, rc.x + 128.0f, rc.y - 22.0f, 14.0f, TextMuted);

			if (signal.empty()) return;

			float maxAbs = 1.0e-12f;
			for (float v : signal) {
				maxAbs = std::max(maxAbs, fabsf(v));
			}

			char peakText[64];
			sprintf(peakText, "peak %.4g", maxAbs);
			DrawRightLabel(peakText, rc.x + rc.width, rc.y - 22.0f, 14.0f, TextMuted);

			const float viewAbs = 2.0f;
			const float yScale = rc.height * 0.5f / viewAbs;
			const float centerY = rc.y + rc.height * 0.5f;

			for (int i = 0; i <= 4; ++i)
			{
				float t = (float)i / 4.0f;
				float x = rc.x + rc.width * t;
				DrawLineEx(Vector2{ x, rc.y }, Vector2{ x, rc.y + rc.height }, 1.0f, Fade(GridColor, 0.17f));

				char txt[32];
				sprintf(txt, "%d", (int)roundf(t * (ResponseSamples - 1)));
				Vector2 m = MeasureLabel(txt, 12.0f);
				DrawLabel(txt, x - m.x * 0.5f, rc.y + rc.height + 4.0f, 12.0f, TextMuted);
			}

			for (int i = -2; i <= 2; ++i)
			{
				float level = (float)i;
				float y = centerY - level * yScale;
				Color grid = (i == 0) ? Fade(DashedRef, 0.50f) : Fade(GridColor, 0.16f);
				DrawLineEx(Vector2{ rc.x, y }, Vector2{ rc.x + rc.width, y }, 1.0f, grid);
				if (i == -1 || i == 0 || i == 1)
				{
					char txt[16];
					if (i == 0) sprintf(txt, "0");
					else sprintf(txt, "%+d", i);
					DrawRightLabel(txt, rc.x - 8.0f, y - 7.0f, 12.0f, TextMuted);
				}
			}

			auto DrawDashedLevel = [&](float level)
				{
					float y = centerY - level * yScale;
					if (y < rc.y || y > rc.y + rc.height) return;
					DrawDashedLine(Vector2{ rc.x, y }, Vector2{ rc.x + rc.width, y }, 9.0f, 6.0f, 1.5f, Fade(DashedRef, 0.72f));
				};

			DrawDashedLevel(1.0f);
			DrawDashedLevel(-1.0f);

			BeginScissorMode((int)rc.x, (int)rc.y, (int)rc.width, (int)rc.height);
			float prevX = rc.x;
			float prevY = centerY - signal[0] * yScale;
			for (size_t i = 1; i < signal.size(); ++i)
			{
				float t = (float)i / (float)(signal.size() - 1);
				float x = rc.x + rc.width * t;
				float y = centerY - signal[i] * yScale;
				DrawLineEx(Vector2{ prevX, prevY }, Vector2{ x, y }, 2.0f, lineColor);
				prevX = x;
				prevY = y;
			}
			EndScissorMode();
		};

	Blep blep;

	std::vector<float> sweepSignal;
	std::vector<float> specDb;
	std::vector<float> blitResponse;
	std::vector<float> blepResponse;
	std::vector<float> blampResponse;
	int specFrames = 0;
	int specBins = 0;

	const float pageMargin = 40.0f;
	const float top = 90.0f;
	const float bottom = 58.0f;
	const float gap = 48.0f;
	const float contentH = (float)ScreenH - top - bottom;
	const float responseW = 520.0f;
	const float responseH = (contentH - gap * 2.0f) / 3.0f;
	const float responseX = (float)ScreenW - pageMargin - responseW;
	const float sweepX = pageMargin + 56.0f;
	Rectangle sweepPlot = {
		sweepX,
		top,
		responseX - gap - 58.0f - sweepX,
		contentH
	};
	Rectangle colorBar = {
		sweepPlot.x + sweepPlot.width + 20.0f,
		sweepPlot.y,
		14.0f,
		220.0f
	};

	Rectangle blitPlot = { responseX, top, responseW, responseH };
	Rectangle blepPlot = { responseX, top + responseH + gap, responseW, responseH };
	Rectangle blampPlot = { responseX, top + (responseH + gap) * 2.0f, responseW, responseH };

	BuildSweepSignal(blep, sweepSignal);
	BuildSpectrogramDb(sweepSignal, specDb, specFrames, specBins);
	BuildResidualResponse(blitResponse, BLIT_MODE);
	BuildResidualResponse(blepResponse, BLEP_MODE);
	BuildResidualResponse(blampResponse, BLAMP_MODE);

	while (!WindowShouldClose())
	{
		int newHopSize = hopSize;
		if (IsKeyPressed(KEY_ONE)) newHopSize = 32;
		if (IsKeyPressed(KEY_TWO)) newHopSize = 64;
		if (IsKeyPressed(KEY_THREE)) newHopSize = 128;
		if (IsKeyPressed(KEY_FOUR)) newHopSize = 256;
		if (newHopSize != hopSize)
		{
			hopSize = newHopSize;
			BuildSpectrogramDb(sweepSignal, specDb, specFrames, specBins);
		}

		BeginDrawing();
		ClearBackground(Bg);

		DrawLabel("Sweep Waterfall / STFT", sweepPlot.x, 26.0f, 24.0f, TextMain);
		DrawLabel("BLIT / BLEP / BLAMP Residual Responses", responseX, 26.0f, 20.0f, TextMain);
		DrawLabel("log sweep validation, single-frame raylib plots", sweepPlot.x, 56.0f, 15.0f, TextMuted);

		DrawSpectrogram(specDb, specFrames, specBins, sweepPlot);
		DrawSpecGrid(sweepPlot, SweepDurSec);
		DrawColorBar(colorBar);

		DrawResponsePlot(blitResponse, blitPlot, "BLIT residual", "single event, stage 0", Color{ 91, 194, 255, 255 });
		DrawResponsePlot(blepResponse, blepPlot, "BLEP residual", "single event, stage 1", Color{ 255, 145, 82, 255 });
		DrawResponsePlot(blampResponse, blampPlot, "BLAMP residual", "single event, stage 2", Color{ 255, 216, 93, 255 });

		char info[512];
		sprintf(info,
			"Sweep %.0f Hz -> %.0f Hz / %.1f s   STFT %d   Hop %d   Keys: [1]32 [2]64 [3]128 [4]256",
			SweepF0, SweepF1, SweepDurSec, STFTSize, hopSize);
		DrawLabel(info, pageMargin, ScreenH - 26.0f, 15.0f, TextMuted);

		EndDrawing();
	}

	if (hasCustomFont)
	{
		UnloadFont(uiFont);
	}
	CloseWindow();
	return 0;
}

int main2()
{
	constexpr int ScreenW = 1500;
	constexpr int ScreenH = 900;
	constexpr float EvalFreqHz = 5000.0f;
	constexpr int EvalSamples = 131072*16;

	std::vector<BlepEvalResult> results;
	results.push_back(Eval<IIRBlep2::IIRBlep>(EvalFreqHz, EvalSamples, "IIRModal"));
	results.push_back(Eval<TableBlep>(EvalFreqHz, EvalSamples, "Table"));
	results.push_back(Eval<Lagrange2thBlep>(EvalFreqHz, EvalSamples, "Lag2"));
	results.push_back(Eval<Lagrange4thBlep>(EvalFreqHz, EvalSamples, "Lag4"));
	results.push_back(Eval<Bspline6thBlep>(EvalFreqHz, EvalSamples, "Bspline6"));

	for (const BlepEvalResult& r : results)
	{
		printf("%-10s  alias suppression = %8.3f dB, generation = %8.4f ms\n",
			r.name, r.aliasSuppressionDb, r.generationMs);
	}

	InitWindow(ScreenW, ScreenH, "BLEP Benchmark");
	SetTargetFPS(60);

	Font uiFont = GetFontDefault();
	bool hasCustomFont = false;
	if (FileExists("C:/Windows/Fonts/segoeui.ttf"))
	{
		uiFont = LoadFontEx("C:/Windows/Fonts/segoeui.ttf", 32, nullptr, 0);
		hasCustomFont = uiFont.texture.id != 0;
	}

	const Color bg = Color{ 18, 18, 22, 255 };
	const Color border = Color{ 150, 150, 156, 255 };
	const Color grid = Color{ 128, 128, 134, 255 };
	const Color text = Color{ 238, 238, 242, 255 };
	const Color muted = Color{ 178, 178, 188, 255 };
	const Color bars[] = {
		Color{ 91, 194, 255, 255 },
		Color{ 255, 145, 82, 255 },
		Color{ 255, 216, 93, 255 },
		Color{ 142, 217, 135, 255 },
		Color{ 218, 130, 255, 255 },
	};

	auto DrawLabel = [&](const char* label, float x, float y, float size, Color color)
		{
			DrawTextEx(uiFont, label, Vector2{ x, y }, size, 1.0f, color);
		};

	auto MeasureLabel = [&](const char* label, float size) -> Vector2
		{
			return MeasureTextEx(uiFont, label, size, 1.0f);
		};

	auto DrawRightLabel = [&](const char* label, float rightX, float y, float size, Color color)
		{
			Vector2 m = MeasureLabel(label, size);
			DrawLabel(label, rightX - m.x, y, size, color);
		};

	auto NiceCeil = [](float x) -> float
		{
			if (x <= 0.0f) return 1.0f;
			const float base = powf(10.0f, floorf(log10f(x)));
			const float n = x / base;
			if (n <= 1.0f) return 1.0f * base;
			if (n <= 2.0f) return 2.0f * base;
			if (n <= 5.0f) return 5.0f * base;
			return 10.0f * base;
		};

	auto DrawBarChart = [&](const Rectangle& rc, const char* title, const char* unit, bool useAliasDb)
		{
			float maxValue = 0.0f;
			for (const BlepEvalResult& r : results)
			{
				const float v = useAliasDb ? r.aliasSuppressionDb : r.generationMs;
				maxValue = std::max(maxValue, v);
			}
			maxValue = NiceCeil(maxValue * 1.08f);

			DrawLabel(title, rc.x, rc.y - 34.0f, 24.0f, text);
			DrawRectangleLinesEx(rc, 1.0f, border);

			for (int i = 0; i <= 5; ++i)
			{
				const float t = (float)i / 5.0f;
				const float value = maxValue * t;
				const float y = rc.y + rc.height * (1.0f - t);
				DrawLineEx(Vector2{ rc.x, y }, Vector2{ rc.x + rc.width, y }, 1.0f, Fade(grid, 0.25f));

				char label[64];
				sprintf(label, "%.3g %s", value, unit);
				DrawRightLabel(label, rc.x - 10.0f, y - 8.0f, 14.0f, muted);
			}

			const int count = (int)results.size();
			const float cellW = rc.width / (float)count;
			const float barW = std::min(110.0f, cellW * 0.45f);

			for (int i = 0; i < count; ++i)
			{
				const BlepEvalResult& r = results[i];
				const float value = useAliasDb ? r.aliasSuppressionDb : r.generationMs;
				const float norm = maxValue > 0.0f ? std::max(0.0f, value / maxValue) : 0.0f;
				const float barH = rc.height * norm;
				const float x = rc.x + cellW * ((float)i + 0.5f) - barW * 0.5f;
				const float y = rc.y + rc.height - barH;
				const Color c = bars[i % (int)(sizeof(bars) / sizeof(bars[0]))];

				DrawRectangleRec(Rectangle{ x, y, barW, barH }, c);
				DrawRectangleLinesEx(Rectangle{ x, y, barW, barH }, 1.0f, Fade(RAYWHITE, 0.6f));

				char valueText[64];
				if (useAliasDb) sprintf(valueText, "%.2f dB", value);
				else sprintf(valueText, "%.3f ms", value);
				Vector2 valueSize = MeasureLabel(valueText, 15.0f);
				DrawLabel(valueText, x + barW * 0.5f - valueSize.x * 0.5f, y - 22.0f, 15.0f, text);

				Vector2 nameSize = MeasureLabel(r.name, 16.0f);
				DrawLabel(r.name, x + barW * 0.5f - nameSize.x * 0.5f, rc.y + rc.height + 10.0f, 16.0f, muted);
			}
		};

	while (!WindowShouldClose())
	{
		BeginDrawing();
		ClearBackground(bg);

		DrawLabel("BLEP Benchmark", 70.0f, 28.0f, 30.0f, text);
		char subtitle[256];
		sprintf(subtitle, "%.0f Hz saw, %d samples, %.0f Hz sample rate", EvalFreqHz, EvalSamples, EvalSampleRate);
		DrawLabel(subtitle, 70.0f, 66.0f, 17.0f, muted);

		DrawBarChart(Rectangle{ 150.0f, 145.0f, 1260.0f, 280.0f },
			"Alias suppression, higher is better", "dB", true);
		DrawBarChart(Rectangle{ 150.0f, 560.0f, 1260.0f, 220.0f },
			"Generation time , lower is better", "ms", false);

		DrawLabel("Metric: Hann-windowed FFT magnitude error vs ideal integrated-Dirac bandlimited saw; phase is discarded, DC/Nyquist bins ignored.",
			70.0f, 840.0f, 15.0f, muted);

		EndDrawing();
	}

	if (hasCustomFont)
	{
		UnloadFont(uiFont);
	}
	CloseWindow();
	return 0;
}

int main()
{
	//main1<IIRBlep2::IIRBlep>();
	return main2();
}
