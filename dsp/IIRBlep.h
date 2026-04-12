#pragma once

#include <cmath>

namespace IIRBlepCoeffs
{
	constexpr static float Ts = 2.08333333333e-05f;
	constexpr static int NumTwoPoles = 6;
	constexpr static int NumOnePoles = 2;

	const float twoPoleParams[NumTwoPoles * 2] =
	{
		-62.8318530718f, 108.827961854f,
		-40189.0272933f, 55793.4533008f,
		-28445.4287912f, 98631.7684816f,
		-16966.8559996f, 125439.49809f,
		-8508.87829805f, 139822.403584f,
		-2542.24706254f, 145902.310941f
	};

	const float onePoleParams[NumOnePoles] =
	{
		-125.663706144f,
		-45507.231287f
	};

	const float twoPoleBlitResidues[NumTwoPoles * 2] =
	{
		-0.0013168263647f, -0.0007512013723f,
		-1.23826777933f, -0.548626744807f,
		0.585097569385f, 0.718909860941f,
		-0.0622758708595f, -0.511084153629f,
		-0.115077178164f, 0.192724755046f,
		0.05683897699f, -0.0233731868357f
	};
	const float onePoleBlitResidues[NumOnePoles] =
	{
		-0.0026337014076f,
		1.55304393381f
	};

	const float twoPoleBlepResidues[NumTwoPoles * 2] =
	{
		0.00299989646564f, 0.00172127986961f,
		0.194464203781f, 0.925225508737f,
		0.247184301023f, -0.356030793163f,
		-0.188889699096f, 0.0493792323782f,
		0.0683120874018f, 0.0353480213983f,
		-0.00801287006699f, -0.0185596799848f
	};
	const float onePoleBlepResidues[NumOnePoles] =
	{
		0.00599983435271f,
		-1.63811567337f
	};

	const float twoPoleBlampResidues[NumTwoPoles * 2] =
	{
		-0.00354322399043f, -1.3210980025f,
		0.444726650264f, -0.487645264048f,
		-0.191989693844f, -0.0649244901803f,
		0.0281565385978f, 0.0684710777029f,
		0.0106680972746f, -0.0241002418088f,
		-0.00605812600442f, 0.00274169074962f
	};
	const float onePoleBlampResidues[NumOnePoles] =
	{
		-2.29176790792f,
		1.72784742332f
	};

	const float blitDirectGain = 0;
	const float blepDirectGain = 0;
	const float blampDirectGain = 0;
}

namespace IIRBlepUtils
{
	constexpr static int TableSize = 64;
	constexpr static int BLIT_MODE = 0;
	constexpr static int BLEP_MODE = 1;
	constexpr static int BLAMP_MODE = 2;

	static bool isTableBuilt = false;

	static float twoPoleBlitG1Table[IIRBlepCoeffs::NumTwoPoles][TableSize];
	static float twoPoleBlitG2Table[IIRBlepCoeffs::NumTwoPoles][TableSize];
	static float onePoleBlitG1Table[IIRBlepCoeffs::NumOnePoles][TableSize];

	static float twoPoleBlepG1Table[IIRBlepCoeffs::NumTwoPoles][TableSize];
	static float twoPoleBlepG2Table[IIRBlepCoeffs::NumTwoPoles][TableSize];
	static float onePoleBlepG1Table[IIRBlepCoeffs::NumOnePoles][TableSize];

	static float twoPoleBlampG1Table[IIRBlepCoeffs::NumTwoPoles][TableSize];
	static float twoPoleBlampG2Table[IIRBlepCoeffs::NumTwoPoles][TableSize];
	static float onePoleBlampG1Table[IIRBlepCoeffs::NumOnePoles][TableSize];

	inline float LerpTable(const float table[TableSize], int index1, int index2, float frac)
	{
		return table[index1] + frac * (table[index2] - table[index1]);
	}

	inline void BuildTwoPoleTable(
		int poleIndex,
		const float residues[],
		float g1Table[IIRBlepCoeffs::NumTwoPoles][TableSize],
		float g2Table[IIRBlepCoeffs::NumTwoPoles][TableSize])
	{
		const float pre = IIRBlepCoeffs::twoPoleParams[poleIndex * 2 + 0];
		const float pim = IIRBlepCoeffs::twoPoleParams[poleIndex * 2 + 1];
		const float rre = residues[poleIndex * 2 + 0];
		const float rim = residues[poleIndex * 2 + 1];
		const float R = expf(pre * IIRBlepCoeffs::Ts);
		const float O = pim * IIRBlepCoeffs::Ts;
		const float a1 = -2.0f * R * cosf(O);
		const float stepRe = R * cosf(O);
		const float stepIm = R * sinf(O);

		for (int i = 0; i < TableSize; ++i) {
			const float tau = (float)i / (float)(TableSize - 1);
			const float dt1 = tau * IIRBlepCoeffs::Ts;
			const float shiftAbs = expf(pre * dt1);
			const float shiftArg = pim * dt1;
			const float shiftRe = shiftAbs * cosf(shiftArg);
			const float shiftIm = shiftAbs * sinf(shiftArg);
			const float A1Re = rre * shiftRe - rim * shiftIm;
			const float A1Im = rre * shiftIm + rim * shiftRe;
			const float A2Re = A1Re * stepRe - A1Im * stepIm;
			const float g1 = 2.0f * A1Re;
			const float g2 = 2.0f * A2Re;
			g1Table[poleIndex][i] = g1;
			g2Table[poleIndex][i] = g2 + a1 * g1;
		}
	}

	inline void BuildOnePoleTable(
		int poleIndex,
		const float residues[],
		float g1Table[IIRBlepCoeffs::NumOnePoles][TableSize])
	{
		const float pre = IIRBlepCoeffs::onePoleParams[poleIndex];
		const float residue = residues[poleIndex];

		for (int i = 0; i < TableSize; ++i) {
			const float tau = (float)i / (float)(TableSize - 1);
			const float dt1 = tau * IIRBlepCoeffs::Ts;
			g1Table[poleIndex][i] = residue * expf(pre * dt1);
		}
	}

	inline void BuildTables()
	{
		if (isTableBuilt) return;

		for (int i = 0; i < IIRBlepCoeffs::NumTwoPoles; ++i) {
			BuildTwoPoleTable(i, IIRBlepCoeffs::twoPoleBlitResidues, twoPoleBlitG1Table, twoPoleBlitG2Table);
			BuildTwoPoleTable(i, IIRBlepCoeffs::twoPoleBlepResidues, twoPoleBlepG1Table, twoPoleBlepG2Table);
			BuildTwoPoleTable(i, IIRBlepCoeffs::twoPoleBlampResidues, twoPoleBlampG1Table, twoPoleBlampG2Table);
		}

		for (int i = 0; i < IIRBlepCoeffs::NumOnePoles; ++i) {
			BuildOnePoleTable(i, IIRBlepCoeffs::onePoleBlitResidues, onePoleBlitG1Table);
			BuildOnePoleTable(i, IIRBlepCoeffs::onePoleBlepResidues, onePoleBlepG1Table);
			BuildOnePoleTable(i, IIRBlepCoeffs::onePoleBlampResidues, onePoleBlampG1Table);
		}

		isTableBuilt = true;
	}
}

namespace IIRBlep2
{
	struct TwoPoleModal
	{
		float a1 = 0.0f, a2 = 0.0f;
		float z1 = 0.0f, z2 = 0.0f;

		inline float ProcessSample()
		{
			float y = z1;
			z1 = -a1 * y + z2;
			z2 = -a2 * y;
			return y;
		}

		void CalcPole(float pre, float pim)
		{
			float R = expf(pre * IIRBlepCoeffs::Ts);
			float O = pim * IIRBlepCoeffs::Ts;
			a1 = -2.0f * R * cosf(O);
			a2 = R * R;
		}

		inline void InjectEvent(float g1, float g2)
		{
			z1 += g1;
			z2 += g2;
		}

		void Reset()
		{
			z1 = 0.0f;
			z2 = 0.0f;
		}
	};

	struct OnePoleModal
	{
		float a1 = 0.0f;
		float z1 = 0.0f;

		inline float ProcessSample()
		{
			float y = z1;
			z1 = -a1 * y;
			return y;
		}

		void CalcPole(float pre)
		{
			a1 = -expf(pre * IIRBlepCoeffs::Ts);
		}

		inline void InjectEvent(float g1)
		{
			z1 += g1;
		}

		void Reset()
		{
			z1 = 0.0f;
		}
	};

	class IIRBlep
	{
	private:
		TwoPoleModal twoPoles[IIRBlepCoeffs::NumTwoPoles];
		OnePoleModal onePoles[IIRBlepCoeffs::NumOnePoles];
		float v = 0.0f;

	public:
		IIRBlep()
		{
			for (int i = 0; i < IIRBlepCoeffs::NumTwoPoles; ++i) {
				const float pre = IIRBlepCoeffs::twoPoleParams[i * 2 + 0];
				const float pim = IIRBlepCoeffs::twoPoleParams[i * 2 + 1];
				twoPoles[i].CalcPole(pre, pim);
			}

			for (int i = 0; i < IIRBlepCoeffs::NumOnePoles; ++i) {
				onePoles[i].CalcPole(IIRBlepCoeffs::onePoleParams[i]);
			}

			IIRBlepUtils::BuildTables();
		}

		void Add(float linear_gain, float tau, int mode = 1)
		{
			if (mode < IIRBlepUtils::BLIT_MODE || mode > IIRBlepUtils::BLAMP_MODE) return;

			if (tau < 0.0f) tau = 0.0f;
			if (tau >= 1.0f) tau = 0.999999999999f;

			const float fpos = tau * (float)(IIRBlepUtils::TableSize - 1);
			const int index1 = (int)fpos;
			const int index2 = index1 + 1;
			const float frac = fpos - (float)index1;

			const float(*twoPoleG1Table)[IIRBlepUtils::TableSize] = IIRBlepUtils::twoPoleBlitG1Table;
			const float(*twoPoleG2Table)[IIRBlepUtils::TableSize] = IIRBlepUtils::twoPoleBlitG2Table;
			const float(*onePoleG1Table)[IIRBlepUtils::TableSize] = IIRBlepUtils::onePoleBlitG1Table;

			if (mode == IIRBlepUtils::BLEP_MODE) {
				twoPoleG1Table = IIRBlepUtils::twoPoleBlepG1Table;
				twoPoleG2Table = IIRBlepUtils::twoPoleBlepG2Table;
				onePoleG1Table = IIRBlepUtils::onePoleBlepG1Table;
			}
			else if (mode == IIRBlepUtils::BLAMP_MODE) {
				twoPoleG1Table = IIRBlepUtils::twoPoleBlampG1Table;
				twoPoleG2Table = IIRBlepUtils::twoPoleBlampG2Table;
				onePoleG1Table = IIRBlepUtils::onePoleBlampG1Table;
			}

			for (int i = 0; i < IIRBlepCoeffs::NumTwoPoles; ++i) {
				const float g1 = IIRBlepUtils::LerpTable(twoPoleG1Table[i], index1, index2, frac) * linear_gain;
				const float g2 = IIRBlepUtils::LerpTable(twoPoleG2Table[i], index1, index2, frac) * linear_gain;
				twoPoles[i].InjectEvent(g1, g2);
			}

			for (int i = 0; i < IIRBlepCoeffs::NumOnePoles; ++i) {
				const float g1 = IIRBlepUtils::LerpTable(onePoleG1Table[i], index1, index2, frac) * linear_gain;
				onePoles[i].InjectEvent(g1);
			}
		}

		void Step()
		{
			float y = 0.0f;
			for (int i = 0; i < IIRBlepCoeffs::NumTwoPoles; ++i) {
				y += twoPoles[i].ProcessSample();
			}
			for (int i = 0; i < IIRBlepCoeffs::NumOnePoles; ++i) {
				y += onePoles[i].ProcessSample();
			}
			v = y;
		}

		float Get()
		{
			return v;
		}

		void Reset()
		{
			for (int i = 0; i < IIRBlepCoeffs::NumTwoPoles; ++i) {
				twoPoles[i].Reset();
			}
			for (int i = 0; i < IIRBlepCoeffs::NumOnePoles; ++i) {
				onePoles[i].Reset();
			}
			v = 0.0f;
		}
	};

}

namespace HybridBlep
{
	namespace HybridBlepFirCoeffs
	{
		constexpr static int FirSize = 16;
		const float firCoeffs[FirSize] = { 1.0, 0, };
	}

	class HybridBlep
	{
	private:
		IIRBlep2::TwoPoleModal twoPoles[IIRBlepCoeffs::NumTwoPoles];
		IIRBlep2::OnePoleModal onePoles[IIRBlepCoeffs::NumOnePoles];
		float v = 0.0f;
		float twoPoleG1States[IIRBlepCoeffs::NumTwoPoles][HybridBlepFirCoeffs::FirSize] = { 0 };
		float twoPoleG2States[IIRBlepCoeffs::NumTwoPoles][HybridBlepFirCoeffs::FirSize] = { 0 };
		float onePoleG1States[IIRBlepCoeffs::NumOnePoles][HybridBlepFirCoeffs::FirSize] = { 0 };
		int firpos = 0;
	public:
		HybridBlep()
		{
			for (int i = 0; i < IIRBlepCoeffs::NumTwoPoles; ++i) {
				const float pre = IIRBlepCoeffs::twoPoleParams[i * 2 + 0];
				const float pim = IIRBlepCoeffs::twoPoleParams[i * 2 + 1];
				twoPoles[i].CalcPole(pre, pim);
			}

			for (int i = 0; i < IIRBlepCoeffs::NumOnePoles; ++i) {
				onePoles[i].CalcPole(IIRBlepCoeffs::onePoleParams[i]);
			}

			IIRBlepUtils::BuildTables();
			Reset();
		}

		void Add(float linear_gain, float tau, int mode = 1)
		{
			if (mode < IIRBlepUtils::BLIT_MODE || mode > IIRBlepUtils::BLAMP_MODE) return;

			if (tau < 0.0f) tau = 0.0f;
			if (tau >= 1.0f) tau = 0.999999999999f;

			const float fpos = tau * (float)(IIRBlepUtils::TableSize - 1);
			const int index1 = (int)fpos;
			const int index2 = index1 + 1;
			const float frac = fpos - (float)index1;

			const float(*twoPoleG1Table)[IIRBlepUtils::TableSize] = IIRBlepUtils::twoPoleBlitG1Table;
			const float(*twoPoleG2Table)[IIRBlepUtils::TableSize] = IIRBlepUtils::twoPoleBlitG2Table;
			const float(*onePoleG1Table)[IIRBlepUtils::TableSize] = IIRBlepUtils::onePoleBlitG1Table;

			if (mode == IIRBlepUtils::BLEP_MODE) {
				twoPoleG1Table = IIRBlepUtils::twoPoleBlepG1Table;
				twoPoleG2Table = IIRBlepUtils::twoPoleBlepG2Table;
				onePoleG1Table = IIRBlepUtils::onePoleBlepG1Table;
			}
			else if (mode == IIRBlepUtils::BLAMP_MODE) {
				twoPoleG1Table = IIRBlepUtils::twoPoleBlampG1Table;
				twoPoleG2Table = IIRBlepUtils::twoPoleBlampG2Table;
				onePoleG1Table = IIRBlepUtils::onePoleBlampG1Table;
			}

			for (int i = 0; i < IIRBlepCoeffs::NumTwoPoles; ++i) {
				const float g1 = IIRBlepUtils::LerpTable(twoPoleG1Table[i], index1, index2, frac) * linear_gain;
				const float g2 = IIRBlepUtils::LerpTable(twoPoleG2Table[i], index1, index2, frac) * linear_gain;
				//twoPoles[i].InjectEvent(g1, g2);
				for (int j = 0; j < HybridBlepFirCoeffs::FirSize; ++j) {
					int index = (firpos + j) % HybridBlepFirCoeffs::FirSize;
					twoPoleG1States[i][index] += g1 * HybridBlepFirCoeffs::firCoeffs[j];
					twoPoleG2States[i][index] += g2 * HybridBlepFirCoeffs::firCoeffs[j];
				}
			}

			for (int i = 0; i < IIRBlepCoeffs::NumOnePoles; ++i) {
				const float g1 = IIRBlepUtils::LerpTable(onePoleG1Table[i], index1, index2, frac) * linear_gain;
				//onePoles[i].InjectEvent(g1);
				for (int j = 0; j < HybridBlepFirCoeffs::FirSize; ++j) {
					int index = (firpos + j) % HybridBlepFirCoeffs::FirSize;
					onePoleG1States[i][index] += g1 * HybridBlepFirCoeffs::firCoeffs[j];
				}
			}
		}

		void Step()
		{
			for (int i = 0; i < IIRBlepCoeffs::NumTwoPoles; ++i) {//inject fir state
				float g1 = twoPoleG1States[i][firpos];
				float g2 = twoPoleG2States[i][firpos];
				twoPoles[i].InjectEvent(g1, g2);
				twoPoleG1States[i][firpos] = 0.0f;
				twoPoleG2States[i][firpos] = 0.0f;
			}
			for (int i = 0; i < IIRBlepCoeffs::NumOnePoles; ++i) {
				float g1 = onePoleG1States[i][firpos];
				onePoles[i].InjectEvent(g1);
				onePoleG1States[i][firpos] = 0.0f;
			}
			firpos++;
			if (firpos >= HybridBlepFirCoeffs::FirSize) firpos = 0;

			float y = 0.0f;
			for (int i = 0; i < IIRBlepCoeffs::NumTwoPoles; ++i) {
				y += twoPoles[i].ProcessSample();
			}
			for (int i = 0; i < IIRBlepCoeffs::NumOnePoles; ++i) {
				y += onePoles[i].ProcessSample();
			}
			v = y;
		}

		float Get()
		{
			return v;
		}

		void Reset()
		{
			for (int i = 0; i < IIRBlepCoeffs::NumTwoPoles; ++i) {
				twoPoles[i].Reset();
			}
			for (int i = 0; i < IIRBlepCoeffs::NumOnePoles; ++i) {
				onePoles[i].Reset();
			}
			v = 0.0f;

			for (int i = 0; i < IIRBlepCoeffs::NumTwoPoles; ++i) {
				for (int j = 0; j < HybridBlepFirCoeffs::FirSize; ++j) {
					twoPoleG1States[i][j] = 0.0f;
					twoPoleG2States[i][j] = 0.0f;
				}
			}
			for (int i = 0; i < IIRBlepCoeffs::NumOnePoles; ++i) {
				for (int j = 0; j < HybridBlepFirCoeffs::FirSize; ++j) {
					onePoleG1States[i][j] = 0.0f;
				}
			}
		}
	};

}
