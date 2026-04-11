#pragma once

#include <cmath>
namespace IIRBlepCoeffs
{
	constexpr static float Ts = 2.08333333333e-05f;
	constexpr static int NumTwoPoles = 6;
	constexpr static int NumOnePoles = 2;

	const float twoPoleParams[NumTwoPoles * 2] =
	{
		-22499.526025f, 25498.4629188f,
		-18558.7771212f, 71663.4117229f,
		-12980.5407875f, 106472.913933f,
		-7916.79591168f, 129087.105763f,
		-4096.74776281f, 141896.234098f,
		-1250.24954596f, 147548.854296f
	};

	const float onePoleParams[NumOnePoles] =
	{
		-282.842712475f,
		-565.685424949f
	};

	const float twoPoleBlitResidues[NumTwoPoles * 2] =
	{
		0.0832021015949f, -0.642350794072f,
		-0.181575534079f, 0.491896025511f,
		0.209339859972f, -0.294667076166f,
		-0.172759144034f, 0.127241233694f,
		0.101732366827f, -0.0257703058395f,
		-0.0308112980905f, -0.00343241885369f
	};
	const float onePoleBlitResidues[NumOnePoles] =
	{
		0.00597529294388f,
		-0.0242351636074f
	};

	const float twoPoleBlepResidues[NumTwoPoles * 2] =
	{
		-0.757563249836f, 0.511839212437f,
		0.338279599672f, 0.0340141486742f,
		-0.142233026f, -0.0770341618383f,
		0.0510612984665f, 0.0611075520522f,
		-0.00970294019698f, -0.034133415449f,
		-0.00103161367176f, 0.0100321489442f
	};
	const float onePoleBlepResidues[NumOnePoles] =
	{
		-0.0140408384458f,
		0.056421823596f
	};

	const float twoPoleBlampResidues[NumTwoPoles * 2] =
	{
		1.24922531497f, 0.32378651737f,
		-0.0336388444697f, -0.217867452746f,
		-0.0265170417089f, 0.0673541328449f,
		0.0214771305525f, -0.0203039054153f,
		-0.0114421920127f, 0.00361261810253f,
		0.00326622753975f, 0.000307924158163f
	};
	const float onePoleBlampResidues[NumOnePoles] =
	{
		2.38280929885f,
		-4.7875504886f
	};

	const float blitDirectGain = 2.33753844646e-11;
	const float blepDirectGain = 0;
	const float blampDirectGain = 0;
}

namespace IIRBlepUtils
{
	constexpr static int TableSize = 32;
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
