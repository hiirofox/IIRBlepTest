#pragma once

#include <cmath>
#include <complex>
#include <tuple>
#include <vector>

struct TwoPoleModal
{
	const float Ts = 1.0f / 48000.0f;

	float a1 = 0, a2 = 0;
	float z1 = 0, z2 = 0;

	float pre = 0, pim = 0, rre = 0, rim = 0;

	std::complex<float> pole = { 0, 0 };
	std::complex<float> residue = { 0, 0 };
	std::complex<float> step1 = { 1, 0 };

	constexpr static int TableSize = 32;
	float g1Table[TableSize];
	float g2Table[TableSize];

	float impNormGain = 1.0;
	float stepNormGain = 1.0;

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

		for (int i = 0; i < TableSize; i++) {
			float tau = (float)i / (TableSize - 2);
			float dt1 = (1.0f - tau) * Ts;
			std::complex<float> shift = std::exp(pole * dt1);
			std::complex<float> A1 = residue * shift;
			std::complex<float> A2 = A1 * step1;
			float g1 = 2.0f * A1.real() * impNormGain;
			float g2 = 2.0f * A2.real() * impNormGain;
			g1Table[i] = g1;
			g2Table[i] = g2 + a1 * g1;
		}
	}
	/*
	void InjectImpulse(float tau, float v)
	{
		if (tau < 0.0f) tau = 0.0f;
		if (tau >= 1.0f) tau = 1.0f;
		float dt1 = (1.0f - tau) * Ts;
		std::complex<float> shift = std::exp(pole * dt1);
		std::complex<float> A1 = v * residue * shift;
		std::complex<float> A2 = A1 * step1;
		float g1 = 2.0f * A1.real() * impNormGain;
		float g2 = 2.0f * A2.real() * impNormGain;
		z1 += g1;
		z2 += g2 + a1 * g1;
	}
	*/
	//void InjectImpulse(float tau, float v)
	inline void InjectImpulse(int index1, float frac, float v)
	{
		//if (tau < 0.0f) tau = 0.0f;
		//if (tau >= 1.0f) tau = 1.0f;
		//float fpos = tau * (TableSize - 2);
		//int index1 = fpos;
		int index2 = index1 + 1;
		//float frac = fpos - index1;
		float g1 = frac * (g1Table[index2] - g1Table[index1]) + g1Table[index1];
		float g2 = frac * (g2Table[index2] - g2Table[index1]) + g2Table[index1];
		z1 += g1 * v;
		z2 += g2 * v;
	}

	void SetNormGain(float impNormGain, float stepNormGain)
	{
		this->impNormGain = impNormGain;
		this->stepNormGain = stepNormGain;
	}

	void Reset()
	{
		z1 = 0;
		z2 = 0;
	}
};

struct OnePoleModal
{
	const float Ts = 1.0f / 48000.0f;

	float a1 = 0;
	float z1 = 0;

	float pre = 0;
	float rre = 0;
	float pole = 0;
	float residue = 0;
	float step1 = 1;
	constexpr static int TableSize = TwoPoleModal::TableSize;
	float g1Table[TableSize];

	float impNormGain = 1.0f;
	float stepNormGain = 1.0f;

	inline float ProcessSample()
	{
		float y = z1;
		z1 = -a1 * y;
		return y;
	}

	void CalcPole(float pre, float rre)
	{
		this->pre = pre;
		this->rre = rre;

		pole = pre;
		residue = rre;
		step1 = expf(pre * Ts);
		a1 = -step1;

		for (int i = 0; i < TableSize; i++) {
			float tau = (float)i / (TableSize - 2);
			float dt1 = (1.0f - tau) * Ts;
			float shift = expf(pole * dt1);
			float g1 = residue * shift * impNormGain;
			g1Table[i] = g1;
		}
	}

	void CalcPole(float pre, float pim, float rre, float rim)
	{
		(void)pim;
		(void)rim;
		CalcPole(pre, rre);
	}

	/*
	void InjectImpulse(float tau, float v)
	{
		if (tau < 0.0f) tau = 0.0f;
		if (tau >= 1.0f) tau = 1.0f;
		float dt1 = (1.0f - tau) * Ts;
		float shift = expf(pole * dt1);
		float g1 = v * residue * shift * impNormGain;
		z1 += g1;
	}*/
	//void InjectImpulse(float tau, float v)
	inline void InjectImpulse(int index1, float frac, float v)
	{
		//if (tau < 0.0f) tau = 0.0f;
		//if (tau >= 1.0f) tau = 1.0f;
		//float fpos = tau * (TableSize - 2);
		//int index1 = fpos;
		int index2 = index1 + 1;
		//float frac = fpos - index1;
		float g1 = frac * (g1Table[index2] - g1Table[index1]) + g1Table[index1];
		z1 += g1 * v;
	}

	void SetNormGain(float impNormGain, float stepNormGain)
	{
		this->impNormGain = impNormGain;
		this->stepNormGain = stepNormGain;
	}

	void Reset()
	{
		z1 = 0;
	}
};

class SystemModal
{
public:
	constexpr static float Ts = 1.0f / 48000.0f;
	constexpr static int TableSize = TwoPoleModal::TableSize;
private:
	std::vector<TwoPoleModal> twoPoles;
	std::vector<OnePoleModal> onePoles;
	int numTwoPoles = 0;
	int numOnePoles = 0;
	bool hasDcPole = false;
	std::vector<float> twoPoleParams, onePoleParams;
	float dcPolePre = 0;
	float dcPoleResidue = 0;
public:
	void CalcPoles(const std::vector<float>& twoPoleParams, const std::vector<float>& onePoleParams)
	{
		this->twoPoleParams = twoPoleParams;
		this->onePoleParams = onePoleParams;
		numTwoPoles = (int)twoPoleParams.size() / 4;
		numOnePoles = (int)onePoleParams.size() / 2;
		twoPoles.resize(numTwoPoles);
		onePoles.resize(numOnePoles);
		for (int i = 0; i < numTwoPoles; i++) {
			float pre = twoPoleParams[i * 4 + 0];
			float pim = twoPoleParams[i * 4 + 1];
			float rre = twoPoleParams[i * 4 + 2];
			float rim = twoPoleParams[i * 4 + 3];
			twoPoles[i].CalcPole(pre, pim, rre, rim);
		}
		for (int i = 0; i < numOnePoles; i++) {
			float pre = onePoleParams[i * 2 + 0];
			float rre = onePoleParams[i * 2 + 1];
			onePoles[i].CalcPole(pre, rre);
		}
	}
	void GetParams(std::vector<float>& twoPoleParams, std::vector<float>& onePoleParams) const
	{
		twoPoleParams = this->twoPoleParams;
		onePoleParams = this->onePoleParams;
	}
	void GetDcPoleParams(float& pre, float& rre) const
	{
		pre = dcPolePre;
		rre = dcPoleResidue;
	}
	void SetNormGain(float impNormGain, float stepNormGain)
	{
		for (int i = 0; i < numTwoPoles; i++) {
			twoPoles[i].SetNormGain(impNormGain, stepNormGain);
		}
		for (int i = 0; i < numOnePoles; i++) {
			onePoles[i].SetNormGain(impNormGain, stepNormGain);
		}
	}
	void InjectImpulse(float tau, float v)
	{
		if (tau < 0.0f) tau = 0.0f;
		if (tau >= 1.0f) tau = 1.0f;
		float fpos = tau * (TableSize - 2);
		int index1 = fpos;
		float frac = fpos - index1;

		for (int i = 0; i < numTwoPoles; i++) {
			twoPoles[i].InjectImpulse(index1, frac, v);
		}
		for (int i = 0; i < numOnePoles; i++) {
			onePoles[i].InjectImpulse(index1, frac, v);
		}
	}
	float ProcessSample()
	{
		float y = 0;
		for (int i = 0; i < numTwoPoles; i++) {
			y += twoPoles[i].ProcessSample();
		}
		for (int i = 0; i < numOnePoles; i++) {
			y += onePoles[i].ProcessSample();
		}
		return y;
	}
	void Reset()
	{
		for (int i = 0; i < numTwoPoles; i++) {
			twoPoles[i].Reset();
		}
		for (int i = 0; i < numOnePoles; i++) {
			onePoles[i].Reset();
		}
	}
};

std::tuple<float, float> NormalizationResidues(std::vector<float>& twoPoleParams, std::vector<float>& onePoleParams)
{
	return std::make_tuple(1.0f * SystemModal::Ts, 1.0f);
	//下面的实现没有考虑到模拟意义。应该直接用模拟的极点和留数就行。

	float magimp = 0;
	float magstep = 0;
	for (size_t i = 0; i < twoPoleParams.size() / 4; i++) {
		std::complex<float> pole = { twoPoleParams[i * 4 + 0], twoPoleParams[i * 4 + 1] };
		std::complex<float> res = { twoPoleParams[i * 4 + 2], twoPoleParams[i * 4 + 3] };
		//magimp += std::abs(res);
		//magstep += std::abs(res / pole);
		magimp += res.real() * res.real();
		res /= pole;
		magstep += res.real() * res.real();
	}
	for (size_t i = 0; i < onePoleParams.size() / 2; i++) {
		float pre = onePoleParams[i * 2 + 0];
		float rre = onePoleParams[i * 2 + 1];
		//magimp += fabsf(rre);
		//magstep += fabsf(rre / pre);
		magimp += rre * rre;
		rre /= pre;
		magstep += rre * rre;
	}
	magimp = sqrtf(magimp);
	magstep = sqrtf(magstep);
	return std::make_tuple(1.0 / magimp, 1.0 / magstep);
}


class IIRBlep
{
private:
	SystemModal modal;
	float v = 0;

	const std::vector<float> twoPoleParams =
	{
		-22499.526025f, 25498.4629188f, -0.757564099837f, 0.511839786731f,
		-18558.7771212f, 71663.4117229f, 0.338279979228f, 0.0340141868387f,
		-12980.5407875f, 106472.913933f, -0.142233185588f, -0.0770342482721f,
		-7916.79591168f, 129087.105763f, 0.0510613557583f, 0.0611076206161f,
		-4096.74776281f, 141896.234098f, -0.00970295108387f, -0.0341334537474f,
		-1250.24954596f, 147548.854296f, -0.00103161482926f, 0.0100321602005f
	};
	const std::vector<float> onePoleParams =
	{
		-282.842712475f, -0.0140408541999f,
		-565.685424949f, 0.0564218869024f
	};
public:
	IIRBlep()
	{
		Setup(twoPoleParams, onePoleParams, 1.0f, 1.0f);
	}
	void Setup(const std::vector<float>& twoPoleParams, const std::vector<float>& onePoleParams, float impNormGain, float stepNormGain)
	{
		modal.CalcPoles(twoPoleParams, onePoleParams);
		modal.SetNormGain(impNormGain, stepNormGain);
	}
	void Add(float tau, float v, int mode)
	{
		if (mode == 1) modal.InjectImpulse(tau, v);
	}
	void Step()
	{
		v = modal.ProcessSample();
	}
	float Get()
	{
		return v;
	}
	void Reset()
	{
		modal.Reset();
		v = 0;
	}
};


class IIRBlamp
{
private:
	SystemModal modal;
	float v = 0;

	std::vector<float> twoPoleParams =
	{
		-22499.526025f, 25498.4629188f, 1.24922671663f, 0.323786880665f,
		-18558.7771212f, 71663.4117229f, -0.0336388822132f, -0.217867697197f,
		-12980.5407875f, 106472.913933f, -0.0265170714616f, 0.0673542084176f,
		-7916.79591168f, 129087.105763f, 0.0214771546503f, -0.0203039281966f,
		-4096.74776281f, 141896.234098f, -0.011442204851f, 0.00361262215595f,
		-1250.24954596f, 147548.854296f, 0.00326623120452f, 0.000307924503661f
	};

	std::vector<float> onePoleParams =
	{
		-282.842712475f, 2.38281197241f,
		-565.685424949f, -4.78755586032f
	};
public:
	IIRBlamp()
	{
		Setup(twoPoleParams, onePoleParams, 1.0f, 1.0f);
	}
	void Setup(const std::vector<float>& twoPoleParams, const std::vector<float>& onePoleParams, float impNormGain, float stepNormGain)
	{
		modal.CalcPoles(twoPoleParams, onePoleParams);
		modal.SetNormGain(impNormGain, stepNormGain);
	}
	void Add(float tau, float v, int mode)
	{
		if (mode == 2) modal.InjectImpulse(tau, v);
	}
	void Step()
	{
		v = modal.ProcessSample();
	}
	float Get()
	{
		return v;
	}
	void Reset()
	{
		modal.Reset();
		v = 0;
	}
};