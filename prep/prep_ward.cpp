#include "main.h"
#include <omp.h>

extern float EMEM_vMF(std::vector<vMF<float>> &result,
	const int k,
	const std::vector<vector3f> &y,
	random_number_generator<float> &rng,
	const int num_init,
	const int init_iter,
	const int max_iter,
	const float epsilon,
	const bool b_update_mu);


template <class Number>
inline Number my_folded_radical_inverse(__int64 n, __int64 base)
{
	Number	val = 0;
	Number	inv_base = (Number)1 / base, inv_bi = inv_base;
	__int64	mod_offset = 0;

	while (val + base * inv_bi != val)
	{
		__int64 digit = (n + mod_offset) % base;

		val += digit * inv_bi;
		n /= base;
		inv_bi *= inv_base;
		mod_offset++;
	}

	return val;
}

template<class Number>
class my_Hammersley_Zaremba_seq
{
private:
	__int64	i, n;

public:
	my_Hammersley_Zaremba_seq(__int64 n, __int64 i = 0) : n(n), i(i) {};

	void get_sample(Number &r)
	{
		r = (Number)i / n;
		i++;
	}

	void get_sample(codex::math::vector::vector3<Number> &r)
	{
		r.x = (Number)i / n;
		r.y = my_folded_radical_inverse<Number>(i, 2);
		r.z = my_folded_radical_inverse<Number>(i, 3);
		i++;
	}
};

void precompute_area(const char *filename, const hemi_octa_frame<float> &env_frame)
{
	int dim = env_frame.get_dim();
	__int64 num_dir = (__int64)dim*(__int64)dim*(__int64)1024;

	codex::math::prob::uniform_hemi_direction_3d<float> prob;

	std::vector<float> area;
	area.resize(dim*dim);

	std::vector<std::vector<int>> counter;
#pragma omp parallel
#pragma omp master
	counter.resize(omp_get_num_threads());

	for (int i = 0; i < counter.size(); i++)
		counter[i].resize(dim*dim);

	codex::utils::timer t;
#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j < 100; j++)
	{
		__int64 i0 = j*num_dir / 100,
			i1 = (j + 1)*num_dir / 100;
		printf_s("%d ", j);

		int idx = omp_get_thread_num();

		my_Hammersley_Zaremba_seq<float> seq(num_dir, i0);
		for (__int64 i = i0; i < i1; i++)
		{
			vector3f rv, p;
			float pdf;
			seq.get_sample(rv);
			prob.sample(p, pdf, rv);

			counter[idx][env_frame.map(p)]++;
		}
	}

	for (int j = 0; j < counter.size(); j++)
		for (int i = 0; i <dim*dim; i++)
			area[i] += counter[j][i] * (float)(2 * PI) / num_dir;

	for (int y = 0; y <dim / 2; y++)
		for (int x = 0; x <dim / 2; x++)
		{
			float avg = (area[x + y*dim] + area[dim - 1 - x + y*dim] + area[x + (dim - 1 - y)*dim] + area[dim - 1 - x + (dim - 1 - y)*dim]) / 4;
			area[x + y*dim] = area[dim - 1 - x + y*dim] = area[x + (dim - 1 - y)*dim] = area[dim - 1 - x + (dim - 1 - y)*dim] = avg;
		}
	FILE *fp;
	FOPEN(fp, filename, "wb");
	fwrite(&area[0], sizeof(float)*area.size(), 1, fp);
	fclose(fp);
}

void load_area(std::vector<float> &area, const char *filename, const hemi_octa_frame<float> &fr)
{
	if (codex::utils::file_exists(filename))
	{
		area.resize(fr.get_dim()*fr.get_dim());
		FILE *fp;
		FOPEN(fp, filename, "rb");
		fread(&area[0], sizeof(float)*area.size(), 1, fp);
		fclose(fp);
	}
	else {
		precompute_area(filename, fr);
	}
}

void fit_BRDF(std::vector<vMF<float>> &result,
	const BRDF<float> &brdf,
	const int num_wo,
	const int num_EM_samples,
	const hemi_octa_frame<float> &fr,
	const std::vector<float> &area)
{
	result.clear();

	int dim = fr.get_dim();
	codex::math::utils::random_number_generator<float> rng;

	std::vector<vector3f> wi;
	for (int i = 0; i < dim*dim; i++)
	{
		vector3f v;
		fr.get_n(v, i);
		wi.push_back(v);
	}

	std::vector<float> pdf, f_brdf, f_ours;
	pdf.resize(dim*dim);
	f_brdf.resize(dim*dim);
	f_ours.resize(dim*dim);

	for (int i = 0; i < num_wo; i++)
	{
		float theta = (float)PI / 2 * (i + 1) / num_wo;
		vector3f wo(cos(theta), 0, sin(theta));

		for (int j = 0; j < wi.size(); j++)
		{
			float f;
			brdf.sample(f, wi[j], wo);
			pdf[j] = area[j] * f*wi[j].z;
			f_brdf[j] = area[j] * f*wi[j].z;
		}

		std::vector<vector3f> brdf_samples;
		codex::math::prob::step_1D_prob<float> cdf(pdf);
		for (int j = 0; j < num_EM_samples; j++)
		{
			int idx;
			float pdf;
			cdf.sample(idx, pdf, (j + 0.5f) / num_EM_samples);

			brdf_samples.push_back(wi[idx]);
		}

		std::vector<vMF<float>> slice;
		slice.resize(1);
		slice[0].w = 1;
		codex::math::vector::reflect(slice[0].mu, wo, vector3f(0, 0, 1));
		EMEM_vMF(slice, 1, brdf_samples, rng, 5, 20, 100, 0, false);

		float xx = 0, xy = 0;
		for (int j = 0; j < wi.size(); j++)
		{
			f_ours[j] = area[j] * slice[0].eval(wi[j]);
			xx += f_ours[j] * f_ours[j];
			xy += f_ours[j] * f_brdf[j];
		}

		slice[0].w = xy / xx;
		result.push_back(slice[0]);
	}
}
