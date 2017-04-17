#include "codex_graphics.h"
#include "advmath.h"
#include "brdflib.h"

#define EM_KMEANS_ITER	20
#define VMF_KAPPA_MAX	1e20f

using codex::math::utils::random_number_generator;
using codex::math::vector3f;
using namespace brdflib;
using namespace advmath;

float EM_vMF_init(std::vector<vMF<float>> &g,
	const int k,
	const std::vector<vector3f> &y,
	const int max_iter,
	std::vector<float> &gamma,
	const bool b_update_mu)
{
	int n = (int)y.size();
	std::vector<float> sum, eta;
	//gamma.resize(n*k, 0);
	sum.resize(n, 0);
	eta.resize(k);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < k; j++)
			sum[i] += gamma[i*k + j];

	int iter = 0;
	do {
		iter++;

		//E-step
		eta.assign(k, 0);

		for (int i = 0; i < n; i++)
			if (sum[i] == 0)
			{
				float best_dot = -1e30f;
				int best_idx = -1;
				for (int j = 0; j < k; j++)
				{
					float temp = y[i] * g[j].mu;
					if (temp > best_dot)
					{
						best_dot = temp;
						best_idx = j;
					}
				}
				gamma[i*k + best_idx] = 1;
				sum[i] = 1;
			}

		for (int i = 0; i < n; i++)
			for (int j = 0; j < k; j++)
			{
				gamma[i*k + j] /= sum[i];
				eta[j] += gamma[i*k + j];
			}

		//M-step
		int i;
		//#pragma omp parallel for private(i)
		for (int j = 0; j < k; j++)
		{
			//w
			g[j].w = eta[j] / n;

			if (eta[j] > 0)
			{
				vector3f mu;
				for (i = 0; i < n; i++)
					mu += gamma[i*k + j] * y[i];
				mu /= eta[j];
				float len_mu = mu.length();

				//kappa
				if (1 - len_mu*len_mu <= 0)
					g[j].kappa = VMF_KAPPA_MAX;
				else
					g[j].kappa = (3 - len_mu*len_mu)*len_mu / (1 - len_mu*len_mu);

				//mu
				if (b_update_mu)
					g[j].mu = mu / len_mu;
			}
		}

		//compute new l
		sum.assign(n, 0);

		int j;
		//#pragma omp parallel for private(j)
		for (int i = 0; i < n; i++)
			for (j = 0; j < k; j++)
			{
				gamma[i*k + j] = g[j].w*g[j].eval(y[i]);
				sum[i] += gamma[i*k + j];
			}
	} while (iter < max_iter);

	int nn = 0;
	float l = 0;
	for (int i = 0; i < n; i++)
		if (sum[i] > 0)
		{
			l += log(sum[i]);
			nn++;
		}
	if (nn > 0) l /= nn;

	std::vector<vMF<float>> result;
	for (int i = 0; i < g.size(); i++)
		if (g[i].w > 0)
			result.push_back(g[i]);
	g = result;

	return l;
}

float EM_random_init(std::vector<vMF<float>> &g,
	std::vector<float> &gamma,
	const int k,
	const int max_iter,
	const std::vector<vector3f> &y,
	const bool b_update_mu,
	random_number_generator<float> &rng)
{
	std::vector<int> perm;
	codex::math::prob::generate_permutation(perm, (int)y.size(), rng, k);

	g.resize(k);
	for (int j = 0; j < k; j++)
	{
		g[j].w = 1.0f / k;
		g[j].kappa = 200;
		if (b_update_mu)
			g[j].mu = y[perm[j]];
	}

	gamma.resize(y.size()*g.size(), 0);
	for (int i = 0; i < y.size(); i++)
	{
		int cluster_id = -1;
		float best_dist = 1e30f;
		for (int j = 0; j < g.size(); j++)
		{
			float dist_sq = (g[j].mu - y[i]).length_sq();
			if (dist_sq < best_dist)
			{
				best_dist = dist_sq;
				cluster_id = j;
			}
		}
		gamma[i*g.size() + cluster_id] = 1;
	}

	return EM_vMF_init(g, (int)g.size(), y, max_iter, gamma, b_update_mu);
}

float EM_kmeans_init(std::vector<vMF<float>> &g,
	std::vector<float> &gamma,
	const int k,
	const int init_iter,
	const std::vector<vector3f> &y,
	random_number_generator<float> &rng)
{
	std::vector<float> y_vec, centers;
	y_vec.resize(3 * y.size());
	for (int i = 0; i < y.size(); i++)
	{
		y_vec[i * 3 + 0] = y[i].x;
		y_vec[i * 3 + 1] = y[i].y;
		y_vec[i * 3 + 2] = y[i].z;
	}
	std::vector<int> cluster_size, cluster_id;
	kmeans_clustering<float>(k, EM_KMEANS_ITER, &y_vec[0], (int)y.size(), 3, &rng, &centers, &cluster_id, &cluster_size);

	g.resize(centers.size() / 3);
	for (int j = 0; j < g.size(); j++)
	{
		g[j].w = (float)cluster_size[j] / y.size();
		g[j].kappa = VMF_KAPPA_MAX;
		g[j].mu = vector3f(centers[j * 3 + 0], centers[j * 3 + 1], centers[j * 3 + 2]);
		g[j].mu.normalize();
	}

	gamma.resize(y.size()*g.size(), 0);
	for (int i = 0; i < y.size(); i++)
		gamma[i*g.size() + cluster_id[i]] = 1;

	return EM_vMF_init(g, (int)g.size(), y, init_iter, gamma, true);
}

float EM_vMF(std::vector<vMF<float>> &g,
	std::vector<float> &gamma,
	const int k,
	const std::vector<vector3f> &y,
	const int max_iter,
	const float epsilon,
	const bool b_update_mu,
	random_number_generator<float> &rng)
{
	int n = (int)y.size();
	double l = -1, previous_l;
	std::vector<float> sum, eta;
	//gamma.resize(n*k, 0);
	sum.resize(n);
	eta.resize(k);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < k; j++)
			sum[i] += gamma[i*k + j];

	int nn;

	int iter = 1;
	do {
		iter++;

		//E-step
		eta.assign(k, 0);
		if (iter > 2)
		{
			for (int i = 0; i < n; i++)
				if (sum[i] == 0)
				{
					float best_dot = -1e30f;
					int best_idx = -1;
					for (int j = 0; j < k; j++)
					{
						float temp = y[i] * g[j].mu;
						if (temp > best_dot)
						{
							best_dot = temp;
							best_idx = j;
						}
					}
					gamma[i*k + best_idx] = 1;
					sum[i] = 1;
				}

			for (int i = 0; i < n; i++)
				for (int j = 0; j < k; j++)
				{
					gamma[i*k + j] /= sum[i];
					eta[j] += gamma[i*k + j];
				}
			//M-step
			//#pragma omp parallel for
			for (int j = 0; j < k; j++)
			{
				//w
				g[j].w = eta[j] / n;

				if (eta[j] > 0)
				{
					vector3f mu;
					for (int i = 0; i < n; i++)
						mu += gamma[i*k + j] * y[i];
					mu /= eta[j];
					float len_mu = mu.length();

					//kappa
					if (1 - len_mu*len_mu <= 0)
						g[j].kappa = VMF_KAPPA_MAX;
					else
						g[j].kappa = (3 - len_mu*len_mu)*len_mu / (1 - len_mu*len_mu);

					//mu
					//if (C <= 0)
					//	g[j].mu = mu / len_mu;
					//else {
					//	g[j].mu = mu + C*neighbor[j].w*neighbor[j].mu;
					//	g[j].mu.normalize();
					//}
					if (b_update_mu)
						g[j].mu = mu / len_mu;
				}
			}
		}

		//compute new l
		previous_l = l;
		l = 0;
		sum.assign(n, 0);

		int j;
		nn = 0;
		//#pragma omp parallel for private(j)
		for (int i = 0; i < n; i++)
		{
			for (j = 0; j < k; j++)
			{
				gamma[i*k + j] = g[j].w*g[j].eval(y[i]);
				sum[i] += gamma[i*k + j];
			}
			if (sum[i] > 0)
			{
				l += log(sum[i]);
				nn++;
			}
		}
		l /= nn;
	} while (abs(l - previous_l) > epsilon && iter - 1 < max_iter);

	std::vector<vMF<float>> result;
	for (int i = 0; i < g.size(); i++)
		if (g[i].w > 0)
			result.push_back(g[i]);
	g = result;
	return (float)l;
}

float EMEM_vMF(std::vector<vMF<float>> &result,
	const int k,
	const std::vector<vector3f> &y,
	random_number_generator<float> &rng,
	const int num_init,
	const int init_iter,
	const int max_iter,
	const float epsilon,
	const bool b_update_mu)
{
	float best = -1e6;
	std::vector<vMF<float>> best_g;
	std::vector<float> best_gamma;

	for (int i = 0; i < num_init; i++)
	{
		std::vector<vMF<float>> g;
		std::vector<float> gamma;

		if (!b_update_mu)
			g = result;
		float p = //EM_kmeans_init(g, gamma, k, init_iter, y, rng);
			EM_random_init(g, gamma, k, init_iter, y, b_update_mu, rng);
		if (p > best)
		{
			best = p;
			best_g = g;
			best_gamma = gamma;
		}
	}

	//printf_s("\n");
	result = best_g;
	if (result.size() == 0)
		return 0;
	else
		return EM_vMF(result, best_gamma, (int)result.size(), y, max_iter, epsilon, b_update_mu, rng);
}
