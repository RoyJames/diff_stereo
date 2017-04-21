#include "main.h"
extern bool debug_flag;

void pref_env::load(const char * filename)
{
	FILE *fp;
	printf("loading env\n");

	FOPEN(fp, filename, "rb");
	fread(&dim, sizeof(int), 1, fp);

	fr.init(dim);
	num_level = 1;
	printf("env loaded\n");
	int d = 1;
	while (d < dim)
	{
		d *= 2;
		num_level++;
	}

	env.resize(3*dim*dim*(num_level + 1));

	fread(&env[0], env.size() * sizeof(float), 1, fp);
	fclose(fp);
}

void pref_env::lookup(vector3f & result, const vector3f & mu, const float kappa_log) const
{
	int kappa_level0 = min(max(0, floor(kappa_log)), num_level - 1);
	int kappa_level1 = min(max(0, kappa_level0 + 1), num_level - 1);

	vector3f c0, c1;

	int num;
	int idx[4];
	float weight[4];
	fr.lerp(num, idx, weight, mu);
	//printf_s("vi=(%g,%g,%g)\n", mu.x, mu.y, mu.z);
	//for (int i = 0; i < num; i++) printf_s("%d\t", idx[i]);
	//printf_s("\n");
	//for (int i = 0; i < num; i++) printf_s("%.3f\t", weight[i]);
	//printf_s("\n");

	/*The following is using Gaussian filter*/
	int w = 6;
	float alpha_fall = 0.01f;
	int *neighbors = new int[w*w];
	fr.filter_prep_index(neighbors, idx[0], w / 2);

	//printf_s("center:%d\n", idx[0]);
	//for (int i = 0; i < w; i++) {
	//	for (int j = 0; j < w; j++) {
	//		printf_s("%8d ", neighbors[i*w + j]);
	//	}
	//	printf_s("\n");
	//}
	//exit(1);

	if (debug_flag) {
		printf_s("passed index:%d\n", idx[0]);
		for (int i = 0; i < w; i++) {
			for (int j = 0; j < w; j++) {
				printf_s("%d ", neighbors[i*w + j]);
			}
			printf_s("\n");
		}
	}
	float p_row = 1.0f - (weight[0] + weight[2]);
	float p_col = 1.0f - (weight[0] + weight[1]);
	float *vals_0 = new float[w * w];
	float *vals_1 = new float[w * w];
	for (int ch = 0; ch < 3; ch++) {
		for (int i = 0; i < w*w; i++) {
			vals_0[i] = env[kappa_level0 * 3 * dim*dim + neighbors[i] * 3 + ch];
			vals_1[i] = env[kappa_level1 * 3 * dim*dim + neighbors[i] * 3 + ch];
		}
		c0.v[ch] = gauss_filter(vals_0, alpha_fall, w, p_row, p_col);
		c1.v[ch] = gauss_filter(vals_1, alpha_fall, w, p_row, p_col);
	}
	delete[]vals_0;
	delete[]vals_1;
	delete[]neighbors;
	/**************************************/

	//for (int i = 0; i < num; i++)
	//	for (int ch = 0; ch < 3; ch++)
	//	{
	//		c0.v[ch] += weight[i] * env[kappa_level0 * 3 * dim*dim + idx[i] * 3 + ch];
	//		c1.v[ch] += weight[i] * env[kappa_level1 * 3 * dim*dim + idx[i] * 3 + ch];
	//	}

	float alpha = min(max(kappa_log - kappa_level0, 0), 1);
	result = c0*(1 - alpha) + c1*alpha;
	if (result.x > 100) printf_s("lookup warning!\n");
}

void pref_env::lookup_cos(vector3f & result, const vector3f & mu) const
{
	vector3f c;
	int num;
	int idx[4];
	float weight[4];
	fr.lerp(num, idx, weight, mu);

	/*The following is using Gaussian filter*/
	//int w = 6;
	//float alpha_fall = 0.01f;
	//int *neighbors = new int[w*w];
	//fr.filter_prep_index(neighbors, idx[0], w / 2);
	//float p_row = 1.0f - (weight[0] + weight[2]);
	//float p_col = 1.0f - (weight[0] + weight[1]);
	//float *vals = new float[w * w];
	//for (int ch = 0; ch < 3; ch++) {
	//	for (int i = 0; i < w*w; i++) {
	//		vals[i] = env[num_level * 3 * dim*dim + neighbors[i] * 3 + ch];
	//	}
	//	c.v[ch] = gauss_filter(vals, alpha_fall, w, p_row, p_col);
	//}
	//delete[]vals;
	//delete[]neighbors;
	/**************************************/

	for (int i = 0; i < num; i++)
		for (int ch = 0; ch < 3; ch++)
			c.v[ch] += weight[i] * env[num_level * 3 * dim*dim + idx[i] * 3 + ch];
	result = c;
}

float pref_env::gauss_filter(const float *data, float alpha, int w, float px, float py) const
{
	float *weight_x = new float[w];
	float *weight_y = new float[w];
	float offset = expf(-alpha * w * w / 4.0f);
	px += (float)w / 2.0f - 1.0f;
	py += (float)w / 2.0f - 1.0f;
	for (int i = 0; i < w; i++) {
		weight_x[i] = expf(-alpha * (px - i) * (px - i)) - offset;
		weight_y[i] = expf(-alpha * (py - i) * (py - i)) - offset;
		//printf_s("%d: %g\t%g\n", i, weight_x[i], weight_y[i]);
	}
	float wt_sum = 0;
	float product_sum = 0;
	int cnt = 0;
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < w; j++) {
			product_sum += data[cnt++] * weight_x[i] * weight_y[j];
			wt_sum += weight_x[i] * weight_y[j];
		}
	}
	//printf_s("weight matrix:\n");
	//for (int i = 0; i < w; i++) {
	//	for (int j = 0; j < w; j++) {
	//		printf_s("%.3f ", weight_x[i] * weight_y[j] / wt_sum);
	//	}
	//	printf("\n");
	//}
	//exit(1);

	delete[]weight_x;
	delete[]weight_y;
	return product_sum / wt_sum;
}


float pref_env::gauss_interp(const int *ax, const float *ay, float p) const
{
	// GAUSS FORWARD INTERPOLATION
	float y = 0;                      // Calculated value for coressponding X.
	float h;                        // Calc. Section
	float diff[5][5];             // to store Y
	float y1, y2, y3, y4;              // Formulae variables.

	int n = 5;	// input 5 pairs, index of interpolation should be 
				// between the 4th and 5th, which means p>0
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			diff[i][j] = 0;
		}
	}
	for (int i = 0; i < n - 1; i++)
		diff[i][1] = ay[i + 1] - ay[i];

	for (int j = 2; j <= 4; j++)
		for (int i = 0; i < n - j; i++)
			diff[i][j] = diff[i + 1][j - 1] - diff[i][j - 1];

	y1 = p*diff[3][1];
	y2 = p*(p - 1.0f)*diff[2][2] / 2.0f;
	y3 = (p + 1.0f)*p*(p - 1.0f)*diff[1][3] / 6.0f;
	y4 = (p + 1.0f)*p*(p - 1.0f)*(p - 2.0f)*diff[0][4] / 24.0f;

	// Taking sum
	y = ay[3] + y1 + y2 + y3 + y4;

	//if (fabs(y) > 0.1) {
	//	for (int k = 0; k < 5; k++) {
	//		printf_s("x=%d:\ty=%g\n", ax[k], ay[k]);
	//	}
	//	printf_s("y1~y4:%g\t%g\t%g\t%g\n", y1, y2, y3, y4);
	//	printf_s("interpolated:%g\n", y);
	//	//exit(1);
	//}
	//exit(1);
	return y;
}
