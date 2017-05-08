#include "main.h"
#include "optimizer.h"

using namespace std;

void prec_brdf()
{
	const int dim = 2048;
	hemi_octa_frame<float> fr;
	fr.init(dim);

	std::vector<float> area;
	load_area(area, "D:/codes/iff/data/area_2048.dat", fr);

	const float alpha = 0.1f;
	WardBRDF<float> brdf(vector3f(1, 1, 1), alpha);
	std::vector<vMF<float>> result;
	fit_BRDF(result, brdf, 90, 1024, fr, area);

	FILE *fp;
	char filename[MAX_PATH];
	sprintf_s(filename, "d:/codes/diff/data/ward_%01.03f.vMF", brdf.a);
	FOPEN(fp, filename, "wb");
	
	fwrite(&brdf.a, sizeof(float), 1, fp);
	int num = result.size();
	fwrite(&num, sizeof(int), 1, fp);
	for (int i = 0; i < result.size(); i++)
	{
		fwrite(&result[i].w, sizeof(float), 1, fp);
		fwrite(&result[i].kappa, sizeof(float), 1, fp);
		fwrite(result[i].mu.v, sizeof(float), 3, fp);
	}
	
	fclose(fp);

	for (int i = 0; i < result.size(); i++)
		printf_s("%d %g %g\n", i, result[i].w, result[i].kappa);
}

double rand_range(double range)
{
	double t = (double)(rand() / (double)RAND_MAX);
	return (t - 0.5) * 2 * range;
}

void test()
{
	la_matrix<float> T(2, 3);
	T.m[0] = 1.0;
	T.m[1] = 1.0;
	T.m[2] = 2.0;
	T.m[3] = 0.0;
	T.m[4] = 0.0;
	T.m[5] = 0.0;
	la_vector<float> v(3);
	v.v[0] = 1.0;
	v.v[1] = 0.1;
	v.v[2] = 0;
	la_vector<float> res;
	smvmul(res, T, v);
	T.save_txt("Tmat");
	v.save_txt("vmat");
	res.save_txt("resmat");
	//std::cout << T.m << std::endl;
	exit(1);
}

extern void compute_rho(vector3f &rho, const vector3f &wo, const vector3f &n, const pref_env &env);


void test_view_diff()
{
	random_number_generator<float> rng;
	float scale_trans = 0.1f;// 0.5f;//0.05f;

	pref_env env;
	env.load("d:/codes/diff/data/uffizi.pref_env");

	std::vector<float> mA, vb;
	for (int i = 0; i < 100; i++)
	{
		vector3f tau(
			(rng.rand_real() * 2 - 1)*scale_trans,
			(rng.rand_real() * 2 - 1)*scale_trans,
			(rng.rand_real() * 2 - 1)*scale_trans
		);

		vector3f n(0.201556f, -0.000390217f, 0.979477f);// (0, 0, 1);
		n.normalize();
		vector3f wo(-0.201933f, 0.00078573f, 4.02294f);// (0, 0, FOCAL_LEN-1);
		wo += tau;
		wo.normalize();

		vector3f c;
		compute_rho(c, wo, n, env);

		mA.push_back(tau.x);
		mA.push_back(tau.y);
		mA.push_back(tau.z);
		vb.push_back((c.x + c.y + c.z) / 3);
	}

	la_matrix<float> A(vb.size(), 4);
	la_vector<float> b(vb);
	for (int i = 0; i < A.row; i++)
	{
		A.m[0 * A.row + i] = mA[i * 3 + 0];
		A.m[1 * A.row + i] = mA[i * 3 + 1];
		A.m[2 * A.row + i] = mA[i * 3 + 2];
		A.m[3 * A.row + i] = 1;
	}

	la_matrix<float> ATA;
	la_vector<float> ATb;
	smmmul(ATA, A.transpose(), A);
	smvmul(ATb, A.transpose(), b);

	la_vector<float> x;
	slsq(x, ATA, ATb, 1);

	printf_s("--org--\n");
	for (int i = 0; i < b.length; i++)
		printf_s("%f\n", b.v[i]);
	
	printf_s("--approx--\n");
	for (int i = 0; i < b.length; i++)
		printf_s("%f\n", 
			mA[i * 3 + 0] * x.v[0] +
			mA[i * 3 + 1] * x.v[1] +
			mA[i * 3 + 2] * x.v[2] +
			x.v[3]
			//x.v[2]
		);

	printf_s("--x--\n");
	for (int i = 0; i < x.length; i++)
		printf_s("%f ", x.v[i]);
	printf_s("\n");

	//{
	//	vector3f wo(0, 0, FOCAL_LEN);
	//	printf_s("----\n");
	//	printf_s("%g %g\n", x.v[2], -wo.x/wo.z*x.v[0] - wo.y / wo.z*x.v[1]);
	//}
}

void compute_brdf(vector3f &rho, const vector3f &wo, const vector3f &n, const float kappa_log, const pref_env &env)
{
	vector3f wi;
	codex::math::vector::reflect(wi, wo, n);
	env.lookup(rho, wi, kappa_log);
}

void test_view_diff2()
{
	random_number_generator<float> rng;
	float scale_trans = 0.1f;// 0.5f;//0.05f;

	pref_env env;
	env.load("d:/codes/diff/data/uffizi.pref_env");

	std::vector<float> angle;
	std::vector<float> grad;

	//vector3f n0(0.201556f, -0.000390217f, 0.979477f);// (0, 0, 1);
	//n0.normalize();
	
	const int dim = 512;
	full_octa_frame<float> fr;
	fr.init(dim);

	const int per_level = 1;
	for (int k = 0; k < per_level*8; k++)
	{
		std::vector<float> img;

		for (int j = 0; j < dim*dim; j++)
		{
			vector3f n0;
			fr.get_n(n0, j);

			vector3f wo0(-0.201933f, 0.00078573f, 4.02294f);// (0, 0, FOCAL_LEN-1);

			float c0;
			{
				vector3f n = n0;
				vector3f wo = wo0;
				wo.normalize();

				vector3f c;
				compute_brdf(c, wo, n, (float)k / per_level, env);
				c0 = (c.x + c.y + c.z) / 3;
				//printf_s("c0 = %g\n", c0);

				img.push_back(c.x);
				img.push_back(c.y);
				img.push_back(c.z);
			}

			std::vector<float> mA, vb;
			for (int i = 0; i < 200; i++)
			{
				vector3f tau(
					(rng.rand_real() * 2 - 1)*scale_trans,
					(rng.rand_real() * 2 - 1)*scale_trans,
					(rng.rand_real() * 2 - 1)*scale_trans
				);

				vector3f n = n0;
				vector3f wo = wo0 + tau;
				wo.normalize();

				vector3f c;
				compute_brdf(c, wo, n, (float)k / per_level, env);

				mA.push_back(tau.x - wo.x / wo.z*tau.z);
				mA.push_back(tau.y - wo.y / wo.z*tau.z);
				vb.push_back((c.x + c.y + c.z) / 3 - c0);
			}

			la_matrix<float> A(vb.size(), 2);
			la_vector<float> b(vb);
			for (int i = 0; i < A.row; i++)
			{
				A.m[0 * A.row + i] = mA[i * 2 + 0];
				A.m[1 * A.row + i] = mA[i * 2 + 1];
			}

			la_matrix<float> ATA;
			la_vector<float> ATb;
			smmmul(ATA, A.transpose(), A);
			smvmul(ATb, A.transpose(), b);

			la_vector<float> x;
			slsq(x, ATA, ATb, 1);

			angle.push_back(atan2(x.v[1], x.v[0]));
			grad.push_back(x.v[0]);
			grad.push_back(x.v[1]);
		}

		char filename[MAX_PATH];
		sprintf_s(filename, "d:/codes/diff/data/ref%02d.raw", k);
		save_img_float(filename, &img[0], dim, dim);
	}

	//std::vector<int> count;
	//count.resize(100);
	//for (int i = 0; i < angle.size(); i++)
	//{
	//	float a = fmodf(angle[i] + 2 * PI, 2 * PI);
	//	count[min(count.size() * a / (2 * PI), count.size() - 1)]++;
	//}

	//printf_s("------\n");
	//for (int i = 0; i < count.size(); i++)
	//	printf_s("%d\n", count[i]);
	
	{
		FILE *fp;
		FOPEN(fp, "d:/codes/diff/data/angle.bat", "wt");

		std::vector<float> img;
		for (int j = 0; j < angle.size() / (dim*dim); j++)
		{
			img.assign(dim*dim * 3, 0);

			float max_len = 0;
			for (int i = 0; i < dim*dim; i++)
			{
				float x = grad[(i + j*dim*dim) * 2 + 0];
				float y = grad[(i + j*dim*dim) * 2 + 1];
				max_len = max(max_len, sqrt(x*x + y*y));
			}

			for (int i = 0; i < dim*dim; i++)
			{
				float x = grad[(i + j*dim*dim) * 2 + 0];
				float y = grad[(i + j*dim*dim) * 2 + 1];
				img[i * 3 + 0] = (x / max_len + 1) / 2;
				img[i * 3 + 1] = (y / max_len + 1) / 2;
				img[i * 3 + 2] = 0;
			}

			char filename[MAX_PATH];
			sprintf_s(filename, "d:/codes/diff/data/angle%02d.raw", j);
			save_img_float(filename, &img[0], dim, dim);

			fprintf_s(fp, "raw2img %s\n", filename);
			sprintf_s(filename, "d:/codes/diff/data/ref%02d.raw", j);
			fprintf_s(fp, "raw2img %s\n", filename);
		}

		fclose(fp);
	}

	FILE *fp;
	FOPEN(fp, "d:/codes/diff/data/angle.txt", "wt");
	for (int i = 0; i < dim*dim; i++)
	{
		for (int j = 0; j < angle.size()/(dim*dim); j++)
			fprintf_s(fp, "%g ", angle[i + j*dim*dim]);
		fprintf_s(fp, "\n");
	}
	fclose(fp);

	//printf_s("--org--\n");
	//for (int i = 0; i < b.length; i++)
	//	printf_s("%f\n", b.v[i]);

	//printf_s("--approx--\n");
	//for (int i = 0; i < b.length; i++)
	//	printf_s("%f\n",
	//		mA[i * 3 + 0] * x.v[0] +
	//		mA[i * 3 + 1] * x.v[1] +
	//		mA[i * 3 + 2] * x.v[2]);

	//printf_s("--x--\n");
	//for (int i = 0; i < x.length; i++)
	//	printf_s("%f ", x.v[i]);
	//printf_s("\n");
}


class diffuse_func : public func<float>
{
public:
	std::vector<float>	coef;
	std::vector<float>	k1;
	float beta;
	int view_num;

	//the function
	virtual void f(la_vector<float> &fx, const la_vector<float> &x)
	{
		fx.init(view_num);
		//printf("calculating f:\n");
		for (int i = 0; i < fx.length; i++)
		{
			fx.v[i] = coef[i] / (k1[i] - beta * x.v[0]);
			//printf("%g ", fx.v[i]);
		}
		//printf("\n");
	}
	//Jacobian
	virtual void j(la_vector<float> &jx, const la_vector<float> &x)
	{
		jx.init(view_num);
		jx.clear();
		//printf("calculating Jacobian:\n");
		for (int i = 0; i < jx.length; i++)
		{
			jx.v[i] = coef[i] * beta / (k1[i] - beta * x.v[0]) / (k1[i] - beta * x.v[0]);
			//printf("%g ", jx.v[i]);
		}
		//printf("\n");
	}
};


//class myfunc : public func<float>
//{
//public:
//	//the function
//	virtual void f(la_vector<float> &fx, const la_vector<float> &x)
//	{
//		fx.init(2);
//		fx.v[0] = x.v[0];
//		fx.v[1] = x.v[1];
//	};
//};

void test_filtered_1st_order()
{
	const int num = 20;
	const int num_gaussian = 5;
	std::vector<float> a;

	for (int i = 0; i < num; i++)
		a.push_back(-1.0f+2.0f*i/num);
	for (int i = 0; i < num; i++)
		a.push_back(i < num / 2 ? -0.5 : 0.5);

	for (int j = 0; j < num_gaussian; j++)
	{
		float c = (j + 1)*3;
		for (int i = 0; i < num; i++)
		{
			float sum_wt = 0, sum = 0;
			for (int k = max(0, i - num / 2); k <= min(i + num / 2, num - 1); k++)
			{
				float wt = exp(-(k - i)*(k - i) / (c*c));
				sum_wt += wt;
				sum += wt*a[num+k];
			}
			a.push_back(sum / sum_wt);
		}
	}

	for (int j = 0; j < num_gaussian+1; j++)
	{
		float xx, xy;
		xx = xy = 0;
		for (int i = 0; i < num; i++)
		{
			xx += a[i] * a[i];
			xy += a[i] * a[(1 + j)*num + i];
		}
		float b = xy / xx;
		for (int i = 0; i < num; i++)
			a.push_back(a[i] * b);
	}

	FILE *fp;
	FOPEN(fp, "d:/codes/diff/data/1st.txt", "wt");
	for (int i = 0; i < num; i++)
	{
		for (int j = 0; j < 2*num_gaussian + 2 + 1; j++)
			fprintf_s(fp, "%g ", a[i + j*num]);
		fprintf_s(fp, "\n");
	}
	fclose(fp);
}


void initialize_samples(object_points &my_points, int x_begin, int y_begin, int x_end, int y_end)
{
	for (int i = x_begin; i <= x_end; i++)
	{
		for (int j = y_begin; j <= y_end; j++)
		{
			my_points.register_point(i, j);
		}
	}
	my_points.match_neighbors();
}

void test_Lambertian_v1(std::vector<std::pair<int, int> > scr, const std::vector<float> groundtruth)
{
	const int nsamples = 1000;
//	const float z0 = 0.1f, z1 = 1.1f;
	const float z0 = -20.0f, z1 = -19.0f;
	std::vector<vector2f> uv;

									 //uv.x = 0.00195312f;
									 //uv.y = -0.00195312f;
									 //DEBUG
	std::vector<float> k1;
	std::vector<std::vector<float> > coeff_A;
	std::vector<std::vector<float> > coeff_b;
	float *recon_err = new float[NUM_IMG_PAIRS * scr.size()];

	printf_s("start preprocessing...\n");
	int width, height;

	std::vector<float> pixel0;
	//std::vector<vector2f> scr;
	{
		char filename[MAX_PATH];
		sprintf_s(filename, PATH "img00-%02d.raw", 0);

		std::vector<float> img;
		load_img_float(filename, img, width, height);

		for (int i_point = 0; i_point < scr.size(); i_point++)
		{
			vector2f scr_i(scr[i_point].first, scr[i_point].second);
			//uv2screen(scr_i, uv[i_point], width, height);
			//printf("%d:%g %g\n", i_point, scr_i.x, scr_i.y);
			//scr.push_back(scr_i);
			//printf("%g %g\n", scr.x, scr.y);
			vector2f uv_i;
			screen2uv(uv_i, scr_i, width, height);
			uv.push_back(uv_i);

			vector3f c;
			lerp_img(c.v, scr_i, img, width, height, 3);
			//pixel0.push_back((c.x + c.y + c.z) / 3);
			pixel0.push_back(c.x);
			//printf_s("pixel0: %g\n", pixel0);
		}
	}

	FILE *prepfile;
	char filename[MAX_PATH];
	sprintf_s(filename, "d:/codes/diff/data/data_lambertian_prep.dat");

#if COMPUTE_PREP == 1
	for (int i = 1; i <= NUM_IMG_PAIRS; i++)
	{
		char filename[MAX_PATH];
		sprintf_s(filename, PATH "img00-%02d-raw.raw", i);

		int width, height;
		std::vector<float> img;
		load_img_float(filename, img, width, height);

		float f;
		matrix4x4f R_tilde, T_tilde;
		vector3f tau_tilde;

		FILE *fp;
		sprintf_s(filename, PATH "img00-%02d-raw.param", i);
		FOPEN(fp, filename, "rb");
		fread(&f, sizeof(float), 1, fp);
		fread(R_tilde.m, sizeof(float) * 16, 1, fp);
		fread(tau_tilde.v, sizeof(float) * 3, 1, fp);
		fclose(fp);
		T_tilde = R_tilde;
		T_tilde._14 = tau_tilde.x;
		T_tilde._24 = tau_tilde.y;
		T_tilde._34 = tau_tilde.z;

		//printf_s("Translation dist = %g; Rot %d\n", tau_tilde.length(), i);
		//for (int mm = 0; mm < 16; mm++)
		//{
		//	printf_s("%g\t", R_tilde.m[mm]);
		//	if (mm % 4 == 0 && mm!=0) printf_s("\n");
		//}
		//printf_s("\n");

		//init
		vector3f vp(0, 0, f);
		matrix4x4f R;
		vector3f tau_pp;
		{
			matrix4x4f T = T_tilde.inversed_matrix();
			R = T;
			R._14 = R._24 = R._34 = 0;
			vector3f tau(T._14, T._24, T._34);
			tau_pp = vp - R.transposed_matrix()*vp + R.transposed_matrix()*tau;
		}

		matrix4x4f T;
		codex::math::matrix::identity(T);
		T._14 = tau_pp.x;
		T._24 = tau_pp.y;
		T._34 = tau_pp.z;

		float beta = 1.0f / f;
		k1.push_back(1 - beta*T._34);
		std::vector<float> coeff_img_A;
		std::vector<float> coeff_img_b;

		for (int i_point = 0; i_point < uv.size(); i_point++)
		{
			float k2 = T._14 + beta*uv[i_point].x*T._34;
			float k3 = T._24 + beta*uv[i_point].y*T._34;
			//printf("k1=%g\n", k1);
			//compute Ia, Ib, Ic
			vector3f gI;

			std::vector<vector3f> A_temp;
			std::vector<float> b_temp;
			float max_du = 0;
			for (int j = 0; j < nsamples; j++)
			{
				float z = z0 + (z1 - z0)*j / nsamples;
				vector2f duv;
				duv.x = k2 / (k1[i - 1] - beta*z);
				duv.y = k3 / (k1[i - 1] - beta*z);

				vector2f uv2 = uv[i_point] + duv;
				vector3f v;
				warp_a_pixel(v, uv2, img, width, height, f, R);
				//max_du = fabs(duv.x) > max_du ? fabs(duv.x) : max_du;

				// check if the point is mapped to empty space
				if ((v - vector3f(1, 0, 1)).length() > 1e-5) {
					A_temp.push_back(vector3f(duv.x, duv.y, 1.0f));
					//b_temp.push_back((v.x + v.y + v.z) / 3);
					b_temp.push_back(v.x);

					// check individual point
/*					{
						if (i_point == 55235) {
							printf("sample %03d:duv=( %.6f %.6f ), pix = %.6f\n",
								j, duv.x, duv.y, (v.x+v.y+v.z)/3);
						}
					}	*/	
				}
				//printf_s("%g %g\n", z, b.v[j]);
			}
			//printf_s("max du:%g\n", max_du);

			int valid_samples = A_temp.size();
			la_matrix<float> A(valid_samples, 3);
			la_vector<float> b(valid_samples);
			for (int j = 0; j < valid_samples; j++) {
				A.m[0 * A.row + j] = A_temp[j].x;
				A.m[1 * A.row + j] = A_temp[j].y;
				A.m[2 * A.row + j] = A_temp[j].z;
				b.v[j] = b_temp[j];
			}
			la_vector<float> x;
			la_matrix<float> A_copy;
			A_copy = A;
			if (valid_samples > 0) {
				slsq(x, A_copy, b, 1.0f);
				gI.x = x.v[0];
				gI.y = x.v[1];
				gI.z = x.v[2];
			}
			else {
				gI = vector3f(0,0,0);
			}
			coeff_img_A.push_back(-k2*gI.x - k3*gI.y);
			coeff_img_b.push_back(gI.z - pixel0[i_point]);

			// this segment lets you check the reconstruction error of fitting (Ia,Ib,Ic)
			//if (i_point == 55235 || i_point == 55695)
			{
				//float z = groundtruth[i_point];
				//vector3f duv;
				//duv.x = k2 / (k1[i - 1] - beta*z);
				//duv.y = k3 / (k1[i - 1] - beta*z);
				//duv.z = 1;
				//float recon = duv*gI;
				//float err = fabs(recon - pixel0[i_point]) / pixel0[i_point] * 100;
				//recon_err[(i - 1)*scr.size() + i_point] = err;
				//if (err > 10.0f) 
				{
					//printf_s("gd=%g,duv=(%g, %g)\n", z, duv.x, duv.y);
					//printf_s("img %d point %d, %d samples: %.6f\t recon: %.6f\t err: %.2f%%\n", i, i_point, valid_samples, pixel0[i_point], recon, err);
					if (valid_samples > 0) {
						la_vector<float> res, diff;
						smvmul(res, A, x);
						svsub(diff, res, b);
						recon_err[(i - 1)*scr.size() + i_point] = svnorm2(diff);
					}
					else {
						recon_err[(i - 1)*scr.size() + i_point] = 10000;
					}

				}
			}					
		}
		coeff_A.push_back(coeff_img_A);
		coeff_b.push_back(coeff_img_b);
		printf("preprocessing done for image pair %d of %d\n", i, NUM_IMG_PAIRS);
	}


	prepfile = fopen(filename, "wb");
	int nk = k1.size();
	int nA = coeff_A.size();
	int nsubA = coeff_A[0].size();	// check possible memory access error later!
	fwrite(&nk, sizeof(int), 1, prepfile);
	fwrite(&nA, sizeof(int), 1, prepfile);
	fwrite(&nsubA, sizeof(int), 1, prepfile);
	for (int ik = 0; ik < k1.size(); ik++) {
		fwrite(&k1[ik], sizeof(float), 1, prepfile);
	}
	for (int iA = 0; iA < coeff_A.size(); iA++) {
		for (int jA = 0; jA < coeff_A[iA].size(); jA++) {
			fwrite(&coeff_A[iA][jA], sizeof(float), 1, prepfile);
		}
	}
	for (int ib = 0; ib < coeff_b.size(); ib++) {
		for (int jb = 0; jb < coeff_b[ib].size(); jb++) {
			fwrite(&coeff_b[ib][jb], sizeof(float), 1, prepfile);
		}
	}
	fwrite(recon_err, sizeof(float), NUM_IMG_PAIRS * scr.size(), prepfile);
	delete[] recon_err;
	fclose(prepfile);
	printf_s("preprocessing complete, written to file\n");
	exit(1);
#else
	prepfile = fopen(filename, "rb");
	int ns[3];
	fread(ns, sizeof(int), 3, prepfile);

	float *k_buffer = new float[ns[0]];
	fread(k_buffer, sizeof(float), ns[0], prepfile);
	for (int ik = 0; ik < ns[0]; ik++) {
		k1.push_back(k_buffer[ik]);
	}
	delete[] k_buffer;

	float *img_buffer = new float[ns[2]];
	for (int i = 0; i < ns[1]; i++) {
		fread(img_buffer, sizeof(float), ns[2], prepfile);
		std::vector<float> vec_buffer(img_buffer, img_buffer+ns[2]);
		coeff_A.push_back(vec_buffer);
	}
	for (int i = 0; i < ns[1]; i++) {
		fread(img_buffer, sizeof(float), ns[2], prepfile);
		std::vector<float> vec_buffer(img_buffer, img_buffer + ns[2]);
		coeff_b.push_back(vec_buffer);
	}
	delete[] img_buffer;
	fread(recon_err, sizeof(float), NUM_IMG_PAIRS * scr.size(), prepfile);

	fclose(prepfile);
	printf_s("parameters read from pre-computed file, solving for depth...\n");

#endif

	FILE *file_obj;
	sprintf_s(filename, "d:/codes/diff/data/lambertian_out.obj");
	if ((file_obj = fopen(filename, "w")) == NULL)
	{
		printf("open file failed!\n");
		exit(1);
	}
	FILE *file_depthmap;
	sprintf_s(filename, "d:/codes/diff/data/lambertian_depthmap.raw");
	if ((file_depthmap = fopen(filename, "wb")) == NULL)
	{
		printf("open file failed!\n");
		exit(1);
	}
	int n_of_point = uv.size();
	float* z_buffer = new float[n_of_point];

	for (int i_point = 0; i_point < n_of_point; i_point++)
	{
		diffuse_func opt_func;
		opt_func.beta = 1.0f / FOCAL_LEN;		
		opt_func.coef.clear();
		opt_func.k1.clear();
		std::vector<float> vref_buffer, A_buffer, b_buffer;
		for (int i_img = 0; i_img < NUM_IMG_PAIRS; i_img++)
		{
			if (recon_err[i_img*scr.size() + i_point] > 1e-1) continue;
			opt_func.coef.push_back(coeff_A[i_img][i_point]);
			opt_func.k1.push_back(k1[i_img]);
			vref_buffer.push_back(coeff_b[i_img][i_point]);
			A_buffer.push_back(coeff_A[i_img][i_point]);
			b_buffer.push_back(coeff_b[i_img][i_point]);
		}
		int n_valid_img = vref_buffer.size();
		float z_est_lm = 0;
		if (n_valid_img == 0) {
			printf_s("point %d has no valid measuerments!\n", i_point);
			z_est_lm = -SHIFT_Z;
		}
		else {
			opt_func_LM<float> opt(1, n_valid_img);
			opt_func.view_num = n_valid_img;
			la_vector<float> vref(n_valid_img);
			la_matrix<float> A_all_pairs(n_valid_img, 1);
			la_vector<float> b_all_pairs(n_valid_img);
			for (int i = 0; i < n_valid_img; i++) {
				A_all_pairs.m[i] = A_buffer[i];
				b_all_pairs.v[i] = b_buffer[i];
				vref.v[i] = vref_buffer[i];
			}

			la_vector<float> z_temp(1);
			z_temp.clear();
			z_temp.v[0] = (z0 + z1) / 2;
			opt.set_f(&opt_func);
			float obj;
			opt_levmar(z_temp, &obj, &opt, vref, 1000);
			z_est_lm = z_temp.v[0];

			float err = 0;
			{
				float cur_z = groundtruth[i_point];
				for (int j = 0; j < n_valid_img; j++)
				{
					float residual = A_all_pairs.m[j] / (opt_func.k1[j] - cur_z / FOCAL_LEN) - b_all_pairs.v[j];
					err += residual * residual;
				}
			}

			//slsq(z_temp, A_all_pairs, b_all_pairs, 1);	// note that slsq will modify A matrix!
			//float z_est_lsq = (1.0 - 1.0 / z_temp.v[0]) * FOCAL_LEN;

			if (z_est_lm < -30)
			{
				printf_s("point %d is bad\n", i_point);
				printf_s("scan=(%d,%d)\tz=%g\tobj=%.6f\tgdz=%g\tgdobj=%.6f\n", scr[i_point].first, scr[i_point].second, z_est_lm, obj, groundtruth[i_point], err);
				printf_s("recon error:\n");
				for (int i = 0; i < NUM_IMG_PAIRS; i++) {
					printf_s("%.3f\t ", recon_err[i*scr.size() + i_point]);
				}
				printf_s("\n");
			}
		}		
		z_buffer[i_point] = z_est_lm;
	}

	// following are the output procedures for .obj and a depthmap
	// 1st pass: building the connectivity
	float* depth_map = new float[IMG_DIM*IMG_DIM];
	int* hash_map = new int[IMG_DIM*IMG_DIM];
	float* z_filtered = new float[n_of_point];
	memset(depth_map, 0, sizeof(depth_map));
	memset(hash_map, 0, sizeof(hash_map));
	for (int i_point = 0; i_point < n_of_point; i_point++)
	{
		int index = scr[i_point].first*IMG_DIM + scr[i_point].second;
		depth_map[index] = z_buffer[i_point];
		hash_map[index] = i_point + 1;
	}
	// 2nd pass: median filtering the depth map
	int kernel_dim = 7;
	for (int i_point = 0; i_point < n_of_point; i_point++)
	{
		int x = scr[i_point].first;
		int y = scr[i_point].second;
		std::vector<float> filter_buffer(kernel_dim*kernel_dim);
		int cnt = 0;
		for (int ix = x - kernel_dim / 2; ix <= x + kernel_dim / 2; ix++)
		{
			for (int iy = y - kernel_dim / 2; iy <= y + kernel_dim / 2; iy++)
			{
				filter_buffer[cnt++] = depth_map[ix*IMG_DIM + iy];
			}
		}
		std::vector<float>::iterator first = filter_buffer.begin();
		std::vector<float>::iterator last = filter_buffer.end();
		std::vector<float>::iterator middle = first + (last - first) / 2;
		std::nth_element(first, middle, last);
		z_filtered[i_point] = *middle;
	}
	delete[] depth_map;

	fwrite(&width, sizeof(int), 1, file_depthmap);
	fwrite(&height, sizeof(int), 1, file_depthmap);
	fwrite(&n_of_point, sizeof(int), 1, file_depthmap);
	for (int i_point = 0; i_point < n_of_point; i_point++)
	{
		{	// output the obj file
			float z = z_filtered[i_point];
			float x = uv[i_point].x * (1 - z / FOCAL_LEN);
			float y = uv[i_point].y * (1 - z / FOCAL_LEN);
			z += SHIFT_Z;
			//if (fabs(x) < 2.2f && fabs(y) < 2.2f && fabs(z) <2.2f)			
			fprintf(file_obj, "v %g %g %g\n", x, y, z);
		}
		{
			// output the depth map
			float z = z_filtered[i_point] + SHIFT_Z;
			fwrite(&scr[i_point].first, sizeof(int), 1, file_depthmap);
			fwrite(&scr[i_point].second, sizeof(int), 1, file_depthmap);
			fwrite(&z, sizeof(float), 1, file_depthmap);
		}
	}
	for (int i_point = 0; i_point < n_of_point; i_point++)
	{
		int x = scr[i_point].first;
		int y = scr[i_point].second;
		if (hash_map[x*IMG_DIM + y] && hash_map[(x + 1)*IMG_DIM + y] && 
			hash_map[x*IMG_DIM + y + 1]) 
		{
			fprintf(file_obj, "f %d %d %d\n", hash_map[x*IMG_DIM + y], 
				hash_map[x*IMG_DIM + y + 1], hash_map[(x + 1)*IMG_DIM + y]);
		}
		if (hash_map[x*IMG_DIM + y] && hash_map[(x - 1)*IMG_DIM + y] &&
			hash_map[x*IMG_DIM + y - 1])
		{
			fprintf(file_obj, "f %d %d %d\n", hash_map[x*IMG_DIM + y],
				hash_map[x*IMG_DIM + y - 1], hash_map[(x - 1)*IMG_DIM + y]);
		}
	}

	delete[] hash_map;
	delete[] z_buffer;
	delete[] z_filtered;
	fclose(file_obj);
	fclose(file_depthmap);
	//la_vector<float> comp1, comp2;
	//opt_func.f(comp1, z_temp);
	//svsub(comp2, comp1, vref);
	//float obj_true = svnorm2(comp2);
	//obj_true *= obj_true;
	//printf("objective: %g\n", obj_true);
}

void test_Specular()
{

	const int nsamples = 1000;
	//	const float z0 = 0.1f, z1 = 1.1f;
	const float z0 = -20.0f, z1 = -19.0f;
	pref_env env;
	env.load("d:/codes/diff/data/uffizi.pref_env");

	// read groundtruth data from file
	char filename_gd[MAX_PATH];
	sprintf_s(filename_gd, "d:/codes/diff/data/groundtruth_depth.dat");
	FILE *groundtruth_file;
	if ((groundtruth_file = fopen(filename_gd, "rb")) == NULL)
	{
		printf("open file failed!\n");
		exit(1);
	}
	int scan_begin, scan_end;
	fread(&scan_begin, sizeof(int), 1, groundtruth_file);
	fread(&scan_end, sizeof(int), 1, groundtruth_file);
	int dim = scan_end - scan_begin + 1;
	double *p = new double[dim * dim];
	fread(p, sizeof(double), dim * dim, groundtruth_file);
	fclose(groundtruth_file);

	// prepare precomputed data structure
	mydata data_prep;
	object_points points_to_solve;
	initialize_samples(points_to_solve, scan_begin, scan_begin, scan_end, scan_end);
	const int n_of_points = points_to_solve.points.size();
	for (int i = 0; i < n_of_points; i++)
	{
		la_vector<float> Ic_Iut(NUM_IMG_PAIRS), k2Ia_k3Ib(NUM_IMG_PAIRS);
		data_prep.Ic_Iut.push_back(Ic_Iut);
		data_prep.k2Ia_k3Ib.push_back(k2Ia_k3Ib);
	}
	data_prep.env = env;
	data_prep.point_set = points_to_solve;

	data_prep.kappa = new float[points_to_solve.points.size()];
	{
		char filename[MAX_PATH];
		sprintf_s(filename, "d:/codes/diff/data/groundtruth_kappa_log.dat");
		FILE *file;
		if ((file = fopen(filename, "rb")) == NULL)
		{
			printf("open file failed!\n");
			exit(1);
		}
		fread(data_prep.kappa, sizeof(float), points_to_solve.points.size(), file);
		fclose(file);
	}
	
#if COMPUTE_PREP == 1

	//data_prep.kappa = 2.5;
	printf_s("start preprocessing...\n");
	std::vector<float> pixel0;
	std::vector<vector2f> scr;
	{
		char filename[MAX_PATH];
		sprintf_s(filename, PATH "img00-%02d.raw", 0);

		int width, height;
		std::vector<float> img;
		load_img_float(filename, img, width, height);

		for (int i_point = 0; i_point < n_of_points; i_point++)
		{
			vector2f scr_i;
			uv2screen(scr_i, points_to_solve.points[i_point].uv, width, height);
			scr.push_back(scr_i);
			//printf("%g %g\n", scr.x, scr.y);

			vector3f c;
			lerp_img(c.v, scr_i, img, width, height, 3);
			pixel0.push_back((c.x + c.y + c.z) / 3);
			//printf_s("pixel0: %g\n", pixel0);
		}
	}

	for (int i = 1; i <= NUM_IMG_PAIRS; i++)
	{
		char filename[MAX_PATH];
		sprintf_s(filename, PATH "img00-%02d-raw.raw", i);

		int width, height;
		std::vector<float> img;
		load_img_float(filename, img, width, height);
		printf_s("Handling %s...\n", filename);

		float f;
		matrix4x4f R_tilde, T_tilde;
		vector3f tau_tilde;
		
		FILE *fp;
		sprintf_s(filename, PATH "img00-%02d-raw.param", i);
		FOPEN(fp, filename, "rb");
		fread(&f, sizeof(float), 1, fp);
		fread(R_tilde.m, sizeof(float) * 16, 1, fp);
		fread(tau_tilde.v, sizeof(float) * 3, 1, fp);
		fclose(fp);
		T_tilde = R_tilde;
		T_tilde._14 = tau_tilde.x;
		T_tilde._24 = tau_tilde.y;
		T_tilde._34 = tau_tilde.z;

		data_prep.Delta_v.push_back(tau_tilde);

		//printf_s("Translation dist = %g; Rot %d\n", tau_tilde.length(), i);
		//for (int mm = 0; mm < 16; mm++)
		//{
		//	printf_s("%g\t", R_tilde.m[mm]);
		//	if (mm % 4 == 0 && mm!=0) printf_s("\n");
		//}
		//printf_s("\n");

		//init
		vector3f vp(0, 0, f);
		matrix4x4f R;
		vector3f tau_pp;
		{
			matrix4x4f T = T_tilde.inversed_matrix();
			R = T;
			R._14 = R._24 = R._34 = 0;
			vector3f tau(T._14, T._24, T._34);
			tau_pp = vp - R.transposed_matrix()*vp + R.transposed_matrix()*tau;
		}

		matrix4x4f T;
		codex::math::matrix::identity(T);
		T._14 = tau_pp.x;
		T._24 = tau_pp.y;
		T._34 = tau_pp.z;

		float beta = 1.0f / f;
		float k1 = 1 - beta * T._34;
		data_prep.k1.push_back(k1);
		//k1.push_back(1 - beta*T._34);

		printf_s("Getting local approximations...\n");
		for (int i_point = 0; i_point < n_of_points; i_point++)
		{
			float k2 = T._14 + beta*points_to_solve.points[i_point].uv.x*T._34;
			float k3 = T._24 + beta*points_to_solve.points[i_point].uv.y*T._34;
			//printf("k1=%g\n", k1);
			//compute Ia, Ib, Ic
			vector3f gI;
			la_matrix<float> A(nsamples, 3);
			la_vector<float> b(nsamples);
			float max_du = 0;
			for (int j = 0; j < nsamples; j++)
			{
				float z = z0 + (z1 - z0)*j / nsamples;
				vector2f duv;
				duv.x = k2 / (k1 - beta*z);
				duv.y = k3 / (k1 - beta*z);
				max_du = fabs(duv.x) > max_du ? fabs(duv.x) : max_du;
				A.m[0 * A.row + j] = duv.x;
				A.m[1 * A.row + j] = duv.y;
				A.m[2 * A.row + j] = 1.0f;

				vector2f uv2 = points_to_solve.points[i_point].uv + duv;
				vector3f v;
				warp_a_pixel(v, uv2, img, width, height, f, R);
				b.v[j] = (v.x + v.y + v.z) / 3;
				//printf_s("%g %g\n", z, b.v[j]);
			}
			//printf_s("max du:%g\n", max_du);
			la_vector<float> x;
			slsq(x, A, b, 1);
			gI.x = x.v[0];
			gI.y = x.v[1];
			gI.z = x.v[2];	
			data_prep.k2Ia_k3Ib[i_point].v[i - 1] = k2*gI.x + k3*gI.y;
			data_prep.Ic_Iut[i_point].v[i - 1] = gI.z - pixel0[i_point];
			if (isnormal(data_prep.k2Ia_k3Ib[i_point].v[i - 1]) == 0 && 
				isnormal(data_prep.Ic_Iut[i_point].v[i - 1]))
			{
				printf_s("point:%d img:%d\nk2=%g\tk3=%g\tgI=(%g,%g,%g)\nk2gI.x+k3gI.y=%g\tgI.z-I=%g\n", 
					i_point, i, k2, k3, gI.x, gI.y, gI.z, data_prep.k2Ia_k3Ib[i_point].v[i - 1],
					data_prep.Ic_Iut[i_point].v[i - 1]);
			}
		}		
		printf_s("\n");
	}
	char filename[MAX_PATH];
	sprintf_s(filename, "d:/codes/diff/data/data_prep.dat");
	save_mydata_file(filename, data_prep);

	printf_s("done preprocessing, solving for depth:\n");
	exit(0);
#else
	char filename[MAX_PATH];
	sprintf_s(filename, "d:/codes/diff/data/data_prep.dat");
	read_mydata_file(filename, data_prep);

	//for (int i = 0; i < 30; i++)
	//{
	//	printf_s("%f\n", data_prep.k2Ia_k3Ib[i].v[0]);
	//}
	printf_s("precomputed data read from file, solving for depth:\n");
#endif

	const int nobs = n_of_points * 2;
	int Jnn = n_of_points * 6;
	double *initial_p = new double[n_of_points];


	double jitters[6] = { 1e-1,1e-2,1e-3,1e-4,1e-5,0 };
	for (int i_test = 0; i_test < 6; i_test++)
	{
		double rms_err_init = 0;
		for (int i = 0; i < n_of_points; i++)
		{
			initial_p[i] = p[i];
			initial_p[i] += rand_range(jitters[i_test]);
			rms_err_init += (initial_p[i] - p[i]) * (initial_p[i] - p[i]);
			//printf_s("%g\n", initial_p[i]);
		}
		rms_err_init = sqrt(rms_err_init / n_of_points);

		double opts[SPLM_OPTS_SZ], info[SPLM_INFO_SZ];
		opts[0] = SPLM_INIT_MU; opts[1] = SPLM_STOP_THRESH; opts[2] = SPLM_STOP_THRESH;
		opts[3] = SPLM_STOP_THRESH;
		opts[4] = 1e-6;
		opts[5] = SPLM_CHOLMOD; // use CHOLMOD

		//double *err = new double[nobs];
		//sparselm_chkjaccrs(ffunc, fjac, p, n_of_points, nobs, Jnn, (void *)&data_prep, err);
		//for (int i = 0; i < nobs; i++)
		//{
		//	printf_s("%g\n", err[i]);
		//}
		//exit(1);

		//calc_pi_ls_err(initial_p, 0, (void *)&data_prep);
		int n_iter = sparselm_dercrs(ffunc, fjac, initial_p, NULL, n_of_points,
			0, nobs, Jnn, -1, 1000, opts, info, (void *)&data_prep);
		//calc_pi_ls_err(initial_p, 0, (void *)&data_prep);

		//int n_iter = sparselm_difcrs(ffunc, fjac_diff, initial_p, NULL, n_of_points,
		//		0, nobs, Jnn, -1, 1000, opts, info, (void *)&data_prep);
		//int n_iter = sparselm_dercrs(ffunc, fjac_mannual_diff, initial_p, NULL, n_of_points,
		//		0, nobs, Jnn, -1, 1000, opts, info, (void *)&data_prep);

		//printf_s("Optimization info:\n");
		//for (int i = 0; i < SPLM_INFO_SZ; i++)
		//{
		//	printf_s("info[%d]: %g\n", i, info[i]);
		//}

		double rms_err = 0;
		for (int i = 0; i < n_of_points; i++)
		{
			rms_err += (initial_p[i] - p[i]) * (initial_p[i] - p[i]);
			//printf_s("%.6f\t%.6f\n", initial_p[i], p[i]);
		}
		rms_err = sqrt(rms_err / n_of_points);
		double init_cost = sqrt(info[0] / 2 / n_of_points);
		double terminal_cost = sqrt(info[1] / 2 / n_of_points);
		//printf_s("Root mean square error of depth:%g\n");
		printf_s("%.6f\t%.6f\t%g\t%.6f\n", init_cost, terminal_cost, rms_err_init, rms_err);
	}
}


void test_gauss_interp() {
	pref_env gauss_interp;
	float y[5] = { 0.0001,   0.0004,   0.0009,   0.0016,   0.0025 };
	int x[5] = { 0.01,0.02,0.03,0.04,0.05 };
	for (int i = 1; i < 10; i++) {
		float delta = 0.01f * i / 10;
		float val = gauss_interp.gauss_interp(x, y, delta/0.01);
		printf_s("%g\t%g\n", delta + 0.04f, val);
	}
	exit(1);
}

void main()
{
	//test();
	//test_gauss_interp();
	//prec_brdf();

	//test_view_diff3();
	//test_filtered_1st_order();
	//synthesize_images();
	//return;


	std::vector<std::pair<int,int> > points_under_test;

	char filename[MAX_PATH];
	sprintf_s(filename, "d:/codes/diff/data/groundtruth_lambertian.dat");
	FILE *file_lambertian;
	if ((file_lambertian = fopen(filename, "rb")) == NULL)
	{
		printf("open file failed!\n");
		exit(1);
	}
	std::vector<float> groundtruth;
	int n_valid;
	fread(&n_valid, sizeof(int), 1, file_lambertian);
	for (int i = 0; i < n_valid; i++) {
		float gd_z;
		int xy[2];
		fread(xy, sizeof(int), 2, file_lambertian);
		fread(&gd_z, sizeof(float), 1, file_lambertian);
		points_under_test.push_back(std::pair<int, int>(xy[0], xy[1]));
		groundtruth.push_back(gd_z);
	}

	test_Lambertian_v1(points_under_test, groundtruth);
	//test_Specular();

	return;
	
	//test_LM_v2();

	////warp_cam();
	////single_diff_analysis();
	//stablize_cam(vector3f(0, 0, 0));
	//solve_diff();
}