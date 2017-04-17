#include "main.h"

const int nsamples = 20;
const float z0 = -0.1f, z1 = 1.1f;

void test_Lambertian_v1()
{
	vector2f uv;
	uv.x = 0.00195312f;
	uv.y = -0.00195312f;

	//DEBUG
	float pixel0;
	{
		char filename[MAX_PATH];
		sprintf_s(filename, PATH "img00-%02d.raw", 0);

		int width, height;
		std::vector<float> img;
		load_img_float(filename, img, width, height);
		
		vector2f scr;
		uv2screen(scr, uv, width, height);

		vector3f c;
		lerp_img(c.v, scr, img, width, height, 3);
		pixel0 = (c.x + c.y + c.z) / 3;
		printf_s("pixel0: %g\n", pixel0);
	}

	float avg = 0;

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
		identity(T);
		T._14 = tau_pp.x;
		T._24 = tau_pp.y;
		T._34 = tau_pp.z;

		float beta = 1.0f / f;
		float k1 = 1 - beta*T._34;
		float k2 = T._14 + beta*uv.x*T._34;
		float k3 = T._24 + beta*uv.y*T._34;

		//compute Ia, Ib, Ic
		vector3f gI;

		la_matrix<float> A(nsamples, 3);
		la_vector<float> b(nsamples);
		for (int j = 0; j < nsamples; j++)
		{
			float z = z0 + (z1 - z0)*j / nsamples;
			vector2f duv;
			duv.x = k2 / (k1 - beta*z);
			duv.y = k3 / (k1 - beta*z);

			A.m[0 * A.row + j] = duv.x;
			A.m[1 * A.row + j] = duv.y;
			A.m[2 * A.row + j] = 1.0f;

			vector2f uv2 = uv + duv;
			vector3f v;
			warp_a_pixel(v, uv2, img, width, height, f, R);
			b.v[j] = (v.x + v.y + v.z) / 3;
			//printf_s("%g %g\n", z, b.v[j]);
		}
		//printf_s("\n");

		la_vector<float> x;
		slsq(x, A, b, 1);
		gI.x = x.v[0];
		gI.y = x.v[1];
		gI.z = x.v[2];

		//{
		//	float z = 1;
		//	vector3f duv;
		//	duv.x = k2 / (k1 - beta*z);
		//	duv.y = k3 / (k1 - beta*z);
		//	duv.z = 1;
		//	printf_s("recon: %g\n", duv*gI);
		//}

		//compute
		//printf_s("Ia Ib Ic = [%g %g %g]\n", gI.x, gI.y, gI.z);

		{
			float p = (gI.z - pixel0) / (-k2*gI.x - k3*gI.y);
			printf_s("z = %g\n", (k1 - 1 / p)/beta);

			avg += p;
		}
	}

	//avg /= NUM_IMG_PAIRS;
	//printf_s("%g\n", FOCAL_LEN*(1 - 1 / avg));
}

class diff_form_func : public func<float>
{
public:
	std::vector<float>	coef;
	std::vector<int>	viewid;
	float beta;

	//the function
	virtual void f(la_vector<float> &fx, const la_vector<float> &x)
	{
		fx.init(viewid.size());
		for (int i = 0; i < fx.length; i++)
		{
			int j = viewid[i];
			fx.v[i] = 
				coef[i * 4 + 0] * x.v[j*3+0] + coef[i * 4 + 1] * x.v[j * 3 + 1] +
				x.v[j * 3 + 2] + coef[i * 4 + 3] / (coef[i * 4 + 2] - beta*x.v[x.length-1]);
		}
	}

	//Jacobian
	virtual void j(la_vector<float> &jx, const la_vector<float> &x)
	{
		jx.init(viewid.size()*x.length);
		jx.clear();
		for (int i = 0; i < viewid.size(); i++)
		{
			int j = viewid[i];

			jx.v[i * x.length + j * 3 + 0] += coef[i * 4 + 0];
			jx.v[i * x.length + j * 3 + 1] += coef[i * 4 + 1];
			jx.v[i * x.length + j * 3 + 2] += 1;
			float temp = (coef[i * 4 + 2] - beta*x.v[x.length-1]);
			jx.v[i * x.length + x.length - 1] = coef[i * 4 + 3] * beta / (temp*temp);
		}
	}
};

class z_solver_data
{
public:
	bool				b_valid;
	diff_form_func		func;
	std::vector<float>	ref;
	vector2f			uv;
};

void test_LM_v1()
{
	vector2f uv;
	uv.x = 0;// 0.687500f;// 0.00195312f;
	uv.y = 0;// -0.00195312f;

	const int dim = 8;
	std::vector<z_solver_data> solver_data;

	//DEBUG
	float pixel0;
	{
		char filename[MAX_PATH];
		sprintf_s(filename, PATH "img00-%02d.raw", 0);

		int width, height;
		std::vector<float> img;
		load_img_float(filename, img, width, height);

		vector2f scr;
		uv2screen(scr, uv, width, height);
		vector3f c;
		lerp_img(c.v, scr, img, width, height, 3);
		pixel0 = (c.x + c.y + c.z) / 3;
		printf_s("pixel0: %g\n", pixel0);

		//solver_data.resize(dim*dim);
		solver_data.resize(dim);
		//for (int y = 0; y < dim; y++)
			for (int x = 0; x < dim; x++)
			{
				//int idx = x +y*dim;
				//solver_data[idx].uv = vector2f((float)x / dim * 2 - 1, (float)y / dim * 2 - 1);
				int idx = x;
				solver_data[idx].uv = vector2f((float)x / dim * 2 - 1, 0);
				printf_s("%f %f\n", solver_data[idx].uv.x, solver_data[idx].uv.y);

				vector2f scr;
				uv2screen(scr, solver_data[idx].uv, width, height);
				vector3f c;
				lerp_img(c.v, scr, img, width, height, 3);

				if (c.x == 1 && c.y == 0 && c.z == 1)
					solver_data[idx].b_valid = false;
				else {
					solver_data[idx].b_valid = true;
					solver_data[idx].func.beta = 1.0f / FOCAL_LEN;
					solver_data[idx].func.coef.clear();
				}
			}

		{
			int idx = 0;

			solver_data.resize(1);
			solver_data[idx].uv = uv;

			solver_data[idx].b_valid = true;
			solver_data[idx].func.beta = 1.0f / FOCAL_LEN;
			solver_data[idx].func.coef.clear();
		}
	}

	//diff_form_func func;
	//func.beta = 1.0f / FOCAL_LEN;
	//func.coef.clear();
	//std::vector<float> ref;

	float gt_err = 0;

	for (int view = 0; view < NUM_IMG_VIEW; view++)
		for (int i = 1; i <= NUM_IMG_PAIRS; i++)
		{
			char filename[MAX_PATH];
			sprintf_s(filename, PATH "img%02d-%02d-raw.raw", view, i);

			int width, height;
			std::vector<float> img;
			load_img_float(filename, img, width, height);

			float f;
			matrix4x4f R_tilde, T_tilde;
			vector3f tau_tilde;

			FILE *fp;
			sprintf_s(filename, PATH "img%02d-%02d-raw.param", view, i);
			FOPEN(fp, filename, "rb");
			fread(&f, sizeof(float), 1, fp);
			fread(R_tilde.m, sizeof(float) * 16, 1, fp);
			fread(tau_tilde.v, sizeof(float) * 3, 1, fp);
			fclose(fp);
			T_tilde = R_tilde;
			T_tilde._14 = tau_tilde.x;
			T_tilde._24 = tau_tilde.y;
			T_tilde._34 = tau_tilde.z;

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
			identity(T);
			T._14 = tau_pp.x;
			T._24 = tau_pp.y;
			T._34 = tau_pp.z;

			for (int k = 0; k < solver_data.size(); k++)
			{
				z_solver_data &sd = solver_data[k];

				if (!sd.b_valid) continue;

				float beta = 1.0f / f;
				float k1 = 1 - beta*T._34;
				float k2 = T._14 + beta*sd.uv.x*T._34;
				float k3 = T._24 + beta*sd.uv.y*T._34;

				//compute Ia, Ib, Ic
				vector3f gI;

				std::vector<float> m_svd, b_svd;
				for (int j = 0; j < nsamples; j++)
				{
					float z = z0 + (z1 - z0)*j / nsamples;
					vector2f duv;
					duv.x = k2 / (k1 - beta*z);
					duv.y = k3 / (k1 - beta*z);

					vector2f uv2 = sd.uv + duv;
					vector3f v;
					bool b_succ = warp_a_pixel(v, uv2, img, width, height, f, R);

					if (b_succ && !(v.x == 1 && v.y == 0 && v.z == 1))
					{
						m_svd.push_back(duv.x);
						m_svd.push_back(duv.y);

						float vv = (v.x + v.y + v.z) / 3;
						b_svd.push_back(vv);
						//printf_s("%f %f %f\n", vv, duv.x, duv.y);
					}
				}

				if (b_svd.size() > nsamples/2)
				{
					la_matrix<float> A(b_svd.size(), 3);
					la_vector<float> b(b_svd);
					for (int j = 0; j < b.length; j++)
					{
						A.m[0 * A.row + j] = m_svd[j * 2 + 0];
						A.m[1 * A.row + j] = m_svd[j * 2 + 1];
						A.m[2 * A.row + j] = 1;
					}

					la_matrix<float> ATA;
					la_vector<float> ATb;
					smmmul(ATA, A.transpose(), A);
					smvmul(ATb, A.transpose(), b);

					la_vector<float> x;
					slsq(x, ATA, ATb, 1);
					gI.x = x.v[0];
					gI.y = x.v[1];
					gI.z = x.v[2];

					//FIX ME: empirical
					if (abs(gI.x) < 10000 && abs(gI.y) < 10000 && abs(gI.z) < 10000)
					{
						sd.func.coef.push_back(tau_tilde.x + beta*sd.uv.x*tau_tilde.z);
						sd.func.coef.push_back(tau_tilde.y + beta*sd.uv.y*tau_tilde.z);
						sd.func.coef.push_back(k1);
						sd.func.coef.push_back(-k2*gI.x - k3*gI.y);
						sd.ref.push_back(gI.z);
						sd.func.viewid.push_back(view);

						//DEBUG
						{
							float z = 1.f;
							vector2f duv;
							duv.x = k2 / (k1 - beta*z);
							duv.y = k3 / (k1 - beta*z);

							float	rho_a = 0.006223f,
									rho_b = 0.012199f,
									rho_d = 0.151765f;
							
							rho_a = -0.002574f;
							rho_b = 0.040749f;
							rho_d = 0.128326f;
							z = 0.97706f;
							//rho_a = 0.0128953f;
							//rho_b = 0.00629066f;
							//rho_d = 0.151853f;
							//z = 1.26f;

							vector2f uv2 = sd.uv + duv;
							vector3f v;
							warp_a_pixel(v, uv2, img, width, height, f, R, true);
							printf_s("#%d:%d gI = [%f %f %f]\n", view, i, gI.x, gI.y, gI.z);

							float gt = (tau_tilde.x + beta*sd.uv.x*tau_tilde.z)*rho_a +
								(tau_tilde.y + beta*sd.uv.y*tau_tilde.z)*rho_b +
								rho_d +
								(-k2*gI.x - k3*gI.y) / (k1 - beta*z);
							printf_s("%f %f %f : %f = %f (%f%%)\n",
								tau_tilde.x + beta*sd.uv.x*tau_tilde.z,
								tau_tilde.y + beta*sd.uv.y*tau_tilde.z,
								-k2*gI.x - k3*gI.y,
								gI.z,
								gt,
								abs(gt-gI.z)/ gI.z*100
							);

							gt_err += (gt - gI.z)*(gt - gI.z);
						}
					}

					////DEBUG
					//printf_s("---\n");
					//for (int j = 0; j < nsamples; j++)
					//{
					//	float z = z0 + (z1 - z0)*j / nsamples;
					//	vector2f duv;
					//	duv.x = k2 / (k1 - beta*z);
					//	duv.y = k3 / (k1 - beta*z);

					//	vector2f uv2 = sd.uv + duv;
					//	vector3f v;
					//	bool b_succ = warp_a_pixel(v, uv2, img, width, height, f, R);

					//	if (b_succ && !(v.x == 1 && v.y == 0 && v.z == 1))
					//	{
					//		printf_s("%f\n", duv.x*gI.x+duv.y*gI.y+gI.z);
					//	}
					//}
					//printf_s("---\n");
				}
			}
		}

		std::vector<vector3f> pos;

		for (int k = 0; k < solver_data.size(); k++)
		{
			z_solver_data &sd = solver_data[k];

			if (sd.b_valid)
			{
				//DEBUG
				{
					la_matrix<float> A(sd.ref.size(), 4);
					for (int i = 0; i < A.row; i++)
					{
						A.m[i + 0 * A.row] = sd.func.coef[i * 4 + 0];
						A.m[i + 1 * A.row] = sd.func.coef[i * 4 + 1];
						A.m[i + 2 * A.row] = 1;
						A.m[i + 3 * A.row] = sd.func.coef[i * 4 + 3];
					}

					la_matrix<float> ATA;
					smmmul(ATA, A.transpose(), A);
					printf_s("cond = %g\n", smcond(ATA));
				}

				la_vector<float> x(3 * NUM_IMG_VIEW + 1);
				x.clear();

				la_vector<float> vref(sd.ref);
				opt_func_LM<float> opt(x.length, vref.length);
				opt.set_f(&sd.func);

				float obj;
				opt_levmar(x, &obj, &opt, vref, 1000);
				printf_s("%f\n", x.v[x.length - 1]);

				printf_s("%g : %g\n", obj, x.v[x.length - 1]);
				for (int i = 0; i < NUM_IMG_VIEW; i++)
					printf_s("%g %g %g\n", x.v[i * 3 + 0], x.v[i * 3 + 1], x.v[i * 3 + 2]);

				vector3f p(sd.uv.x, sd.uv.y, x.v[x.length - 1]);
				if (!_isnanf(p.x) && !_isnanf(p.y) && !_isnanf(p.z))
					pos.push_back(p);
			} else {
				printf_s("%f\n", -1.0f);
			}
		}

	printf_s("gt err = %f\n", gt_err);

	export_pointcloud_ply("d:/temp/diffstereo/result.ply", pos);
}