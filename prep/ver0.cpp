#include "main.h"

class eq_coef
{
public:
	vector3f	d_v;
	vector2f	I_u;
	float		m1, m2;
	float		k1, k2, k3, k4;
	float		alpha1, alpha2;
	float		d_I;
};

void stablize_cam(const vector3f &center)
{
	for (int i = 1; i <= NUM_IMG_PAIRS; i++)
	{
		char filename[MAX_PATH];
		sprintf_s(filename, PATH "img00-%02d-raw.raw", i);

		int width, height;
		std::vector<float> img0;
		load_img_float(filename, img0, width, height);

		float f;
		matrix4x4f R_tilde;
		vector3f tau_tilde;

		FILE *fp;
		sprintf_s(filename, PATH "img00-%02d-raw.param", i);
		FOPEN(fp, filename, "rb");
		fread(&f, sizeof(float), 1, fp);
		fread(R_tilde.m, sizeof(float) * 16, 1, fp);
		fread(tau_tilde.v, sizeof(float) * 3, 1, fp);
		fclose(fp);

		vector3f p(0, 0, f);
		vector3f pp, v;
		vector3f dx(1, 0, 0), dy(0, 1, 0), dz(0, 0, 1);
		v = (center - p).normalized_vector();
		vector3f ray0(v*dx, v*dy, v*dz);

		vector3f dx2, dy2, dz2;
		dx2 = R_tilde*dx;
		dy2 = R_tilde*dy;
		dz2 = R_tilde*dz;

		pp = R_tilde*p + tau_tilde;
		v = (center - pp).normalized_vector();
		vector3f ray1(v*dx2, v*dy2, v*dz2);

		vector3f axis = (ray1^ray0).normalized_vector();
		matrix4x4f rot;
		rotation_axis(rot, axis, acos(ray1*ray0));
		rot.transpose();

		matrix4x4f new_R_tilde = rot*R_tilde;

		vector3f new_tau_tilde;
		std::vector<float> img;
		warp_view(img, new_R_tilde, new_tau_tilde, img0, width, height, f, R_tilde, tau_tilde);

		sprintf_s(filename, PATH "img00-%02d.raw", i);
		save_img_float(filename, &img[0], width, height);

		sprintf_s(filename, PATH "img00-%02d.param", i);
		FOPEN(fp, filename, "wb");
		fwrite(&f, sizeof(float), 1, fp);
		fwrite(new_R_tilde.m, sizeof(float) * 16, 1, fp);
		fwrite(new_tau_tilde.v, sizeof(float) * 3, 1, fp);
		fclose(fp);
	}
}

void compute_eq(eq_coef &coef, const vector2f &uv, const int id)
{
	float f;
	matrix4x4f R_tilde;
	vector3f tau_tilde;

	char filename[MAX_PATH];

	FILE *fp;
	sprintf_s(filename, PATH "img00-%02d.param", id);
	FOPEN(fp, filename, "rb");
	fread(&f, sizeof(float), 1, fp);
	fread(R_tilde.m, sizeof(float) * 16, 1, fp);
	fread(tau_tilde.v, sizeof(float) * 3, 1, fp);
	fclose(fp);

	int width, height;
	std::vector<float> img0, img1;
	load_img_float(PATH "img00-00.raw", img0, width, height);
	sprintf_s(filename, PATH "img00-%02d.raw", id);
	load_img_float(filename, img1, width, height);

	float beta = 1.0f / f;
	vector3f p(0, 0, f);

	matrix4x4f R;
	vector3f tau;
	R = R_tilde.transposed_matrix();
	tau = -(R*tau_tilde);

	matrix4x4f E;
	identity(E);
	matrix4x4f RmE = R - E;

	float m1 = beta*uv.x;
	float m2 = beta*uv.y;
	float l1 = RmE._11*uv.x + RmE._12*uv.y + tau.x;
	float l3 = RmE._21*uv.x + RmE._22*uv.y + tau.y;
	float l5 = RmE._31*uv.x + RmE._32*uv.y + tau.z;
	float l2 = RmE._13 - RmE._11*m1 - RmE._12*m2;
	float l4 = RmE._23 - RmE._21*m1 - RmE._22*m2;
	float l6 = RmE._33 - RmE._31*m1 - RmE._32*m2;

	float alpha1 = 1 - beta*l5;
	float alpha2 = -beta*(1 + l6);
	float k1 = (l2 + m1*l6) / alpha2;
	float k3 = (l4 + m2*l6) / alpha2;
	float k2 = l1 + m1*l5 - alpha1 / alpha2*(l2 + m1*l6);
	float k4 = l3 + m2*l5 - alpha1 / alpha2*(l4 + m2*l6);

	coef.alpha1 = alpha1;
	coef.alpha2 = alpha2;
	coef.k1 = k1;
	coef.k2 = k2;
	coef.k3 = k3;
	coef.k4 = k4;
	coef.m1 = m1;
	coef.m2 = m2;

	vector2f scr;
	uv2screen(scr, uv, width, height);
	vector3f I0, I1;
	lerp_img(I0.v, scr, img0, width, height, 3);
	lerp_img(I1.v, scr, img1, width, height, 3);
	coef.d_I = I1.y - I0.y;

	vector3f grad_I;
	lerp_img(I0.v, scr + vector2f(-1, 0), img1, width, height, 3);
	lerp_img(I1.v, scr + vector2f(1, 0), img1, width, height, 3);
	grad_I.x = (I1.y - I0.y) / 2;
	lerp_img(I0.v, scr + vector2f(0, -1), img1, width, height, 3);
	lerp_img(I1.v, scr + vector2f(0, 1), img1, width, height, 3);
	grad_I.y = (I1.y - I0.y) / 2;
	coef.I_u.x = grad_I.x * (width / 2);
	coef.I_u.y = -grad_I.y * (height / 2);

	coef.d_v = R_tilde*p - p + tau_tilde;
	//printf_s("%g %g\n", alpha1, alpha2);

	//
	float z0 = -0.5f, z1 = 1.5f;
	//printf_s("%g [%g %g] ", -c.alpha1 / c.alpha2, c.k1 + c.k2 / (c.alpha1 + c.alpha2*z1), c.k1 + c.k2 / (c.alpha1 + c.alpha2*z0));
	//printf_s("[%g %g]\n", c.k3 + c.k4 / (c.alpha1 + c.alpha2*z1), c.k3 + c.k4 / (c.alpha1 + c.alpha2*z0));

	{
		vector2f scr;
		uv2screen(scr, uv, width, height);
		lerp_img(I0.v, scr, img1, width, height, 3);
	}

	const int nsamples = 20;
	la_matrix<float> A(nsamples, 2);
	la_vector<float> b(nsamples);
	for (int j = 0; j < nsamples; j++)
	{
		float z = z0 + (z1 - z0)*j / nsamples;
		vector2f duv;
		duv.x = k1 + k2 / (alpha1 + alpha2*z);
		duv.y = k3 + k4 / (alpha1 + alpha2*z);

		vector2f scr;
		uv2screen(scr, uv + duv, width, height);
		vector3f v;
		lerp_img(v.v, scr, img1, width, height, 3);
		b.v[j] = v.y - I0.y;
		//printf_s("%g\n", v.x + v.y + v.z);
		A.m[j + 0 * A.row] = duv.x;
		A.m[j + 1 * A.row] = duv.y;
		//printf_s("%g\n", v.y-I0.y);
		//printf_s("%g\n", coef.I_u.x*duv.x + coef.I_u.y*duv.y);
	}

	la_vector<float> x;
	slsq(x, A, b, 1);

	coef.I_u.x = x.v[0];
	coef.I_u.y = x.v[1];
	//float err = 0;
	//for (int j = 0; j < nsamples; j++)
	//{
	//	float z = z0 + (z1 - z0)*j / nsamples;
	//	vector2f duv;
	//	duv.x = k1 + k2 / (alpha1 + alpha2*z);
	//	duv.y = k3 + k4 / (alpha1 + alpha2*z);

	//	vector2f scr;
	//	uv2screen(scr, uv + duv, width, height);
	//	vector3f v;
	//	lerp_img(v.v, scr, img1, width, height, 3);

	//	printf_s("%g %g\n", x.v[0] * duv.x + x.v[1] * duv.y, v.y - I0.y);

	//	//coef.I_u.x*duv.x + coef.I_u.y*duv.y
	//	float e = v.y - I0.y - (coef.I_u.x*duv.x + coef.I_u.y*duv.y);
	//	//float e = v.y - I0.y - (x.v[0] * duv.x + x.v[1] * duv.y);
	//	err += e*e;
	//}
	//printf_s("e = %g\n", err);

	//throw 1;
}

void solve_diff()
{
	const int num = NUM_IMG_PAIRS;
	vector2f uv;
	uv.x = 0.00195312f;
	uv.y = -0.00195312f;

	la_matrix<float> A(num, 3);
	la_vector<float> b(num);

	std::vector<eq_coef> coef;
	for (int i = 1; i <= num; i++)
	{
		eq_coef c;
		compute_eq(c, uv, i);
		coef.push_back(c);
	}

	float zz = 0;
	for (int i = 1; i <= num; i++)
	{
		const eq_coef& c = coef[i - 1];
		A.m[(i - 1) + 0 * A.row] = c.d_v.x + c.m1*c.d_v.z;
		A.m[(i - 1) + 1 * A.row] = c.d_v.y + c.m2*c.d_v.z;
		A.m[(i - 1) + 2 * A.row] = -c.I_u.x*c.k2 - c.I_u.y*c.k4;
		b.v[i - 1] = c.I_u.x*c.k1 + c.I_u.y*c.k3 + c.d_I;

		float temp = b.v[i - 1] / A.m[(i - 1) + 2 * A.row];
		printf_s("%g\n", temp);
		zz += (1 / temp - coef[i - 1].alpha1) / coef[i - 1].alpha2;
	}

	zz /= num;
	printf_s("Lambertian z = %g\n", zz);

	la_vector<float> x;
	la_matrix<float> T;
	T = A;
	slsq(x, T, b);
	printf_s("\n");
	for (int i = 0; i < x.length; i++)
		printf_s("%g ", x.v[i]);
	printf_s("\n");

	float res = 0;
	la_vector<float> Ax, Axb;
	smvmul(Ax, A, x);
	svsub(Axb, Ax, b);
	res = svnorm2(Axb);
	printf_s("res : %g\n", res);

	for (int i = 0; i < coef.size(); i++)
		printf_s("z = %g\n", (1 / x.v[2] - coef[i].alpha1) / coef[i].alpha2);

	la_vector<float> xx(3);
	res = 0;
	xx.v[0] = 0.0f;
	xx.v[1] = 0.0f;
	xx.v[2] = 1.0f / (coef[0].alpha1 + coef[0].alpha2 * 1);

	printf_s("\n");
	for (int i = 0; i < xx.length; i++)
		printf_s("%g ", xx.v[i]);
	printf_s("\n");
	printf_s("%g %g\n", coef[0].alpha1, coef[0].alpha2);

	smvmul(Ax, A, xx);
	svsub(Axb, Ax, b);
	res = svnorm2(Axb);
	printf_s("res : %g\n", res);


	la_matrix<float> U, Vt;
	la_vector<float> sigma;
	T = A;
	ssvd(U, Vt, sigma, T, 1.0f);
	printf_s("\nsigma : ");
	for (int i = 0; i < sigma.length; i++)
		printf_s("%g ", sigma.v[i]);
	printf_s("\n");

	for (int i = 0; i < 3; i++)
		printf_s("[%f %f %f]\n", Vt.m[0 + i * 3], Vt.m[1 + i * 3], Vt.m[2 + i * 3]);

	{
		la_vector<float> temp;
		svsub(temp, xx, x);
		for (int i = 0; i < 3; i++)
			printf_s("%g %g\n", svdot(temp, Vt.column_vec(i)), svdot(temp, Vt.column_vec(i)) / sigma.v[i]);

		svaxpy(svdot(temp, Vt.column_vec(2)), temp, x);
		for (int i = 0; i < xx.length; i++)
			printf_s("%g ", x.v[i]);
	}

	//printf_s("******\n");
	//print_matrix(A);
	//printf_s("******\n");
	//print_vector(b);
	//printf_s("******\n");

	//for (int i = 0; i < 10; i++)
	//{
	//	la_vector<float> xx;
	//	xx = x;
	//	float t = float(i) / 10 * 2 - 1;

	//	for (int j = 0; j < 3; j++)
	//		xx.v[j] += Vt.m[j + 2 * 3]*t*0.1f;

	//	la_vector<float> Ax, Axb;
	//	smvmul(Ax, A, xx);
	//	svsub(Axb, Ax, b);
	//	res = svnorm2(Axb);

	//	printf_s("%f %f %f : %f\n", xx.v[0], xx.v[1], xx.v[2], res);
	//}
}
