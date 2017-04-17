//void test_Lambertian_v1()
//{
//	vector2f uv;
//	uv.x = 0.00195312f;
//	uv.y = -0.00195312f;
//
//	//DEBUG
//	float pixel0;
//	{
//		char filename[MAX_PATH];
//		sprintf_s(filename, PATH "img00-%02d.raw", 0);
//
//		int width, height;
//		std::vector<float> img;
//		load_img_float(filename, img, width, height);
//
//		vector2f scr;
//		uv2screen(scr, uv, width, height);
//
//		vector3f c;
//		lerp_img(c.v, scr, img, width, height, 3);
//		pixel0 = (c.x + c.y + c.z) / 3;
//		printf_s("pixel0: %g\n", pixel0);
//	}
//
//	float avg = 0;
//
//	for (int i = 1; i <= NUM_IMG_PAIRS; i++)
//	{
//		char filename[MAX_PATH];
//		sprintf_s(filename, PATH "img00-%02d-raw.raw", i);
//
//		int width, height;
//		std::vector<float> img;
//		load_img_float(filename, img, width, height);
//
//		float f;
//		matrix4x4f R_tilde, T_tilde;
//		vector3f tau_tilde;
//
//		FILE *fp;
//		sprintf_s(filename, PATH "img00-%02d-raw.param", i);
//		FOPEN(fp, filename, "rb");
//		fread(&f, sizeof(float), 1, fp);
//		fread(R_tilde.m, sizeof(float) * 16, 1, fp);
//		fread(tau_tilde.v, sizeof(float) * 3, 1, fp);
//		fclose(fp);
//		T_tilde = R_tilde;
//		T_tilde._14 = tau_tilde.x;
//		T_tilde._24 = tau_tilde.y;
//		T_tilde._34 = tau_tilde.z;
//
//		//init
//		vector3f vp(0, 0, f);
//		matrix4x4f R;
//		vector3f tau_pp;
//		{
//			matrix4x4f T = T_tilde.inversed_matrix();
//			R = T;
//			R._14 = R._24 = R._34 = 0;
//			vector3f tau(T._14, T._24, T._34);
//			tau_pp = vp - R.transposed_matrix()*vp + R.transposed_matrix()*tau;
//		}
//
//		matrix4x4f T;
//		identity(T);
//		T._14 = tau_pp.x;
//		T._24 = tau_pp.y;
//		T._34 = tau_pp.z;
//
//		float beta = 1.0f / f;
//		float k1 = 1 - beta*T._34;
//		float k2 = T._14 + beta*uv.x*T._34;
//		float k3 = T._24 + beta*uv.y*T._34;
//
//		//compute Ia, Ib, Ic
//		vector3f gI;
//
//		la_matrix<float> A(nsamples, 3);
//		la_vector<float> b(nsamples);
//		for (int j = 0; j < nsamples; j++)
//		{
//			float z = z0 + (z1 - z0)*j / nsamples;
//			vector2f duv;
//			duv.x = k2 / (k1 - beta*z);
//			duv.y = k3 / (k1 - beta*z);
//
//			A.m[0 * A.row + j] = duv.x;
//			A.m[1 * A.row + j] = duv.y;
//			A.m[2 * A.row + j] = 1.0f;
//
//			vector2f uv2 = uv + duv;
//			vector3f v;
//			warp_a_pixel(v, uv2, img, width, height, f, R);
//			b.v[j] = (v.x + v.y + v.z) / 3;
//			//printf_s("%g %g\n", z, b.v[j]);
//		}
//		//printf_s("\n");
//
//		la_vector<float> x;
//		slsq(x, A, b, 1);
//		gI.x = x.v[0];
//		gI.y = x.v[1];
//		gI.z = x.v[2];
//
//		//{
//		//	float z = 1;
//		//	vector3f duv;
//		//	duv.x = k2 / (k1 - beta*z);
//		//	duv.y = k3 / (k1 - beta*z);
//		//	duv.z = 1;
//		//	printf_s("recon: %g\n", duv*gI);
//		//}
//
//		//compute
//		//printf_s("Ia Ib Ic = [%g %g %g]\n", gI.x, gI.y, gI.z);
//
//		{
//			float p = (gI.z - pixel0) / (-k2*gI.x - k3*gI.y);
//			printf_s("z = %g\n", (k1 - 1 / p) / beta);
//
//			avg += p;
//		}
//	}
//
//	//avg /= NUM_IMG_PAIRS;
//	//printf_s("%g\n", FOCAL_LEN*(1 - 1 / avg));
//}


//{
		//	for (int j = 0; j < nsamples; j++)
		//	{
		//		float z = z0 + (z1 - z0)*j / nsamples;
		//		vector2f duv;
		//		duv.x = k2 / (k1 - beta*z);
		//		duv.y = k3 / (k1 - beta*z);

		//		//printf_s("%g\n", duv.x*gI.x+duv.y*gI.y+gI.z);

		//		vector2f scr;
		//		uv2screen(scr, uv + duv, width, height);
		//		vector3f v;
		//		lerp_img(v.v, scr, img, width, height, 3);

		//		printf_s("%g %g\n", z, (v.x + v.y + v.z) / 3);
		//	}
		//	//printf_s("\n");
		//}

//{
		//	matrix4x4f R = trueT;
		//	R._14 = R._24 = R._34 = 0;

		//	std::vector<float> result;
		//	result.resize(width*height * 3, 0);
		//	for (int y = 0; y < height; y++)
		//		for (int x = 0; x < width; x++)
		//		{
		//			vector2f uv;
		//			screen2uv(uv, vector2f(x, y), width, height);

		//			vector3f v;
		//			warp_a_pixel(v, uv, img, width, height, f, R);
		//			result[(x + y*width) * 3 + 0] = v.x;
		//			result[(x + y*width) * 3 + 1] = v.y;
		//			result[(x + y*width) * 3 + 2] = v.z;
		//		}

		//	sprintf_s(filename, PATH "img00-%02d-warp.raw", i);
		//	save_img_float(filename, &result[0], width, height);
		//	return;
		//}

////DEBUG
//		{
//			vector3f tau_pp_tilde = -tau_pp;
//			img.clear();
//			render(img, width, height,
//				vp + tau_pp_tilde,
//				vector3f(1, 0, 0), vector3f(0, 1, 0),
//				f, m, device);
//
//			sprintf_s(filename, PATH "img00-%02d-verify.raw", i);
//			save_img_float(filename, &img[0], width, height);
//		}


////DEBUG
//T_tilde = R_tilde;
//T_tilde._14 = tau_tilde.x;
//T_tilde._24 = tau_tilde.y;
//T_tilde._34 = tau_tilde.z;
//T = T_tilde.inversed_matrix();
//
//matrix4x4f R;
//R = T;
//R._14 = R._24 = R._34 = 0;
//vector3f tau(T._14, T._24, T._34);
//
//vector3f tau_pp;
//tau_pp = vp - R.transposed_matrix()*vp + R.transposed_matrix()*tau;
//
//{
//	vector3f p(0, 0, 1), tp;
//	tp = T*p;
//
//	float beta = 1.0f / f;
//	vector2f uv;
//	uv.x = tp.x / (1 - beta*tp.z);
//	uv.y = tp.y / (1 - beta*tp.z);
//
//	//{
//	//	vector3f dir = tp - vector3f(0, 0, f);
//	//	dir.normalize();
//	//	printf_s("dir1 %g %g %g\n", dir.x, dir.y, dir.z);
//	//}
//
//	vector2f scr;
//	uv2screen(scr, uv, width, height);
//	printf_s("scr (raw): %g %g\n", scr.x, scr.y);
//}
//
//{
//	matrix4x4f newT;
//	identity(newT);
//	newT._14 = tau_pp.x;
//	newT._24 = tau_pp.y;
//	newT._34 = tau_pp.z;
//
//	vector3f p(0, 0, 1), tp;
//	tp = newT*p;
//
//	float beta = 1.0f / f;
//	vector2f uv;
//	uv.x = tp.x / (1 - beta*tp.z);
//	uv.y = tp.y / (1 - beta*tp.z);
//
//	//{
//	//	vector3f dir = tp - vector3f(0, 0, f);
//	//	dir.normalize();
//	//	printf_s("dir2 %g %g %g\n", dir.x, dir.y, dir.z);
//	//}
//
//	vector2f scr;
//	uv2screen(scr, uv, width, height);
//	printf_s("scr (verify): %g %g uv(%gf, %gf)\n", scr.x, scr.y, uv.x, uv.y);
//}
//
//{
//	matrix4x4f newT;
//	identity(newT);
//	newT._14 = tau_pp.x;
//	newT._24 = tau_pp.y;
//	newT._34 = tau_pp.z;
//
//	vector3f x(0, 0, 1);
//	vector3f xpmp, xppmp;
//	xpmp = T*x - vp;
//	xppmp = newT*x - vp;
//	printf_s("x'-p     [%g %g %g]\n", xpmp.x, xpmp.y, xpmp.z);
//	printf_s("x''-p    [%g %g %g]\n", xppmp.x, xppmp.y, xppmp.z);
//	printf_s("R(x''-p) [%g %g %g]\n", (R*xppmp).x, (R*xppmp).y, (R*xppmp).z);
//}
//
////{
////	vector2f uv(-0.419542f, -0.744388f);
//
////	vector3f ray(uv.x, uv.y, -f), oldray;
////	oldray = R*ray;
//
////	vector2f uv2, scr;
////	uv2.x = -f / oldray.z * oldray.x;
////	uv2.y = -f / oldray.z * oldray.y;
////	uv2screen(scr, uv2, width, height);
////	printf_s("scr (warp): %g %g\n", scr.x, scr.y);
////}

//void single_diff_analysis()
//{
//	float f;
//	matrix4x4f R_tilde;
//	vector3f tau_tilde;
//
//	FILE *fp;
//	FOPEN(fp, PATH "img00-03.param", "rb");
//	fread(&f, sizeof(float), 1, fp);
//	fread(R_tilde.m, sizeof(float) * 16, 1, fp);
//	fread(tau_tilde.v, sizeof(float) * 3, 1, fp);
//	fclose(fp);
//
//	float beta = 1.0f / f;
//	vector3f p(0, 0, f);
//
//	int width, height;
//	std::vector<float> img0, img1;
//	load_img_float(PATH "img00-00.raw", img0, width, height);
//	load_img_float(PATH "img00-03.raw", img1, width, height);
//
//	pref_env env;
//	env.load(PATH "uffizi.pref_env");
//
//	matrix4x4f R;
//	vector3f tau;
//	R = R_tilde.transposed_matrix();
//	tau = -(R*tau_tilde);
//
//	//vector3f	x(0.573567f, -0.000929606f, -9.51916f + f), 
//	//			n(0.766875f, -0.00136332f, 0.641796f);
//	vector3f	x(0.00180695f, -0.00180695f, -9.2516f + f),
//		n(0.0026028f, -0.00210045f, 0.999994f);
//	n.normalize();
//
//	vector3f x2;
//	x2 = R*x + tau;
//
//	vector3f delta_v, v;
//	delta_v = R_tilde*p - p + tau_tilde;
//	v = p - x;
//
//	{
//		printf_s("v : %g %g %g\n", v.x, v.y, v.z);
//		vector3f v2;
//		v2 = R_tilde*p + tau_tilde - x;
//		printf_s("v': %g %g %g\n", v2.x, v2.y, v2.z);
//	}
//
//	vector2f uv, uv2, delta_uv;
//	p2uv(uv, x, beta);
//	printf_s("uv: (%g %g)\n", uv.x, uv.y);
//	p2uv(uv2, x2, beta);
//	delta_uv = uv2 - uv;
//
//	float m1 = beta*uv.x;
//	float m2 = beta*uv.y;
//
//	const float delta = 5e-3f;
//	vector3f nabla_v_rho;
//	{
//		float dlen = delta*f * 2;
//		vector3f rho_p, rho_m;
//
//		compute_rho(rho_p, (v + vector3f(dlen, 0, 0)).normalized_vector(), n, env);
//		compute_rho(rho_m, (v - vector3f(dlen, 0, 0)).normalized_vector(), n, env);
//		nabla_v_rho.x = (rho_p.y - rho_m.y) / (dlen * 2);
//
//		compute_rho(rho_p, (v + vector3f(0, dlen, 0)).normalized_vector(), n, env);
//		compute_rho(rho_m, (v - vector3f(0, dlen, 0)).normalized_vector(), n, env);
//		nabla_v_rho.y = (rho_p.y - rho_m.y) / (dlen * 2);
//
//		compute_rho(rho_p, (v + vector3f(0, 0, dlen)).normalized_vector(), n, env);
//		compute_rho(rho_m, (v - vector3f(0, 0, dlen)).normalized_vector(), n, env);
//		nabla_v_rho.z = (rho_p.y - rho_m.y) / (dlen * 2);
//	}
//
//	printf_s("delta_v    : %g %g %g\n", delta_v.x, delta_v.y, delta_v.z);
//	printf_s("nabla_v_rho: %g %g %g\n", nabla_v_rho.x, nabla_v_rho.y, nabla_v_rho.z);
//	printf_s("nabla_v_rho*v = %g\n", nabla_v_rho*v);
//	nabla_v_rho.z = m1*nabla_v_rho.x + m2*nabla_v_rho.y;
//	printf_s("nabla_v_rho*v = %g [our formula]\n", nabla_v_rho*v);
//
//	printf_s("m1*nabla_v_rho_x + m2*nabla_v_rho_y: %g\n", m1*nabla_v_rho.x + m2*nabla_v_rho.y);
//
//	//std::vector<std::vector<float>> img_grad;
//	//compute_img_grad(img_grad, img1, width, height);
//
//	vector2f scr;
//	uv2screen(scr, uv, width, height);
//
//	vector2f grad_I;
//	//lerp_img(grad_I.v, scr, img_grad[1], width, height, 2);
//
//	float I_u, I_v;
//	//I_u = grad_I.x * (width / 2);
//	//I_v = -grad_I.y * (height / 2);
//	//printf_s("I_u: (%g %g)\n", I_u, I_v);
//
//	{
//		vector3f I0, I1;
//		lerp_img(I0.v, scr + vector2f(-1, 0), img1, width, height, 3);
//		lerp_img(I1.v, scr + vector2f(1, 0), img1, width, height, 3);
//		//printf_s("x-/x+ : %g %g (%g)\n", I0.y, I1.y, I1.y - I0.y);
//		grad_I.x = (I1.y - I0.y) / 2;
//
//		lerp_img(I0.v, scr + vector2f(0, -1), img1, width, height, 3);
//		lerp_img(I1.v, scr + vector2f(0, 1), img1, width, height, 3);
//		//printf_s("y-/y+ : %g %g (%g)\n", I0.y, I1.y, I1.y - I0.y);
//		grad_I.y = (I1.y - I0.y) / 2;
//	}
//	I_u = grad_I.x * (width / 2);
//	I_v = -grad_I.y * (height / 2);
//	//printf_s("I_u: (%g %g)\n", I_u, I_v);
//
//	float delta_I;
//	vector3f I0, I1;
//	lerp_img(I0.v, scr, img0, width, height, 3);
//	lerp_img(I1.v, scr, img1, width, height, 3);
//	delta_I = I1.y - I0.y;
//
//	{
//		printf_s("\n");
//
//		vector2f scrscr;
//		uv2screen(scrscr, uv2, width, height);
//
//		printf_s("delta_uv: (%g %g)\n", delta_uv.x, delta_uv.y);
//		printf_s("delta_scr: (%g %g)\n", (scrscr - scr).x, (scrscr - scr).y);
//
//		vector3f rho0, rho1;
//		compute_rho(rho0, v.normalized_vector(), n, env);
//		compute_rho(rho1, (v + delta_v).normalized_vector(), n, env);
//
//		vector3f I0, I1;
//		lerp_img(I0.v, scr, img0, width, height, 3);
//		lerp_img(I1.v, scrscr, img1, width, height, 3);
//		printf_s("scr: (%g %g) (%g %g)\n", scr.x, scr.y, scrscr.x, scrscr.y);
//
//		//{
//		//	printf_s("rho': %g\n", rho1.y);
//		//	vector3f n(0.00249023, -0.00101858, 0.999996);
//		//	vector3f v(-0.00615035, -0.00206208, 0.999979);
//		//	vector3f rho;
//		//	compute_rho(rho, v.normalized_vector(), n, env);
//		//	printf_s("imaged rho': %g\n", rho.y);
//		//	printf_s("lerped rho': %g\n", I1.y);
//		//	printf_s("lerped rho: %g\n", I0.y);
//		//}
//
//		printf_s("%g = %g [dI : diff-rho]\n",
//			I1.y - I0.y,
//			(rho1 - rho0).y);
//
//		//printf_s("%g (dI) = %g (drho) = %g\n",
//		//	I1.y-I0.y,
//		//	(rho1 - rho0).y,
//		//	nabla_v_rho.x*delta_v.x + nabla_v_rho.y*delta_v.y + (m1*nabla_v_rho.x + m2*nabla_v_rho.y)*delta_v.z);
//
//		//DEBUG
//		//{
//		//	vector3f rho_p, rho_m, vp, vm, rho1, rho0;
//		//	vp = (v + delta_v.normalized_vector()*delta*f*2).normalized_vector();
//		//	compute_rho(rho_p, vp, n, env);
//		//	vm = (v - delta_v.normalized_vector()*delta*f*2).normalized_vector();
//		//	compute_rho(rho_m, vm, n, env);
//
//		//	for (int i = 0; i < 20; i++)
//		//	{
//		//		float len = delta_v.length()*(i - 19.5f) / 10;
//		//		vector3f d = delta_v.normalized_vector()*len;
//		//		compute_rho(rho0, v.normalized_vector(), n, env);
//		//		compute_rho(rho1, (v + d).normalized_vector(), n, env);
//
//		//		printf("%g\t%g\t%g\n", (rho1 - rho0).y, (rho_p.y - rho_m.y) / (delta*f*4)*len, d.length());
//		//	}
//		//}
//
//		lerp_img(I0.v, scr, img1, width, height, 3);
//		lerp_img(I1.v, scrscr, img1, width, height, 3);
//
//		printf_s("%g = %g [dI : diff-approx]\n",
//			I1.y - I0.y,
//			I_u*delta_uv.x + I_v*delta_uv.y);
//
//		////DEBUG
//		//for (int i = 0; i < 20; i++)
//		//{
//		//	float len = delta_uv.length()*(i - 9.5f) / 10;
//		//	vector2f d = delta_uv.normalized_vector()*len;
//
//		//	vector2f uvuv = uv+d, s;
//		//	uv2screen(s, uvuv, width, height);
//		//	lerp_img(I1.v, s, img1, width, height, 3);
//
//		//	printf("%g\t%g\t%g\n", I1.y - I0.y, I_u*d.x + I_v*d.y, I1.y);
//		//}
//
//		printf_s("ground-truth: %g = %g+%g = [%g]\n", (rho1 - rho0).y, I1.y - I0.y, delta_I, I1.y - I0.y + delta_I);
//		printf_s("ours: %g = %g+%g = [%g]\n",
//			nabla_v_rho.x*delta_v.x + nabla_v_rho.y*delta_v.y + (m1*nabla_v_rho.x + m2*nabla_v_rho.y)*delta_v.z,
//			I_u*delta_uv.x + I_v*delta_uv.y, delta_I,
//			I_u*delta_uv.x + I_v*delta_uv.y + delta_I);
//
//		printf_s("\n");
//		//printf_s("dx = (%g %g %g)\n", (x2 - x).x, (x2 - x).y, (x2 - x).z);
//
//		matrix4x4f E;
//		identity(E);
//		matrix4x4f RmE = R - E;
//
//		float m1 = beta*uv.x;
//		float m2 = beta*uv.y;
//		float l1 = RmE._11*uv.x + RmE._12*uv.y + tau.x;
//		float l3 = RmE._21*uv.x + RmE._22*uv.y + tau.y;
//		float l5 = RmE._31*uv.x + RmE._32*uv.y + tau.z;
//		float l2 = RmE._13 - RmE._11*m1 - RmE._12*m2;
//		float l4 = RmE._23 - RmE._21*m1 - RmE._22*m2;
//		float l6 = RmE._33 - RmE._31*m1 - RmE._32*m2;
//		//printf_s("dx = (%g %g %g)\n", l1+l2*x.z, l3 + l4*x.z, l5 + l6*x.z);
//
//		float alpha1 = 1 - beta*l5;
//		float alpha2 = -beta*(1 + l6);
//		float k1 = (l2 + m1*l6) / alpha2;
//		float k3 = (l4 + m2*l6) / alpha2;
//		float k2 = l1 + m1*l5 - alpha1 / alpha2*(l2 + m1*l6);
//		float k4 = l3 + m2*l5 - alpha1 / alpha2*(l4 + m2*l6);
//		//printf_s("delta_uv = (%g %g)\n", k1 + k2 / (alpha1 + alpha2*x.z), k3 + k4 / (alpha1 + alpha2*x.z));
//	}
//	//printf_s("\n%g = %g + %g?\n",
//	//nabla_v_rho.x*delta_v.x + nabla_v_rho.y*delta_v.y + (m1*nabla_v_rho.x + m2*nabla_v_rho.y)*delta_v.z,
//	//I_u*delta_uv.x + I_v*delta_uv.y, delta_I);
//}


//void warp_cam()//const vector3f &center)
//{
//	for (int i = 1; i <= NUM_IMG_PAIRS; i++)
//	{
//		char filename[MAX_PATH];
//		sprintf_s(filename, PATH "img00-%02d-raw.raw", i);
//
//		int width, height;
//		std::vector<float> img0;
//		load_img_float(filename, img0, width, height);
//
//		float f;
//		matrix4x4f R_tilde;
//		vector3f tau_tilde;
//
//		FILE *fp;
//		sprintf_s(filename, PATH "img00-%02d-raw.param", i);
//		FOPEN(fp, filename, "rb");
//		fread(&f, sizeof(float), 1, fp);
//		fread(R_tilde.m, sizeof(float) * 16, 1, fp);
//		fread(tau_tilde.v, sizeof(float) * 3, 1, fp);
//		fclose(fp);
//
//		//vector3f p(0, 0, f);
//		//vector3f pp, v;
//		//vector3f dx(1, 0, 0), dy(0, 1, 0), dz(0, 0, 1);
//		//v = (center - p).normalized_vector();
//		//vector3f ray0(v*dx, v*dy, v*dz);
//
//		//vector3f dx2, dy2, dz2;
//		//dx2 = R_tilde*dx;
//		//dy2 = R_tilde*dy;
//		//dz2 = R_tilde*dz;
//
//		//pp = R_tilde*p + tau_tilde;
//		//v = (center - pp).normalized_vector();
//		//vector3f ray1(v*dx2, v*dy2, v*dz2);
//
//		//vector3f axis = (ray1^ray0).normalized_vector();
//		//matrix4x4f rot;
//		//rotation_axis(rot, axis, acos(ray1*ray0));
//		//rot.transpose();
//
//		//matrix4x4f new_R_tilde = rot*R_tilde;
//
//		matrix4x4f new_R_tilde;
//		identity(new_R_tilde);
//		vector3f new_tau_tilde;
//		std::vector<float> img;
//		warp_view(img, new_R_tilde, new_tau_tilde, img0, width, height, f, R_tilde, tau_tilde);
//
//		sprintf_s(filename, PATH "img00-%02d.raw", i);
//		save_img_float(filename, &img[0], width, height);
//
//		sprintf_s(filename, PATH "img00-%02d.param", i);
//		FOPEN(fp, filename, "wb");
//		fwrite(&f, sizeof(float), 1, fp);
//		fwrite(new_R_tilde.m, sizeof(float) * 16, 1, fp);
//		fwrite(new_tau_tilde.v, sizeof(float) * 3, 1, fp);
//		fclose(fp);
//	}
//}

//void compute_img_grad(std::vector<std::vector<float>> &grad, const std::vector<float> &img, const int width, const int height)
//{
//	grad.resize(3);
//	for (int i = 0; i < grad.size(); i++)
//		grad[i].resize(width*height*2);
//	for (int y = 0; y < height; y++)
//		for (int x = 0; x < width; x++)
//		{
//			float step = 0;
//			
//			int idx;
//			vector3f m, p, gx, gy;
//
//			if (x - 1 >= 0)
//			{
//				idx = (x-1 + y * width) * 3;
//				step += 1;
//			} else {
//				idx = (x + y * width) * 3;
//			}
//			m = vector3f(img[idx + 0], img[idx + 1], img[idx + 2]);
//
//			if (x + 1 <= width-1)
//			{
//				idx = (x + 1 + y * width) * 3;
//				step += 1;
//			} else {
//				idx = (x + y * width) * 3;
//			}
//			p = vector3f(img[idx + 0], img[idx + 1], img[idx + 2]);
//			gx = (p - m) / step;
//
//			if (y - 1 >= 0)
//			{
//				idx = (x + (y-1) * width) * 3;
//				step += 1;
//			} else {
//				idx = (x + y * width) * 3;
//			}
//			m = vector3f(img[idx + 0], img[idx + 1], img[idx + 2]);
//
//			if (y + 1 <= height - 1)
//			{
//				idx = (x + (y+1) * width) * 3;
//				step += 1;
//			} else {
//				idx = (x + y * width) * 3;
//			}
//			p = vector3f(img[idx + 0], img[idx + 1], img[idx + 2]);
//			gy = (p - m) / step;
//
//			grad[0][(x + y*width) * 2 + 0] = gx.x;
//			grad[0][(x + y*width) * 2 + 1] = gy.x;
//			grad[1][(x + y*width) * 2 + 0] = gx.y;
//			grad[1][(x + y*width) * 2 + 1] = gy.y;
//			grad[2][(x + y*width) * 2 + 0] = gx.z;
//			grad[2][(x + y*width) * 2 + 1] = gy.z;
//		}
//}