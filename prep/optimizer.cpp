#include "optimizer.h"

template <typename T>
void swap(T &x, T &y)
{
	T z;
	z = x; x = y; y = z;
}

void save_mydata_file(const char *filename, const mydata &dat)
{
	FILE *file;
	if ((file = fopen(filename, "wb")) == NULL)
	{
		printf("open file failed!\n");
		exit(1);
	}

	int n_of_points = dat.point_set.points.size();
	fwrite(&n_of_points, sizeof(int), 1, file);
	//fwrite(&dat.kappa, sizeof(float), 1, file);

	for (int i = 0; i < NUM_IMG_PAIRS; i++)
	{
		fwrite(&dat.k1[i], sizeof(float), 1, file);
		fwrite(dat.Delta_v[i].v, sizeof(float), 3, file);
		//fwrite(&dat.Delta_v[i].x, sizeof(float), 1, file);
		//fwrite(&dat.Delta_v[i].y, sizeof(float), 1, file);
		//fwrite(&dat.Delta_v[i].z, sizeof(float), 1, file);
	}
	for (int i = 0; i < n_of_points; i++)
	{
		fwrite(dat.Ic_Iut[i].v, sizeof(float), NUM_IMG_PAIRS, file);
		fwrite(dat.k2Ia_k3Ib[i].v, sizeof(float), NUM_IMG_PAIRS, file);
	}
	fclose(file);
}

void read_mydata_file(const char *filename, mydata &dat)
{
	FILE *file;
	if ((file = fopen(filename, "rb")) == NULL)
	{
		printf("open file failed!\n");
		exit(1);
	}

	int n_of_points;
	fread(&n_of_points, sizeof(int), 1, file);
	//fread(&dat.kappa, sizeof(float), 1, file);

	for (int i = 0; i < NUM_IMG_PAIRS; i++)
	{
		//float tmp[4];
		//fread(tmp, sizeof(float), 4, file);
		//dat.k1.push_back(tmp[0]);
		//dat.Delta_v.push_back(vector3f(tmp[1], tmp[2], tmp[3]));
		float k1;
		fread(&k1, sizeof(float), 1, file);
		dat.k1.push_back(k1);
		vector3f Dv;
		fread(Dv.v, sizeof(float), 3, file);
		dat.Delta_v.push_back(Dv);
	}
	for (int i = 0; i < n_of_points; i++)
	{
		fread(dat.Ic_Iut[i].v, sizeof(float), NUM_IMG_PAIRS, file);
		fread(dat.k2Ia_k3Ib[i].v, sizeof(float), NUM_IMG_PAIRS, file);
	}
	fclose(file);
}

void get_pi_map(const vector3f &ref_n, const vector3f &wo, const float &kappa, const pref_env &env)
{
	vector3f n0;
	hemi_octa_frame<float> pi_map_by_n;
	int dim = 8192;
	pi_map_by_n.init(dim);
	printf_s("roughness=%g\n", kappa);
	n0 = ref_n;
	n0.normalize();
	int idx = pi_map_by_n.map(n0);
	int row = idx / dim;
	int col = idx % dim;
	printf_s("mapped index is (%d, %d)\n", row, col);
	float *pi_val = new float[dim * 3];
	int cnt = 0;
	//vector3f pre_n;
	for (int i = row - 1; i <= row + 1; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			int cur = i * dim + j;
			vector3f n;
			la_vector<float> pi;
			pi_map_by_n.get_n(n, cur);
			//if (cnt > 0) {
			//	printf_s("%.6f\t%.6f\t%.6f\n", n.x - pre_n.x, n.y - pre_n.y, n.z - pre_n.z);
			//}
			//pre_n = n;
			calc_pi_hat(pi, wo, n, kappa, env);
			pi_val[cnt++] = pi.v[0];
			if (j == col && i == row)
			{
				printf_s("mapped pi_a is %g\n", pi.v[0]);
			}
		}
	}
	char filename[MAX_PATH];
	sprintf_s(filename, "d:/codes/diff/data/pi_map_by_n.dat");
	FILE *file;
	if ((file = fopen(filename, "wb")) == NULL)
	{
		printf("open file failed!\n");
		exit(1);
	}
	fwrite(pi_val, sizeof(float), dim * 3, file);
	fclose(file);

	//{
	//	float delta = 1e-2;
	//	vector3f rho;
	//	int range = 100;
	//	float *rho_val = new float[range * 2 + 1];
	//	int cnt = 0;
	//	for (int i = -range; i <= range; i++)
	//	{
	//		vector3f wo_d = wo + vector3f(i*delta, 0, 0);
	//		wo_d.normalize();
	//		compute_rho(rho, wo_d, ref_n, kappa, env);
	//		rho_val[cnt++] = (rho.x + rho.y + rho.z) / 3;
	//		printf_s("%g\n", (rho.x + rho.y + rho.z) / 3);
	//	}		
	//}
	exit(1);
}


void calc_pi_ls_err(double *p, int p_id, void *adata)
{
	struct mydata *dptr = (struct mydata *)adata;
	float beta = 1 / FOCAL_LEN;
	// update world coordinates of all points
	dptr->point_set.update_xyz(p);

	la_matrix<float> A(NUM_IMG_PAIRS, 2);
	float betau = beta * dptr->point_set.points[p_id].uv.x;
	float betav = beta * dptr->point_set.points[p_id].uv.y;
	for (int j = 0; j < NUM_IMG_PAIRS; j++)
	{
		A.m[0 * A.row + j] = dptr->Delta_v[j].x + betau * dptr->Delta_v[j].z;
		A.m[1 * A.row + j] = dptr->Delta_v[j].y + betav * dptr->Delta_v[j].z;
	}

	const la_vector<float> b1 = dptr->Ic_Iut[p_id];
	const la_vector<float> b2 = dptr->k2Ia_k3Ib[p_id];
	const la_vector<float> k1 = dptr->k1;
	float z = p[p_id];

	assert(A.row == b1.length);
	int dim = b1.length;
	la_vector<float> b0(dim);
	for (int i = 0; i < dim; i++)
	{
		b0.v[i] = b1.v[i] + b2.v[i] / (k1.v[i] - beta*z);
	}

	vector3f point_normal;
	calc_normal(point_normal, dptr->point_set, p_id);
	la_vector<float> pi;
	vector3f view_vec =
		vector3f(0, 0, FOCAL_LEN) - dptr->point_set.points[p_id].world_coord;
	calc_pi_hat(pi, view_vec, point_normal, dptr->kappa[p_id], dptr->env);

	la_vector<float> temp1, temp2;
	smvmul(temp1, A, pi);
	svsub(temp2, temp1, b0);
	printf_s("%g\t%g\t%g\n", svnorm2(b0), svnorm2(temp1), svnorm2(temp2));
	for (int i = 0; i < b0.length; i++)
	{
		//printf_s("%.6f\t%.6f\t\t%.6f\t%.6f%%\n",
		printf_s("%g\t%g\t\t%g\t%g%%\n",
			b0.v[i], temp1.v[i], temp2.v[i], fabs(temp2.v[i] / b0.v[i]) * 100);
	}
	//exit(1);
}

void ffunc(double *p, double *hx, int nvars, int nobs, void *adata)
/* functional relation describing measurements. Given a parameter vector p,
* computes a prediction of the measurements \hat{x}. p is nvarsx1,
* \hat{x} is nobsx1, maximum
*/
{
	static int cnt = 0;
	//printf_s("counter for accessing func: %d\n", ++cnt);
	struct mydata *dptr = (struct mydata *)adata;
	float beta = 1 / FOCAL_LEN;
	// update world coordinates of all points
	dptr->point_set.update_xyz(p);
	// compute normalized pi by least squares
	for (int i = 0; i < nvars; i++)
	{
		la_vector<float> pi;
		la_matrix<float> A(NUM_IMG_PAIRS, 2);
		float betau = beta * dptr->point_set.points[i].uv.x;
		float betav = beta * dptr->point_set.points[i].uv.y;
		for (int j = 0; j < NUM_IMG_PAIRS; j++)
		{
			A.m[0 * A.row + j] = dptr->Delta_v[j].x + betau * dptr->Delta_v[j].z;
			A.m[1 * A.row + j] = dptr->Delta_v[j].y + betav * dptr->Delta_v[j].z;
		}		
		calc_pi_ls(pi, A, dptr->Ic_Iut[i], dptr->k2Ia_k3Ib[i],
			dptr->k1, beta, p[i]);
		//calc_pi_ls_err(pi, A, dptr->Ic_Iut[i], dptr->k2Ia_k3Ib[i],
		//	dptr->k1, beta, p[i]);

		la_vector<float> pi_norm;
		calc_vecnorm(pi_norm, pi);
		//pi_norm = pi;
		hx[i * 2] = pi_norm.v[0];
		hx[i * 2 + 1] = pi_norm.v[1];
	}
	//{
		char filename[MAX_PATH];
		sprintf_s(filename, "d:/codes/diff/data/groundtruth_normal.dat");
		FILE *file_normal;
		if ((file_normal = fopen(filename, "rb")) == NULL)
		{
			printf("open file failed!\n");
			exit(1);
		}
		//sprintf_s(filename, "d:/codes/diff/data/dot_of_normals.txt");
		//FILE *file_dot;
		//if ((file_dot = fopen(filename, "w")) == NULL)
		//{
		//	printf("open file failed!\n");
		//	exit(1);
		//}
	//}
	// compute normalized pi by lookup table
	for (int i = 0; i < nvars; i++)
	{
		vector3f point_normal;
		calc_normal(point_normal, dptr->point_set, i);
		//{	// check how well we approximate normals
			//vector3f true_normal;
			//fread(true_normal.v, sizeof(float), 3, file_normal);
			//point_normal = true_normal;
			//vector3f temp;
			//temp = point_normal;
			//temp.normalize();
			//float dot_pro = temp * true_normal;
			//fprintf_s(file_dot, "%g\n", dot_pro);
			//printf_s("%d:\t%g\n", i, dot_pro);
		//}
		//point_normal.normalize();
		//std::cout << point_normal << std::endl;
		la_vector<float> pi;
		// strictly, view_vec is affected by z but we ignore such minor difference
		vector3f view_vec = 
			vector3f(0, 0, FOCAL_LEN) - dptr->point_set.points[i].world_coord;
		//view_vec.normalize();	

		//int x = dptr->point_set.points[i].scr_u;
		//int y = dptr->point_set.points[i].scr_v;
		//if (x == IMG_DIM / 2 - 24 && y == IMG_DIM / 2 + 24)
		//{
		//	vector3f c;
		//	compute_rho(c, view_vec, point_normal, dptr->env);
		//	printf_s("n(%gf, %gf, %gf)\n",
		//		point_normal.x, point_normal.y, point_normal.z);
		//	printf_s("v(%gf, %gf, %gf)\n", view_vec.x, view_vec.y, view_vec.z);
		//	printf_s("c(%gf, %gf, %gf)\n", c.x, c.y, c.z);
		//	//exit(1);
		//}
		//get_pi_map(point_normal, view_vec, dptr->kappa[i], dptr->env);
		calc_pi_hat(pi, view_vec, point_normal, dptr->kappa[i], dptr->env);

		{	// check the error of substituting differential pi into ls formula

		}

		la_vector<float> pi_norm;
		calc_vecnorm(pi_norm, pi);
		//pi_norm = pi;
		//printf_s("point %d:(%g,%g) - (%g,%g) flag:%d\n", 
		//	i, hx[i * 2], hx[i*2+1], pi_norm.v[0], pi_norm.v[1], dptr->point_set.points[i].inv_cross);
		hx[i * 2] = pi_norm.v[0] - hx[i * 2];
		hx[i * 2 + 1] = pi_norm.v[1] - hx[i * 2 + 1];
	}
	//fclose(file_dot);
	//exit(1);
}

void func_point(double *p, la_vector<float> &pi_diff, int nvars, int nobs, void *adata, int p_id)
{
	struct mydata *dptr = (struct mydata *)adata;
	float beta = 1 / FOCAL_LEN;
	// update world coordinates of all points
	dptr->point_set.update_xyz(p);
	// compute normalized pi by least squares

	la_vector<float> pi;
	la_matrix<float> A(NUM_IMG_PAIRS, 2);
	float betau = beta * dptr->point_set.points[p_id].uv.x;
	float betav = beta * dptr->point_set.points[p_id].uv.y;
	for (int j = 0; j < NUM_IMG_PAIRS; j++)
	{
		A.m[0 * A.row + j] = dptr->Delta_v[j].x + betau * dptr->Delta_v[j].z;
		A.m[1 * A.row + j] = dptr->Delta_v[j].y + betav * dptr->Delta_v[j].z;
	}
	calc_pi_ls(pi, A, dptr->Ic_Iut[p_id], dptr->k2Ia_k3Ib[p_id],
		dptr->k1, beta, p[p_id]);

	la_vector<float> pi_norm;
	calc_vecnorm(pi_norm, pi);
	pi_diff.v[0] = pi_norm.v[0];
	pi_diff.v[1] = pi_norm.v[1];

	// compute normalized pi by lookup table

	vector3f point_normal;
	calc_normal(point_normal, dptr->point_set, p_id);
	//point_normal.normalize();

	//std::cout << point_normal << std::endl;
	//la_vector<float> pi;
	// strictly, view_vec is affected by z but we ignore such minor difference
	vector3f view_vec =
		vector3f(0, 0, FOCAL_LEN) - dptr->point_set.points[p_id].world_coord;
	//view_vec.normalize();

	//int x = dptr->point_set.points[p_id].scr_u;
	//int y = dptr->point_set.points[p_id].scr_v;
	//if (x == IMG_DIM / 2 - 24 && y == IMG_DIM / 2 + 24)
	//{
	//	vector3f c;
	//	compute_rho(c, view_vec, point_normal, dptr->env);
	//	printf_s("n(%gf, %gf, %gf)\n", 
	//		point_normal.x, point_normal.y, point_normal.z);
	//	printf_s("v(%gf, %gf, %gf)\n", view_vec.x, view_vec.y, view_vec.z);
	//	printf_s("c(%gf, %gf, %gf)\n", c.x, c.y, c.z);
	//	exit(1);
	//}
	//printf_s("n(%gf, %gf, %gf)\n", 
	//	point_normal.x, point_normal.y, point_normal.z);
	//printf_s("v(%gf, %gf, %gf)\n", view_vec.x, view_vec.y, view_vec.z);
	printf_s("(%g, %g) - ", pi.v[0], pi.v[1]);
	calc_pi_hat(pi, view_vec, point_normal, dptr->kappa[p_id], dptr->env);
	printf_s("(%g, %g)\n", pi.v[0], pi.v[1]);
	//la_vector<float> pi_norm;
	calc_vecnorm(pi_norm, pi);
	pi_diff.v[0] -= pi_norm.v[0];
	pi_diff.v[1] -= pi_norm.v[1];
}

void fjac(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata)
/* function to supply the nonzero pattern of the sparse Jacobian of func and
* evaluate it at p in CRS format. Non-zero elements are to be stored in jac
* which has been preallocated with a capacity of Jnnz
*/
{
	static int cnt = 0;
	//printf_s("counter for accessing fjac: %d\n", ++cnt);
	struct mydata *dptr = (struct mydata *)adata;
	float beta = 1 / FOCAL_LEN;	
	//printf_s("Parameters: nvars=%d nobs=%d\n", nvars, nobs);
	//printf_s("Jacobian info: Dim:%dx%d nnz#=%d\n", jac->nr, jac->nc, jac->nnz);
	//exit(1);

	// update world coordinates of all points
	dptr->point_set.update_xyz(p);

	// compute gradient of pi from least squares
	float *pi_ls_grad = new float[nobs];
	for (int i = 0; i < nvars; i++)
	{
		la_vector<float> pi_grad, pi;
		la_matrix<float> A(NUM_IMG_PAIRS, 2);
		float betau = beta * dptr->point_set.points[i].uv.x;
		float betav = beta * dptr->point_set.points[i].uv.y;
		for (int j = 0; j < NUM_IMG_PAIRS; j++)
		{
			A.m[0 * A.row + j] = dptr->Delta_v[j].x + betau * dptr->Delta_v[j].z;
			A.m[1 * A.row + j] = dptr->Delta_v[j].y + betav * dptr->Delta_v[j].z;
		}
		calc_pi_ls(pi, A, dptr->Ic_Iut[i], dptr->k2Ia_k3Ib[i],
			dptr->k1, beta, p[i]);
		grad_pi_ls(pi_grad, A, dptr->Ic_Iut[i], dptr->k2Ia_k3Ib[i],
			dptr->k1, beta, p[i]);
		//printf_s("point %d:pi_a=%.6f, pi_b=%.6f\ngradient:(%g,%g)\n", 
		//	i, pi.v[0], pi.v[1], pi_grad.v[0], pi_grad.v[1]);

		la_matrix<float> pi_norm_grad;
		grad_vec2norm(pi_norm_grad, pi);
		la_vector<float> mul_grad;
		smvmul(mul_grad, pi_norm_grad, pi_grad);
		pi_ls_grad[i * 2] = mul_grad.v[0];
		pi_ls_grad[i * 2 + 1] = mul_grad.v[1];
	}
	// compute normalized pi by lookup table
	for (int i = 0; i < nvars; i++)
	{
		vector3f point_normal;
		calc_normal(point_normal, dptr->point_set, i);
		la_matrix<float> grad_point_normal;
		grad_normal(grad_point_normal, dptr->point_set, i);
		//la_matrix<float> grad_unit_normal;
		//grad_vec3norm(grad_unit_normal, point_normal);
		//point_normal.normalize();
		la_vector<float> pi;
		// strictly, view_vec is affected by z but we ignore such minor difference
		vector3f view_vec =
			vector3f(0, 0, FOCAL_LEN) - dptr->point_set.points[i].world_coord;
		//view_vec.normalize();
		calc_pi_hat(pi, view_vec, point_normal, dptr->kappa[i], dptr->env);
		la_matrix<float> grad_pi;
		grad_pi_hat(grad_pi, view_vec, point_normal, dptr->kappa[i], dptr->env);
		la_matrix<float> grad_pi_norm;
		grad_vec2norm(grad_pi_norm, pi);
		la_matrix<float> mul_temp1, mul_temp2, mul_grad;
		//smmmul(mul_temp1, grad_unit_normal, grad_point_normal);
		//smmmul(mul_temp2, grad_pi, mul_temp1);
		smmmul(mul_temp2, grad_pi, grad_point_normal);
		smmmul(mul_grad, grad_pi_norm, mul_temp2);
		mul_grad.m[0 * mul_grad.row + 0] -= pi_ls_grad[i * 2];
		mul_grad.m[0 * mul_grad.row + 1] -= pi_ls_grad[i * 2 + 1];

		//printf_s("point %d:pi_a=%.6f, pi_b=%.6f\ngradient:\n", i, pi.v[0], pi.v[1]);
		//for (int j = 0; j < 6; j++)
		//{
		//	printf_s("%g ", grad_pi.m[j]);
		//}
		//printf_s("\n");

		// first sort then fill
		int base = i * 6;
		int p1 = dptr->point_set.points[i].neighbor1;
		int p2 = dptr->point_set.points[i].neighbor2;
		if (dptr->point_set.points[i].inv_cross)
		{
			swap<int>(p1, p2);
		}
		int index[3] = {i, p1, p2};
		for (int j = 0; j < 3; j++)
		{
			jac->val[base + j] = mul_grad.m[j * mul_grad.row + 0];
			jac->val[base + j + 3] = mul_grad.m[j * mul_grad.row + 1];
		}
		for (int rep = 0; rep < 2; rep++)
		{
			for (int j = 0; j < 2 - rep; j++)
			{
				if (index[j] > index[j + 1])
				{
					swap<int>(index[j], index[j + 1]);
					swap<double>(jac->val[base + j], jac->val[base + j + 1]);
					swap<double>(jac->val[base + 3 + j], jac->val[base + 3 + j + 1]);
				}
			}
		}
		assert(index[0] < index[1]);
		assert(index[1] < index[2]);
		for (int j = 0; j < 3; j++)
		{
			jac->colidx[base + j] = index[j];
			jac->colidx[base + j + 3] = index[j];
		}
		jac->rowptr[i * 2] = base;
		jac->rowptr[i * 2 + 1] = base + 3;
		//printf_s("Jacobian line %d/%d:\n", i + 1, nvars);
		//for (int j = 0; j < 3; j++)
		//{
		//	printf_s("%.5f\t\t", jac->val[base + j]);
		//}
		//printf_s("\n");
		//for (int j = 0; j < 3; j++)
		//{
		//	printf_s("%.5f\t\t", jac->val[base + j + 3]);
		//}
		//printf_s("\n");
	}
	jac->rowptr[nobs] = nvars * 6;
	//exit(1);
}

void fjac_diff(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata)
{
	struct mydata *dptr = (struct mydata *)adata;
	for (int i = 0; i < nvars; i++)
	{
		// first sort then fill
		int base = i * 6;
		int index[3] = { i, dptr->point_set.points[i].neighbor1 , dptr->point_set.points[i].neighbor2 };
		for (int rep = 0; rep < 2; rep++)
		{
			for (int j = 0; j < 2 - rep; j++)
			{
				if (index[j] > index[j + 1])
				{
					swap<int>(index[j], index[j + 1]);
				}
			}
		}
		assert(index[0] < index[1]);
		assert(index[1] < index[2]);
		for (int j = 0; j < 3; j++)
		{
			jac->colidx[base + j] = index[j];
			jac->colidx[base + j + 3] = index[j];
		}
		jac->rowptr[i * 2] = base;
		jac->rowptr[i * 2 + 1] = base + 3;
	}
	jac->rowptr[nobs] = nvars * 6;
}

void fjac_mannual_diff(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata)
{
	struct mydata *dptr = (struct mydata *)adata;
	for (int i = 0; i < nvars; i++)
	{
		// first sort then fill
		int base = i * 6;
		int index[3] = { i, dptr->point_set.points[i].neighbor1 , dptr->point_set.points[i].neighbor2 };

		for (int rep = 0; rep < 2; rep++)
		{
			for (int j = 0; j < 2 - rep; j++)
			{
				if (index[j] > index[j + 1])
				{
					swap<int>(index[j], index[j + 1]);
				}
			}
		}
		assert(index[0] < index[1]);
		assert(index[1] < index[2]);

		la_vector<float> pi_0(2), pi_alter(2);
		func_point(p, pi_0, nvars, nobs, adata, i);
		float recover, delta = 1e-4;

		printf_s("original pi:(%g, %g), delta=%g\n", pi_0.v[0], pi_0.v[1], delta);

		for (int j = 0; j < 3; j++)
		{
			jac->colidx[base + j] = index[j];
			jac->colidx[base + j + 3] = index[j];
			recover = p[index[j]];
			p[index[j]] += delta;
			func_point(p, pi_alter, nvars, nobs, adata, i);
			printf_s("pi alter%d:(%g, %g)\n", j, pi_alter.v[0], pi_alter.v[1]);
			p[index[j]] = recover;
			jac->val[base + j] = (pi_alter.v[0] - pi_0.v[0]) / delta;
			jac->val[base + j + 3] = (pi_alter.v[1] - pi_0.v[1]) / delta;
		}
		jac->rowptr[i * 2] = base;
		jac->rowptr[i * 2 + 1] = base + 3;
		printf_s("Jacobian line %d/%d:\n", i + 1, nvars);
		for (int j = 0; j < 3; j++)
		{
			printf_s("%.5f\t\t", jac->val[base + j]);
		}
		printf_s("\n");
		for (int j = 0; j < 3; j++)
		{
			printf_s("%.5f\t\t", jac->val[base + j + 3]);
		}
		printf_s("\n");
	}
	jac->rowptr[nobs] = nvars * 6;
	exit(1);
}

void calc_pi_ls(la_vector<float> &pi, const la_matrix<float> &A,
	const la_vector<float> &b1, const la_vector<float> &b2,
	const la_vector<float> &k1, float beta, float z)
{
	assert(A.row == b1.length);

	int dim = b1.length;
	la_vector<float> b0(dim);
	for (int i = 0; i < dim; i++)
	{
		b0.v[i] = b1.v[i] + b2.v[i] / (k1.v[i] - beta*z);
	}

	la_matrix<float> A_copy;
	A_copy = A;	// create a deep copy since slsq() will change it
	slsq(pi, A_copy, b0, 1);
}

void grad_pi_ls(la_vector<float> &grad, const la_matrix<float> &A,
	const la_vector<float> &b1, const la_vector<float> &b2,
	const la_vector<float> &k1, float beta, float z, 
	boolean do_verify)
{
	assert(A.row == b1.length);
	//la_matrix<float> ATA, invATA, C;
	//smmmul(ATA, A.transpose(), A);
	//sminv(invATA, ATA);
	//smmmul(C, invATA, A.transpose());

	int dim = b1.length;
	la_vector<float> b0(dim);
	for (int i = 0; i < dim; i++)
	{
		b0.v[i] = beta*b2.v[i] / ((k1.v[i] - beta*z)*(k1.v[i] - beta*z));
	}
	//smvmul(grad, C, b0);
	la_matrix<float> A_copy;
	A_copy = A;	// create a deep copy since slsq() will change it
	slsq(grad, A_copy, b0, 1);

	assert(grad.length == 2);
	if (do_verify)
	{
		la_vector<float> numeric_grad(2), pi0, pi1;
		float delta = 1e-3;
		calc_pi_ls(pi0, A, b1, b2, k1, beta, z);
		calc_pi_ls(pi1, A, b1, b2, k1, beta, z + delta);
		numeric_grad.v[0] = (pi1.v[0] - pi0.v[0]) / delta;
		numeric_grad.v[1] = (pi1.v[1] - pi0.v[1]) / delta;
		printf("Analytical\tNumerical (grad_pi_ls)\n");
		for (int i = 0; i < grad.length; i++)
		{
			printf("%.6f\t%.6f\n", grad.v[i], numeric_grad.v[i]);
		}
	}
}

void calc_vecnorm(la_vector<float> &vecnorm, const la_vector<float> &pi)
{
	vecnorm = pi;
	svnormalize(vecnorm);
}

void grad_vec2norm(la_matrix<float> &grad, const la_vector<float> &pi, 
	boolean do_verify)
{
	float len = svnorm2(pi);
	float len_cube = len*len*len;
	grad.init(2, 2);	
	grad.m[0] = 1.0 / len - pi.v[0] * pi.v[0] / len_cube;
	grad.m[1] = -pi.v[0] * pi.v[1] / len_cube;
	grad.m[2] = grad.m[1];
	grad.m[3] = 1.0 / len - pi.v[1] * pi.v[1] / len_cube;
	
	if (do_verify)
	{		
		la_matrix<float> numeric_grad(2, 2);
		la_vector<float> pi1, norm0, norm1;
		float delta = 1e-6;
		pi1 = pi;
		pi1.v[0] += delta;
		calc_vecnorm(norm0, pi);
		calc_vecnorm(norm1, pi1);
		numeric_grad.m[0] = (norm1.v[0] - norm0.v[0]) / delta;
		numeric_grad.m[2] = (norm1.v[1] - norm0.v[1]) / delta;
		pi1 = pi;
		pi1.v[1] += delta;
		//calc_vecnorm(norm0, pi);
		calc_vecnorm(norm1, pi1);
		numeric_grad.m[1] = (norm1.v[0] - norm0.v[0]) / delta;
		numeric_grad.m[3] = (norm1.v[1] - norm0.v[1]) / delta;
		printf("Analytical\tNumerical (grad_vec2norm)\n");
		for (int i = 0; i < grad.row * grad.column; i++)
		{
			printf("%.6f\t%.6f\n", grad.m[i], numeric_grad.m[i]);
		}
	}
}

void grad_vec3norm(la_matrix<float> &grad, const vector3f &n,
	boolean do_verify)
{
	const float len = n.length();
	assert((_fpclassf(len) == FP_ZERO) || (_fpclassf(len) == FP_SUBNORMAL));
	float len_cube = len*len*len;
	grad.init(3, 3);
	grad.m[0] = 1.0 / len - n.v[0] * n.v[0] / len_cube;
	grad.m[1] = -n.v[0] * n.v[1] / len_cube;
	grad.m[2] = -n.v[0] * n.v[2] / len_cube;
	grad.m[3] = grad.m[1];
	grad.m[4] = 1.0 / len - n.v[1] * n.v[1] / len_cube;
	grad.m[5] = -n.v[1] * n.v[2] / len_cube;
	grad.m[6] = grad.m[2];
	grad.m[7] = grad.m[5];
	grad.m[8] = 1.0 / len - n.v[2] * n.v[2] / len_cube;
	if (do_verify)
	{
		la_matrix<float> numeric_grad(3, 3);
		vector3f n0, n1, n2, n3;
		float delta = 1e-7;
		n0 = n;
		n1 = n + vector3f(delta, 0, 0);
		n2 = n + vector3f(0, delta, 0);
		n3 = n + vector3f(0, 0, delta);
		n0.normalize();
		n1.normalize();
		n2.normalize();
		n3.normalize();
		for (int i = 0; i < 3; i++)
		{
			numeric_grad.m[i] = (n1.v[i] - n0.v[i]) / delta;
			numeric_grad.m[i + 3] = (n2.v[i] - n0.v[i]) / delta;
			numeric_grad.m[i + 6] = (n3.v[i] - n0.v[i]) / delta;
		}
		printf("Analytical\tNumerical (grad_vec3norm)\n");
		for (int i = 0; i < grad.row * grad.column; i++)
		{
			printf("%.6f\t%.6f\n", grad.m[i], numeric_grad.m[i]);
		}
	}
}

void calc_normal(vector3f &point_normal, const object_points& point_set,
	int p_id)
{
	// remember to update world coordinates of the point set externally
	point_normal = 
		(point_set.points[point_set.points[p_id].neighbor1].world_coord - 
		point_set.points[p_id].world_coord) ^
		(point_set.points[point_set.points[p_id].neighbor2].world_coord - 
			point_set.points[p_id].world_coord);
	//point_normal.normalize();
}

void grad_normal(la_matrix<float> &grad, const object_points& point_set,
	int p_id, boolean do_verify)
{
	//int n_points = point_set.points.size();
	//grad.init(3, n_points);
	// actually only need 3 columns: self, neighbor1, neighbor2
	int p1 = point_set.points[p_id].neighbor1;
	int p2 = point_set.points[p_id].neighbor2;

	if (point_set.points[p_id].inv_cross)
	{
		swap<int>(p1, p2);
		//printf_s("neighbors inversed\n");
	}

	float u = point_set.points[p_id].uv.x;
	float v = point_set.points[p_id].uv.y;
	float du = point_set.points[p1].uv.x - u;
	float dv = point_set.points[p2].uv.y - v;
	float z = point_set.points[p_id].world_coord.z;
	float zdu = point_set.points[p1].world_coord.z;
	float zdv = point_set.points[p2].world_coord.z;
	//printf_s("u=%g, v=%g, du=%g, dv=%g\nz=%g, zdu=%g, zdv=%g\n",
	//	u, v, du, dv, z, zdu, zdv);

	grad.init(3, 3);	
	grad.m[0] = (1 - zdv / FOCAL_LEN) * dv;
	grad.m[3] = -grad.m[0];
	grad.m[6] = -dv / FOCAL_LEN * (z - zdu);

	grad.m[1] = (1 - zdu / FOCAL_LEN) * du;
	grad.m[4] = -du / FOCAL_LEN * (z - zdv);
	grad.m[7] = -grad.m[1];

	grad.m[2] = (1 - zdv / FOCAL_LEN) * u * dv / FOCAL_LEN +
		(1 - zdu / FOCAL_LEN) * du * v / FOCAL_LEN;
	grad.m[5] = (zdv / FOCAL_LEN - 1) * u * dv / FOCAL_LEN + 
		((zdv - z) / FOCAL_LEN) * du * v / FOCAL_LEN +
		(zdv / FOCAL_LEN - 1) / FOCAL_LEN * du * dv;
	grad.m[8] = ((zdu - z) / FOCAL_LEN) * u * dv / FOCAL_LEN + 
		(zdu / FOCAL_LEN - 1) * du * v / FOCAL_LEN + 
		(zdu / FOCAL_LEN - 1) / FOCAL_LEN * du * dv;

	if (point_set.points[p_id].inv_cross)
	{
		for (int i = 0; i < 9; i++)
		{
			grad.m[i] = -grad.m[i];
		}
	}
	if (do_verify)
	{
		la_matrix<float> numeric_grad(3, 3);
		vector3f n, n0, n1, n2;
		calc_normal(n, point_set, p_id);
		float delta = 1e-3;
		object_points op_z0, op_z1, op_z2;
		op_z0 = point_set;
		op_z0.points[p_id].world_coord.z += delta;
		op_z0.update_xyz(NULL);
		calc_normal(n0, op_z0, p_id);
		op_z1 = point_set;
		op_z1.points[p1].world_coord.z += delta;
		op_z1.update_xyz(NULL);
		calc_normal(n1, op_z1, p_id);
		op_z2 = point_set;
		op_z2.points[p2].world_coord.z += delta;
		op_z2.update_xyz(NULL);
		calc_normal(n2, op_z2, p_id);
		for (int i = 0; i < 3; i++)
		{
			numeric_grad.m[i] = (n0.v[i] - n.v[i]) / delta;
			numeric_grad.m[i + 3] = (n1.v[i] - n.v[i]) / delta;
			numeric_grad.m[i + 6] = (n2.v[i] - n.v[i]) / delta;
		}
		printf("Analytical\tNumerical (grad_normal)\n");
		for (int i = 0; i < grad.row * grad.column; i++)
		{
			printf("%.6f\t%.6f\n", grad.m[i], numeric_grad.m[i]);
		}
	}	
}


float opt_compute_w(const float dot)
{
	return 0.0181197371f*dot*dot + 0.9212926391f*dot + 0.01913864113f;
}

void compute_rho(vector3f &rho, const vector3f &wo, const vector3f &n, 
	const float &kappa, const pref_env &env)
{
	vector3f wi;	
	codex::math::vector::reflect(wi, wo, n);
	float wt = opt_compute_w(wi*n);
	vector3f c_s, c_d, c;
	env.lookup(c_s, wi, kappa);
	env.lookup_cos(c_d, n);
	rho = 0.2f*c_d + 0.8f*wt*c_s;

	//rho = c_d;
	//rho = wt*c_s;
}

void vectord2f(vector3f &tar, const vector3d &src)
{
	for (int i = 0; i < 3; i++)
	{
		tar.v[i] = src.v[i];
	}
}
void vectorf2d(vector3d &tar, const vector3f &src)
{
	for (int i = 0; i < 3; i++)
	{
		tar.v[i] = src.v[i];
	}
}

void calc_pi_hat(la_vector<float> &pi, const vector3f &wo, const vector3f &n,
	const float &kappa, const pref_env &env)
{
	// let's use numerical method for now
	// pi_hat itself is a gradient of rho
	pi.init(2);
	double delta_vals[9] = { 5,2,1,0.5,0.1,0.05,1e-3,1e-4,1e-5 };
	double delta = 1;
	
	//printf_s("wo=%g,%g,%g\n", wo.x, wo.y, wo.z);

	//for (int i = 0; i < 20; i++)
	//{
		//delta *= 0.5;
		//double delta = delta_vals[i];
		vector3f rho, rho0, rho1;
		vector3f wo_norm, wo0, wo1, n_norm;
		vector3d wo_norm_3d, wo0_3d, wo1_3d, n_norm_3d;

		vectorf2d(n_norm_3d, n);
		n_norm_3d.normalize();
		vectord2f(n_norm, n_norm_3d);

		vectorf2d(wo_norm_3d, wo);
		wo_norm_3d.normalize();
		vectord2f(wo_norm, wo_norm_3d);
		//wo_norm = wo;
		//wo_norm.normalize();
		compute_rho(rho, wo_norm, n_norm, kappa, env);

		vectorf2d(wo0_3d, wo);
		wo0_3d.x += delta;
		wo0_3d.normalize();
		vectord2f(wo0, wo0_3d);
		//wo0 = wo + vector3f(delta, 0, 0);
		//wo0.normalize();
		compute_rho(rho0, wo0, n_norm, kappa, env);

		vectorf2d(wo1_3d, wo);
		wo1_3d.y += delta;
		wo1_3d.normalize();
		vectord2f(wo1, wo1_3d);
		//wo1 = wo + vector3f(0, delta, 0);
		//wo1.normalize();
		compute_rho(rho1, wo1, n_norm, kappa, env);

		float rho_mean, rho0_mean, rho1_mean;
		rho_mean = (rho.x + rho.y + rho.z) / 3;
		rho0_mean = (rho0.x + rho0.y + rho0.z) / 3;
		rho1_mean = (rho1.x + rho1.y + rho1.z) / 3;
		pi.v[0] = (rho0_mean - rho_mean) / delta;
		pi.v[1] = (rho1_mean - rho_mean) / delta;
		//printf_s("%g\t%g\t%g\npi=(%g, %g)\n", rho0_mean, rho1_mean, rho_mean, pi.v[0], pi.v[1]);

	//	printf_s("%.6f\t%.6f\n", delta, pi.v[1]);
	//}
	//printf_s("\n");
	//exit(1);
	
	//printf_s("delta=%g\n", delta);
	//for (int i = 0; i < 3; i++) printf_s("%g\t", wo0_3d.v[i] - wo_norm_3d.v[i]);
	//printf_s("\n");
	//for (int i = 0; i < 3; i++) printf_s("%g\t", wo1_3d.v[i] - wo_norm_3d.v[i]);
	//printf_s("\n");
	//exit(0);
}

void grad_pi_hat(la_matrix<float> &grad, const vector3f &wo, const vector3f &n,
	const float &kappa, const pref_env &env, boolean do_verify)
{
	// let's use numerical method for now
	// pay attention to view direction
	grad.init(2, 3);
	float delta_x = 1e-6;
	float delta_y = 1e-6; // or 1e-6 but not 1e-3
	float delta_z = 1e-7;
	//printf_s("n=%g,%g,%g\n", n.x, n.y, n.z);

	delta_x = delta_y = delta_z = 1e-6;

	//delta_x = 1e-2;
	//for (int i = 0; i < 20; i++)
	//{
	//	delta_x *= 0.5;
	//	delta_x = delta_y = delta_z;
		la_vector<float> pi, pi0, pi1, pi2;
		vector3f n0, n1, n2;
		//printf_s("reference point:\n");
		calc_pi_hat(pi, wo, n, kappa, env);
		n0 = n + vector3f(delta_x, 0, 0);
		//printf_s("n + ndx:\n");
		calc_pi_hat(pi0, wo, n0, kappa, env);
		n1 = n + vector3f(0, delta_y, 0);
		//printf_s("n + ndy:\n");
		calc_pi_hat(pi1, wo, n1, kappa, env);
		n2 = n + vector3f(0, 0, delta_z);
		//printf_s("n + ndz:\n");
		calc_pi_hat(pi2, wo, n2, kappa, env);
		grad.m[0] = (pi0.v[0] - pi.v[0]) / delta_x;
		grad.m[1] = (pi0.v[1] - pi.v[1]) / delta_x;

		grad.m[2] = (pi1.v[0] - pi.v[0]) / delta_y;
		grad.m[3] = (pi1.v[1] - pi.v[1]) / delta_y;

		grad.m[4] = (pi2.v[0] - pi.v[0]) / delta_z;
		grad.m[5] = (pi2.v[1] - pi.v[1]) / delta_z;

	//	printf_s("%g\t\t%g\n", delta_x, grad.m[0]);
	//	//printf_s("n_dx=(%g,%g,%g)\n", n0.x, n0.y, n0.z);
	//	//n0.normalize();
	//	//printf_s("n_dx_normalized=(%g,%g,%g)\n", n0.x, n0.y, n0.z);
	//}
	//exit(1);
	
	//printf_s("pi=(%g, %g)\n", pi.v[0], pi.v[1]);
	//printf_s("pi_dx=(%g, %g)\n", pi0.v[0], pi0.v[1]);
	//printf_s("pi_dy=(%g, %g)\n", pi1.v[0], pi1.v[1]);
	//printf_s("pi_dz=(%g, %g)\n", pi2.v[0], pi2.v[1]);

	if (do_verify)
	{
		printf_s("Numerical gradient for pi_hat:\n");
		for (int i = 0; i < grad.row; i++)
		{
			for (int j = 0; j < grad.column; j++)
			{
				printf_s("%.6f\t",grad.m[j * grad.row + i]);
			}
			printf_s("\n");
		}
	}
}

/*

void calc_wo(vector3f &wo, const vector3f &n,
	const vector3f &wi)
{
	// pay attention to the direction of wi (the view vector)
	codex::math::vector::reflect(wo, wi, n); 
}

void grad_wo(la_matrix<float> &grad, const vector3f &n,
	const vector3f &wi, boolean do_verify = false)
{
	// normalization of wo is not considered here
	grad.init(3, 3);
	
	grad.m[0] = -2.0 * (n * wi + n.x * wi.x);
	grad.m[1] = -2.0 * wi.y;
	grad.m[2] = -2.0 * wi.z;

	grad.m[3] = -2.0 * wi.x;
	grad.m[4] = -2.0 * (n * wi + n.y * wi.y);
	grad.m[5] = -2.0 * wi.z;

	grad.m[6] = -2.0 * wi.x;
	grad.m[7] = -2.0 * wi.y;
	grad.m[8] = -2.0 * (n * wi + n.z * wi.z);

	if (do_verify)
	{
		std::cout << "gradient of vector normalization" << std::endl << grad.m;
		la_matrix<float> numeric_grad(3, 3);
		vector3f wo, wi0, wi1, wi2, wo0, wo1, wo2;
		calc_wo(wo, n, wi);
		float delta = 1e-6;
		wi0 = wi + vector3f(delta, 0, 0);
		calc_wo(wo0, n, wi0);
		wi1 = wi + vector3f(0, delta, 0);
		calc_wo(wo1, n, wi1);
		wi2 = wi + vector3f(0, 0, delta);
		calc_wo(wo2, n, wi2);
		numeric_grad.m[0] = (wi0.x - wi.x) / delta;
		numeric_grad.m[3] = (wi0.y - wi.y) / delta;
		numeric_grad.m[6] = (wi0.z - wi.z) / delta;

		numeric_grad.m[1] = (wi1.x - wi.x) / delta;
		numeric_grad.m[4] = (wi1.y - wi.y) / delta;
		numeric_grad.m[7] = (wi1.z - wi.z) / delta;

		numeric_grad.m[2] = (wi2.x - wi.x) / delta;
		numeric_grad.m[5] = (wi2.y - wi.y) / delta;
		numeric_grad.m[8] = (wi2.z - wi.z) / delta;
		std::cout << "numeric gradient:" << std::endl << numeric_grad.m;
	}
}

void grad_uv(la_matrix<float> &grad, const la_vector<float> &wo, 
	boolean do_verify = false)
{

}
*/
