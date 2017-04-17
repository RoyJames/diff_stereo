#pragma once
#include "main.h"
#include <embree2\rtcore_ray.h>
#include <splm.h>

typedef struct {
	int scr_u, scr_v;
	vector2f uv;
	vector3f world_coord;
	int up = -1;
	int down = -1;
	int left = -1;
	int right = -1;
	int neighbor1 = -1;	// the first neighbor to appear in cross product
	int neighbor2 = -1;	// the second neighbor to appear in cross product
	bool inv_cross = false;	// true if we use cross(v,u) instead of cross(u,v)
}point_info;

class object_points
{
public:
	std::vector<point_info> points;	// neighbors stored in order: right, down, left, up
public:
	void register_point(int scr_u, int scr_v, float depth = 0.0)
	{
		point_info new_point;
		new_point.scr_u = scr_u;	new_point.scr_v = scr_v;	// (u,v) should range from (-1,1)
		screen2uv(new_point.uv, vector2f(scr_u, scr_v), IMG_DIM, IMG_DIM);
		new_point.world_coord = vector3f(0, 0, depth);
		points.push_back(new_point);
	};
	void update_xyz(const double *z)
	{
		if (z == NULL)
		{
			for (int i = 0; i < points.size(); i++)
			{
				points[i].world_coord.x = points[i].uv.x *
					(1 - points[i].world_coord.z / FOCAL_LEN);
				points[i].world_coord.y = points[i].uv.y *
					(1 - points[i].world_coord.z / FOCAL_LEN);
			}
		}
		else
		{
			for (int i = 0; i < points.size(); i++)
			{
				points[i].world_coord.z = z[i];
				points[i].world_coord.x = points[i].uv.x *
					(1 - points[i].world_coord.z / FOCAL_LEN);
				points[i].world_coord.y = points[i].uv.y *
					(1 - points[i].world_coord.z / FOCAL_LEN);
			}
		}
	}
	void match_neighbors()
	{
		int n_points = points.size();
		for (int i = 0; i < n_points; i++)
		{
			// find all its neighbors
			for (int j = 0; j < n_points; j++)
			{
				if (points[i].scr_u + 1 == points[j].scr_u &&
					points[i].scr_v == points[j].scr_v)
					points[i].right = j;
				else if (points[i].scr_u == points[j].scr_u &&
					points[i].scr_v + 1 == points[j].scr_v)
					points[i].down = j;
				else if (points[i].scr_u - 1 == points[j].scr_u &&
					points[i].scr_v == points[j].scr_v)
					points[i].left = j;
				else if (points[i].scr_u == points[j].scr_u &&
					points[i].scr_v - 1 == points[j].scr_v)
					points[i].up = j;
			}
			// tag the representative neighbors
			if (points[i].right != -1)
			{
				if (points[i].down != -1)
				{
					points[i].neighbor1 = points[i].down;
					points[i].neighbor2 = points[i].right;
					points[i].inv_cross = true;
					continue;
				}
				else if (points[i].up != -1)
				{
					points[i].neighbor1 = points[i].right;
					points[i].neighbor2 = points[i].up;
					continue;
				}
			}
			else if (points[i].down != -1)
			{
				if (points[i].left != -1)
				{
					points[i].neighbor1 = points[i].left;
					points[i].neighbor2 = points[i].down;
					continue;
				}
			}
			else if (points[i].left != -1)
			{
				if (points[i].up != -1)
				{
					points[i].neighbor1 = points[i].up;
					points[i].neighbor2 = points[i].left;
					points[i].inv_cross = true;
					continue;
				}
			}
		}
	};
};


struct mydata {
	object_points point_set;
	pref_env env;
	float *kappa;
	// the following vectors should all have dimension = NUM_IMG_PAIRS
	std::vector<vector3f> Delta_v;
	std::vector<float> k1;
	// the following vectors should all have 1st dimension = # of unknowns
	// and 2nd dimension = NUM_IMG_PAIRS
	std::vector<la_vector<float> > Ic_Iut;
	std::vector<la_vector<float> > k2Ia_k3Ib;
};

void save_mydata_file(const char *filename, const mydata &dat);
void read_mydata_file(const char *filename, mydata &dat);

void ffunc(double *p, double *hx, int nvars, int nobs, void *adata);
void fjac(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata);
void fjac_diff(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata);
void func_point(double *p, la_vector<float> &pi_diff, int nvars, int nobs, void *adata, int p_id);
void fjac_mannual_diff(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata);

void calc_pi_ls_err(double *p, int p_id, void *adata);

void calc_pi_ls(la_vector<float> &pi, const la_matrix<float> &A,
	const la_vector<float> &b1, const la_vector<float> &b2,
	const la_vector<float> &k1, float beta, float z);
void grad_pi_ls(la_vector<float> &grad, const la_matrix<float> &A,
	const la_vector<float> &b1, const la_vector<float> &b2,
	const la_vector<float> &k1, float beta, float z,
	boolean do_verify = false);
void calc_vecnorm(la_vector<float> &vecnorm, const la_vector<float> &pi);
void grad_vec2norm(la_matrix<float> &grad, const la_vector<float> &pi,
	boolean do_verify = false);
void grad_vec3norm(la_matrix<float> &grad, const vector3f &n,
	boolean do_verify = false);
void calc_normal(vector3f &point_normal, const object_points& point_set,
	int p_id);
void grad_normal(la_matrix<float> &grad, const object_points& point_set,
	int p_id, boolean do_verify = false);
void calc_pi_hat(la_vector<float> &pi, const vector3f &wo, const vector3f &n,
	const float &kappa, const pref_env &env);
void grad_pi_hat(la_matrix<float> &grad, const vector3f &wo, const vector3f &n,
	const float &kappa, const pref_env &env, boolean do_verify = false);
