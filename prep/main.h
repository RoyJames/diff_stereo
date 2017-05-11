#pragma once

#include <embree2\rtcore.h>
#include "BRDFLib.h"
#include "advmath.h"

#define PATH			"D:/codes/diff/data/"
#define NUM_IMG_PAIRS	32
#define NUM_IMG_VIEW	2
#define IMG_DIM			2048//512///4096//2048
#define FOCAL_LEN		5.0f
#define SHIFT_Z			20.0f
#define COMPUTE_PREP	0

using namespace codex::math;
using codex::math::utils::random_number_generator;
//using codex::math::vector2f;
//using codex::math::vector3f;
//using codex::math::vector3d;
using namespace brdflib;
using namespace advmath;

class pref_env
{
private:
	full_octa_frame<float> fr;
	int num_level, dim;
	std::vector<float> env;

public:
	void load(const char *filename);
	void lookup(vector3f &result, const vector3f &mu, const float kappa_log) const;
	void lookup_cos(vector3f & result, const vector3f & mu) const;
	float gauss_interp(const int *ax, const float *ay, float p) const;
	float gauss_filter(const float *data, float alpha, int w, float px, float py) const;
};

extern void load_area(std::vector<float> &area, const char *filename, const hemi_octa_frame<float> &fr);
extern void fit_BRDF(std::vector<vMF<float>> &result,
	const BRDF<float> &brdf,
	const int num_wo,
	const int num_EM_samples,
	const hemi_octa_frame<float> &fr,
	const std::vector<float> &area);

extern void render(
	std::vector<float> &img,
	const int width,
	const int height,

	const vector3f &center,
	const vector3f &dx,
	const vector3f &dy,
	const float f,

	const base_tri_mesh &m,
	RTCDevice device);

extern void save_img_float(const char *outputname, const float *img, const int new_width, const int new_height);
extern void load_img_float(const char *inputname, std::vector<float> &img_input, int &width, int &height);
extern void export_pointcloud_ply(const char *filename, const std::vector<codex::math::vector3f> &pos,
	const std::vector<codex::math::vector3f> *p_normal = nullptr,
	const std::vector<codex::math::vector3f> *p_color = nullptr);

extern void compute_rho(vector3f &rho, const vector3f &wo, const vector3f &n, const pref_env &env);
extern void lerp_img(float *c, const vector2f &scr, const std::vector<float> &img, const int width, const int height, const int num_channel);
extern void warp_view(std::vector<float> &img,
	const matrix4x4f &new_R_tilde,
	vector3f &new_tau_tilde,
	const std::vector<float> &img0, const int width, const int height,
	const float f, const matrix4x4f &R_tilde, const vector3f &tau_tilde);
extern bool warp_a_pixel(vector3f &result,
	const vector2f &uv,
	const std::vector<float> &img,
	const int width, const int height,
	const float f,
	const matrix4x4f &R,
	const bool b_debug = false);

extern void screen2uv(vector2f &uv, const vector2f &scr, const int width, const int height);
extern void uv2screen(vector2f &scr, const vector2f &uv, const int width, const int height);
extern void p2uv(vector2f &uv, const vector3f &xyz, const float beta);

extern void synthesize_images();

extern void test_Lambertian_v1();
extern void test_LM_v1();
extern void test_LM_v2();

