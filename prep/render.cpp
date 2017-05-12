#include "main.h"
#include <embree2\rtcore_ray.h>

bool debug_flag;

class rt_result
{
public:
	vector2f	bary;
	int			tri_idx;
};

class embree_cam
{
public:
	embree_cam(
		const vector3f&     cam_center,
		const vector3f&     cam_up,
		const vector3f&     cam_lookat,
		const float         fx,
		const int           width,
		const int           height
	);

	void gen_ray(RTCRay& ray, const int x, const int y) const;

private:
	const int       w, h;
	const float     cx, cy;
	const vector3f  aim, right, up;
	const vector3f  org;
	const float     px;
	const vector3f  dx, dy;
};

embree_cam::embree_cam(
	const vector3f&     cam_center,
	const vector3f&     cam_up,
	const vector3f&     cam_lookat,
	const float         fx,
	const int           width,
	const int           height
) :
	w(width), h(height),
	cx(width / 2.f), cy(height / 2.f),
	aim((cam_lookat - cam_center).normalized_vector()),
	right((aim^cam_up).normalized_vector()),
	up(right^aim),
	org(cam_center),
	//px(tanf(fovx*float(PI) / 2 / 180.f) / (width / 2.f)),
	px(1.0f / fx / (width / 2.f)),
	dx(right*px),
	dy(-up*px)
{
	//printf_s("%g %g %g\n", cam_center.x, cam_center.y, cam_center.z);
	//printf_s("%g %g %g\n", org.x, org.y, org.z);
	//printf_s("%g %g %g\n", cam_lookat.x, cam_lookat.y, cam_lookat.z);
	//printf_s("aim %g %g %g\n", aim.x, aim.y, aim.z);
	//printf_s("up %g %g %g\n", up.x, up.y, up.z);
	//printf_s("%g %g %g\n", right.x, right.y, right.z);
}

void embree_cam::gen_ray(RTCRay& ray, const int x, const int y) const
{
	ray.org[0] = org.x;
	ray.org[1] = org.y;
	ray.org[2] = org.z;

	vector3f dir = (x + 0.5f - cx)*dx + (y + 0.5f - cy)*dy + aim;
	dir.normalize();

	ray.dir[0] = dir.x;
	ray.dir[1] = dir.y;
	ray.dir[2] = dir.z;

	ray.tnear = 0.f;
	ray.tfar = FLT_MAX;
	ray.geomID = RTC_INVALID_GEOMETRY_ID;
	ray.primID = RTC_INVALID_GEOMETRY_ID;
	ray.mask = -1;
	ray.time = 0;
}

void setup_scene(
	RTCScene&                           scene,
	const std::vector<vector3f>&        tri_points,
	const std::vector<int>&             tri_index,
	RTCDevice							device)
{
	// Initialize the scene
	scene = rtcDeviceNewScene(device, RTC_SCENE_STATIC, RTC_INTERSECT1);

	int triNum = tri_index.size() / 3;
	int vertNum = tri_points.size();
	//printf_s("%d tri %d vertex\n", triNum, vertNum);

	unsigned geomID = rtcNewTriangleMesh(scene, RTC_GEOMETRY_STATIC, triNum, vertNum, 1);

	// Set vertices
	float *vertex = (float *)rtcMapBuffer(scene, geomID, RTC_VERTEX_BUFFER);

	for (int vertId = 0; vertId<vertNum; vertId++) {
		int local_index = vertId * 4;
		vertex[local_index] = tri_points[vertId].x;
		vertex[local_index + 1] = tri_points[vertId].y;
		vertex[local_index + 2] = tri_points[vertId].z;
	}

	rtcUnmapBuffer(scene, geomID, RTC_VERTEX_BUFFER);

	// Set triangles
	int* index = (int *)rtcMapBuffer(scene, geomID, RTC_INDEX_BUFFER);
	memcpy_s(index, sizeof(int) * 3 * triNum, tri_index.size() ? &tri_index[0] : nullptr, sizeof(int) * 3 * triNum);
	rtcUnmapBuffer(scene, geomID, RTC_INDEX_BUFFER);

	// Commit the scene
	rtcCommit(scene);
	//printf_s("Mesh init done.\n");
}

void intersect(rt_result &result, const RTCScene& scene, RTCRay& ray)
{
	rtcIntersect(scene, ray);

	// Get the correspondence
	result.tri_idx = -1;
	if (ray.primID != RTC_INVALID_GEOMETRY_ID || ray.tfar != FLT_MAX)
	{
		result.tri_idx = ray.primID;
		result.bary = vector2f(ray.u, ray.v);
	}
}

float compute_w(const float dot)
{
	return 0.0181197371f*dot*dot + 0.9212926391f*dot + 0.01913864113f;
}

float compute_kappa_log(const float dot)
{
	return 0.5173923471f*dot*dot - 1.45329596f*dot + 3.321729293f;
}

void compute_rho(vector3f &rho, const vector3f &wo, const vector3f &n, const pref_env &env)
{
	vector3f wi;
	codex::math::vector::reflect(wi, wo, n);

	//float wt = compute_w(wi*n);
	//float kappa_log = compute_kappa_log(wi*n);
	vector3f c_s, c_d, c;
	//env.lookup(c_s, wi, kappa_log);
	env.lookup_cos(c_d, n);
	//rho = 0.2f*c_d + 0.8f*wt*c_s;
	rho = c_d;
	if (debug_flag) {
		printf_s("(%g,%g,%g)\n", rho.x, rho.y, rho.z);
	}

	{
		//char filename[MAX_PATH];
		//sprintf_s(filename, "d:/codes/diff/data/groundtruth_kappa_log.dat");
		//FILE *file;
		//if ((file = fopen(filename, "ab")) == NULL)
		//{
		//	printf("open file failed!\n");
		//	exit(1);
		//}
		//printf_s("%g\n", kappa_log);
		//fwrite(&kappa_log, sizeof(float), 1, file);
		//fclose(file);
	}
	//rho = c_d;
	//rho = wt*c_s;
}

void render(
	std::vector<float> &img,
	const int width,
	const int height,

	const vector3f &center,
	const vector3f &dx,
	const vector3f &dy,
	const float f,

	const base_tri_mesh &m,
	RTCDevice device,
	int gd_idx)
{
	std::vector<int> face;
	for (int i = 0; i < m.num_faces(); i++)
	{
		const triangle_index f = m.face(i);
		face.push_back(f.i0);
		face.push_back(f.i1);
		face.push_back(f.i2);
	}

	std::vector<vector3f> pos;
	for (int i = 0; i < m.num_vertices(); i++)
	{
		vector3 p = m.position(i);
		vector3f v(float(p.x), float(p.y), float(p.z));
		pos.push_back(v);		
	}

	img.resize(3 * width*height);
	for (int i = 0; i < width*height; i++)
	{
		img[i * 3 + 0] = 1;
		img[i * 3 + 1] = 0;
		img[i * 3 + 2] = 1;
	}

	RTCScene scene;
	setup_scene(scene, pos, face, device);

	pref_env env;
	env.load("d:/codes/diff/data/uffizi.pref_env");

	embree_cam camera(center, dy, center-(dx^dy)*f, f, width, height);	
	//printf_s("up   [%g %g %g]\n", dy.x, dy.y, dy.z);
	//printf_s("view [%g %g %g]\n", (- (dx^dy)).x, (-(dx^dy)).y, (-(dx^dy)).z);

	// compute groundtruth depth and kappa_log
	// you MUST comment this segment out when not used!
	if(gd_idx >=0)
	{
		//printf("%g %g %g\n", center.x, center.y, center.z);
		RTCRay ray;
		char filename[MAX_PATH];
		//sprintf_s(filename, "d:/codes/diff/data/groundtruth_depth.dat");
		//FILE *file;
		//if ((file = fopen(filename, "wb")) == NULL)
		//{
		//	printf("open file failed!\n");
		//	exit(1);
		//}
		//sprintf_s(filename, "d:/codes/diff/data/groundtruth_normal.dat");
		//FILE *file_normal;
		//if ((file_normal = fopen(filename, "wb")) == NULL)
		//{
		//	printf("open file failed!\n");
		//	exit(1);
		//}
		sprintf_s(filename, "d:/codes/diff/data/groundtruth_lambertian_view%02d.dat", gd_idx);
		FILE *file_lambertian;
		if ((file_lambertian = fopen(filename, "wb")) == NULL)
		{
			printf("open file failed!\n");
			exit(1);
		}
		//int scan_begin = 1000, scan_end = 1048;
		//fwrite(&scan_begin, sizeof(int), 1, file);
		//fwrite(&scan_end, sizeof(int), 1, file);
		//remove("d:/codes/diff/data/groundtruth_kappa_log.dat");
		//for (int i = scan_begin; i <= scan_end; i++)
		std::vector<float> gd_lambertian;
		std::vector<std::pair<int, int> > coord_lambertian;
		int cnt = 0;
		//int sample_dist = 10;
		for (int i = 0; i < height; i++)
		{
			//for (int j = scan_begin; j <= scan_end; j++)
			for (int j = 0; j < width; j++)
			{
				camera.gen_ray(ray, i, j);
				vector2f src(i, j);
				vector2f uv;
				screen2uv(uv, src, IMG_DIM, IMG_DIM);
				rt_result result;
				intersect(result, scene, ray);
				if (result.tri_idx != -1)
				{
					cnt++;
					const triangle_index t = m.face(result.tri_idx);
					vector3d inter_pos = m.position(t.i1)*result.bary.x +
						m.position(t.i2)*result.bary.y +
						m.position(t.i0)*(1 - result.bary.x - result.bary.y);
					//printf("scan=%.6f\tz=%g\n", uv.x, inter_pos.z);
					vector3d nd = m.normal(t.i1)*result.bary.x +
						m.normal(t.i2)*result.bary.y +
						m.normal(t.i0)*(1 - result.bary.x - result.bary.y);
					nd.normalize();

					vector3f n(nd.x, nd.y, nd.z);
					vector3f wo(-ray.dir[0], -ray.dir[1], -ray.dir[2]);

					//vector3f c;
					//compute_rho(c, wo, n, env);
					//fwrite(&inter_pos.z, sizeof(double), 1, file);
					//fwrite(n.v, sizeof(float), 3, file_normal);
					//printf_s("%g\n", inter_pos.z);
					//if(cnt % sample_dist == 0)
					gd_lambertian.push_back(inter_pos.z);
					coord_lambertian.push_back(std::pair<int, int>(i, j));
				}
				else {
					//printf("No intersection.\n");
				}
			}
		}
		//fclose(file);
		//fclose(file_normal);
		int n_valid = gd_lambertian.size();
		printf("%d points are valid\n", n_valid);
		fwrite(&n_valid, sizeof(int), 1, file_lambertian);
		for (int i = 0; i < gd_lambertian.size(); i++) {
			fwrite(&coord_lambertian[i].first, sizeof(int), 1, file_lambertian);
			fwrite(&coord_lambertian[i].second, sizeof(int), 1, file_lambertian);
			fwrite(&gd_lambertian[i], sizeof(float), 1, file_lambertian);
		}
		fclose(file_lambertian);
		printf("writing groundtruth file for view %d done\n", gd_idx);
		//exit(1);
	}

#pragma omp parallel for
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			int idx = x + y*width;

			RTCRay ray;
			camera.gen_ray(ray, x, y);

			//if (x == 637 && y == 259)
			//{
			//	printf_s("ray(%gf, %gf, %gf)\n", ray.dir[0], ray.dir[1], ray.dir[2]);
			//	printf_s("org(%gf, %gf, %gf)\n", ray.org[0], ray.org[1], ray.org[2]);
			//}
			if (x == 1023 && y == 1171) {
				debug_flag = true;
			}
			else {
				debug_flag = false;
			}
			debug_flag = false;

			rt_result result;
			intersect(result, scene, ray);
			if (result.tri_idx != -1)
			{
				const triangle_index t = m.face(result.tri_idx);
				vector3d nd = m.normal(t.i1)*result.bary.x + 
					m.normal(t.i2)*result.bary.y + 
					m.normal(t.i0)*(1 - result.bary.x - result.bary.y);
				nd.normalize();

				vector3f n(nd.x, nd.y, nd.z);
				vector3f wo(-ray.dir[0], -ray.dir[1], -ray.dir[2]);

				vector3f c;
				compute_rho(c, wo, n, env);

				//DEBUG
				//c = (n+vector3f(1, 1, 1))/2;
				img[idx * 3 + 0] = c.x;
				img[idx * 3 + 1] = c.y;
				img[idx * 3 + 2] = c.z;

				//////DEBUG
				//////if (x == 2 && y == 512)
				//if (x == IMG_DIM / 2 - 24 && y == IMG_DIM / 2 + 24)
				//{
				//	vector3f pos = vector3f(ray.dir[0], ray.dir[1], ray.dir[2])*ray.tfar + vector3f(ray.org[0], ray.org[1], ray.org[2]);
				//	printf_s("ray(%gf, %gf, %gf)\n", ray.dir[0], ray.dir[1], ray.dir[2]);
				//	printf_s("org(%gf, %gf, %gf)\n", ray.org[0], ray.org[1], ray.org[2]);
				//	printf_s("n(%gf, %gf, %gf)\n", nd.x, nd.y, nd.z);
				//	printf_s("v(%gf, %gf, %gf)\n", wo.x, wo.y, wo.z);
				//	printf_s("normalized v(%gf, %gf, %gf)\n", ray.tfar*wo.x, ray.tfar*wo.y, ray.tfar*wo.z);
				//	//	//printf_s("t(%g %g)\n", ray.tnear, ray.tfar);
				//	printf_s("pos(%gf, %gf, %gf)\n", pos.x, pos.y, pos.z);
				//	printf_s("c(%gf, %gf, %gf)\n", c.x, c.y, c.z);
				//	printf_s("\n");
				//	exit(1);
				//}
			}
		}
	rtcDeleteScene(scene);
}