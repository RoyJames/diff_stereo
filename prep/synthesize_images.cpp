#include "main.h"

void synthesize_images()
{
	const float f = FOCAL_LEN;
	const int width = IMG_DIM;
	const int height = IMG_DIM;
	vector3f vp(0, 0, f);

	RTCDevice device = rtcNewDevice(NULL);
	base_tri_mesh m;
	m.load_obj("d:/codes/diff/data/bunny.obj");

	for (int i = 0; i < m.num_vertices(); i++)
		m.position(i).z -= SHIFT_Z;// f / 2;

	m.compute_normal();

	std::vector<float> img;

	render(img, width, height,
		vp, vector3f(1, 0, 0), vector3f(0, 1, 0),
		f, m, device);
	
	save_img_float(PATH "img00-00.raw", &img[0], width, height);

	random_number_generator<float> rng;
	for (int j = 0; j < NUM_IMG_VIEW; j++)
	{
		matrix4x4f R_abs;
		vector3f T_abs;

		//rotation_x(R_abs, PI / 6);
		identity(R_abs);
		if (j == 0)
		{
			identity(R_abs);
		} else {
			if (j == 1)
				rotation_x(R_abs, PI / 6);
			else if (j == 2)
				rotation_y(R_abs, PI / 6);
			else if (j == 3)
				rotation_x(R_abs, -PI / 6);
			else if (j == 4)
				rotation_y(R_abs, -PI / 6);
		}

		for (int i = 1; i <= NUM_IMG_PAIRS; i++)
		{
			float	scale_rot = 0.0f;//0.03f;//3.0f,
			float	scale_trans = 0.5f;// 0.5f;//0.05f;

			char filename[MAX_PATH];
			matrix4x4f rotx, roty, rotz, R_tilde, T_tilde, T;
			identity(rotx);
			identity(roty);
			identity(rotz);

			rotation_x(rotx, (rng.rand_real() * 2 - 1)*scale_rot);
			rotation_y(roty, (rng.rand_real() * 2 - 1)*scale_rot);
			rotation_z(rotz, (rng.rand_real() * 2 - 1)*scale_rot);
			//R_tilde = R_abs;
			R_tilde = rotx*roty*rotz*R_abs;
			vector3f tau_tilde(
				(rng.rand_real() * 2 - 1)*scale_trans,
				(rng.rand_real() * 2 - 1)*scale_trans,
				(rng.rand_real() * 2 - 1)*scale_trans
			);
			tau_tilde += T_abs;

			img.clear();
			render(img, width, height,
				R_tilde*vp + tau_tilde,
				R_tilde*vector3f(1, 0, 0), R_tilde*vector3f(0, 1, 0),
				f, m, device);

			sprintf_s(filename, PATH "img%02d-%02d-raw.raw", j, i);
			save_img_float(filename, &img[0], width, height);

			FILE *fp;
			sprintf_s(filename, PATH "img%02d-%02d-raw.param", j, i);
			FOPEN(fp, filename, "wb");
			fwrite(&f, sizeof(float), 1, fp);
			fwrite(R_tilde.m, sizeof(float) * 16, 1, fp);
			fwrite(tau_tilde.v, sizeof(float) * 3, 1, fp);
			fclose(fp);
		}
	}

	rtcDeleteDevice(device);
}
