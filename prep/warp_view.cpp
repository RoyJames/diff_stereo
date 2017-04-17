#include "main.h"

void warp_view(std::vector<float> &img, 
	const matrix4x4f &new_R_tilde,
	vector3f &new_tau_tilde,
	const std::vector<float> &img0, const int width, const int height,
	const float f, const matrix4x4f &R_tilde, const vector3f &tau_tilde)
{
	vector3f p(0, 0, f);
	//vector3f pp;
	//pp = R_tilde*p + tau_tilde;
	vector3f dx(1, 0, 0), dy(0, 1, 0), dz(0, 0, 1);

	vector3f dx2, dy2, dz2, ndx2, ndy2, ndz2;
	dx2 = R_tilde*dx;
	dy2 = R_tilde*dy;
	dz2 = R_tilde*dz;
	ndx2 = new_R_tilde*dx;
	ndy2 = new_R_tilde*dy;
	ndz2 = new_R_tilde*dz;

	img.resize(width*height * 3, 0);// 0.5f);
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			vector2f uv;
			screen2uv(uv, vector2f(x, y), width, height);

			vector3f ray = (uv.x*ndx2 + uv.y*ndy2 + (-f)*ndz2);
			vector3f oldray(ray*dx2, ray*dy2, ray*dz2);

			vector2f uv2, scr;
			uv2.x = -f / oldray.z * oldray.x;
			uv2.y = -f / oldray.z * oldray.y;
			uv2screen(scr, uv2, width, height);

			if (scr.x >= 0 && scr.x <= width - 1 && scr.y >= 0 && scr.y <= height - 1)
			{
				vector3f c;
				lerp_img(c.v, scr, img0, width, height, 3);

				int idx = (x + y*width) * 3;
				img[idx + 0] = c.x;
				img[idx + 1] = c.y;
				img[idx + 2] = c.z;
			}
		}
	new_tau_tilde = tau_tilde + (R_tilde*p) - (new_R_tilde*p);
}

bool warp_a_pixel(vector3f &result, 
	const vector2f &uv,
	const std::vector<float> &img,
	const int width, const int height, 
	const float f, 
	const matrix4x4f &R,
	const bool b_debug)
{
	vector3f ray(uv.x, uv.y, -f), oldray;
	oldray = R*ray;

	vector2f uv2, scr;
	uv2.x = -f / oldray.z * oldray.x;
	uv2.y = -f / oldray.z * oldray.y;

	uv2screen(scr, uv2, width, height);
	if (b_debug)
		printf_s("[%f %f]\n", scr.x, scr.y);
	if (scr.x >= 0 && scr.x < width && scr.y >= 0 && scr.y < height)
	{
		lerp_img(result.v, scr, img, width, height, 3);
		return true;
	}

	return false;
}