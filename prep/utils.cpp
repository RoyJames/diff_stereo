#include "main.h"

void screen2uv(vector2f &uv, const vector2f &scr, const int width, const int height)
{
	uv.x = (scr.x + 0.5f) / (width / 2) - 1;
	uv.y = -(scr.y + 0.5f) / (height / 2) + 1;
}

void uv2screen(vector2f &scr, const vector2f &uv, const int width, const int height)
{
	scr.x = (uv.x + 1)*width / 2 - 0.5f;
	scr.y = (-uv.y + 1)*height / 2 - 0.5f;
}

void p2uv(vector2f &uv, const vector3f &xyz, const float beta)
{
	uv.x = xyz.x / (1.0f - beta*xyz.z);
	uv.y = xyz.y / (1.0f - beta*xyz.z);
}

void lerp_img(float *c, const vector2f &scr, const std::vector<float> &img, const int width, const int height, const int num_channel)
{
	int x[2], y[2];
	float wx[2], wy[2];

	x[0] = floor(scr.x);
	y[0] = floor(scr.y);
	x[1] = x[0] + 1;
	y[1] = y[0] + 1;

	wx[1] = scr.x - x[0];
	wx[0] = 1.0f - wx[1];
	wy[1] = scr.y - y[0];
	wy[0] = 1.0f - wy[1];

	x[0] = min(max(x[0], 0), width - 1);
	x[1] = min(max(x[1], 0), width - 1);
	y[0] = min(max(y[0], 0), height - 1);
	y[1] = min(max(y[1], 0), height - 1);

	for (int ch = 0; ch < num_channel; ch++)
		c[ch] = 0;

	for (int j = 0; j < 2; j++)
		for (int i = 0; i < 2; i++)
		{
			int idx = (x[i] + y[j] * width) * num_channel;

			for (int ch = 0; ch < num_channel; ch++)
				c[ch] += img[idx + ch] * wx[i] * wy[j];
		}
}

void save_img_float(const char *outputname, const float *img, const int new_width, const int new_height)
{
	FILE *fp;
	printf("saving image to %s\n", outputname);
	FOPEN(fp, outputname, "wb");

	codex::graphics::raw_file_header header;

	header.version = RAW_VERSION;
	header.flag = RAW_HEADER_FLOAT_STORAGE;

	header.width = new_width;
	header.height = new_height;
	header.bytes_per_channel
		= 4;
	header.num_channels
		= 3;

	fwrite(&header, sizeof(header), 1, fp);

	for (int i = 0; i < new_width*new_height; i++)
	{
		float r, g, b;

		r = img[i * 3 + 0];
		g = img[i * 3 + 1];
		b = img[i * 3 + 2];

		fwrite(&r, sizeof(r), 1, fp);
		fwrite(&g, sizeof(g), 1, fp);
		fwrite(&b, sizeof(b), 1, fp);
	}
	fclose(fp);
}

void load_img_float(const char *inputname, std::vector<float> &img_input, int &width, int &height)
{
	FILE *fp;

	FOPEN(fp, inputname, "rb");

	codex::graphics::raw_file_header header;

	fread(&header, sizeof(header), 1, fp);
	width = header.width;
	height = header.height;

	img_input.resize(width*height*header.num_channels);

	for (int i = 0; i < width*height; i++)
	{
		float fr, fg, fb, fa;

		if ((header.flag & RAW_HEADER_FLOAT_STORAGE) == 0)
		{
			BYTE	r, g, b, a;

			fread(&r, sizeof(r), 1, fp);
			fread(&g, sizeof(g), 1, fp);
			fread(&b, sizeof(b), 1, fp);
			if (header.num_channels == 4)
				fread(&a, sizeof(a), 1, fp);
			else
				a = 255;

			fr = r / 255.0f;
			fg = g / 255.0f;
			fb = b / 255.0f;
			fa = a / 255.0f;
		}
		else {
			fread(&fr, sizeof(fr), 1, fp);
			fread(&fg, sizeof(fg), 1, fp);
			fread(&fb, sizeof(fb), 1, fp);
			if (header.num_channels == 4)
				fread(&fa, sizeof(fa), 1, fp);
			else
				fa = 1.0f;
		}

		img_input[i*header.num_channels + 0] = fr;
		img_input[i*header.num_channels + 1] = fg;
		img_input[i*header.num_channels + 2] = fb;
		if (header.num_channels == 4)
			img_input[i*header.num_channels + 3] = fa;
	}

	fclose(fp);
}

void export_pointcloud_ply(const char *filename, const std::vector<codex::math::vector3f> &pos,
	const std::vector<codex::math::vector3f> *p_normal,
	const std::vector<codex::math::vector3f> *p_color)
{
	FILE *fp;

	FOPEN(fp, filename, "wt");

	fprintf_s(fp, "ply\n");
	fprintf_s(fp, "format ascii 1.0\n");
	fprintf_s(fp, "comment (C) Hongzhi Wu, Sep 2013.\n");

	fprintf_s(fp, "element vertex %d\n", pos.size());
	fprintf_s(fp, "property float x\n");
	fprintf_s(fp, "property float y\n");
	fprintf_s(fp, "property float z\n");
	if (p_normal)
	{
		fprintf_s(fp, "property float nx\n");
		fprintf_s(fp, "property float ny\n");
		fprintf_s(fp, "property float nz\n");
	}
	if (p_color)
	{
		fprintf_s(fp, "property uchar red\n");
		fprintf_s(fp, "property uchar green\n");
		fprintf_s(fp, "property uchar blue\n");
		fprintf_s(fp, "property uchar alpha\n");
	}

	fprintf_s(fp, "end_header\n");

	for (int i = 0; i < pos.size(); i++)
	{
		fprintf_s(fp, "%g %g %g ", pos[i].x, pos[i].y, pos[i].z);

		if (p_normal)
		{
			fprintf_s(fp, "%g %g %g ", (*p_normal)[i].x, (*p_normal)[i].y, (*p_normal)[i].z);
		}

		if (p_color)
		{
			int r = max(min(int((*p_color)[i].x * 255), 255), 0),
				g = max(min(int((*p_color)[i].y * 255), 255), 0),
				b = max(min(int((*p_color)[i].z * 255), 255), 0);
			fprintf_s(fp, "%d %d %d %d ", r, g, b, 255);
		}

		fprintf_s(fp, "\n");
	}

	fclose(fp);
}
