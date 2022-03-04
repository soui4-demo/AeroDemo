#include "StdAfx.h"
#include "GaussianBlur.h"
#include <math.h>
#pragma warning(disable:4244)

class Vec3 {
public:
	Vec3(float _r, float _g, float _b)
		:r(_r), g(_g), b(_b)
	{
	}

	Vec3(float v) : r(v), g(v), b(v) {}


	float r, g, b;
};

static const double PI = 3.141592653589793;

void CGaussianBlur::BoxesForGauss(int sigma, int *pBox, int n)  // standard deviation, number of boxes
{
	float wIdeal = sqrt((12.0 * sigma*sigma / n) + 2);  // Ideal averaging filter w 
	int wl = floor(wIdeal);
	if (wl % 2 == 0) wl--;
	int wu = wl + 2;

	float mIdeal = (12 * sigma*sigma - n*wl*wl - 4 * n*wl - 3 * n) / (-4 * wl - 4);
	int m = sround(mIdeal);
	// var sigmaActual = Math.sqrt( (m*wl*wl + (n-m)*wu*wu - n)/12 );

	for (int i = 0; i < n; i++)
		pBox[i] = i < m ? wl : wu;
}
// standard gaussian

void CGaussianBlur::GaussianBlur(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r)
{
	int radius = static_cast<int>(ceil(r * 2.57));
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w*ch; j += ch)
		{
			Vec3 color = Vec3(0.0f);
			float allWeights = 0.0f;
			for (int ix = i - radius; ix < i + radius + 1; ix++)
			{
				for (int iy = j - radius*ch; iy < j + (radius + 1)*ch; iy += ch)
				{
					int dsq = (iy / ch - j / ch)*(iy / ch - j / ch) + (ix - i)*(ix - i);// x^2 + y^2
					float weight = exp(float(-dsq / (2 * r * r))) / (PI * 2 * r * r);// gaussian function: 1/(2*pi*sgima^2) * e^(-(x^2+y^2)/(2*sigma^2))

					int x = smin((w - 1) * ch, smax(0, iy));
					int y = smin(h - 1, smax(0, ix));

					color.r += scl[y*w*ch + x] * weight;
					color.g += scl[y*w*ch + x + 1] * weight;
					color.b += scl[y*w*ch + x + 2] * weight;

					allWeights += weight;
				}
			}
			tcl[i*w*ch + j] = sround(color.r / allWeights);
			tcl[i*w*ch + j + 1] = sround(color.g / allWeights);
			tcl[i*w*ch + j + 2] = sround(color.b / allWeights);
		}
	}
}
// algorithm2
void CGaussianBlur::BoxBlur_2(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r)
{
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w*ch; j += ch)
		{
			Vec3 val = Vec3(0.0f);
			for (int ix = i - r; ix < i + r + 1; ix++)
			{
				for (int iy = j - r*ch; iy < j + (r + 1)*ch; iy += ch)
				{
					int x = smin((w - 1) * ch, smax(0, iy));
					int y = smin(h - 1, smax(0, ix));

					val.r += scl[y*w*ch + x];
					val.g += scl[y*w*ch + x + 1];
					val.b += scl[y*w*ch + x + 2];
				}
			}
			tcl[i*w*ch + j] = val.r / ((r + r + 1)*(r + r + 1));
			tcl[i*w*ch + j + 1] = val.g / ((r + r + 1)*(r + r + 1));
			tcl[i*w*ch + j + 2] = val.b / ((r + r + 1)*(r + r + 1));
		}
	}
}

void CGaussianBlur::GaussianBlur2(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r) {
	int bxs[3];
	BoxesForGauss(r,bxs, 3);
	BoxBlur_2(scl, tcl, w, h, ch, (bxs[0] - 1) / 2);
	BoxBlur_2(tcl, scl, w, h, ch, (bxs[1] - 1) / 2);
	BoxBlur_2(scl, tcl, w, h, ch, (bxs[2] - 1) / 2);
}

// algorithm3
void CGaussianBlur::BoxBlurH_3(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r)
{
	for (int i = 0; i < h; i++)
		for (int j = 0; j < w * ch; j += ch)
		{
			Vec3 color = Vec3(0.0f);
			for (int ix = j - r*ch; ix < j + r*ch + ch; ix += ch)
			{
				int x = smin(w * ch - ch, smax(0, ix));
				color.r += scl[i*w*ch + x];
				color.g += scl[i*w*ch + x + 1];
				color.b += scl[i*w*ch + x + 2];
			}
			tcl[i*w*ch + j] = color.r / (r + r + 1);
			tcl[i*w*ch + j + 1] = color.g / (r + r + 1);
			tcl[i*w*ch + j + 2] = color.b / (r + r + 1);
		}
}

void CGaussianBlur::BoxBlurT_3(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r)
{
	for (int i = 0; i < h; i++)
		for (int j = 0; j < w*ch; j += ch)
		{
			Vec3 color = Vec3(0.0f);
			for (int iy = i - r; iy < i + r + 1; iy++)
			{
				int y = smin(h - 1, smax(0, iy));
				color.r += scl[y*w*ch + j];
				color.g += scl[y*w*ch + j + 1];
				color.b += scl[y*w*ch + j + 2];
			}
			tcl[i*w*ch + j] = color.r / (r + r + 1);
			tcl[i*w*ch + j + 1] = color.g / (r + r + 1);
			tcl[i*w*ch + j + 2] = color.b / (r + r + 1);
		}
}

void CGaussianBlur::BoxBlur_3(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r)
{
	memcpy(tcl, scl, w*h*ch);
	BoxBlurH_3(tcl, scl, w, h, ch, r);
	BoxBlurT_3(scl, tcl, w, h, ch, r);
}

void CGaussianBlur::GaussianBlur3(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r) {
	int bxs[3];
	BoxesForGauss(r, bxs, 3);
	BoxBlur_3(scl, tcl, w, h, ch, (bxs[0] - 1) / 2);
	BoxBlur_3(tcl, scl, w, h, ch, (bxs[1] - 1) / 2);
	BoxBlur_3(scl, tcl, w, h, ch, (bxs[2] - 1) / 2);
}

// algorithm4
void CGaussianBlur::BoxBlurH_4(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r)
{
	float iarr = 1.0f / (r + r + 1.0f);
	for (int i = 0; i < h; i++) {
		int ti = i*LINEWID(w*ch);// middle index
		int li = ti;// left index
		int ri = ti + r*ch;// right index
		Vec3 fv = Vec3(scl[ti], scl[ti + 1], scl[ti + 2]);// first value
		Vec3 lv = Vec3(scl[ti + (w - 1)*ch], scl[ti + (w - 1)*ch + 1], scl[ti + (w - 1)*ch + 2]);// last value
		Vec3 val = Vec3(fv.r*(r + 1), fv.g*(r + 1), fv.b*(r + 1));// (r+1)/(2r+1)
		for (int j = 0; j < r*ch; j += ch)
		{
			val.r += scl[ti + j];
			val.g += scl[ti + j + 1];
			val.b += scl[ti + j + 2];
		}
		for (int j = 0; j <= r*ch; j += ch)
		{
			val.r += scl[ri] - fv.r;
			val.g += scl[ri + 1] - fv.g;
			val.b += scl[ri + 2] - fv.b;

			tcl[ti] = val.r*iarr;
			tcl[ti + 1] = val.g*iarr;
			tcl[ti + 2] = val.b*iarr;

			ri += ch;
			ti += ch;
		}
		for (int j = (r + 1)*ch; j < (w - r)*ch; j += ch)
		{
			val.r += scl[ri] - scl[li];
			val.g += scl[ri + 1] - scl[li + 1];
			val.b += scl[ri + 2] - scl[li + 2];

			tcl[ti] = val.r*iarr;
			tcl[ti + 1] = val.g*iarr;
			tcl[ti + 2] = val.b*iarr;

			ri += ch;
			li += ch;
			ti += ch;
		}
		for (int j = (w - r)*ch; j < w*ch; j += ch)
		{
			val.r += lv.r - scl[li];
			val.g += lv.g - scl[li + 1];
			val.b += lv.b - scl[li + 2];

			tcl[ti] = val.r*iarr;
			tcl[ti + 1] = val.g*iarr;
			tcl[ti + 2] = val.b*iarr;

			li += ch;
			ti += ch;
		}
	}
}

void CGaussianBlur::BoxBlurT_4(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r)
{
	int nLineWid = LINEWID(w*ch);
	float iarr = 1.0f / (r + r + 1.0f);
	for (int i = 0; i < w*ch; i += ch) {
		int ti = i;
		int li = ti;
		int ri = ti + r*nLineWid;
		Vec3 fv = Vec3(scl[ti], scl[ti + 1], scl[ti + 2]);
		Vec3 lv = Vec3(scl[ti + (h - 1)*nLineWid], scl[ti + (h - 1)*nLineWid + 1], scl[ti + (h - 1)*nLineWid + 2]);
		Vec3 val = Vec3((r + 1)*fv.r, (r + 1)*fv.g, (r + 1)*fv.b);
		for (int j = 0; j < r; j++)
		{
			val.r += scl[ti + j*nLineWid];
			val.g += scl[ti + j*nLineWid + 1];
			val.b += scl[ti + j*nLineWid + 2];
		}
		for (int j = 0; j <= r; j++)
		{
			val.r += scl[ri] - fv.r;
			val.g += scl[ri + 1] - fv.g;
			val.b += scl[ri + 2] - fv.b;

			tcl[ti] = val.r*iarr;
			tcl[ti + 1] = val.g*iarr;
			tcl[ti + 2] = val.b*iarr;

			ri += nLineWid;
			ti += nLineWid;
		}
		for (int j = r + 1; j < h - r; j++)
		{
			val.r += scl[ri] - scl[li];
			val.g += scl[ri + 1] - scl[li + 1];
			val.b += scl[ri + 2] - scl[li + 2];

			tcl[ti] = val.r*iarr;
			tcl[ti + 1] = val.g*iarr;
			tcl[ti + 2] = val.b*iarr;

			li += nLineWid;
			ri += nLineWid;
			ti += nLineWid;
		}
		for (int j = h - r; j < h; j++)
		{
			val.r += lv.r - scl[li];
			val.g += lv.g - scl[li + 1];
			val.b += lv.b - scl[li + 2];

			tcl[ti] = val.r*iarr;
			tcl[ti + 1] = val.g*iarr;
			tcl[ti + 2] = val.b*iarr;

			li += nLineWid;
			ti += nLineWid;
		}
	}
}

void CGaussianBlur::BoxBlur_4(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r)
{
	memcpy(tcl, scl, w*h*ch);
	BoxBlurH_4(tcl, scl, w, h, ch, r);
	BoxBlurT_4(scl, tcl, w, h, ch, r);
}

void CGaussianBlur::GaussianBlur4(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r)
{
	int bxs[3];
	BoxesForGauss(r, bxs, 3);
	BoxBlur_4(scl, tcl, w, h, ch, (bxs[0] - 1) / 2);
	BoxBlur_4(tcl, scl, w, h, ch, (bxs[1] - 1) / 2);
	BoxBlur_4(scl, tcl, w, h, ch, (bxs[2] - 1) / 2);
}
