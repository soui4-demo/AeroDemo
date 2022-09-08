#pragma once

/*
 copy from https://github.com/bigt1234/FastGaussianBlur.git
*/

class CGaussianBlur
{
public:
	static void GaussianBlur4(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void GaussianBlur3(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void GaussianBlur2(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void GaussianBlur(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);

private:
	static void BoxBlur_4(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void BoxBlurV_4(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void BoxBlurH_4(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void BoxBlur_3(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void BoxBlurV_3(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void BoxBlurH_3(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void BoxBlur_2(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void BoxesForGauss(int sigma, int pBox[3]); // standard deviation, number of boxes;
	static void BoxBlur_4_SSE(float *dst, float *src, int w, int h, int ch, int r);
	static void BoxBlurH_4_SSE(float *src, float *dst, int w, int h, int ch, int r);
	static void BoxBlurV_4_SSE(float *dst, float *src, int w, int h, int ch, int r);
};
