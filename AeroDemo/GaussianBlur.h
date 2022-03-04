#pragma once

/*
 copy from https://github.com/bigt1234/FastGaussianBlur.git
*/
#define BYTE_ALIGN 4
#define LINEWID(x) (((x)+BYTE_ALIGN-1)/BYTE_ALIGN*BYTE_ALIGN)

class CGaussianBlur
{
public:
	static void GaussianBlur4(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void GaussianBlur3(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void GaussianBlur2(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void GaussianBlur(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);

private:
	static void BoxBlur_4(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void BoxBlurT_4(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void BoxBlurH_4(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void BoxBlur_3(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void BoxBlurT_3(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void BoxBlurH_3(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void BoxBlur_2(unsigned char *scl, unsigned char *tcl, int w, int h, int ch, int r);
	static void BoxesForGauss(int sigma, int *pBox, int n); // standard deviation, number of boxes;
};
