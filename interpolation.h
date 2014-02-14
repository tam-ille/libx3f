#ifndef __LIB_X3F_INTERPOLATION__
#define __LIB_X3F_INTERPOLATION__

#include "raw_x3f.h"

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#define SQR(x) ((x)*(x))
#define LIM(x,min,max) MAX(min,MIN(x,max))
#define ULIM(x,y,z) ((y) < (z) ? LIM(x,y,z) : LIM(x,z,y))
#define CLIP(x) LIM(x,0,65535)
#define SWAP(a,b) { a=a+b; b=a-b; a=a-b; }

void x3f_interpolate(INTERPOLATED_IMG *i_img, X3F *x3f);
void x3f_simple_coeff (INTERPOLATED_IMG *img, int index);
float x3f_foveon_avg (short *pix, int range[2], float cfilt);
short *x3f_foveon_make_curve (double max, double mul, double filt);
void x3f_foveon_make_curves(short **curvep, float dq[3], float div[3], float filt);
int x3f_foveon_apply_curve (short *curve, int i);
void x3f_apply_gamma(INTERPOLATED_IMG *interpolated, float bright);
void x3f_output_ppm(INTERPOLATED_IMG *interpolated, int output_bps, int flip, char *ifname);
void x3f_trueII_interpolate(INTERPOLATED_IMG *i_img, X3F *x3f);
#endif
