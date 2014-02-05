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

void interpolate(INTERPOLATED_IMG *i_img, X3F *x3f);
void simple_coeff (INTERPOLATED_IMG *img, int index);
float foveon_avg (short *pix, int range[2], float cfilt);
short *foveon_make_curve (double max, double mul, double filt);
void foveon_make_curves(short **curvep, float dq[3], float div[3], float filt);
int foveon_apply_curve (short *curve, int i);
void apply_gamma(INTERPOLATED_IMG *interpolated, int output_bps, float bright, char *ifname);
void foveon_f20_interpolate(INTERPOLATED_IMG *i_img, X3F *x3f);
#endif
