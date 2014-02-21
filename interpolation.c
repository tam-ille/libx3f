
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <glib.h>
#include <float.h>
#include <limits.h>
#include <sys/types.h>
#include <math.h>
#include <wchar.h>
#include <locale.h>
#include "interpolation.h"
#include "raw_x3f.h"

float x3f_foveon_avg (short *pix, int range[2], float cfilt)
{
  int i;
  float val, min=FLT_MAX, max=-FLT_MAX, sum=0;

  for (i=range[0]; i <= range[1]; i++) {
    sum += val = pix[i*4] + (pix[i*4]-pix[(i-1)*4]) * cfilt;
    if (min > val) min = val;
    if (max < val) max = val;
  }
  if (range[1] - range[0] == 1) return sum/2;
  return (sum - min - max) / (range[1] - range[0] - 1);
}

short *x3f_foveon_make_curve (double max, double mul, double filt)
{
  short *curve;
  unsigned i, size;
  double x;

  if (!filt) filt = 0.8;
  size = 4*M_PI*max / filt;
  if (size == UINT_MAX) size--;
  curve = (short *) calloc (size+1, sizeof *curve);
/*   merror (curve, "foveon_make_curve()"); */
  curve[0] = size;
  for (i=0; i < size; i++) {
    x = i*filt/max/4;
    curve[i+1] = (cos(x)+1)/2 * tanh(i*filt/mul) * mul + 0.5;
  }
  return curve;
}

void x3f_foveon_make_curves
	(short **curvep, float dq[3], float div[3], float filt)
{
  double mul[3], max=0;
  int c;

  FORC3 mul[c] = dq[c]/div[c];
  FORC3 if (max < mul[c]) max = mul[c];
  FORC3 curvep[c] = x3f_foveon_make_curve (max, mul[c], filt);
}

int x3f_foveon_apply_curve (short *curve, int i)
{
  if (abs(i) >= curve[0]) return 0;
  return i < 0 ? -curve[1-i] : curve[1+i];
}

void x3f_gamma_curve(INTERPOLATED_IMG *interpolated, double pwr, double ts, int mode, int imax){
  int i, c;
  double g[6], bnd[2]={0,0}, r;

  printf("=> gamma_curve\n");
  g[0] = pwr;
  g[1] = ts;
  g[2] = g[3] = g[4] = 0;
  bnd[g[1] >= 1] = 1;
  if (g[1] && (g[1]-1)*(g[0]-1) <= 0) {
    for (i=0; i < 48; i++) {
      g[2] = (bnd[0] + bnd[1])/2;
      if (g[0]) bnd[(pow(g[2]/g[1],-g[0]) - 1)/g[0] - 1/g[2] > -1] = g[2];
      else	bnd[g[2]/exp(1-1/g[2]) < g[1]] = g[2];
    }
    g[3] = g[2] / g[1];
    if (g[0]) g[4] = g[2] * (1/g[0] - 1);
  }
  if (g[0]) g[5] = 1 / (g[1]*SQR(g[3])/2 - g[4]*(1 - g[3]) +
		(1 - pow(g[3],1+g[0]))*(1 + g[4])/(1 + g[0])) - 1;
  else      g[5] = 1 / (g[1]*SQR(g[3])/2 + 1
		- g[2] - g[3] -	g[2]*g[3]*(log(g[3]) - 1)) - 1;
  if (!mode--) {
    memcpy (interpolated->gamm, g, sizeof interpolated->gamm);
    return;
  }
  for (i=0; i < 0x10000; i++) {
    interpolated->curve[i] = 0xffff;
    if ((r = (double) i / imax) < 1)
      interpolated->curve[i] = 0x10000 * ( mode
	? (r < g[3] ? r*g[1] : (g[0] ? pow( r,g[0])*(1+g[4])-g[4]    : log(r)*g[2]+1))
	: (r < g[2] ? r/g[1] : (g[0] ? pow((r+g[4])/(1+g[4]),1/g[0]) : exp((r-1)/g[2]))));
  }
}

void x3f_convert_to_rgb(INTERPOLATED_IMG *interpolated, int output_color)
{
  uint row, col, c, i, j, k;
/*   unsigned *oprof; */
  uint16_t *img, (*tmp)[4];
  float out[3], out_cam[3][3];
  static const double xyzd50_srgb[3][3] =
  { { 0.436083, 0.385083, 0.143055 },
    { 0.222507, 0.716888, 0.060608 },
    { 0.013930, 0.097097, 0.714022 } };
  static const double rgb_rgb[3][3] =
  { { 1,0,0 }, { 0,1,0 }, { 0,0,1 } };
  static const double adobe_rgb[3][3] =
  { { 0.715146, 0.284856, 0.000000 },
    { 0.000000, 1.000000, 0.000000 },
    { 0.000000, 0.041166, 0.958839 } };
  static const double wide_rgb[3][3] =
  { { 0.593087, 0.404710, 0.002206 },
    { 0.095413, 0.843149, 0.061439 },
    { 0.011621, 0.069091, 0.919288 } };
  static const double prophoto_rgb[3][3] =
  { { 0.529317, 0.330092, 0.140588 },
    { 0.098368, 0.873465, 0.028169 },
    { 0.016879, 0.117663, 0.865457 } };
  static const double xyz_rgb[3][3] = {			/* XYZ from RGB */
    { 0.412453, 0.357580, 0.180423 },
    { 0.212671, 0.715160, 0.072169 },
    { 0.019334, 0.119193, 0.950227 } };
  static const double (*out_rgb[])[3] =
  { rgb_rgb, adobe_rgb, wide_rgb, prophoto_rgb, xyz_rgb };
  static const char *name[] =
  { "sRGB", "Adobe RGB (1998)", "WideGamut D65", "ProPhoto D65", "XYZ" };

  printf("Converting to %s colorspace...\n", name[output_color-1]);

  x3f_gamma_curve (interpolated, interpolated->gamm[0], interpolated->gamm[1], 0, 0);
  memcpy (out_cam, interpolated->rgb_cam, sizeof out_cam);

  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      for (out_cam[i][j] = k=0; k < 3; k++)
	out_cam[i][j] += out_rgb[output_color-1][i][k] * interpolated->rgb_cam[k][j];


  memset (interpolated->histogram, 0, sizeof interpolated->histogram);
  tmp=(uint16_t (*)[4])interpolated->img;
  for (img=tmp[0], row=0; row < interpolated->height; row++)
    for (col=0; col < interpolated->width; col++, img+=4) {
      out[0] = out[1] = out[2] = 0;
      FORC3 {
		out[0] += out_cam[0][c] * img[c];
		out[1] += out_cam[1][c] * img[c];
		out[2] += out_cam[2][c] * img[c];
      }
      FORC3 img[c] = CLIP((uint) out[c]);
      FORC3 interpolated->histogram[c][img[c] >> 3]++;
    }
}

void x3f_simple_coeff (INTERPOLATED_IMG *img, int index)
{
  static const float table[][12] = {
  /* index 0 -- all Foveon cameras */
  {   1.4032, -0.2231, -0.1016,
	-0.5263, 1.4816, 0.017,
	-0.0112, 0.0183, 0.9113 },

  {    1.0832,-0.1931,-0.1016,
	-0.5763,1.0616,0.017,
	-0.0612,0.0183,0.8613 },

  {   1.2032, -0.1231, -0.1016,
	-0.6263, 1.4816, 0.017,
	-0.1112, 0.0183, 0.9113 },

  {   1.068924, -0.098145, 0.029205,
	-0.103180, 0.923837, 0.072983,
	-0.037827, -0.023300, 1.085474 },

  { 1.1032,-0.5231,-0.3016,-0.7263,1.1816,-0.017,-0.0112,0.0183,0.9113 },
  };
  int i, c;

  for (i=0; i < 3; i++)
    FORC3 img->rgb_cam[i][c] = table[index][i*3+c];
}


int flip_index (int row, int col, int height, int width, int flip)
{
  if (flip & 4) SWAP(row,col);
  if (flip & 2) row = height - 1 - row;
  if (flip & 1) col = width  - 1 - col;
  return row * width + col;
}

void x3f_stretch(){
  /* for resolution=MED, we will have to:
for each channel, average the value of pixel[-1]+pixel[0]+pixel[+1] and move 2 pixels forward. (ok for width)
then row[0] do not change, for each row, average the value of pix[0-width]+pix[0]+pix[0+width]

   }
  */
}

void x3f_apply_gamma(INTERPOLATED_IMG *interpolated, float bright){
  int perc, val, total, white=0x2000;
  int flip=0;
  uint c, row, col;
  uint16_t *corrected;
  uint16_t (*image)[4];

  image=(uint16_t (*)[4])interpolated->img;

  perc = interpolated->width * interpolated->height * 0.01;		/* 99th percentile white level */
  for (white=c=0; c < 3; c++) {
  	for (val=0x2000, total=0; --val > 32; )
  	  if ((total += interpolated->histogram[c][val]) > perc) break;
  	if (white < val) white = val;
  }

  x3f_gamma_curve (interpolated, interpolated->gamm[0], interpolated->gamm[1], 2, (white << 3)/bright);
  for (row=0; row < interpolated->height; row++) {
	for (col=0; col < interpolated->width; col++) {
	  FORC3 image[row*interpolated->width+col][c] = interpolated->curve[image[row*interpolated->width+col][c]];
	}
  }
  printf("Applied gamma curve\n");
}

void x3f_output_ppm(INTERPOLATED_IMG *interpolated, int output_bps, int flip, char *ifname){
  FILE *ofp;
  char *ofname, *cp;
  unsigned char *ppm;
  ushort *ppm2;
  uint c, row, col;
  int soff, rstep, cstep;
  int perc, val, total, white=0x2000;
  uint16_t (*image)[4];

  image=(uint16_t (*)[4])interpolated->img;

  ofname = (char *) malloc (strlen(ifname) + 64);

  strcpy (ofname, ifname);
  if ((cp = strrchr (ofname, '.'))) *cp = 0;
  strcat (ofname, ".ppm");
  ofp = fopen (ofname, "wb");
  if (!ofp) {
	return;
  }
  printf ("Writing data to %s ...\n", ofname);
  ppm = (unsigned char *) calloc (interpolated->width, 3*output_bps/8);
  ppm2 = (ushort *) ppm;
  /*     if (oprof) */
  /*       fwrite (oprof, ntohl(oprof[0]), 1, ofp); */
  fprintf (ofp, "P6\n%d %d\n%d\n",
		   interpolated->width, interpolated->height, (1 << output_bps)-1);

  soff  = flip_index (0, 0, interpolated->height, interpolated->width, flip);
  cstep = flip_index (0, 1, interpolated->height, interpolated->width, flip) - soff;
  rstep = flip_index (1, 0, interpolated->height, interpolated->width, flip) - flip_index (0, interpolated->width, interpolated->height, interpolated->width, flip);
  for (row=0; row < interpolated->height; row++, soff += rstep) {
	for (col=0; col < interpolated->width; col++, soff += cstep) {
	  if (output_bps==8)
		FORC3 ppm[col*3+c] =image[soff][c] >> 8;
	  else {
		FORC3 ppm2[col*3+c] = image[soff][c];
	  }
	}
/* 	memcpy(image[row*col], ppm, sizeof ppm); */
	if (output_bps == 16 && htons(0x55aa) != 0x55aa)
	  swab (ppm2, ppm2, interpolated->width*6);
	fwrite (ppm, 3*output_bps/8, interpolated->width, ofp);
  }
  free (ppm);
  fclose(ofp);
  free(ofname);

}


int x3f_foveon_fixed (CAMF *camf, void *ptr, int size, const char *name)
{
  void *dp;
  unsigned dim[3];

  if (!name) return 0;
  dp = X3F_foveon_camf_matrix (camf, dim, name);
  if (!dp) return 0;
  memcpy (ptr, dp, size*4);
  free (dp);
  return 1;
}

void x3f_trueII_interpolate(INTERPOLATED_IMG *i_img, X3F *x3f)
{
  static const short hood[] = { -1,-1, -1,0, -1,1, 0,-1, 0,1, 1,-1, 1,0, 1,1 };
  IMA *ima=(IMA *)x3f->raw->datas;
  DIR_ENTRY *section;
  CAMF *camf;
  ushort (*img)[4], *ipix, iipix[3], ival;
  double (*image)[4], *pix, val, fpix[3], xyz_pix[3];;
  int width=(int)i_img->width, height=(int)i_img->height;
  short  prev[3], *curve[8], (*shrink)[3];
  float cfilt=0, ddft[3][3][2], ppm[3][3][3], cp1[2][3][3], cp2[3][3], wbxy[2], x, y;
  float cam_xyz[3][3], correct[3][3]={{1,0,0},{0,1,0},{0,0,1}}, last[3][3], last2[3][3], trans[3][3];
  float chroma_dq[3]={4,4,4}, color_dq[3], diag[3][3], div[3], tempgainfact[3];
  float (*black)[3], (*sgain)[3], (*sgrow)[3];
  float fsum[3], frow, num, pgain[4];;
  int row, col, c, i, j, diff, sgx, irow, sum, min, max, limit;
  int dscr[2][2], dstb[4], (*smrow[7])[3], total[4];
  int work[3][3], smlast, smred, smred_p=0, dev[3];
  int satlev[3], keep[4], active[4], version[2];
  unsigned dim[3], *badpix, rgb_max[3], rgb_min[3], rgb_maxi[3], rgb_mini[3];
  double dsum=0, trsum[3];
  char str[128];
  int dp1=0;
  const char* cp;

  section=X3F_get_section(x3f, X3F_CAMF);
  camf=section->datas;

  x3f_foveon_fixed(camf, version, 2, "ContentVersionNumber");
  printf("%d\n", version[1]);
  if (version[1]==70){ /* make DP1 correctly processed */
	dp1=1;
  }
  printf ("TRUEII interpolation...\n");
/*   free(i_img->img); */
  img=(ushort (*)[4])i_img->img;


  x3f_foveon_fixed (camf, dscr, 4, "DarkShieldColRange");
  printf("dscr: %d %d %d %d\n", dscr[0][0], dscr[0][1], dscr[1][0],
  dscr[1][1]);


  x3f_foveon_fixed (camf, keep, 4, "KeepImageArea");
  printf("keep: %d %d %d\n", keep[0], keep[1], keep[2], keep[3]);
  x3f_foveon_fixed (camf,active, 4, "ActiveImageArea");
  printf("active: %d %d %d\n", active[0], active[1], active[2], active[3]);
  if (X3F_foveon_camf_param (camf,"IncludeBlocks", "ChromaDQ"))
	x3f_foveon_fixed (camf,chroma_dq, 3, "ChromaDQ");
  printf("chroma_dq: %f %f %f\n", chroma_dq[0], chroma_dq[1], chroma_dq[2]);
  x3f_foveon_fixed (camf,color_dq, 3,
	X3F_foveon_camf_param (camf,"IncludeBlocks", "ColorDQ") ?
		"ColorDQ" : "ColorDQCamRGB");
  printf("color_dq: %f %f %f\n", color_dq[0], color_dq[1], color_dq[2]);
  if (X3F_foveon_camf_param (camf,"IncludeBlocks", "ColumnFilter"))
		 x3f_foveon_fixed (camf,&cfilt, 1, "ColumnFilter");
  printf("cfilt: %f\n", cfilt);

  memset (ddft, 0, sizeof ddft);
  if (!X3F_foveon_camf_param (camf,"IncludeBlocks", "DarkDrift")
	 || !x3f_foveon_fixed (camf,ddft[1][0], 12, "DarkDrift"))
    for (i=0; i < 2; i++) {
      x3f_foveon_fixed (camf,dstb, 4, i ? "DarkShieldBottom":"DarkShieldTop");
      for (row = dstb[1]; row <= dstb[3]; row++)
	for (col = dstb[0]; col <= dstb[2]; col++)
	  FORC3 ddft[i+1][c][1] += (short) img[row*width+col][c];
      FORC3 ddft[i+1][c][1] /= (dstb[3]-dstb[1]+1) * (dstb[2]-dstb[0]+1);
    }
  printf("ddft: %f %f %f\n", ddft[1][0][1], ddft[1][0][2], ddft[1][1][1]);

  memset (last2, 0, sizeof last2);
  if (dp1){
  if (X3F_foveon_camf_param (camf,"IncludeBlocks", "DP1_WhiteBalanceColorCorrections"))
	x3f_foveon_fixed (camf,cam_xyz, 9,
				  X3F_foveon_camf_param (camf,"DP1_WhiteBalanceColorCorrections", (char *)x3f->header->whiteBalanceString));
  printf("cam_xyz: %f %f %f\n", cam_xyz[0][0], cam_xyz[0][1],
  cam_xyz[0][2]);
  if (X3F_foveon_camf_param (camf,"IncludeBlocks", "ColorModeCompensations"))
	x3f_foveon_fixed (camf,correct, 9,
				  X3F_foveon_camf_param (camf,"ColorModeCompensations", X3F_foveon_get_property(x3f->property, "CM_DESC")));
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
	  last2[i][j] = cam_xyz[i][j];
  } else {
  if (X3F_foveon_camf_param (camf,"IncludeBlocks", "WhiteBalanceColorCorrections"))
	x3f_foveon_fixed (camf,cam_xyz, 9,
				  X3F_foveon_camf_param (camf,"WhiteBalanceColorCorrections", (char *)x3f->header->whiteBalanceString));
  printf("cam_xyz: %f %f %f\n", cam_xyz[0][0], cam_xyz[0][1], cam_xyz[0][2]);
  
  if (X3F_foveon_camf_param (camf,"IncludeBlocks", "ColorModeCompensations"))
	x3f_foveon_fixed (camf,correct, 9,
				  X3F_foveon_camf_param (camf,"ColorModeCompensations", X3F_foveon_get_property(x3f->property, "CM_DESC")));
  printf("correct: %f %f %f\n", correct[0][0], correct[0][1], correct[0][2]);

  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
	  FORC3 last2[i][j] += cam_xyz[i][c]* correct[c][j];
  }

  #define LAST(x,y) last2[(i+x)%3][(c+y)%3]
  for (i=0; i < 3; i++)
    FORC3 diag[c][i] = LAST(1,1)*LAST(2,2) - LAST(1,2)*LAST(2,1);
  #undef LAST
  FORC3 div[c] = diag[c][0]*0.3127 + diag[c][1]*0.329 + diag[c][2]*0.3583;
  FORC3 printf("div[%d]=%f\n", c, div[c]);

/*   sprintf(str, "%sWBGains", (char
	 *)x3f->header->whiteBalanceString); */
  if (dp1){
  if (X3F_foveon_camf_param (camf,"IncludeBlocks", "DP1_WhiteBalanceGains"))
	x3f_foveon_fixed (camf,div, 3, X3F_foveon_camf_param (camf,"DP1_WhiteBalanceGains", (char *)x3f->header->whiteBalanceString));
  } else {
  if (X3F_foveon_camf_param (camf,"IncludeBlocks", "WhiteBalanceGains"))
	x3f_foveon_fixed (camf,div, 3, X3F_foveon_camf_param (camf,"WhiteBalanceGains", (char *)x3f->header->whiteBalanceString));
  }
  if (X3F_foveon_camf_param(camf,"IncludeBlocks", "TempGainFact")){
	x3f_foveon_fixed(camf,tempgainfact, 3, "TempGainFact");
	FORC3 div[c]*=tempgainfact[c];
  }
  if (X3F_foveon_camf_param(camf,"IncludeBlocks",
							"SensorAdjustmentGainFact")){
	x3f_foveon_fixed(camf,tempgainfact, 3, "SensorAdjustmentGainFact");
	FORC3 div[c]*=tempgainfact[c];
  }
  if (X3F_foveon_camf_param(camf,"IncludeBlocks", "CorrectColorGain_BR")){
	x3f_foveon_fixed(camf,tempgainfact, 3, "CorrectColorGain_BR");
	FORC3 div[c]*=tempgainfact[c];
  }
  if (X3F_foveon_camf_param(camf,"IncludeBlocks", "CorrectColorGain_GR")){
	x3f_foveon_fixed(camf,tempgainfact, 3, "CorrectColorGain_GR");
	FORC3 div[c]*=tempgainfact[c];
  }
  if (X3F_foveon_camf_param(camf,"IncludeBlocks", "CorrectColorGain_RR")){
	x3f_foveon_fixed(camf,tempgainfact, 3, "CorrectColorGain_RR");
	FORC3 div[c]*=tempgainfact[c];
  }
  if (X3F_foveon_camf_param(camf,"IncludeBlocks", "despAdjust")){
	x3f_foveon_fixed(camf,tempgainfact, 3, "DespAdjust");
	FORC3 div[c]*=tempgainfact[c];
  }
  num = 0;
  FORC3 if (num < div[c]) num = div[c];
  FORC3 div[c] /= num;
  FORC3 printf("div[%d]=%f\n", c, div[c]);

/*   FORC3 color_dq[c]/=16; */
  x3f_foveon_make_curves (curve, color_dq, div, cfilt);
  FORC3 chroma_dq[c] /= 3;
  x3f_foveon_make_curves (curve+3, chroma_dq, div, cfilt);
  FORC3 dsum += chroma_dq[c] / div[c];
  curve[6] = x3f_foveon_make_curve (dsum, dsum, cfilt);
  curve[7] = x3f_foveon_make_curve (dsum*2, dsum*2, cfilt);


  if ((badpix = (unsigned int *) X3F_foveon_camf_matrix (camf, dim, "BadPixelsF20"))) {
	for (i=0; i < (int)dim[0]*(int)dim[1]; i+=3) {
	  col = (badpix[i+1]) -keep[0];
	  row = (badpix[i] ) -keep[1];
	  if ((row-1) > height-3 || (col-1) > width-3)
		continue;
	  memset (fsum, 0, sizeof fsum);
	  for (sum=j=0; j < 8; j++){
		/* 	if (badpix[i] & (1 << j)) { */
		FORC3 fsum[c] += (short)
		  img[(row+hood[j*2])*width+col+hood[j*2+1]][c];
		sum++;
	  }
	  if (sum) FORC3 img[row*width+col][c] = fsum[c]/sum;
	}
	free (badpix);
  } else   if ((badpix = (unsigned int *) X3F_foveon_camf_matrix (camf, dim, "BadPixels"))) {
	for (i=0; i < (int)dim[0]; i++) {
	  col = (badpix[i] >> 8 & 0xfff) - keep[0];
	  row = (badpix[i] >> 20       ) - keep[1];
	  if ((row-1) > height-3 ||(col-1) > width-3)
		continue;
	  memset (fsum, 0, sizeof fsum);
	  for (sum=j=0; j < 8; j++)
		if (badpix[i] & (1 << j)) {
		  FORC3 fsum[c] += (short)
			img[(row+hood[j*2])*width+col+hood[j*2+1]][c];
		  sum++;
		}
	  if (sum) FORC3 img[row*width+col][c] = fsum[c]/sum;
	}
	free (badpix);
  }

  /* Array for 5x5 Gaussian averaging of red values */
  smrow[6] = (int (*)[3]) calloc (width*5, sizeof **smrow);
  /*   merror (smrow[6], "x3f_foveon_interpolate()"); */
  for (i=0; i < 5; i++)
	smrow[i] = smrow[6] + i*width;

  /* Sharpen the reds against these Gaussian averages */
  for (smlast=-1, row=2; row < height-2; row++) {
	while (smlast < row+2) {
	  for (i=0; i < 6; i++)
		smrow[(i+5) % 6] = smrow[i];
	  ipix = img[++smlast*width+2];
	  for (col=2; col < width-2; col++) {
		smrow[4][col][0] =
		  (ipix[0]*6 + (ipix[-4]+ipix[4])*4 + ipix[-8]+ipix[8] + 8) >> 4;
		ipix += 4;
	  }
	}
	ipix = img[row*width+2];
	for (col=2; col < width-2; col++) {
	  smred = ( 6 *  smrow[2][col][0]
				+ 4 * (smrow[1][col][0] + smrow[3][col][0])
				+      smrow[0][col][0] + smrow[4][col][0] + 8 ) >> 4;
	  if (col == 2)
		smred_p = smred;
	  i = ipix[0] + ((ipix[0] - ((smred*7 + smred_p) >> 3)) >> 3);
	  if (i > 32000) i = 32000;
	  ipix[0] = i;
	  smred_p = smred;
	  ipix += 4;
	}
  }


  /*   Adjust the brighter pixels for better linearity */
  if (X3F_foveon_camf_param (camf,"IncludeBlocks", "SaturationLevel")) {
	x3f_foveon_fixed (camf,satlev, 3, "SaturationLevel");
	min = 0xffff;
	FORC3 {
	  i = satlev[c] / div[c];
	  if (min > i) min = i;
	}
	limit = min * 9 >> 4;
	for (ipix=img[0]; ipix < img[height*width]; ipix+=4) {
	  if (ipix[0] <= limit || ipix[1] <= limit || ipix[2] <= limit)
		continue;
	  min = max = ipix[0];
	  for (c=1; c < 3; c++) {
		if (min > ipix[c]) min = ipix[c];
		if (max < ipix[c]) max = ipix[c];
	  }
	  if (min >= limit*2) {
		ipix[0] = ipix[1] = ipix[2] = max;
	  } else {
		i = 0x4000 - ((min - limit) << 14) / limit;
		i = 0x4000 - (i*i >> 14);
		i = i*i >> 14;
		FORC3 ipix[c] += (max - ipix[c]) * i >> 14;
	  }
	}
  }

  /*    Because photons that miss one detector often hit another, */
  /*    the sum R+G+B is much less noisy than the individual colors. */
  /*    So smooth the hues without smoothing the total. */
  for (smlast=-1, row=2; row < height-2; row++) {
	while (smlast < row+2) {
	  for (i=0; i < 6; i++)
		smrow[(i+5) % 6] = smrow[i];
	  ipix = img[++smlast*width+2];
	  for (col=2; col < width-2; col++) {
		FORC3 smrow[4][col][c] = (ipix[c-4]+2*ipix[c]+ipix[c+4]+2) >> 2;
		ipix += 4;
	  }
	}
	ipix = img[row*width+2];
	for (col=2; col < width-2; col++) {
	  FORC3 dev[c] = -x3f_foveon_apply_curve (curve[7], ipix[c] -
											  ((smrow[1][col][c] + 2*smrow[2][col][c] + smrow[3][col][c]) >> 2));
	  sum = (dev[0] + dev[1] + dev[2]) >> 3;
	  FORC3 ipix[c] += dev[c] - sum;
	  ipix += 4;
	}
  }
  for (smlast=-1, row=2; row < height-2; row++) {
	while (smlast < row+2) {
	  for (i=0; i < 6; i++)
		smrow[(i+5) % 6] = smrow[i];
	  ipix = img[++smlast*width+2];
	  for (col=2; col < width-2; col++) {
		FORC3 smrow[4][col][c] =
		  (ipix[c-8]+ipix[c-4]+ipix[c]+ipix[c+4]+ipix[c+8]+2) >> 2;
		ipix += 4;
	  }
	}
	ipix = img[row*width+2];
	for (col=2; col < width-2; col++) {
	  for (total[3]=375, sum=60, c=0; c < 3; c++) {
		for (total[c]=i=0; i < 5; i++)
		  total[c] += smrow[i][col][c];
		total[3] += total[c];
		sum += ipix[c];
	  }
	  if (sum < 0) sum = 0;
	  j = total[3] > 375 ? (sum << 16) / total[3] : sum * 174;
	  FORC3 ipix[c] += x3f_foveon_apply_curve (curve[6],
											   ((j*total[c] + 0x8000) >> 16) - ipix[c]);
	  ipix += 4;
	}
  }

/*   if (dp1)  x3f_foveon_fixed (camf, cp2, 9, "DP1_CP2_Matrix"); */
/*   else  x3f_foveon_fixed (camf, cp2, 9, "CP2_Matrix"); */
/*   for (i=0; i<3;i++) */
/* 	FORC3 i_img->rgb_cam[c][i]=cp2[c][i]; */
  min=0xffff;
  max=0;
  for (ipix=img[0]; ipix < img[height*width]; ipix+=4) {
    FORC3 ipix[c] -= x3f_foveon_apply_curve (curve[c], ipix[c]);
    sum = (ipix[0]+ipix[1]+ipix[1]+ipix[2]) >> 2;
    FORC3 ipix[c] -= x3f_foveon_apply_curve (curve[c], ipix[c]-sum);
    FORC3 {
      for (dsum=i=0; i < 3; i++)
	dsum += i_img->rgb_cam[c][i] * last2[c][i] * div[i] * ipix[i];
      if (dsum < 0)  dsum = 0;
      if (dsum > 24000) dsum = 24000;
      iipix[c] = dsum + 0.5;
    }
    FORC3 {
	  ipix[c] = iipix[c];
	  if (ipix[c]<min)min=ipix[c];
	  if (ipix[c]>max)max=ipix[c];
	}
  }

/*   for (row=0; row<height; row++){ */
/* 	ipix = img[row*width]; */
/* 	for (col=0; col<width; col++){ */
/* 	  for (c=0; c < 3; c++) iipix[c]=ipix[c]; */
/* 	  for (c=0; c < 3; c++) { */
/* 		float tmp; */

/* 		tmp=i_img->rgb_cam[c][0]*last2[c][0]* div[0]*iipix[0]+i_img->rgb_cam[c][1]*last2[c][1]*div[1]*iipix[1]+i_img->rgb_cam[c][2]*last2[c][2]*div[2]*iipix[2]; */
/* 		if (tmp<0) ipix[c]=0; */
/* /\* 		else if (tmp >3500) ipix[c]=3500; *\/ */
/* 		else ipix[c]=floor(tmp); */
/* 		if (ipix[c]<min)min=ipix[c]; */
/* 		if (ipix[c]>max)max=ipix[c]; */
/* 	  } */
/* 	  ipix+=4; */
/* 	} */

/*   } */
/*   printf("min: %d\t max: %d\n", min, max); */


/* 	/\*   Smooth the image bottom-to-top and save at 1/4 scale *\/ */
/* 	shrink = (short (*)[3]) calloc ((height/4), (width/4)*sizeof *shrink); */
/* 	/\*   merror (shrink, "x3f_foveon_interpolate()"); *\/ */
/* 	for (row = height/4; row--; ) */
/* 	  for (col=0; col < width/4; col++) { */
/* 		iipix[0] = iipix[1] = iipix[2] = 0; */
/* 		for (i=0; i < 4; i++) */
/* 		  for (j=0; j < 4; j++) */
/* 			FORC3 iipix[c] += img[(row*4+i)*width+col*4+j][c]; */
/* 		FORC3 */
/* 		  if (row+2 > height/4) */
/* 			shrink[row*(width/4)+col][c] = iipix[c] >> 4; */
/* 		  else */
/* 			shrink[row*(width/4)+col][c] = */
/* 			  (shrink[(row+1)*(width/4)+col][c]*1840 + iipix[c]*141 + 2048) >> 12; */
/* 	  } */

/* 	/\* From the 1/4-scale image, smooth right-to-left *\/ */
/* 	for (row=0; row < (height & ~3); row++) { */
/* 	  iipix[0] = iipix[1] = iipix[2] = 0; */
/* /\* 	  printf("Here %d\n", row); *\/ */
/* 	  if ((row & 3) == 0) */
/* 		for (col = width & ~3 ; col--; ) */
/* 		  FORC3 smrow[0][col][c] = iipix[c] = */
/* 			(shrink[(row/4)*(width/4)+col/4][c]*1485 + iipix[c]*6707 + 4096) >> 13; */

/* 	  /\* Then smooth left-to-right *\/ */
/* 	  iipix[0] = iipix[1] = iipix[2] = 0; */
/* 	  for (col=0; col < (width & ~3); col++) */
/* 		FORC3 smrow[1][col][c] = iipix[c] = */
/* 		  (smrow[0][col][c]*1485 + iipix[c]*6707 + 4096) >> 13; */

/* 	  /\* Smooth top-to-bottom *\/ */
/* 	  if (row == 0) */
/* 		memcpy (smrow[2], smrow[1], sizeof **smrow * width); */
/* 	  else */
/* 		for (col=0; col < (width & ~3); col++) */
/* 		  FORC3 smrow[2][col][c] = */
/* 			(smrow[2][col][c]*6707 + smrow[1][col][c]*1485 + 4096) >> 13; */

/* 	  /\* Adjust the chroma toward the smooth values *\/ */
/* 	  for (col=0; col < (width & ~3); col++) { */
/* 		for (i=j=30, c=0; c < 3; c++) { */
/* 		  i += smrow[2][col][c]; */
/* 		  j += img[row*width+col][c]; */
/* 		} */
/* 		j = (j << 16) / i; */
/* 		for (sum=c=0; c < 3; c++) { */
/* 		  iipix[c] = x3f_foveon_apply_curve (curve[c+3], */
/* 										((smrow[2][col][c] * j + 0x8000) >> 16) - img[row*width+col][c]); */
/* 		  sum += iipix[c]; */
/* 		} */
/* 		sum >>= 3; */
/* 		FORC3 { */
/* 		  i = img[row*width+col][c] + iipix[c] - sum; */
/* 		  if (i < 0) i = 0; */
/* 		  img[row*width+col][c] = i; */
/* 		} */
/* 	  } */
/* 	} */
/* 	free (shrink); */
	free (smrow[6]);


 for (i=0; i < 8; i++)
   free (curve[i]);

/*  for (ipix=img[0]; ipix<img[height*width]; ipix+=4){ */
/*    FORC3 { */
/* 	 ipix[c]=((ipix[c]-min)*65535/(max-min)); */
/*    } */
/*  } */

 /* adjust Saturation */
 /* adjust black */

  /* Trim off the black border */
  active[1] -= keep[1];
  active[3] -= 2;
  i = active[2] - active[0];
  for (row=0; row < active[3]-active[1]; row++)
    memcpy (img[row*i], img[(row+active[1])*width+active[0]],
	 i * sizeof *img);
  i_img->width = i;
  i_img->height = row;
/*   i_img->img=image; */
  printf("Done\n");
}

void spatial_gain(CAMF *camf, INTERPOLATED_IMG *i_img, float last[3][3], float div[3], float cfilt, char* wbstring){
  int16_t (*image)[4], *pix;
  int width=(int)i_img->width, height=(int)i_img->height;
  int work[3][3];
  short prev[3];
  int row, col, c, i, j, diff, sgx, irow, ipix[3] ;
  float (*black)[3], (*sgain)[3], (*sgrow)[3];
  float ddft[3][3][2], ppm[3][3][3], frow, fsum[3];
  float val;
  int dscr[2][2], dstb[4], total[4];
  int dim[3];

  printf("SpatialGain\n");
  image=(int16_t(*)[4])i_img->img;

  memset (ddft, 0, sizeof ddft);
  if (!X3F_foveon_camf_param (camf, "IncludeBlocks", "DarkDrift")
	  || !x3f_foveon_fixed (camf, ddft[1][0], 12, "DarkDrift"))
	for (i=0; i < 2; i++) {
	  x3f_foveon_fixed (camf, dstb, 4, i ? "DarkShieldBottom":"DarkShieldTop");
	  for (row = dstb[1]; row <= dstb[3]; row++)
		for (col = dstb[0]; col <= dstb[2]; col++)
		  FORC3 ddft[i+1][c][1] += (short) image[row*width+col][c];
	  FORC3 ddft[i+1][c][1] /= (dstb[3]-dstb[1]+1) * (dstb[2]-dstb[0]+1);
	}

  black = (float (*)[3]) calloc (height, sizeof *black);
  for (row=0; row < height; row++) {
	for (i=0; i < 6; i++)
	  ddft[0][0][i] = ddft[1][0][i] +
		row / (height-1.0) * (ddft[2][0][i] - ddft[1][0][i]);
	FORC3 black[row][c] =
	  ( x3f_foveon_avg (image[row*width]+c, dscr[0], cfilt) +
		x3f_foveon_avg (image[row*width]+c, dscr[1], cfilt) * 3
		- ddft[0][c][0] ) / 4 - ddft[0][c][1];
  }
  memcpy (black, black+8, sizeof *black*8);
  memcpy (black+height-11, black+height-22, 11*sizeof *black);
  memcpy (last, black, sizeof last);

  for (row=1; row < height-1; row++) {
	FORC3 if (last[1][c] > last[0][c]) {
	  if (last[1][c] > last[2][c])
		black[row][c] = (last[0][c] > last[2][c]) ? last[0][c]:last[2][c];
	} else
	  if (last[1][c] < last[2][c])
		black[row][c] = (last[0][c] < last[2][c]) ? last[0][c]:last[2][c];
	memmove (last, last+1, 2*sizeof last[0]);
	memcpy (last[2], black[row+1], sizeof last[2]);
  }
  FORC3 black[row][c] = (last[0][c] + last[1][c])/2;
  FORC3 black[0][c] = (black[1][c] + black[3][c])/2;

  val = 1 - exp(-1/24.0);
  memcpy (fsum, black, sizeof fsum);
  for (row=1; row < height; row++)
	FORC3 fsum[c] += black[row][c] =
	  (black[row][c] - black[row-1][c])*val + black[row-1][c];
  memcpy (last[0], black[height-1], sizeof last[0]);
  FORC3 fsum[c] /= height;
  for (row = height; row--; )
	FORC3 last[0][c] = black[row][c] =
	  (black[row][c] - fsum[c] - last[0][c])*val + last[0][c];

  memset (total, 0, sizeof total);
  for (row=2; row < height; row+=4)
	for (col=2; col < width; col+=4) {
	  FORC3 total[c] += (short) image[row*width+col][c];
	  total[3]++;
	}
  for (row=0; row < height; row++)
	FORC3 black[row][c] += fsum[c]/2 + total[c]/(total[3]*100.0);

  if (X3F_foveon_camf_param (camf, "IncludeBlocks", "SpatialGainTables"))
	sgain=(float (*)[3])X3F_foveon_camf_matrix (camf, dim, X3F_foveon_camf_param (camf, "SpatialGainTables", wbstring));
  else
	sgain = (float (*)[3]) X3F_foveon_camf_matrix (camf, dim, "SpatialGain");
  if (!sgain) return;
  sgrow = (float (*)[3]) calloc (dim[1], sizeof *sgrow);
  sgx = (width + dim[1]-2) / (dim[1]-1);

  for (row=0; row < height; row++) {
	for (i=0; i < 6; i++)
	  ddft[0][0][i] = ddft[1][0][i] +
		row / (height-1.0) * (ddft[2][0][i] - ddft[1][0][i]);
	pix = image[row*width];
	memcpy (prev, pix, sizeof prev);
	frow = row / (height-1.0) * (dim[2]-1);
	if ((irow = frow) == dim[2]-1) irow--;
	frow -= irow;
	for (i=0; i < dim[1]; i++)
	  FORC3 sgrow[i][c] = sgain[ irow   *dim[1]+i][c] * (1-frow) +
		sgain[(irow+1)*dim[1]+i][c] *    frow;
	for (col=0; col < width; col++) {
	  FORC3 {
		diff = pix[c] - prev[c];
		prev[c] = pix[c];
		ipix[c] = pix[c] + floor ((diff + (diff*diff >> 14)) * cfilt
								  - ddft[0][c][1] - ddft[0][c][0] * ((float) col/width - 0.5)
								  - black[row][c] );
	  }
	  FORC3 {
		work[0][c] = ipix[c] * ipix[c] >> 14;
		work[2][c] = ipix[c] * work[0][c] >> 14;
		work[1][2-c] = ipix[(c+1) % 3] * ipix[(c+2) % 3] >> 14;
	  }
	  FORC3 {
		for (val=i=0; i < 3; i++)
		  for (  j=0; j < 3; j++)
			val += ppm[c][i][j] * work[i][j];
		ipix[c] = floor ((ipix[c] + floor(val)) *
						 ( sgrow[col/sgx  ][c] * (sgx - col%sgx) +
						   sgrow[col/sgx+1][c] * (col%sgx) ) / sgx / div[c]);
		if (ipix[c] > 32000) ipix[c] = 32000;
		pix[c] = ipix[c];
	  }
	  pix += 4;
	}
  }
  free (sgrow);
  free (black);
  free (sgain);
  return;
}


void x3f_interpolate(INTERPOLATED_IMG *i_img, X3F *x3f) {
  static const short hood[] = { -1,-1, -1,0, -1,1, 0,-1, 0,1, 1,-1, 1,0, 1,1 };
  IMA *ima=(IMA *)x3f->raw->datas;
  int16_t (*image)[4];
  int width=(int)i_img->width, height=(int)i_img->height;
  DIR_ENTRY *section;
  CAMF *camf;
/*   uint16_t *pix; */
/*   short  prev[3], *curve[8], (*shrink)[3]; */
/*   float cfilt=0, ddft[3][3][2], ppm[3][3][3]; */
/*   float cam_xyz[3][3], correct[3][3]={{1,0,0},{0,1,0},{0,0,1}}; */
/*   float last[3][3], trans[3][3]; */
/*   float chroma_dq[3]={10,10,10}; */
/*   float color_dq[3], diag[3][3], div[3]; */
/*   float wbgain[3]; */
/*   float (*black)[3], (*sgrow)[3], (*sgain)[3]; */
/*   float fsum[3], val, frow, num; */
/*   float *mat, rgb_cam[3][3]; */
  int row, col;
  int valueCount;
  int16_t *pix;
  short  prev[3], *curve[8], (*shrink)[3];
  float cfilt=0, ddft[3][3][2], ppm[3][3][3];
  float cam_xyz[3][3], correct[3][3]={{1,0,0},{0,1,0},{0,0,1}}, last[3][3], trans[3][3];
  float chroma_dq[3], color_dq[3], diag[3][3], div[3];
  float (*black)[3], (*sgain)[3], (*sgrow)[3];
  float fsum[3], val, frow, num;
  int  c, i, j, diff, sgx, irow, sum, min, max, limit;
  int dscr[2][2], dstb[4], (*smrow[7])[3], total[4], ipix[3];
  int work[3][3], smlast, smred, smred_p=0, dev[3];
  int satlev[3], keep[4], active[4];
  unsigned dim[3], *badpix;
  double dsum=0, trsum[3];
  char str[128];
  const char* cp;

 /*  image=(int16_t (*)[4])i_img->img; */
  printf ("Foveon interpolation...\n");

  section=X3F_get_section(x3f, X3F_CAMF);
  camf=(CAMF *)section->datas;

  x3f_foveon_fixed (camf, dscr, 4, "DarkShieldColRange");
  memset(ppm, 0, sizeof ppm);
  if (X3F_foveon_camf_param (camf, "IncludeBlocks", "PostPolyMatrix"))
	x3f_foveon_fixed (camf, ppm[0][0], 27, "PostPolyMatrix");
  x3f_foveon_fixed (camf, satlev, 3, "SaturationLevel");
  x3f_foveon_fixed (camf, keep, 4, "KeepImageArea");
  x3f_foveon_fixed (camf, active, 4, "ActiveImageArea");
  x3f_foveon_fixed (camf, chroma_dq, 3, "ChromaDQ");
  x3f_foveon_fixed (camf, color_dq, 3,
				X3F_foveon_camf_param (camf, "IncludeBlocks", "ColorDQ") ?
				"ColorDQ" : "ColorDQCamRGB");
  if (X3F_foveon_camf_param (camf, "IncludeBlocks", "ColumnFilter"))
	x3f_foveon_fixed (camf, &cfilt, 1, "ColumnFilter");


  if (!(cp = X3F_foveon_camf_param (camf, "WhiteBalanceIlluminants", (char *)x3f->header->whiteBalanceString)))
	{ fprintf (stderr,_("Invalid white balance \"%s\"\n"), (char *)x3f->header->whiteBalanceString);
	  return; }
  x3f_foveon_fixed (camf, cam_xyz, 9, cp);
  if (X3F_foveon_camf_param (camf, "IncludeBlocks", "WhiteBalanceCorrections"))
	x3f_foveon_fixed (camf, correct, 9,
				  X3F_foveon_camf_param (camf, "WhiteBalanceCorrections", (char *)x3f->header->whiteBalanceString));
/*   else  */
/* 	x3f_foveon_fixed (camf, correct, 9, "CP2_Matrix"); */

  memset (last, 0, sizeof last);
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      FORC3 last[i][j] += correct[i][c] * cam_xyz[c][j];

#define LAST(x,y) last[(i+x)%3][(c+y)%3]
  for (i=0; i < 3; i++)
    FORC3 diag[c][i] = LAST(1,1)*LAST(2,2) - LAST(1,2)*LAST(2,1);
#undef LAST
  FORC3 div[c] = diag[c][0]*0.3127 + diag[c][1]*0.329 + diag[c][2]*0.3583;

  sprintf (str, "%sRGBNeutral", (char *)x3f->header->whiteBalanceString);
  if (X3F_foveon_camf_param (camf, "IncludeBlocks", str))
    x3f_foveon_fixed (camf, div, 3, str);
  else {
/* 	sprintf(str, "%sWBGains", (char *)x3f->header->whiteBalanceString); */
	if (X3F_foveon_camf_param (camf, "IncludeBlocks", "WhiteBalanceGains"))
	  x3f_foveon_fixed (camf, div, 3, X3F_foveon_camf_param (camf, "WhiteBalanceGains", (char *)x3f->header->whiteBalanceString));
  }

  num = 0;
  FORC3 if (num < div[c]) num = div[c];
  FORC3 div[c] /= num;

  memset (trans, 0, sizeof trans);
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      FORC3 trans[i][j] += i_img->rgb_cam[i][c] * last[c][j] * div[j];
  FORC3 trsum[c] = trans[c][0] + trans[c][1] + trans[c][2];
  dsum = (6*trsum[0] + 11*trsum[1] + 3*trsum[2]) / 20;
  for (i=0; i < 3; i++)
    FORC3 last[i][c] = trans[i][c] * dsum / trsum[i];
  memset (trans, 0, sizeof trans);
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      FORC3 trans[i][j] += (i==c ? 32 : -1) * last[c][j] / 30;

  x3f_foveon_make_curves (curve, color_dq, div, cfilt);
  FORC3 chroma_dq[c] /= 3;
  x3f_foveon_make_curves (curve+3, chroma_dq, div, cfilt);
  FORC3 dsum += chroma_dq[c] / div[c];
  curve[6] = x3f_foveon_make_curve (dsum, dsum, cfilt);
  curve[7] = x3f_foveon_make_curve (dsum*2, dsum*2, cfilt);

/*   spatial_gain(camf, i_img, last, div, cfilt, (char *)x3f->header->whiteBalanceString); */

  image=(int16_t(*)[4])calloc(width*height*4, sizeof (*image)[4]);
  memcpy(image,i_img->img,width*height*8);
  free(i_img->img);

  memset (ddft, 0, sizeof ddft);
  if (!X3F_foveon_camf_param (camf, "IncludeBlocks", "DarkDrift")
	  || !x3f_foveon_fixed (camf, ddft[1][0], 12, "DarkDrift"))
	for (i=0; i < 2; i++) {
	  x3f_foveon_fixed (camf, dstb, 4, i ? "DarkShieldBottom":"DarkShieldTop");
	  for (row = dstb[1]; row <= dstb[3]; row++)
		for (col = dstb[0]; col <= dstb[2]; col++)
		  FORC3 ddft[i+1][c][1] += (short) image[row*width+col][c];
	  FORC3 ddft[i+1][c][1] /= (dstb[3]-dstb[1]+1) * (dstb[2]-dstb[0]+1);
	}

  black = (float (*)[3]) calloc (height, sizeof *black);
  for (row=0; row < height; row++) {
	for (i=0; i < 6; i++)
	  ddft[0][0][i] = ddft[1][0][i] +
		row / (height-1.0) * (ddft[2][0][i] - ddft[1][0][i]);
	FORC3 black[row][c] =
	  ( x3f_foveon_avg (image[row*width]+c, dscr[0], cfilt) +
		x3f_foveon_avg (image[row*width]+c, dscr[1], cfilt) * 3
		- ddft[0][c][0] ) / 4 - ddft[0][c][1];
  }
  memcpy (black, black+8, sizeof *black*8);
  memcpy (black+height-11, black+height-22, 11*sizeof *black);
  memcpy (last, black, sizeof last);

  for (row=1; row < height-1; row++) {
	FORC3 if (last[1][c] > last[0][c]) {
	  if (last[1][c] > last[2][c])
		black[row][c] = (last[0][c] > last[2][c]) ? last[0][c]:last[2][c];
	} else
	  if (last[1][c] < last[2][c])
		black[row][c] = (last[0][c] < last[2][c]) ? last[0][c]:last[2][c];
	memmove (last, last+1, 2*sizeof last[0]);
	memcpy (last[2], black[row+1], sizeof last[2]);
  }
  FORC3 black[row][c] = (last[0][c] + last[1][c])/2;
  FORC3 black[0][c] = (black[1][c] + black[3][c])/2;

  val = 1 - exp(-1/24.0);
  memcpy (fsum, black, sizeof fsum);
  for (row=1; row < height; row++)
	FORC3 fsum[c] += black[row][c] =
	  (black[row][c] - black[row-1][c])*val + black[row-1][c];
  memcpy (last[0], black[height-1], sizeof last[0]);
  FORC3 fsum[c] /= height;
  for (row = height; row--; )
	FORC3 last[0][c] = black[row][c] =
	  (black[row][c] - fsum[c] - last[0][c])*val + last[0][c];

  memset (total, 0, sizeof total);
  for (row=2; row < height; row+=4)
	for (col=2; col < width; col+=4) {
	  FORC3 total[c] += (short) image[row*width+col][c];
	  total[3]++;
	}
  for (row=0; row < height; row++)
	FORC3 black[row][c] += fsum[c]/2 + total[c]/(total[3]*100.0);

  if (X3F_foveon_camf_param (camf, "IncludeBlocks", "SpatialGainTables"))
	sgain=(float (*)[3])X3F_foveon_camf_matrix (camf, dim, X3F_foveon_camf_param (camf, "SpatialGainTables", (char *)x3f->header->whiteBalanceString));
  else
	sgain = (float (*)[3]) X3F_foveon_camf_matrix (camf, dim, "SpatialGain");
  if (!sgain) return;
  sgrow = (float (*)[3]) calloc (dim[1], sizeof *sgrow);
  sgx = (width + dim[1]-2) / (dim[1]-1);

  for (row=0; row < height; row++) {
	for (i=0; i < 6; i++)
	  ddft[0][0][i] = ddft[1][0][i] +
		row / (height-1.0) * (ddft[2][0][i] - ddft[1][0][i]);
	pix = image[row*width];
	memcpy (prev, pix, sizeof prev);
	frow = row / (height-1.0) * (dim[2]-1);
	if ((irow = frow) == dim[2]-1) irow--;
	frow -= irow;
	for (i=0; i < dim[1]; i++)
	  FORC3 sgrow[i][c] = sgain[ irow   *dim[1]+i][c] * (1-frow) +
		sgain[(irow+1)*dim[1]+i][c] *    frow;
	for (col=0; col < width; col++) {
	  FORC3 {
		diff = pix[c] - prev[c];
		prev[c] = pix[c];
		ipix[c] = pix[c] + floor ((diff + (diff*diff >> 14)) * cfilt
								  - ddft[0][c][1] - ddft[0][c][0] * ((float) col/width - 0.5)
								  - black[row][c] );
	  }
	  FORC3 {
		work[0][c] = ipix[c] * ipix[c] >> 14;
		work[2][c] = ipix[c] * work[0][c] >> 14;
		work[1][2-c] = ipix[(c+1) % 3] * ipix[(c+2) % 3] >> 14;
	  }
	  FORC3 {
		for (val=i=0; i < 3; i++)
		  for (  j=0; j < 3; j++)
			val += ppm[c][i][j] * work[i][j];
		ipix[c] = floor ((ipix[c] + floor(val)) *
						 ( sgrow[col/sgx  ][c] * (sgx - col%sgx) +
						   sgrow[col/sgx+1][c] * (col%sgx) ) / sgx / div[c]);
		if (ipix[c] > 32000) ipix[c] = 32000;
		pix[c] = ipix[c];
	  }
	  pix += 4;
	}
  }
  free (sgrow);
  free (black);
  free (sgain);


  if ((badpix = (unsigned int *) X3F_foveon_camf_matrix (camf, dim, "BadPixels"))) {
    for (i=0; i < dim[0]; i++) {
      col = (badpix[i] >> 8 & 0xfff) - keep[0];
      row = (badpix[i] >> 20       ) - keep[1];
      if ((unsigned)(row-1) > height-3 || (unsigned)(col-1) > width-3)
		continue;
      memset (fsum, 0, sizeof fsum);
      for (sum=j=0; j < 8; j++)
		if (badpix[i] & (1 << j)) {
		  FORC3 fsum[c] += (short)
			image[(row+hood[j*2])*width+col+hood[j*2+1]][c];
		  sum++;
		}
      if (sum) FORC3 image[row*width+col][c] = fsum[c]/sum;
    }
    free (badpix);
  }

  /* Array for 5x5 Gaussian averaging of red values */
  smrow[6] = (int (*)[3]) calloc (width*5, sizeof **smrow);
  /*   merror (smrow[6], "x3f_foveon_interpolate()"); */
  for (i=0; i < 5; i++)
    smrow[i] = smrow[6] + i*width;

  /* Sharpen the reds against these Gaussian averages */
  for (smlast=-1, row=2; row < height-2; row++) {
    while (smlast < (int)row+2) {
      for (i=0; i < 6; i++)
		smrow[(i+5) % 6] = smrow[i];
      pix = image[++smlast*width+2];
      for (col=2; col < width-2; col++) {
		smrow[4][col][0] =
		  (pix[0]*6 + (pix[-4]+pix[4])*4 + pix[-8]+pix[8] + 8) >> 4;
		pix += 4;
      }
    }
    pix = image[row*width+2];
    for (col=2; col < width-2; col++) {
      smred = ( 6 *  smrow[2][col][0]
				+ 4 * (smrow[1][col][0] + smrow[3][col][0])
				+      smrow[0][col][0] + smrow[4][col][0] + 8 ) >> 4;
      if (col == 2)
		smred_p = smred;
      i = pix[0] + ((pix[0] - ((smred*7 + smred_p) >> 3)) >> 3);
      if (i > 32000) i = 32000;
      pix[0] = i;
      smred_p = smred;
      pix += 4;
    }
  }


  /*   Adjust the brighter pixels for better linearity */
  min = 0xffff;
  FORC3 {
	i = satlev[c] / div[c];
	if (min > i) min = i;
  }
  limit = min * 9 >> 4;
  for (pix=image[0]; pix < image[height*width]; pix+=4) {
	if (pix[0] <= limit || pix[1] <= limit || pix[2] <= limit)
	  continue;
	min = max = pix[0];
	for (c=1; c < 3; c++) {
	  if (min > pix[c]) min = pix[c];
	  if (max < pix[c]) max = pix[c];
	}
	if (min >= limit*2) {
	  pix[0] = pix[1] = pix[2] = max;
	} else {
	  i = 0x4000 - ((min - limit) << 14) / limit;
	  i = 0x4000 - (i*i >> 14);
	  i = i*i >> 14;
	  FORC3 pix[c] += (max - pix[c]) * i >> 14;
	}
  }

  /*
    Because photons that miss one detector often hit another,
    the sum R+G+B is much less noisy than the individual colors.
    So smooth the hues without smoothing the total.
  */
  for (smlast=-1, row=2; row < height-2; row++) {
    while (smlast < (int)row+2) {
      for (i=0; i < 6; i++)
		smrow[(i+5) % 6] = smrow[i];
      pix = image[++smlast*width+2];
      for (col=2; col < width-2; col++) {
		FORC3 smrow[4][col][c] = (pix[c-4]+2*pix[c]+pix[c+4]+2) >> 2;
		pix += 4;
      }
    }
    pix = image[row*width+2];
    for (col=2; col < width-2; col++) {
      FORC3 dev[c] = -x3f_foveon_apply_curve (curve[7], pix[c] -
										  ((smrow[1][col][c] + 2*smrow[2][col][c] + smrow[3][col][c]) >> 2));
      sum = (dev[0] + dev[1] + dev[2]) >> 3;
      FORC3 pix[c] += dev[c] - sum;
      pix += 4;
    }
  }
  for (smlast=-1, row=2; row < height-2; row++) {
    while (smlast < (int)row+2) {
      for (i=0; i < 6; i++)
		smrow[(i+5) % 6] = smrow[i];
      pix = image[++smlast*width+2];
      for (col=2; col < width-2; col++) {
		FORC3 smrow[4][col][c] =
		  (pix[c-8]+pix[c-4]+pix[c]+pix[c+4]+pix[c+8]+2) >> 2;
		pix += 4;
      }
    }
    pix = image[row*width+2];
    for (col=2; col < width-2; col++) {
      for (total[3]=375, sum=60, c=0; c < 3; c++) {
		for (total[c]=i=0; i < 5; i++)
		  total[c] += smrow[i][col][c];
		total[3] += total[c];
		sum += pix[c];
      }
      if (sum < 0) sum = 0;
      j = total[3] > 375 ? (sum << 16) / total[3] : sum * 174;
      FORC3 pix[c] += x3f_foveon_apply_curve (curve[6],
										  ((j*total[c] + 0x8000) >> 16) - pix[c]);
      pix += 4;
    }
  }

  /*   Transform the image to a different colorspace */
  for (pix=image[0]; pix < image[height*width]; pix+=4) {
    FORC3 pix[c] -= x3f_foveon_apply_curve (curve[c], pix[c]);
    sum = (pix[0]+pix[1]+pix[1]+pix[2]) >> 2;
    FORC3 pix[c] -= x3f_foveon_apply_curve (curve[c], pix[c]-sum);
    FORC3 {
      for (dsum=i=0; i < 3; i++)
		dsum += trans[c][i] * pix[i];
      if (dsum < 0)  dsum = 0;
      if (dsum > 24000) dsum = 24000;
      ipix[c] = dsum + 0.5;
    }
    FORC3 pix[c] = ipix[c];
  }

  /* Smooth the image bottom-to-top and save at 1/4 scale */
  shrink = (short (*)[3]) calloc ((height/4), (width/4)*sizeof *shrink);
  /*   merror (shrink, "x3f_foveon_interpolate()"); */
  for (row = height/4; row--; )
    for (col=0; col < width/4; col++) {
      ipix[0] = ipix[1] = ipix[2] = 0;
      for (i=0; i < 4; i++)
		for (j=0; j < 4; j++)
		  FORC3 ipix[c] += image[(row*4+i)*width+col*4+j][c];
      FORC3
		if (row+2 > height/4)
		  shrink[row*(width/4)+col][c] = ipix[c] >> 4;
		else
		  shrink[row*(width/4)+col][c] =
			(shrink[(row+1)*(width/4)+col][c]*1840 + ipix[c]*141 + 2048) >> 12;
    }
  /* From the 1/4-scale image, smooth right-to-left */
  for (row=0; row < (height & ~3); row++) {
    ipix[0] = ipix[1] = ipix[2] = 0;
    if ((row & 3) == 0)
      for (col = width & ~3 ; col--; )
		FORC3 smrow[0][col][c] = ipix[c] =
		  (shrink[(row/4)*(width/4)+col/4][c]*1485 + ipix[c]*6707 + 4096) >> 13;

    /* Then smooth left-to-right */
    ipix[0] = ipix[1] = ipix[2] = 0;
    for (col=0; col < (width & ~3); col++)
      FORC3 smrow[1][col][c] = ipix[c] =
		(smrow[0][col][c]*1485 + ipix[c]*6707 + 4096) >> 13;

    /* Smooth top-to-bottom */
    if (row == 0)
      memcpy (smrow[2], smrow[1], sizeof **smrow * width);
    else
      for (col=0; col < (width & ~3); col++)
		FORC3 smrow[2][col][c] =
		  (smrow[2][col][c]*6707 + smrow[1][col][c]*1485 + 4096) >> 13;

    /* Adjust the chroma toward the smooth values */
    for (col=0; col < (width & ~3); col++) {
      for (i=j=30, c=0; c < 3; c++) {
		i += smrow[2][col][c];
		j += image[row*width+col][c];
      }
      j = (j << 16) / i;
      for (sum=c=0; c < 3; c++) {
		ipix[c] = x3f_foveon_apply_curve (curve[c+3],
									  ((smrow[2][col][c] * j + 0x8000) >> 16) - image[row*width+col][c]);
		sum += ipix[c];
      }
      sum >>= 3;
      FORC3 {
		i = image[row*width+col][c] + ipix[c] - sum;
		if (i < 0) i = 0;
		image[row*width+col][c] = i;
      }
    }
  }
  free (shrink);
  free (smrow[6]);
  for (i=0; i < 8; i++)
    free (curve[i]);

  /* Trim off the black border */
  active[1] -= keep[1];
  active[3] -= 2;
  i = active[2] - active[0];
  for (row=0; row < active[3]-active[1]; row++)
	memcpy (image[row*i], image[(row+active[1])*width+active[0]],
			i * sizeof *image);
  i_img->img=image;
  i_img->width = i;
  i_img->height = row;
  printf("End of interpolation process\n");
/*   convert_to_rgb(i_img, 1); */

}
