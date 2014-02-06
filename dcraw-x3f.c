/*
   dcraw.c -- Dave Coffin's raw photo decoder
   Copyright 1997-2014 by Dave Coffin, dcoffin a cybercom o net

   This is a command-line ANSI C program to convert raw photos from
   any digital camera on any computer running any operating system.

   No license is required to download and use dcraw.c.  However,
   to lawfully redistribute dcraw, you must either (a) offer, at
   no extra charge, full source code* for all executable files
   containing RESTRICTED functions, (b) distribute this code under
   the GPL Version 2 or later, (c) remove all RESTRICTED functions,
   re-implement them, or copy them from an earlier, unrestricted
   Revision of dcraw.c, or (d) purchase a license from the author.

   The functions that process Foveon images have been RESTRICTED
   since Revision 1.237.  All other code remains free for all uses.

   *If you have not modified dcraw.c in any way, a link to my
   homepage qualifies as "full source code".

   $Revision: 1.459 $
   $Date: 2014/01/16 21:27:02 $
 */

#define DCRAW_VERSION "9.20"

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define _USE_MATH_DEFINES
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>

#if defined(DJGPP) || defined(__MINGW32__)
#define fseeko fseek
#define ftello ftell
#else
#define fgetc getc_unlocked
#endif
#ifdef __CYGWIN__
#include <io.h>
#endif
#ifdef WIN32
#include <sys/utime.h>
#include <winsock2.h>
#pragma comment(lib, "ws2_32.lib")
#define snprintf _snprintf
#define strcasecmp stricmp
#define strncasecmp strnicmp
typedef __int64 INT64;
typedef unsigned __int64 UINT64;
#else
#include <unistd.h>
#include <utime.h>
#include <netinet/in.h>
typedef long long INT64;
typedef unsigned long long UINT64;
#endif

#ifdef NODEPS
#define NO_JASPER
#define NO_JPEG
#define NO_LCMS
#endif
#ifndef NO_JASPER
#include <jasper/jasper.h>	/* Decode Red camera movies */
#endif
#ifndef NO_JPEG
#include <jpeglib.h>		/* Decode compressed Kodak DC120 photos */
#endif				/* and Adobe Lossy DNGs */
#ifndef NO_LCMS
#include <lcms2.h>		/* Support color profiles */
#endif
#ifdef LOCALEDIR
#include <libintl.h>
#define _(String) gettext(String)
#else
#define _(String) (String)
#endif

#ifdef LJPEG_DECODE
#error Please compile dcraw.c by itself.
#error Do not link it with ljpeg_decode.
#endif

#ifndef LONG_BIT
#define LONG_BIT (8 * sizeof (long))
#endif

#if !defined(uchar)
#define uchar unsigned char
#endif
#if !defined(ushort)
#define ushort unsigned short
#endif

/*
   All global variables are defined here, and all functions that
   access them are prefixed with "CLASS".  Note that a thread-safe
   C++ class cannot have non-const static local variables.
 */
FILE *ifp, *ofp;
short order;
const char *ifname;
char *meta_data, xtrans[6][6];
char cdesc[5], desc[512], make[64], model[64], model2[64], artist[64], cm_desc[64];
float flash_used, canon_ev, iso_speed, shutter, aperture, focal_len;
time_t timestamp;
off_t strip_offset, data_offset;
off_t thumb_offset, meta_offset, profile_offset;
unsigned shot_order, kodak_cbpp, exif_cfa, unique_id;
unsigned thumb_length, meta_length, profile_length;
unsigned thumb_misc, *oprof, fuji_layout, shot_select=0, multi_out=0;
unsigned tiff_nifds, tiff_samples, tiff_bps, tiff_compress;
unsigned black, maximum, mix_green, raw_color, zero_is_bad;
unsigned zero_after_ff, is_raw, dng_version, is_foveon, data_error;
unsigned tile_width, tile_length, gpsdata[32], load_flags;
unsigned flip, tiff_flip, filters, colors;
ushort raw_height, raw_width, height, width, top_margin, left_margin;
ushort shrink, iheight, iwidth, fuji_width, thumb_width, thumb_height;
ushort *raw_image, (*image)[4], cblack[4102];
ushort white[8][8], curve[0x10000], cr2_slice[3], sraw_mul[4];
double pixel_aspect, aber[4]={1,1,1,1}, gamm[6]={ 0.45,4.5,0,0,0,0 };
float bright=1, user_mul[4]={0,0,0,0}, threshold=0;
int mask[8][4];
int half_size=0, four_color_rgb=0, document_mode=0, highlight=0;
int verbose=0, use_auto_wb=0, use_camera_wb=0, use_camera_matrix=1;
int output_color=1, output_bps=8, output_tiff=0, med_passes=0;
int no_auto_bright=0;
unsigned greybox[4] = { 0, 0, UINT_MAX, UINT_MAX };
float cam_mul[4], pre_mul[4], cmatrix[3][4], rgb_cam[3][4];
const double xyz_rgb[3][3] = {			/* XYZ from RGB */
  { 0.412453, 0.357580, 0.180423 },
  { 0.212671, 0.715160, 0.072169 },
  { 0.019334, 0.119193, 0.950227 } };
const float d65_white[3] = { 0.950456, 1, 1.088754 };
int histogram[4][0x2000];
void (*write_thumb)(), (*write_fun)();
void (*load_raw)(), (*thumb_load_raw)();
jmp_buf failure;

struct decode {
  struct decode *branch[2];
  int leaf;
} first_decode[2048], *second_decode, *free_decode;

struct tiff_ifd {
  int width, height, bps, comp, phint, offset, flip, samples, bytes;
  int tile_width, tile_length;
} tiff_ifd[10];

struct ph1 {
  int format, key_off, black, black_off, split_col, tag_21a;
  float tag_210;
} ph1;

#define CLASS

#define FORC(cnt) for (c=0; c < cnt; c++)
#define FORC3 FORC(3)
#define FORC4 FORC(4)
#define FORCC FORC(colors)

#define SQR(x) ((x)*(x))
#define ABS(x) (((int)(x) ^ ((int)(x) >> 31)) - ((int)(x) >> 31))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define LIM(x,min,max) MAX(min,MIN(x,max))
#define ULIM(x,y,z) ((y) < (z) ? LIM(x,y,z) : LIM(x,z,y))
#define CLIP(x) LIM(x,0,65535)
#define SWAP(a,b) { a=a+b; b=a-b; a=a-b; }

/* #define RAW(row,col) \ */
/* 	raw_image[(row)*raw_width+(col)] */

#ifndef __GLIBC__
char *my_memmem (char *haystack, size_t haystacklen,
	      char *needle, size_t needlelen)
{
  char *c;
  for (c = haystack; c <= haystack + haystacklen - needlelen; c++)
    if (!memcmp (c, needle, needlelen))
      return c;
  return 0;
}
#define memmem my_memmem
char *my_strcasestr (char *haystack, const char *needle)
{
  char *c;
  for (c = haystack; *c; c++)
    if (!strncasecmp(c, needle, strlen(needle)))
      return c;
  return 0;
}
#define strcasestr my_strcasestr
#endif

void CLASS merror (void *ptr, const char *where)
{
  if (ptr) return;
  fprintf (stderr,_("%s: Out of memory in %s\n"), ifname, where);
  longjmp (failure, 1);
}

void CLASS derror()
{
  if (!data_error) {
    fprintf (stderr, "%s: ", ifname);
    if (feof(ifp))
      fprintf (stderr,_("Unexpected end of file\n"));
    else
      fprintf (stderr,_("Corrupt data near 0x%llx\n"), (INT64) ftello(ifp));
  }
  data_error++;
}

ushort CLASS sget2 (uchar *s)
{
  if (order == 0x4949)		/* "II" means little-endian */
    return s[0] | s[1] << 8;
  else				/* "MM" means big-endian */
    return s[0] << 8 | s[1];
}

ushort CLASS get2()
{
  uchar str[2] = { 0xff,0xff };
  fread (str, 1, 2, ifp);
  return sget2(str);
}

unsigned CLASS sget4 (uchar *s)
{
  if (order == 0x4949)
    return s[0] | s[1] << 8 | s[2] << 16 | s[3] << 24;
  else
    return s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3];
}
#define sget4(s) sget4((uchar *)s)

unsigned CLASS get4()
{
  uchar str[4] = { 0xff,0xff,0xff,0xff };
  fread (str, 1, 4, ifp);
  return sget4(str);
}

unsigned CLASS getint (int type)
{
  return type == 3 ? get2() : get4();
}

float CLASS int_to_float (int i)
{
  union { int i; float f; } u;
  u.i = i;
  return u.f;
}

double CLASS getreal (int type)
{
  union { char c[8]; double d; } u;
  int i, rev;

  switch (type) {
    case 3: return (unsigned short) get2();
    case 4: return (unsigned int) get4();
    case 5:  u.d = (unsigned int) get4();
      return u.d / (unsigned int) get4();
    case 8: return (signed short) get2();
    case 9: return (signed int) get4();
    case 10: u.d = (signed int) get4();
      return u.d / (signed int) get4();
    case 11: return int_to_float (get4());
    case 12:
      rev = 7 * ((order == 0x4949) == (ntohs(0x1234) == 0x1234));
      for (i=0; i < 8; i++)
	u.c[i ^ rev] = fgetc(ifp);
      return u.d;
    default: return fgetc(ifp);
  }
}

void CLASS read_shorts (ushort *pixel, int count)
{
  if (fread (pixel, 2, count, ifp) < count) derror();
  if ((order == 0x4949) == (ntohs(0x1234) == 0x1234))
    swab (pixel, pixel, count*2);
}

unsigned CLASS getbithuff (int nbits, ushort *huff)
{
  static unsigned bitbuf=0;
  static int vbits=0, reset=0;
  unsigned c;

  if (nbits > 25) return 0;
  if (nbits < 0)
    return bitbuf = vbits = reset = 0;
  if (nbits == 0 || vbits < 0) return 0;
  while (!reset && vbits < nbits && (c = fgetc(ifp)) != EOF &&
    !(reset = zero_after_ff && c == 0xff && fgetc(ifp))) {
    bitbuf = (bitbuf << 8) + (uchar) c;
    vbits += 8;
  }
  c = bitbuf << (32-vbits) >> (32-nbits);
  if (huff) {
    vbits -= huff[c] >> 8;
    c = (uchar) huff[c];
  } else
    vbits -= nbits;
  if (vbits < 0) derror();
  return c;
}

#define getbits(n) getbithuff(n,0)
#define gethuff(h) getbithuff(*h,h+1)

ushort * CLASS make_decoder_ref (const uchar **source)
{
  int max, len, h, i, j;
  const uchar *count;
  ushort *huff;

  count = (*source += 16) - 17;
  for (max=16; max && !count[max]; max--);
  huff = (ushort *) calloc (1 + (1 << max), sizeof *huff);
  merror (huff, "make_decoder()");
  huff[0] = max;
  for (h=len=1; len <= max; len++)
    for (i=0; i < count[len]; i++, ++*source)
      for (j=0; j < 1 << (max-len); j++)
	if (h <= 1 << max)
	  huff[h++] = len << 8 | **source;
  return huff;
}

/*
   Not a full implementation of Lossless JPEG, just
   enough to decode Canon, Kodak and Adobe DNG images.
 */
struct jhead {
  int bits, high, wide, clrs, sraw, psv, restart, vpred[6];
  ushort *huff[6], *free[4], *row;
};

int CLASS ljpeg_start (struct jhead *jh, int info_only)
{
  int c, tag, len;
  uchar data[0x10000];
  const uchar *dp;

  memset (jh, 0, sizeof *jh);
  jh->restart = INT_MAX;
  fread (data, 2, 1, ifp);
  if (data[1] != 0xd8) return 0;
  do {
    fread (data, 2, 2, ifp);
    tag =  data[0] << 8 | data[1];
    len = (data[2] << 8 | data[3]) - 2;
    if (tag <= 0xff00) return 0;
    fread (data, 1, len, ifp);
    switch (tag) {
      case 0xffc3:
	jh->sraw = ((data[7] >> 4) * (data[7] & 15) - 1) & 3;
      case 0xffc0:
	jh->bits = data[0];
	jh->high = data[1] << 8 | data[2];
	jh->wide = data[3] << 8 | data[4];
	jh->clrs = data[5] + jh->sraw;
	if (len == 9 && !dng_version) getc(ifp);
	break;
      case 0xffc4:
	if (info_only) break;
	for (dp = data; dp < data+len && (c = *dp++) < 4; )
	  jh->free[c] = jh->huff[c] = make_decoder_ref (&dp);
	break;
      case 0xffda:
	jh->psv = data[1+data[0]*2];
	jh->bits -= data[3+data[0]*2] & 15;
	break;
      case 0xffdd:
	jh->restart = data[0] << 8 | data[1];
    }
  } while (tag != 0xffda);
  if (info_only) return 1;
  if (jh->clrs > 6 || !jh->huff[0]) return 0;
  FORC(5) if (!jh->huff[c+1]) jh->huff[c+1] = jh->huff[c];
  if (jh->sraw) {
    FORC(4)        jh->huff[2+c] = jh->huff[1];
    FORC(jh->sraw) jh->huff[1+c] = jh->huff[0];
  }
  jh->row = (ushort *) calloc (jh->wide*jh->clrs, 4);
  merror (jh->row, "ljpeg_start()");
  return zero_after_ff = 1;
}

void CLASS ljpeg_end (struct jhead *jh)
{
  int c;
  FORC4 if (jh->free[c]) free (jh->free[c]);
  free (jh->row);
}

int CLASS ljpeg_diff (ushort *huff)
{
  int len, diff;

  len = gethuff(huff);
  if (len == 16 && (!dng_version || dng_version >= 0x1010000))
    return -32768;
  diff = getbits(len);
  if ((diff & (1 << (len-1))) == 0)
    diff -= (1 << len) - 1;
  return diff;
}

ushort * CLASS ljpeg_row (int jrow, struct jhead *jh)
{
  int col, c, diff, pred, spred=0;
  ushort mark=0, *row[3];

  if (jrow * jh->wide % jh->restart == 0) {
    FORC(6) jh->vpred[c] = 1 << (jh->bits-1);
    if (jrow) {
      fseek (ifp, -2, SEEK_CUR);
      do mark = (mark << 8) + (c = fgetc(ifp));
      while (c != EOF && mark >> 4 != 0xffd);
    }
    getbits(-1);
  }
  FORC3 row[c] = jh->row + jh->wide*jh->clrs*((jrow+c) & 1);
  for (col=0; col < jh->wide; col++)
    FORC(jh->clrs) {
      diff = ljpeg_diff (jh->huff[c]);
      if (jh->sraw && c <= jh->sraw && (col | c))
		    pred = spred;
      else if (col) pred = row[0][-jh->clrs];
      else	    pred = (jh->vpred[c] += diff) - diff;
      if (jrow && col) switch (jh->psv) {
	case 1:	break;
	case 2: pred = row[1][0];					break;
	case 3: pred = row[1][-jh->clrs];				break;
	case 4: pred = pred +   row[1][0] - row[1][-jh->clrs];		break;
	case 5: pred = pred + ((row[1][0] - row[1][-jh->clrs]) >> 1);	break;
	case 6: pred = row[1][0] + ((pred - row[1][-jh->clrs]) >> 1);	break;
	case 7: pred = (pred + row[1][0]) >> 1;				break;
	default: pred = 0;
      }
      if ((**row = pred + diff) >> jh->bits) derror();
      if (c <= jh->sraw) spred = **row;
      row[0]++; row[1]++;
    }
  return row[2];
}

void CLASS jpeg_thumb();

void CLASS ppm_thumb()
{
  char *thumb;
  thumb_length = thumb_width*thumb_height*3;
  thumb = (char *) malloc (thumb_length);
  merror (thumb, "ppm_thumb()");
  fprintf (ofp, "P6\n%d %d\n255\n", thumb_width, thumb_height);
  fread  (thumb, 1, thumb_length, ifp);
  fwrite (thumb, 1, thumb_length, ofp);
  free (thumb);
}

void CLASS ppm16_thumb()
{
  int i;
  char *thumb;
  thumb_length = thumb_width*thumb_height*3;
  thumb = (char *) calloc (thumb_length, 2);
  merror (thumb, "ppm16_thumb()");
  read_shorts ((ushort *) thumb, thumb_length);
  for (i=0; i < thumb_length; i++)
    thumb[i] = ((ushort *) thumb)[i] >> 8;
  fprintf (ofp, "P6\n%d %d\n255\n", thumb_width, thumb_height);
  fwrite (thumb, 1, thumb_length, ofp);
  free (thumb);
}

void CLASS layer_thumb()
{
  int i, c;
  char *thumb, map[][4] = { "012","102" };

  colors = thumb_misc >> 5 & 7;
  thumb_length = thumb_width*thumb_height;
  thumb = (char *) calloc (colors, thumb_length);
  merror (thumb, "layer_thumb()");
  fprintf (ofp, "P%d\n%d %d\n255\n",
	5 + (colors >> 1), thumb_width, thumb_height);
  fread (thumb, thumb_length, colors, ifp);
  for (i=0; i < thumb_length; i++)
    FORCC putc (thumb[i+thumb_length*(map[thumb_misc >> 8][c]-'0')], ofp);
  free (thumb);
}

/* RESTRICTED code starts here */

void CLASS foveon_decoder (unsigned size, unsigned code)
{
  static unsigned huff[1024];
  struct decode *cur;
  int i, len;

  if (!code) {
    for (i=0; i < size; i++)
      huff[i] = get4();
    memset (first_decode, 0, sizeof first_decode);
    free_decode = first_decode;
  }
  cur = free_decode++;
  if (free_decode > first_decode+2048) {
    fprintf (stderr,_("%s: decoder table overflow\n"), ifname);
    longjmp (failure, 2);
  }
  if (code)
    for (i=0; i < size; i++)
      if (huff[i] == code) {
	cur->leaf = i;
	return;
      }
  if ((len = code >> 27) > 26) return;
  code = (len+1) << 27 | (code & 0x3ffffff) << 1;

  cur->branch[0] = free_decode;
  foveon_decoder (size, code);
  cur->branch[1] = free_decode;
  foveon_decoder (size, code+1);
}

void CLASS foveon_thumb()
{
  unsigned bwide, row, col, bitbuf=0, bit=1, c, i;
  char *buf;
  struct decode *dindex;
  short pred[3];

  bwide = get4();
  fprintf (ofp, "P6\n%d %d\n255\n", thumb_width, thumb_height);
  if (bwide > 0) {
    if (bwide < thumb_width*3) return;
    buf = (char *) malloc (bwide);
    merror (buf, "foveon_thumb()");
    for (row=0; row < thumb_height; row++) {
      fread  (buf, 1, bwide, ifp);
      fwrite (buf, 3, thumb_width, ofp);
    }
    free (buf);
    return;
  }
  foveon_decoder (256, 0);

  for (row=0; row < thumb_height; row++) {
    memset (pred, 0, sizeof pred);
    if (!bit) get4();
    for (bit=col=0; col < thumb_width; col++)
      FORC3 {
	for (dindex=first_decode; dindex->branch[0]; ) {
	  if ((bit = (bit-1) & 31) == 31)
	    for (i=0; i < 4; i++)
	      bitbuf = (bitbuf << 8) + fgetc(ifp);
	  dindex = dindex->branch[bitbuf >> bit & 1];
	}
	pred[c] += dindex->leaf;
	fputc (pred[c], ofp);
      }
  }
}

void CLASS foveon_sd_load_raw()
{
  struct decode *dindex;
  short diff[1024];
  unsigned bitbuf=0;
  int pred[3], row, col, bit=-1, c, i;

  read_shorts ((ushort *) diff, 1024);
  if (!load_flags) foveon_decoder (1024, 0);

  for (row=0; row < height; row++) {
    memset (pred, 0, sizeof pred);
    if (!bit && !load_flags && atoi(model+2) < 14) get4();
    for (col=bit=0; col < width; col++) {
      if (load_flags) {
	bitbuf = get4();
	FORC3 pred[2-c] += diff[bitbuf >> c*10 & 0x3ff];
      }
      else FORC3 {
	for (dindex=first_decode; dindex->branch[0]; ) {
	  if ((bit = (bit-1) & 31) == 31)
	    for (i=0; i < 4; i++)
	      bitbuf = (bitbuf << 8) + fgetc(ifp);
	  dindex = dindex->branch[bitbuf >> bit & 1];
	}
	pred[c] += diff[dindex->leaf];
	if (pred[c] >> 16 && ~pred[c] >> 16) derror();
      }
      FORC3 image[row*width+col][c] = pred[c];
    }
  }
}

void CLASS foveon_huff (ushort *huff)
{
  int i, j, clen, code;

  huff[0] = 8;
  for (i=0; i < 13; i++) {
    clen = getc(ifp);
    code = getc(ifp);
    for (j=0; j < 256 >> clen; )
      huff[code+ ++j] = clen << 8 | i;
  }
  get2();
}

void CLASS foveon_dp_load_raw()
{
  unsigned c, roff[4], row, col, diff;
  ushort huff[512], vpred[2][2], hpred[2];

  fseek (ifp, 8, SEEK_CUR);
  foveon_huff (huff);
  roff[0] = 48;
  FORC3 roff[c+1] = -(-(roff[c] + get4()) & -16);
  FORC3 {
    fseek (ifp, data_offset+roff[c], SEEK_SET);
    getbits(-1);
    vpred[0][0] = vpred[0][1] = vpred[1][0] = vpred[1][1] = 512;
    for (row=0; row < height; row++) {
      for (col=0; col < width; col++) {
	diff = ljpeg_diff(huff);
	if (col < 2) hpred[col] = vpred[row & 1][col] += diff;
	else hpred[col & 1] += diff;
	image[row*width+col][c] = hpred[col & 1];
      }
    }
  }
}

void CLASS foveon_load_camf()
{
  unsigned type, wide, high, i, j, row, col, diff;
  ushort huff[258], vpred[2][2] = {{512,512},{512,512}}, hpred[2];

  fseek (ifp, meta_offset, SEEK_SET);
  type = get4();  get4();  get4();
  wide = get4();
  high = get4();
  if (type == 2) {
    fread (meta_data, 1, meta_length, ifp);
    for (i=0; i < meta_length; i++) {
      high = (high * 1597 + 51749) % 244944;
      wide = high * (INT64) 301593171 >> 24;
      meta_data[i] ^= ((((high << 8) - wide) >> 1) + wide) >> 17;
    }
  } else if (type == 4) {
    free (meta_data);
    meta_data = (char *) malloc (meta_length = wide*high*3/2);
    merror (meta_data, "foveon_load_camf()");
    foveon_huff (huff);
    get4();
    getbits(-1);
    for (j=row=0; row < high; row++) {
      for (col=0; col < wide; col++) {
	diff = ljpeg_diff(huff);
	if (col < 2) hpred[col] = vpred[row & 1][col] += diff;
	else         hpred[col & 1] += diff;
	if (col & 1) {
	  meta_data[j++] = hpred[0] >> 4;
	  meta_data[j++] = hpred[0] << 4 | hpred[1] >> 8;
	  meta_data[j++] = hpred[1];
	}
      }
    }
  } else
    fprintf (stderr,_("%s has unknown CAMF type %d.\n"), ifname, type);
}

const char * CLASS foveon_camf_param (const char *block, const char *param)
{
  unsigned idx, num;
  char *pos, *cp, *dp;

  for (idx=0; idx < meta_length; idx += sget4(pos+8)) {
    pos = meta_data + idx;
    if (strncmp (pos, "CMb", 3)) break;
    if (pos[3] != 'P') continue;
    if (strcmp (block, pos+sget4(pos+12))) continue;
    cp = pos + sget4(pos+16);
    num = sget4(cp);
    dp = pos + sget4(cp+4);
    while (num--) {
      cp += 8;
      if (!strcmp (param, dp+sget4(cp)))
	return dp+sget4(cp+4);
    }
  }
  return 0;
}

void * CLASS foveon_camf_matrix (unsigned dim[3], const char *name)
{
  unsigned i, idx, type, ndim, size, *mat;
  char *pos, *cp, *dp;
  double dsize;

  for (idx=0; idx < meta_length; idx += sget4(pos+8)) {
    pos = meta_data + idx;
    if (strncmp (pos, "CMb", 3)) break;
    if (pos[3] != 'M') continue;
    if (strcmp (name, pos+sget4(pos+12))) continue;
    dim[0] = dim[1] = dim[2] = 1;
    cp = pos + sget4(pos+16);
    type = sget4(cp);
    if ((ndim = sget4(cp+4)) > 3) break;
    dp = pos + sget4(cp+8);
    for (i=ndim; i--; ) {
      cp += 12;
      dim[i] = sget4(cp);
    }
    if ((dsize = (double) dim[0]*dim[1]*dim[2]) > meta_length/4) break;
    mat = (unsigned *) malloc ((size = dsize) * 4);
    merror (mat, "foveon_camf_matrix()");
    for (i=0; i < size; i++)
      if (type && type != 6)
	mat[i] = sget4(dp + i*4);
      else
	mat[i] = sget4(dp + i*2) & 0xffff;
    return mat;
  }
  fprintf (stderr,_("%s: \"%s\" matrix not found!\n"), ifname, name);
  return 0;
}

int CLASS foveon_fixed (void *ptr, int size, const char *name)
{
  void *dp;
  unsigned dim[3];

  if (!name) return 0;
  dp = foveon_camf_matrix (dim, name);
  if (!dp) return 0;
  memcpy (ptr, dp, size*4);
  free (dp);
  return 1;
}

float CLASS foveon_avg (short *pix, int range[2], float cfilt)
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

short * CLASS foveon_make_curve (double max, double mul, double filt)
{
  short *curve;
  unsigned i, size;
  double x;

  if (!filt) filt = 0.8;
  size = 4*M_PI*max / filt;
  if (size == UINT_MAX) size--;
  curve = (short *) calloc (size+1, sizeof *curve);
  merror (curve, "foveon_make_curve()");
  curve[0] = size;
  for (i=0; i < size; i++) {
    x = i*filt/max/4;
    curve[i+1] = (cos(x)+1)/2 * tanh(i*filt/mul) * mul + 0.5;
  }
  return curve;
}

void CLASS foveon_make_curves
	(short **curvep, float dq[3], float div[3], float filt)
{
  double mul[3], max=0;
  int c;

  FORC3 mul[c] = dq[c]/div[c];
  FORC3 if (max < mul[c]) max = mul[c];
  FORC3 curvep[c] = foveon_make_curve (max, mul[c], filt);
}

int CLASS foveon_apply_curve (short *curve, int i)
{
  if (abs(i) >= curve[0]) return 0;
  return i < 0 ? -curve[1-i] : curve[1+i];
}

#define image ((short (*)[4]) image)

void CLASS foveon_f20_interpolate()
{
  static const short hood[] = { -1,-1, -1,0, -1,1, 0,-1, 0,1, 1,-1, 1,0, 1,1 };
  short *pix, prev[3], *curve[8], (*shrink)[3];
  float cfilt=0, ddft[3][3][2], ppm[3][3][3];
  float cam_xyz[3][3], correct[3][3]={{1,0,0},{0,1,0},{0,0,1}},cp2[3][3], last[3][3], last2[3][3], trans[3][3];
  float chroma_dq[3]={4,4,4}, color_dq[3], diag[3][3], div[3], tempgainfact[3];
  float (*black)[3], (*sgain)[3], (*sgrow)[3];
  float fsum[3], val, frow, num;
  int row, col, c, i, j, diff, sgx, irow, sum, min, max, limit;
  int dscr[2][2], dstb[4], (*smrow[7])[3], total[4], ipix[3];
  int work[3][3], smlast, smred, smred_p=0, dev[3];
  int satlev[3], keep[4], active[4];
  unsigned dim[3], *badpix;
  double dsum=0, trsum[3];
  char str[128];
  const char* cp;

  if (verbose)
    fprintf (stderr,_("Foveon F20 interpolation...\n"));

  foveon_load_camf();
  foveon_fixed (dscr, 4, "DarkShieldColRange");
  printf("dscr: %d %d %d %d\n", dscr[0][0], dscr[0][1], dscr[1][0],
  dscr[1][1]);
  memset(ppm, 0, sizeof ppm);
  foveon_fixed (ppm, 9, "CP2_Matrix");

  foveon_fixed (keep, 4, "KeepImageArea");
  printf("keep: %d %d %d\n", keep[0], keep[1], keep[2], keep[3]);
  foveon_fixed (active, 4, "ActiveImageArea");
  printf("active: %d %d %d\n", active[0], active[1], active[2], active[3]);
  if (foveon_camf_param ("IncludeBlocks", "ChromaDQ"))
	foveon_fixed (chroma_dq, 3, "ChromaDQ");
  printf("chroma_dq: %f %f %f\n", chroma_dq[0], chroma_dq[1], chroma_dq[2]);
  foveon_fixed (color_dq, 3,
	foveon_camf_param ("IncludeBlocks", "ColorDQ") ?
		"ColorDQ" : "ColorDQCamRGB");
  printf("color_dq: %f %f %f\n", color_dq[0], color_dq[1], color_dq[2]);
  if (foveon_camf_param ("IncludeBlocks", "ColumnFilter"))
		 foveon_fixed (&cfilt, 1, "ColumnFilter");
  printf("cfilt: %f\n", cfilt);

  memset (ddft, 0, sizeof ddft);
  if (!foveon_camf_param ("IncludeBlocks", "DarkDrift")
	 || !foveon_fixed (ddft[1][0], 12, "DarkDrift"))
    for (i=0; i < 2; i++) {
      foveon_fixed (dstb, 4, i ? "DarkShieldBottom":"DarkShieldTop");
      for (row = dstb[1]; row <= dstb[3]; row++)
	for (col = dstb[0]; col <= dstb[2]; col++)
	  FORC3 ddft[i+1][c][1] += (short) image[row*width+col][c];
      FORC3 ddft[i+1][c][1] /= (dstb[3]-dstb[1]+1) * (dstb[2]-dstb[0]+1);
    }
  printf("ddft: %f %f %f\n", ddft[1][0][1], ddft[1][0][2], ddft[1][1][1]);


  if (foveon_camf_param ("IncludeBlocks", "WhiteBalanceColorCorrections"))
	foveon_fixed (cam_xyz, 9,
				  foveon_camf_param ("WhiteBalanceColorCorrections", model2));
  printf("cam_xyz: %f %f %f\n", cam_xyz[0][0], cam_xyz[0][1], cam_xyz[0][2]);
  if (foveon_camf_param ("IncludeBlocks", "ColorModeCompensations"))
	foveon_fixed (correct, 9,
				  foveon_camf_param ("ColorModeCompensations", cm_desc));
  printf("correct: %f %f %f\n", correct[0][0], correct[0][1], correct[0][2]);

  memset (last2, 0, sizeof last2);
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
/* 	  last[i][j]=cam_xyz[i][j]; */
	  FORC3 last2[i][j] += cam_xyz[i][c]* correct[c][j];
/*   if (foveon_camf_param ("IncludeBlocks", "CP2_Matrix")) */
/* 	foveon_fixed (cp2, 9, */
/* 				  foveon_camf_param ("CP2_Matrix", cm_desc)); */
/*   printf("correct: %f %f %f\n", correct[0][0], correct[0][1], correct[0][2]); */

/*   memset (last2, 0, sizeof last2); */
/*   for (i=0; i < 3; i++) */
/*     for (j=0; j < 3; j++) */
/* /\* 	  last[i][j]=cam_xyz[i][j]; *\/ */
/* 	  FORC3 last2[i][j] += last[i][c]* cp2[c][j]; */

  #define LAST(x,y) last2[(i+x)%3][(c+y)%3]
  for (i=0; i < 3; i++)
    FORC3 diag[c][i] = LAST(1,1)*LAST(2,2) - LAST(1,2)*LAST(2,1);
  #undef LAST
  FORC3 div[c] = diag[c][0]*0.3127 + diag[c][1]*0.329 + diag[c][2]*0.3583;
  FORC3 printf("div[%d]=%f\n", c, div[c]);


  sprintf(str, "%sWBGains", model2);
  if (foveon_camf_param ("IncludeBlocks", "WhiteBalanceGains"))
	foveon_fixed (div, 3, foveon_camf_param ("WhiteBalanceGains", model2));
  if (foveon_camf_param("IncludeBlocks", "TempGainFact"))
	foveon_fixed(tempgainfact, 3, "TempGainFact");
  FORC3 div[c]*=tempgainfact[c];
  if (foveon_camf_param("IncludeBlocks", "SensorAdjustmentGainFact"))
	foveon_fixed(tempgainfact, 3, "SensorAdjustmentGainFact");
  FORC3 div[c]*=tempgainfact[c];
  num = 0;
  FORC3 if (num < div[c]) num = div[c];
  FORC3 div[c] /= num;
  FORC3 printf("div[%d]=%f\n", c, div[c]);

  memset (trans, 0, sizeof trans);
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      FORC3 { trans[i][j] += rgb_cam[i][c] * last2[c][j] * div[j];
	}
  FORC3 trsum[c] = trans[c][0] + trans[c][1] + trans[c][2];
  dsum = (6*trsum[0] + 11*trsum[1] + 3*trsum[2]) / 20;
  for (i=0; i < 3; i++)
    FORC3 last[i][c] = trans[i][c] * dsum / trsum[i];
/*   memset (trans, 0, sizeof trans); */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
	  printf("trans[%d][%d]: %f\n", i, j, trans[i][j]);


/*   FORC3 color_dq[c]/=8; */
  foveon_make_curves (curve, color_dq, div, cfilt);
  FORC3 chroma_dq[c] /= 3;
  foveon_make_curves (curve+3, chroma_dq, div, cfilt);
  FORC3 dsum += chroma_dq[c] / div[c];
  curve[6] = foveon_make_curve (dsum, dsum, cfilt);
  curve[7] = foveon_make_curve (dsum*2, dsum*2, cfilt);

/*   black = (float (*)[3]) calloc (height, sizeof *black); */
/*   for (row=0; row < height; row++) { */
/*     for (i=0; i < 6; i++) */
/*       ddft[0][0][i] = ddft[1][0][i] + */
/* 	row / (height-1.0) * (ddft[2][0][i] - ddft[1][0][i]); */
/*     FORC3 black[row][c] = */
/* 	( foveon_avg (image[row*width]+c, dscr[0], cfilt) + */
/* 	  foveon_avg (image[row*width]+c, dscr[1], cfilt) * 3 */
/* 	  - ddft[0][c][0] ) / 4 - ddft[0][c][1]; */
/*   } */
/*   memcpy (black, black+8, sizeof *black*8); */
/*   memcpy (black+height-11, black+height-22, 11*sizeof *black); */
/*   memcpy (last, black, sizeof last); */

/*   for (row=1; row < height-1; row++) { */
/*     FORC3 if (last[1][c] > last[0][c]) { */
/* 	if (last[1][c] > last[2][c]) */
/* 	  black[row][c] = (last[0][c] > last[2][c]) ? last[0][c]:last[2][c]; */
/*       } else */
/* 	if (last[1][c] < last[2][c]) */
/* 	  black[row][c] = (last[0][c] < last[2][c]) ? last[0][c]:last[2][c]; */
/*     memmove (last, last+1, 2*sizeof last[0]); */
/*     memcpy (last[2], black[row+1], sizeof last[2]); */
/*   } */
/*   FORC3 black[row][c] = (last[0][c] + last[1][c])/2; */
/*   FORC3 black[0][c] = (black[1][c] + black[3][c])/2; */

/*   val = 1 - exp(-1/24.0); */
/*   memcpy (fsum, black, sizeof fsum); */
/*   for (row=1; row < height; row++) */
/*     FORC3 fsum[c] += black[row][c] = */
/* 	(black[row][c] - black[row-1][c])*val + black[row-1][c]; */
/*   memcpy (last[0], black[height-1], sizeof last[0]); */
/*   FORC3 fsum[c] /= height; */
/*   for (row = height; row--; ) */
/*     FORC3 last[0][c] = black[row][c] = */
/* 	(black[row][c] - fsum[c] - last[0][c])*val + last[0][c]; */

/*   memset (total, 0, sizeof total); */
/*   for (row=2; row < height; row+=4) */
/*     for (col=2; col < width; col+=4) { */
/*       FORC3 total[c] += (short) image[row*width+col][c]; */
/*       total[3]++; */
/*     } */
/*   for (row=0; row < height; row++) */
/*     FORC3 black[row][c] += fsum[c]/2 + total[c]/(total[3]*100.0); */

/* /\*   sgain = (float (*)[3]) foveon_camf_matrix (dim, "SpatialGain"); *\/ */
/* /\*   if (!sgain) return; *\/ */
/* /\*   sgrow = (float (*)[3]) calloc (dim[1], sizeof *sgrow); *\/ */
/* /\*   sgx = (width + dim[1]-2) / (dim[1]-1); *\/ */

/*   for (row=0; row < height; row++) { */
/*     for (i=0; i < 6; i++) */
/*       ddft[0][0][i] = ddft[1][0][i] + */
/* 	row / (height-1.0) * (ddft[2][0][i] - ddft[1][0][i]); */
/*     pix = image[row*width]; */
/*     memcpy (prev, pix, sizeof prev); */
/* /\*     frow = row / (height-1.0) * (dim[2]-1); *\/ */
/* /\*     if ((irow = frow) == dim[2]-1) irow--; *\/ */
/* /\*     frow -= irow; *\/ */
/* /\*     for (i=0; i < dim[1]; i++) *\/ */
/* /\*       FORC3 sgrow[i][c] = sgain[ irow   *dim[1]+i][c] * (1-frow) + *\/ */
/* /\* 			  sgain[(irow+1)*dim[1]+i][c] *    frow; *\/ */
/*     for (col=0; col < width; col++) { */
/* /\*       FORC3 { *\/ */
/* /\* 	diff = pix[c] - prev[c]; *\/ */
/* /\* 	prev[c] = pix[c]; *\/ */
/* /\* 	ipix[c] = pix[c] + floor ((diff + (diff*diff >> 14)) * cfilt *\/ */
/* /\* 		- ddft[0][c][1] - ddft[0][c][0] * ((float) col/width - 0.5) *\/ */
/* /\* 		- black[row][c] ); *\/ */
/* /\*       } *\/ */
/* /\*       FORC3 { *\/ */
/* /\* 	work[0][c] = ipix[c] * ipix[c] >> 14; *\/ */
/* /\* 	work[2][c] = ipix[c] * work[0][c] >> 14; *\/ */
/* /\* 	work[1][2-c] = ipix[(c+1) % 3] * ipix[(c+2) % 3] >> 14; *\/ */
/* /\*       } *\/ */
/*       FORC3 { */
/* /\* 	for (val=i=0; i < 3; i++) *\/ */
/* /\* 	  for (  j=0; j < 3; j++) *\/ */
/* /\* 	    val += ppm[c][i][j] * work[i][j]; *\/ */
/* 	ipix[c] = floor ((ipix[c]/\*  + floor(val) *\/) / div[c]); */
/* 	if (ipix[c] > 32000) ipix[c] = 32000; */
/* 	pix[c] = ipix[c]; */
/*       } */
/*       pix += 4; */
/*     } */
/*   } */
/*  /\*  free (sgrow); *\/ */
/*   free (black); */
/* /\*   free (sgain); *\/ */

/*   if ((badpix = (unsigned int *) foveon_camf_matrix (dim, "BadPixels"))) { */
/*     for (i=0; i < dim[0]; i++) { */
/*       col = (badpix[i] >> 8 & 0xfff) - keep[0]; */
/*       row = (badpix[i] >> 20       ) - keep[1]; */
/*       if ((unsigned)(row-1) > height-3 || (unsigned)(col-1) > width-3) */
/* 	continue; */
/*       memset (fsum, 0, sizeof fsum); */
/*       for (sum=j=0; j < 8; j++) */
/* 	if (badpix[i] & (1 << j)) { */
/* 	  FORC3 fsum[c] += (short) */
/* 		image[(row+hood[j*2])*width+col+hood[j*2+1]][c]; */
/* 	  sum++; */
/* 	} */
/*       if (sum) FORC3 image[row*width+col][c] = fsum[c]/sum; */
/*     } */
/*     free (badpix); */
/*   } */

  /* Array for 5x5 Gaussian averaging of red values */
  smrow[6] = (int (*)[3]) calloc (width*5, sizeof **smrow);
  merror (smrow[6], "foveon_interpolate()");
  for (i=0; i < 5; i++)
    smrow[i] = smrow[6] + i*width;

  /* Sharpen the reds against these Gaussian averages */
  for (smlast=-1, row=2; row < height-2; row++) {
    while (smlast < row+2) {
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

  /* Adjust the brighter pixels for better linearity */
  min = 0xffff;
  FORC3 {
    i = 6000 / div[c];
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
    while (smlast < row+2) {
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
      FORC3 dev[c] = -foveon_apply_curve (curve[7], pix[c] -
	((smrow[1][col][c] + 2*smrow[2][col][c] + smrow[3][col][c]) >> 2));
      sum = (dev[0] + dev[1] + dev[2]) >> 3;
      FORC3 pix[c] += dev[c] - sum;
      pix += 4;
    }
  }
  for (smlast=-1, row=2; row < height-2; row++) {
    while (smlast < row+2) {
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
      FORC3 pix[c] += foveon_apply_curve (curve[6],
		((j*total[c] + 0x8000) >> 16) - pix[c]);
      pix += 4;
    }
  }

/*   /\* Transform the image to a different colorspace *\/ */
/*   for (pix=image[0]; pix < image[height*width]; pix+=4) { */
/*     FORC3 pix[c] -= foveon_apply_curve (curve[c], pix[c]); */
/*     sum = (pix[0]+pix[1]+pix[1]+pix[2]) >> 2; */
/*     FORC3 pix[c] -= foveon_apply_curve (curve[c], pix[c]-sum); */
/*     FORC3 { */
/*       for (dsum=i=0; i < 3; i++) */
/* 	dsum += trans[c][i] * pix[i]; */
/*       if (dsum < 0)  dsum = 0; */
/*       if (dsum > 24000) dsum = 24000; */
/*       ipix[c] = dsum + 0.5; */
/*     } */
/*     FORC3 pix[c] = ipix[c]; */
/*   } */


  for (row=0; row<height; row++){
	pix = image[row*width];
	for (col=0; col<width; col++){
	  for (c=0; c < 3; c++) ipix[c]=pix[c];
	  for (c=0; c < 3; c++) {
		float tmp;
		tmp=rgb_cam[c][0]*last2[c][0]* div[0]*ipix[0]+rgb_cam[c][1]*last2[c][1]*div[1]*ipix[1]+rgb_cam[c][2]*last2[c][2]*div[2]*ipix[2];
/* 		tmp=trans[c][0]*1.2/\* * div[0] *\/\*ipix[0]+trans[c][1]/\* *div[1] *\/\*ipix[1]+trans[c][2]*1.1/\* *div[2] *\/\*ipix[2]; */
/* 		if (c==0) tmp*=1.96875; */
/* 		if (c==2) tmp*=.625; */
		if (tmp<0) pix[c]=0;
/* 		else if (tmp >24000) pix[c]=24000; */
		else pix[c]=floor(tmp);
	  }
	  pix+=4;
	}

  }

  /* Smooth the image bottom-to-top and save at 1/4 scale */
  shrink = (short (*)[3]) calloc ((height/4), (width/4)*sizeof *shrink);
  merror (shrink, "foveon_interpolate()");
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
	ipix[c] = foveon_apply_curve (curve[c+3],
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
  width = i;
  height = row;
}

void CLASS foveon_interpolate()
{
  static const short hood[] = { -1,-1, -1,0, -1,1, 0,-1, 0,1, 1,-1, 1,0, 1,1 };
  short *pix, prev[3], *curve[8], (*shrink)[3];
  float cfilt=0, ddft[3][3][2], ppm[3][3][3];
  float cam_xyz[3][3], correct[3][3]={{1,0,0},{0,1,0},{0,0,1}}, last[3][3], trans[3][3];
  float chroma_dq[3], color_dq[3], diag[3][3], div[3];
  float (*black)[3], (*sgain)[3], (*sgrow)[3];
  float fsum[3], val, frow, num;
  int row, col, c, i, j, diff, sgx, irow, sum, min, max, limit;
  int dscr[2][2], dstb[4], (*smrow[7])[3], total[4], ipix[3];
  int work[3][3], smlast, smred, smred_p=0, dev[3];
  int satlev[3], keep[4], active[4];
  unsigned dim[3], *badpix;
  double dsum=0, trsum[3];
  char str[128];
  const char* cp;

  if (verbose)
    fprintf (stderr,_("Foveon interpolation...\n"));

  foveon_load_camf();
  foveon_fixed (dscr, 4, "DarkShieldColRange");
  printf("dscr: %d %d %d %d\n", dscr[0][0], dscr[0][1], dscr[1][0],
  dscr[1][1]);
  memset(ppm, 0, sizeof ppm);
  if (foveon_camf_param ("IncludeBlocks", "PostPolyMatrix"))
	foveon_fixed (ppm[0][0], 27, "PostPolyMatrix");
  printf("ppm: %f %f %f\n", ppm[0][0][0], ppm[1][0][0], ppm[2][0][0]);
  foveon_fixed (satlev, 3, "SaturationLevel");
  printf("satlev: %d %d %d\n", satlev[0], satlev[1], satlev[2]);
  foveon_fixed (keep, 4, "KeepImageArea");
  printf("keep: %d %d %d\n", keep[0], keep[1], keep[2], keep[3]);
  foveon_fixed (active, 4, "ActiveImageArea");
  printf("active: %d %d %d\n", active[0], active[1], active[2], active[3]);
  foveon_fixed (chroma_dq, 3, "ChromaDQ");
  printf("chroma_dq: %f %f %f\n", chroma_dq[0], chroma_dq[1], chroma_dq[2]);
  foveon_fixed (color_dq, 3,
	foveon_camf_param ("IncludeBlocks", "ColorDQ") ?
		"ColorDQ" : "ColorDQCamRGB");
  printf("color_dq: %f %f %f\n", color_dq[0], color_dq[1], color_dq[2]);
  if (foveon_camf_param ("IncludeBlocks", "ColumnFilter"))
		 foveon_fixed (&cfilt, 1, "ColumnFilter");
  printf("cfilt: %f\n", cfilt);

  memset (ddft, 0, sizeof ddft);
  if (!foveon_camf_param ("IncludeBlocks", "DarkDrift")
	 || !foveon_fixed (ddft[1][0], 12, "DarkDrift"))
    for (i=0; i < 2; i++) {
      foveon_fixed (dstb, 4, i ? "DarkShieldBottom":"DarkShieldTop");
      for (row = dstb[1]; row <= dstb[3]; row++)
	for (col = dstb[0]; col <= dstb[2]; col++)
	  FORC3 ddft[i+1][c][1] += (short) image[row*width+col][c];
      FORC3 ddft[i+1][c][1] /= (dstb[3]-dstb[1]+1) * (dstb[2]-dstb[0]+1);
    }
  printf("ddft: %f %f %f\n", ddft[1][0][1], ddft[1][0][2], ddft[1][1][1]);

  if (!(cp = foveon_camf_param ("WhiteBalanceIlluminants", model2)))
  { fprintf (stderr,_("%s: Invalid white balance \"%s\"\n"), ifname, model2);
    return; }
  foveon_fixed (cam_xyz, 9, cp);
  printf("cam_xyz: %f %f %f\n", cam_xyz[0][0], cam_xyz[0][1], cam_xyz[0][2]);
  if (foveon_camf_param ("IncludeBlocks", "WhiteBalanceCorrections"))
	foveon_fixed (correct, 9,
				  foveon_camf_param ("WhiteBalanceCorrections", model2));
  printf("correct: %f %f %f\n", correct[0][0], correct[0][1],
  correct[0][2]);

  memset (last, 0, sizeof last);
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      FORC3 last[i][j] += correct[i][c] * cam_xyz[c][j];

  #define LAST(x,y) last[(i+x)%3][(c+y)%3]
  for (i=0; i < 3; i++)
    FORC3 diag[c][i] = LAST(1,1)*LAST(2,2) - LAST(1,2)*LAST(2,1);
  #undef LAST
  FORC3 div[c] = diag[c][0]*0.3127 + diag[c][1]*0.329 + diag[c][2]*0.3583;
  FORC3 printf("div[%d]=%f\n", c, div[c]);

  sprintf (str, "%sRGBNeutral", model2);
  if (foveon_camf_param ("IncludeBlocks", str))
    foveon_fixed (div, 3, str);
  else {
	sprintf(str, "%sWBGains", model2);
	if (foveon_camf_param ("IncludeBlocks", "WhiteBalanceGains"))
	  foveon_fixed (div, 3, foveon_camf_param ("WhiteBalanceGains", model2));
  }
  FORC3 printf("div[%d]=%f\n", c, div[c]);

  num = 0;
  FORC3 if (num < div[c]) num = div[c];
  FORC3 div[c] /= num;

  memset (trans, 0, sizeof trans);
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      FORC3 { trans[i][j] += rgb_cam[i][c] * last[c][j] * div[j];
	  printf("rgb_cam[%d][%d]: %f\n", i, c, rgb_cam[i][c]);
	}
  FORC3 trsum[c] = trans[c][0] + trans[c][1] + trans[c][2];
  dsum = (6*trsum[0] + 11*trsum[1] + 3*trsum[2]) / 20;
  for (i=0; i < 3; i++)
    FORC3 last[i][c] = trans[i][c] * dsum / trsum[i];
  memset (trans, 0, sizeof trans);
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++){
      FORC3 trans[i][j] += (i==c ? 32 : -1) * last[c][j] / 30;
	  printf("trans[%d][%d]: %f\n", i, j, trans[i][j]);
	}

  foveon_make_curves (curve, color_dq, div, cfilt);
  FORC3 chroma_dq[c] /= 3;
  foveon_make_curves (curve+3, chroma_dq, div, cfilt);
  FORC3 dsum += chroma_dq[c] / div[c];
  curve[6] = foveon_make_curve (dsum, dsum, cfilt);
  curve[7] = foveon_make_curve (dsum*2, dsum*2, cfilt);

  black = (float (*)[3]) calloc (height, sizeof *black);
  for (row=0; row < height; row++) {
    for (i=0; i < 6; i++)
      ddft[0][0][i] = ddft[1][0][i] +
	row / (height-1.0) * (ddft[2][0][i] - ddft[1][0][i]);
    FORC3 black[row][c] =
	( foveon_avg (image[row*width]+c, dscr[0], cfilt) +
	  foveon_avg (image[row*width]+c, dscr[1], cfilt) * 3
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

  sgain = (float (*)[3]) foveon_camf_matrix (dim, "SpatialGain");
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

  if ((badpix = (unsigned int *) foveon_camf_matrix (dim, "BadPixels"))) {
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
  merror (smrow[6], "foveon_interpolate()");
  for (i=0; i < 5; i++)
    smrow[i] = smrow[6] + i*width;

  /* Sharpen the reds against these Gaussian averages */
  for (smlast=-1, row=2; row < height-2; row++) {
    while (smlast < row+2) {
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

  /* Adjust the brighter pixels for better linearity */
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
    while (smlast < row+2) {
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
      FORC3 dev[c] = -foveon_apply_curve (curve[7], pix[c] -
	((smrow[1][col][c] + 2*smrow[2][col][c] + smrow[3][col][c]) >> 2));
      sum = (dev[0] + dev[1] + dev[2]) >> 3;
      FORC3 pix[c] += dev[c] - sum;
      pix += 4;
    }
  }
  for (smlast=-1, row=2; row < height-2; row++) {
    while (smlast < row+2) {
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
      FORC3 pix[c] += foveon_apply_curve (curve[6],
		((j*total[c] + 0x8000) >> 16) - pix[c]);
      pix += 4;
    }
  }

  /* Transform the image to a different colorspace */
  for (pix=image[0]; pix < image[height*width]; pix+=4) {
    FORC3 pix[c] -= foveon_apply_curve (curve[c], pix[c]);
    sum = (pix[0]+pix[1]+pix[1]+pix[2]) >> 2;
    FORC3 pix[c] -= foveon_apply_curve (curve[c], pix[c]-sum);
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
  merror (shrink, "foveon_interpolate()");
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
	ipix[c] = foveon_apply_curve (curve[c+3],
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
  width = i;
  height = row;
}
#undef image

/* RESTRICTED code ends here */

void CLASS gamma_curve (double pwr, double ts, int mode, int imax)
{
  int i, c;
  double g[6], bnd[2]={0,0}, r;

  printf("gamma_curve\n");
  FORC(6) printf("gamm[%d]=%f\n",c, gamm[c]);
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
  printf("mode: %d\n", mode);
  if (!mode--) {
    memcpy (gamm, g, sizeof gamm);
    return;
  }
  for (i=0; i < 0x10000; i++) {
    curve[i] = 0xffff;
    if ((r = (double) i / imax) < 1)
      curve[i] = 0x10000 * ( mode
	? (r < g[3] ? r*g[1] : (g[0] ? pow( r,g[0])*(1+g[4])-g[4]    : log(r)*g[2]+1))
	: (r < g[2] ? r/g[1] : (g[0] ? pow((r+g[4])/(1+g[4]),1/g[0]) : exp((r-1)/g[2]))));
  }
}

void CLASS pseudoinverse (double (*in)[3], double (*out)[3], int size)
{
  double work[3][6], num;
  int i, j, k;

  for (i=0; i < 3; i++) {
    for (j=0; j < 6; j++)
      work[i][j] = j == i+3;
    for (j=0; j < 3; j++)
      for (k=0; k < size; k++)
	work[i][j] += in[k][i] * in[k][j];
  }
  for (i=0; i < 3; i++) {
    num = work[i][i];
    for (j=0; j < 6; j++)
      work[i][j] /= num;
    for (k=0; k < 3; k++) {
      if (k==i) continue;
      num = work[k][i];
      for (j=0; j < 6; j++)
	work[k][j] -= work[i][j] * num;
    }
  }
  for (i=0; i < size; i++)
    for (j=0; j < 3; j++)
      for (out[i][j]=k=0; k < 3; k++)
	out[i][j] += work[j][k+3] * in[i][k];
}

void CLASS cam_xyz_coeff (float rgb_cam[3][4], double cam_xyz[4][3])
{
  double cam_rgb[4][3], inverse[4][3], num;
  int i, j, k;

  for (i=0; i < colors; i++)		/* Multiply out XYZ colorspace */
    for (j=0; j < 3; j++)
      for (cam_rgb[i][j] = k=0; k < 3; k++)
	cam_rgb[i][j] += cam_xyz[i][k] * xyz_rgb[k][j];

  for (i=0; i < colors; i++) {		/* Normalize cam_rgb so that */
    for (num=j=0; j < 3; j++)		/* cam_rgb * (1,1,1) is (1,1,1,1) */
      num += cam_rgb[i][j];
    for (j=0; j < 3; j++)
      cam_rgb[i][j] /= num;
    pre_mul[i] = 1 / num;
  }
  pseudoinverse (cam_rgb, inverse, colors);
  for (i=0; i < 3; i++)
    for (j=0; j < colors; j++)
      rgb_cam[i][j] = inverse[j][i];
}

int CLASS parse_jpeg (int offset)
{
  int len, save, hlen, mark;

  fseek (ifp, offset, SEEK_SET);
  if (fgetc(ifp) != 0xff || fgetc(ifp) != 0xd8) return 0;

  while (fgetc(ifp) == 0xff && (mark = fgetc(ifp)) != 0xda) {
    order = 0x4d4d;
    len   = get2() - 2;
    save  = ftell(ifp);
    if (mark == 0xc0 || mark == 0xc3) {
      fgetc(ifp);
      raw_height = get2();
      raw_width  = get2();
    }
    order = get2();
    hlen  = get4();
    fseek (ifp, save+len, SEEK_SET);
  }
  return 1;
}

char * CLASS foveon_gets (int offset, char *str, int len)
{
  int i;
  fseek (ifp, offset, SEEK_SET);
  for (i=0; i < len-1; i++)
    if ((str[i] = get2()) == 0) break;
  str[i] = 0;
  return str;
}

void CLASS parse_foveon()
{
  int entries, img=0, off, len, tag, save, i, wide, high, pent, poff[256][2];
  char name[64], value[64];

  order = 0x4949;			/* Little-endian */
  fseek (ifp, 36, SEEK_SET);
  flip = get4();
  fseek (ifp, -4, SEEK_END);
  fseek (ifp, get4(), SEEK_SET);
  if (get4() != 0x64434553) return;	/* SECd */
  entries = (get4(),get4());
  while (entries--) {
    off = get4();
    len = get4();
    tag = get4();
    save = ftell(ifp);
    fseek (ifp, off, SEEK_SET);
    if (get4() != (0x20434553 | (tag << 24))) return;
    switch (tag) {
      case 0x47414d49:			/* IMAG */
      case 0x32414d49:			/* IMA2 */
	fseek (ifp, 8, SEEK_CUR);
	pent = get4();
	wide = get4();
	high = get4();
	if (wide > raw_width && high > raw_height) {
	  switch (pent) {
	    case  5:  load_flags = 1;
	    case  6:  load_raw = &CLASS foveon_sd_load_raw;  break;
	    case 30:  load_raw = &CLASS foveon_dp_load_raw;  break;
	    default:  load_raw = 0;
	  }
	  raw_width  = wide;
	  raw_height = high;
	  data_offset = off+28;
	  is_foveon = 1;
	}
	fseek (ifp, off+28, SEEK_SET);
	if (fgetc(ifp) == 0xff && fgetc(ifp) == 0xd8
		&& thumb_length < len-28) {
	  thumb_offset = off+28;
	  thumb_length = len-28;
	  write_thumb = &CLASS jpeg_thumb;
	}
	if (++img == 2 && !thumb_length) {
	  thumb_offset = off+24;
	  thumb_width = wide;
	  thumb_height = high;
	  write_thumb = &CLASS foveon_thumb;
	}
	break;
      case 0x464d4143:			/* CAMF */
	meta_offset = off+8;
	meta_length = len-28;
	break;
      case 0x504f5250:			/* PROP */
	pent = (get4(),get4());
	fseek (ifp, 12, SEEK_CUR);
	off += pent*8 + 24;
	if ((unsigned) pent > 256) pent=256;
	for (i=0; i < pent*2; i++)
	  poff[0][i] = off + get4()*2;
	for (i=0; i < pent; i++) {
	  foveon_gets (poff[i][0], name, 64);
	  foveon_gets (poff[i][1], value, 64);
	  if (!strcmp (name, "ISO"))
	    iso_speed = atoi(value);
	  if (!strcmp (name, "CAMMANUF"))
	    strcpy (make, value);
	  if (!strcmp (name, "CAMMODEL"))
	    strcpy (model, value);
	  if (!strcmp (name, "WB_DESC"))
	    strcpy (model2, value);
	  if (!strcmp (name, "TIME"))
	    timestamp = atoi(value);
	  if (!strcmp (name, "EXPTIME"))
	    shutter = atoi(value) / 1000000.0;
	  if (!strcmp (name, "APERTURE"))
	    aperture = atof(value);
	  if (!strcmp (name, "FLENGTH"))
	    focal_len = atof(value);
	  if (!strcmp (name, "CM_DESC"))
		strcpy (cm_desc, value);
	}
#ifdef LOCALTIME
	timestamp = mktime (gmtime (&timestamp));
#endif
    }
    fseek (ifp, save, SEEK_SET);
  }
}


void CLASS simple_coeff (int index)
{
  static const float table[][12] = {
  /* index 0 -- all Foveon cameras */
  { 1.4032,-0.2231,-0.1016,-0.5263,1.4816,0.017,-0.0112,0.0183,0.9113 },
  { 1.5032,-0.2231,-0.1016,-0.0263,1.5816,0.017,-0.0112,0.0183,0.9113 },
  };
  int i, c;

  for (raw_color = i=0; i < 3; i++)
    FORCC rgb_cam[i][c] = table[index][i*colors+c];
}

short CLASS guess_byte_order (int words)
{
  uchar test[4][2];
  int t=2, msb;
  double diff, sum[2] = {0,0};

  fread (test[0], 2, 2, ifp);
  for (words-=2; words--; ) {
    fread (test[t], 2, 1, ifp);
    for (msb=0; msb < 2; msb++) {
      diff = (test[t^2][msb] << 8 | test[t^2][!msb])
	   - (test[t  ][msb] << 8 | test[t  ][!msb]);
      sum[msb] += diff*diff;
    }
    t = (t+1) & 3;
  }
  return sum[0] < sum[1] ? 0x4d4d : 0x4949;
}

void CLASS stretch()
{
  ushort newdim, (*img)[4], *pix0, *pix1;
  int row, col, c;
  double rc, frac;

  if (pixel_aspect == 1) return;
  if (verbose) fprintf (stderr,_("Stretching the image...\n"));
  if (pixel_aspect < 1) {
    newdim = height / pixel_aspect + 0.5;
    img = (ushort (*)[4]) calloc (width, newdim*sizeof *img);
    merror (img, "stretch()");
    for (rc=row=0; row < newdim; row++, rc+=pixel_aspect) {
      frac = rc - (c = rc);
      pix0 = pix1 = image[c*width];
      if (c+1 < height) pix1 += width*4;
      for (col=0; col < width; col++, pix0+=4, pix1+=4)
	FORCC img[row*width+col][c] = pix0[c]*(1-frac) + pix1[c]*frac + 0.5;
    }
    height = newdim;
  } else {
    newdim = width * pixel_aspect + 0.5;
    img = (ushort (*)[4]) calloc (height, newdim*sizeof *img);
    merror (img, "stretch()");
    for (rc=col=0; col < newdim; col++, rc+=1/pixel_aspect) {
      frac = rc - (c = rc);
      pix0 = pix1 = image[c];
      if (c+1 < width) pix1 += 4;
      for (row=0; row < height; row++, pix0+=width*4, pix1+=width*4)
	FORCC img[row*newdim+col][c] = pix0[c]*(1-frac) + pix1[c]*frac + 0.5;
    }
    width = newdim;
  }
  free (image);
  image = img;
}

/*
   Identify which camera created this file, and set global variables
   accordingly.
 */
void CLASS identify()
{

  static const char *corp[] =
    {  "Sigma" };
  char head[32], *cp;
  int hlen, flen, fsize, zero_fsize=1, i, c;
  struct jhead jh;

  tiff_flip = flip = filters = UINT_MAX;	/* unknown */
  raw_height = raw_width = fuji_width = fuji_layout = cr2_slice[0] = 0;
  maximum = height = width = top_margin = left_margin = 0;
  cdesc[0] = desc[0] = artist[0] = make[0] = model[0] = model2[0] = 0;
  iso_speed = shutter = aperture = focal_len = unique_id = 0;
  tiff_nifds = 0;
  memset (tiff_ifd, 0, sizeof tiff_ifd);
  memset (gpsdata, 0, sizeof gpsdata);
  memset (cblack, 0, sizeof cblack);
  memset (white, 0, sizeof white);
  memset (mask, 0, sizeof mask);
  thumb_offset = thumb_length = thumb_width = thumb_height = 0;
  load_raw = thumb_load_raw = 0;
  write_thumb = &CLASS jpeg_thumb;
  data_offset = meta_length = tiff_bps = tiff_compress = 0;
  kodak_cbpp = zero_after_ff = dng_version = load_flags = 0;
  timestamp = shot_order = tiff_samples = black = is_foveon = 0;
  mix_green = profile_length = data_error = zero_is_bad = 0;
  pixel_aspect = is_raw = raw_color = 1;
  tile_width = tile_length = 0;
  for (i=0; i < 4; i++) {
    cam_mul[i] = i == 1;
    pre_mul[i] = i < 3;
    FORC3 cmatrix[c][i] = 0;
    FORC3 rgb_cam[c][i] = c == i;
  }
  colors = 3;
  for (i=0; i < 0x10000; i++) curve[i] = i;

  order = get2();
  hlen = get4();
  fseek (ifp, 0, SEEK_SET);
  fread (head, 1, 32, ifp);
  fseek (ifp, 0, SEEK_END);
  flen = fsize = ftell(ifp);
  if (!memcmp (head,"FOVb",4))
    parse_foveon();

  if (zero_fsize) fsize = 0;

  for (i=0; i < sizeof corp / sizeof *corp; i++)
    if (strcasestr (make, corp[i]))	/* Simplify company names */
	    strcpy (make, corp[i]);
  cp = make + strlen(make);		/* Remove trailing spaces */
  while (*--cp == ' ') *cp = 0;
  cp = model + strlen(model);
  while (*--cp == ' ') *cp = 0;
  i = strlen(make);			/* Remove make from model */
  if (!strncasecmp (model, make, i) && model[i++] == ' ')
    memmove (model, model+i, 64-i);
  desc[511] = artist[63] = make[63] = model[63] = model2[63] = 0;
  if (!is_raw) goto notraw;

  if (!height) height = raw_height;
  if (!width)  width  = raw_width;

/* Set parameters based on camera name (for non-DNG files). */

  if (is_foveon) {
    if (height*2 < width) pixel_aspect = 0.5;
    if (height   > width) pixel_aspect = 2;
    filters = 0;
/* 	if (load_raw == &CLASS foveon_dp_load_raw) */
/* 	  simple_coeff(0); */
/* 	else */
	  simple_coeff(0);
	stretch();
  }
  if (!model[0])
    sprintf (model, "%dx%d", width, height);
  if (filters == UINT_MAX) filters = 0x94949494;
  if (thumb_offset && !thumb_height) {
    fseek (ifp, thumb_offset, SEEK_SET);
    if (ljpeg_start (&jh, 1)) {
      thumb_width  = jh.wide;
      thumb_height = jh.high;
    }
  }

notraw:
  if (flip == UINT_MAX) flip = tiff_flip;
  if (flip == UINT_MAX) flip = 0;
}

#ifndef NO_LCMS
void CLASS apply_profile (const char *input, const char *output)
{
  char *prof;
  cmsHPROFILE hInProfile=0, hOutProfile=0;
  cmsHTRANSFORM hTransform;
  FILE *fp;
  unsigned size;

  if (strcmp (input, "embed"))
    hInProfile = cmsOpenProfileFromFile (input, "r");
  else if (profile_length) {
    prof = (char *) malloc (profile_length);
    merror (prof, "apply_profile()");
    fseek (ifp, profile_offset, SEEK_SET);
    fread (prof, 1, profile_length, ifp);
    hInProfile = cmsOpenProfileFromMem (prof, profile_length);
    free (prof);
  } else
    fprintf (stderr,_("%s has no embedded profile.\n"), ifname);
  if (!hInProfile) return;
  if (!output)
    hOutProfile = cmsCreate_sRGBProfile();
  else if ((fp = fopen (output, "rb"))) {
    fread (&size, 4, 1, fp);
    fseek (fp, 0, SEEK_SET);
    oprof = (unsigned *) malloc (size = ntohl(size));
    merror (oprof, "apply_profile()");
    fread (oprof, 1, size, fp);
    fclose (fp);
    if (!(hOutProfile = cmsOpenProfileFromMem (oprof, size))) {
      free (oprof);
      oprof = 0;
    }
  } else
    fprintf (stderr,_("Cannot open file %s!\n"), output);
  if (!hOutProfile) goto quit;
  if (verbose)
    fprintf (stderr,_("Applying color profile...\n"));
  hTransform = cmsCreateTransform (hInProfile, TYPE_RGBA_16,
	hOutProfile, TYPE_RGBA_16, INTENT_PERCEPTUAL, 0);
  cmsDoTransform (hTransform, image, image, width*height);
  raw_color = 1;		/* Don't use rgb_cam with a profile */
  cmsDeleteTransform (hTransform);
  cmsCloseProfile (hOutProfile);
quit:
  cmsCloseProfile (hInProfile);
}
#endif

void CLASS convert_to_rgb()
{
  int row, col, c, i, j, k;
  ushort *img;
  float out[3], out_cam[3][4];
  double num, inverse[3][3];
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
  static const double (*out_rgb[])[3] =
  { rgb_rgb, adobe_rgb, wide_rgb, prophoto_rgb, xyz_rgb };
  static const char *name[] =
  { "sRGB", "Adobe RGB (1998)", "WideGamut D65", "ProPhoto D65", "XYZ" };
  static const unsigned phead[] =
  { 1024, 0, 0x2100000, 0x6d6e7472, 0x52474220, 0x58595a20, 0, 0, 0,
    0x61637370, 0, 0, 0x6e6f6e65, 0, 0, 0, 0, 0xf6d6, 0x10000, 0xd32d };
  unsigned pbody[] =
  { 10, 0x63707274, 0, 36,	/* cprt */
	0x64657363, 0, 40,	/* desc */
	0x77747074, 0, 20,	/* wtpt */
	0x626b7074, 0, 20,	/* bkpt */
	0x72545243, 0, 14,	/* rTRC */
	0x67545243, 0, 14,	/* gTRC */
	0x62545243, 0, 14,	/* bTRC */
	0x7258595a, 0, 20,	/* rXYZ */
	0x6758595a, 0, 20,	/* gXYZ */
	0x6258595a, 0, 20 };	/* bXYZ */
  static const unsigned pwhite[] = { 0xf351, 0x10000, 0x116cc };
  unsigned pcurve[] = { 0x63757276, 0, 1, 0x1000000 };

  gamma_curve (gamm[0], gamm[1], 0, 0);
  memcpy (out_cam, rgb_cam, sizeof out_cam);
  raw_color |= colors == 1 || document_mode ||
		output_color < 1 || output_color > 5;
  if (!raw_color) {
    oprof = (unsigned *) calloc (phead[0], 1);
    merror (oprof, "convert_to_rgb()");
    memcpy (oprof, phead, sizeof phead);
    if (output_color == 5) oprof[4] = oprof[5];
    oprof[0] = 132 + 12*pbody[0];
    for (i=0; i < pbody[0]; i++) {
      oprof[oprof[0]/4] = i ? (i > 1 ? 0x58595a20 : 0x64657363) : 0x74657874;
      pbody[i*3+2] = oprof[0];
      oprof[0] += (pbody[i*3+3] + 3) & -4;
    }
    memcpy (oprof+32, pbody, sizeof pbody);
    oprof[pbody[5]/4+2] = strlen(name[output_color-1]) + 1;
    memcpy ((char *)oprof+pbody[8]+8, pwhite, sizeof pwhite);
    pcurve[3] = (short)(256/gamm[5]+0.5) << 16;
    for (i=4; i < 7; i++)
      memcpy ((char *)oprof+pbody[i*3+2], pcurve, sizeof pcurve);
    pseudoinverse ((double (*)[3]) out_rgb[output_color-1], inverse, 3);
    for (i=0; i < 3; i++)
      for (j=0; j < 3; j++) {
	for (num = k=0; k < 3; k++)
	  num += xyzd50_srgb[i][k] * inverse[j][k];
	oprof[pbody[j*3+23]/4+i+2] = num * 0x10000 + 0.5;
      }
    for (i=0; i < phead[0]/4; i++)
      oprof[i] = htonl(oprof[i]);
    strcpy ((char *)oprof+pbody[2]+8, "auto-generated by dcraw");
    strcpy ((char *)oprof+pbody[5]+12, name[output_color-1]);
    for (i=0; i < 3; i++)
      for (j=0; j < colors; j++)
	for (out_cam[i][j] = k=0; k < 3; k++)
	  out_cam[i][j] += out_rgb[output_color-1][i][k] * rgb_cam[k][j];
  }
  if (verbose)
    fprintf (stderr, raw_color ? _("Building histograms...\n") :
	_("Converting to %s colorspace...\n"), name[output_color-1]);

  memset (histogram, 0, sizeof histogram);
  for (img=image[0], row=0; row < height; row++)
    for (col=0; col < width; col++, img+=4) {
      if (!raw_color) {
	out[0] = out[1] = out[2] = 0;
	FORCC {
	  out[0] += out_cam[0][c] * img[c];
	  out[1] += out_cam[1][c] * img[c];
	  out[2] += out_cam[2][c] * img[c];
	}
	FORC3 img[c] = CLIP((int) out[c]);
      }
/*       else if (document_mode) */
/* 	img[0] = img[fcol(row,col)]; */
      FORCC histogram[c][img[c] >> 3]++;
    }
  if (colors == 4 && output_color) colors = 3;
  if (document_mode && filters) colors = 1;
}


int CLASS flip_index (int row, int col)
{
  if (flip & 4) SWAP(row,col);
  if (flip & 2) row = iheight - 1 - row;
  if (flip & 1) col = iwidth  - 1 - col;
  return row * iwidth + col;
}

struct tiff_tag {
  ushort tag, type;
  int count;
  union { char c[4]; short s[2]; int i; } val;
};

struct tiff_hdr {
  ushort order, magic;
  int ifd;
  ushort pad, ntag;
  struct tiff_tag tag[23];
  int nextifd;
  ushort pad2, nexif;
  struct tiff_tag exif[4];
  ushort pad3, ngps;
  struct tiff_tag gpst[10];
  short bps[4];
  int rat[10];
  unsigned gps[26];
  char desc[512], make[64], model[64], soft[32], date[20], artist[64];
};

void CLASS tiff_set (ushort *ntag,
	ushort tag, ushort type, int count, int val)
{
  struct tiff_tag *tt;
  int c;

  tt = (struct tiff_tag *)(ntag+1) + (*ntag)++;
  tt->tag = tag;
  tt->type = type;
  tt->count = count;
  if (type < 3 && count <= 4)
    FORC(4) tt->val.c[c] = val >> (c << 3);
  else if (type == 3 && count <= 2)
    FORC(2) tt->val.s[c] = val >> (c << 4);
  else tt->val.i = val;
}

#define TOFF(ptr) ((char *)(&(ptr)) - (char *)th)

void CLASS tiff_head (struct tiff_hdr *th, int full)
{
  int c, psize=0;
  struct tm *t;

  memset (th, 0, sizeof *th);
  th->order = htonl(0x4d4d4949) >> 16;
  th->magic = 42;
  th->ifd = 10;
  if (full) {
    tiff_set (&th->ntag, 254, 4, 1, 0);
    tiff_set (&th->ntag, 256, 4, 1, width);
    tiff_set (&th->ntag, 257, 4, 1, height);
    tiff_set (&th->ntag, 258, 3, colors, output_bps);
    if (colors > 2)
      th->tag[th->ntag-1].val.i = TOFF(th->bps);
    FORC4 th->bps[c] = output_bps;
    tiff_set (&th->ntag, 259, 3, 1, 1);
    tiff_set (&th->ntag, 262, 3, 1, 1 + (colors > 1));
  }
  tiff_set (&th->ntag, 270, 2, 512, TOFF(th->desc));
  tiff_set (&th->ntag, 271, 2, 64, TOFF(th->make));
  tiff_set (&th->ntag, 272, 2, 64, TOFF(th->model));
  if (full) {
    if (oprof) psize = ntohl(oprof[0]);
    tiff_set (&th->ntag, 273, 4, 1, sizeof *th + psize);
    tiff_set (&th->ntag, 277, 3, 1, colors);
    tiff_set (&th->ntag, 278, 4, 1, height);
    tiff_set (&th->ntag, 279, 4, 1, height*width*colors*output_bps/8);
  } else
    tiff_set (&th->ntag, 274, 3, 1, "12435867"[flip]-'0');
  tiff_set (&th->ntag, 282, 5, 1, TOFF(th->rat[0]));
  tiff_set (&th->ntag, 283, 5, 1, TOFF(th->rat[2]));
  tiff_set (&th->ntag, 284, 3, 1, 1);
  tiff_set (&th->ntag, 296, 3, 1, 2);
  tiff_set (&th->ntag, 305, 2, 32, TOFF(th->soft));
  tiff_set (&th->ntag, 306, 2, 20, TOFF(th->date));
  tiff_set (&th->ntag, 315, 2, 64, TOFF(th->artist));
  tiff_set (&th->ntag, 34665, 4, 1, TOFF(th->nexif));
  if (psize) tiff_set (&th->ntag, 34675, 7, psize, sizeof *th);
  tiff_set (&th->nexif, 33434, 5, 1, TOFF(th->rat[4]));
  tiff_set (&th->nexif, 33437, 5, 1, TOFF(th->rat[6]));
  tiff_set (&th->nexif, 34855, 3, 1, iso_speed);
  tiff_set (&th->nexif, 37386, 5, 1, TOFF(th->rat[8]));
  if (gpsdata[1]) {
    tiff_set (&th->ntag, 34853, 4, 1, TOFF(th->ngps));
    tiff_set (&th->ngps,  0, 1,  4, 0x202);
    tiff_set (&th->ngps,  1, 2,  2, gpsdata[29]);
    tiff_set (&th->ngps,  2, 5,  3, TOFF(th->gps[0]));
    tiff_set (&th->ngps,  3, 2,  2, gpsdata[30]);
    tiff_set (&th->ngps,  4, 5,  3, TOFF(th->gps[6]));
    tiff_set (&th->ngps,  5, 1,  1, gpsdata[31]);
    tiff_set (&th->ngps,  6, 5,  1, TOFF(th->gps[18]));
    tiff_set (&th->ngps,  7, 5,  3, TOFF(th->gps[12]));
    tiff_set (&th->ngps, 18, 2, 12, TOFF(th->gps[20]));
    tiff_set (&th->ngps, 29, 2, 12, TOFF(th->gps[23]));
    memcpy (th->gps, gpsdata, sizeof th->gps);
  }
  th->rat[0] = th->rat[2] = 300;
  th->rat[1] = th->rat[3] = 1;
  FORC(6) th->rat[4+c] = 1000000;
  th->rat[4] *= shutter;
  th->rat[6] *= aperture;
  th->rat[8] *= focal_len;
  strncpy (th->desc, desc, 512);
  strncpy (th->make, make, 64);
  strncpy (th->model, model, 64);
  strcpy (th->soft, "dcraw v"DCRAW_VERSION);
  t = localtime (&timestamp);
  sprintf (th->date, "%04d:%02d:%02d %02d:%02d:%02d",
      t->tm_year+1900,t->tm_mon+1,t->tm_mday,t->tm_hour,t->tm_min,t->tm_sec);
  strncpy (th->artist, artist, 64);
}

void CLASS jpeg_thumb()
{
  char *thumb;
  ushort exif[5];
  struct tiff_hdr th;

  thumb = (char *) malloc (thumb_length);
  merror (thumb, "jpeg_thumb()");
  fread (thumb, 1, thumb_length, ifp);
  fputc (0xff, ofp);
  fputc (0xd8, ofp);
  if (strcmp (thumb+6, "Exif")) {
    memcpy (exif, "\xff\xe1  Exif\0\0", 10);
    exif[1] = htons (8 + sizeof th);
    fwrite (exif, 1, sizeof exif, ofp);
    tiff_head (&th, 0);
    fwrite (&th, 1, sizeof th, ofp);
  }
  fwrite (thumb+2, 1, thumb_length-2, ofp);
  free (thumb);
}

void CLASS write_ppm_tiff()
{
  struct tiff_hdr th;
  uchar *ppm;
  ushort *ppm2;
  int c, row, col, soff, rstep, cstep;
  int perc, val, total, white=0x2000;

  perc = width * height * 0.01;		/* 99th percentile white level */
/*   if (fuji_width) perc /= 2; */
  if (!((highlight & ~2) || no_auto_bright))
    for (white=c=0; c < colors; c++) {
      for (val=0x2000, total=0; --val > 32; )
	if ((total += histogram[c][val]) > perc) break;
      if (white < val) white = val;
    }
  gamma_curve (gamm[0], gamm[1], 2, (white << 3)/bright);
  iheight = height;
  iwidth  = width;
  if (flip & 4) SWAP(height,width);
  ppm = (uchar *) calloc (width, colors*output_bps/8);
  ppm2 = (ushort *) ppm;
  merror (ppm, "write_ppm_tiff()");
  printf("colors: %d\n", colors);
  if (output_tiff) {
    tiff_head (&th, 1);
    fwrite (&th, sizeof th, 1, ofp);
    if (oprof)
      fwrite (oprof, ntohl(oprof[0]), 1, ofp);
  } else if (colors > 3)
    fprintf (ofp,
      "P7\nWIDTH %d\nHEIGHT %d\nDEPTH %d\nMAXVAL %d\nTUPLTYPE %s\nENDHDR\n",
	width, height, colors, (1 << output_bps)-1, cdesc);
  else
    fprintf (ofp, "P%d\n%d %d\n%d\n",
	colors/2+5, width, height, (1 << output_bps)-1);
  printf("flip: %d\n", flip);
  soff  = flip_index (0, 0);
  cstep = flip_index (0, 1) - soff;
  rstep = flip_index (1, 0) - flip_index (0, width);
  for (row=0; row < height; row++, soff += rstep) {
    for (col=0; col < width; col++, soff += cstep)
      if (output_bps == 8)
		FORCC ppm [col*colors+c] =/*  image[row*width+col][c] */curve[image[soff][c]]>> 8;
      else FORCC ppm2[col*colors+c] = /* image[row*width+col][c] */curve[image[soff][c]];
    if (output_bps == 16 && !output_tiff && htons(0x55aa) != 0x55aa)
      swab (ppm2, ppm2, width*colors*2);
    fwrite (ppm, colors*output_bps/8, width, ofp);
  }
  free (ppm);
}

int CLASS main (int argc, const char **argv)
{
  int arg, status=0, quality, i, c;
  int timestamp_only=0, thumbnail_only=0, identify_only=0;
  int user_qual=-1, user_black=-1, user_sat=-1, user_flip=-1;
  int write_to_stdout=0, read_from_stdin=0;
  const char *sp, *bpfile=0, *dark_frame=0, *write_ext;
  char opm, opt, *ofname, *cp;
  struct utimbuf ut;
#ifndef NO_LCMS
  const char *cam_profile=0, *out_profile=0;
#endif

#ifndef LOCALTIME
  putenv ((char *) "TZ=UTC");
#endif
#ifdef LOCALEDIR
  setlocale (LC_CTYPE, "");
  setlocale (LC_MESSAGES, "");
  bindtextdomain ("dcraw", LOCALEDIR);
  textdomain ("dcraw");
#endif

  if (argc == 1) {
    printf(_("\nRaw photo decoder \"dcraw\" v%s"), DCRAW_VERSION);
    printf(_("\nby Dave Coffin, dcoffin a cybercom o net\n"));
    printf(_("\nUsage:  %s [OPTION]... [FILE]...\n\n"), argv[0]);
    puts(_("-v        Print verbose messages"));
    puts(_("-c        Write image data to standard output"));
    puts(_("-e        Extract embedded thumbnail image"));
    puts(_("-i        Identify files without decoding them"));
    puts(_("-i -v     Identify files and show metadata"));
    puts(_("-z        Change file dates to camera timestamp"));
    puts(_("-w        Use camera white balance, if possible"));
    puts(_("-a        Average the whole image for white balance"));
    puts(_("-A <x y w h> Average a grey box for white balance"));
    puts(_("-r <r g b g> Set custom white balance"));
    puts(_("+M/-M     Use/don't use an embedded color matrix"));
    puts(_("-C <r b>  Correct chromatic aberration"));
    puts(_("-P <file> Fix the dead pixels listed in this file"));
    puts(_("-K <file> Subtract dark frame (16-bit raw PGM)"));
    puts(_("-k <num>  Set the darkness level"));
    puts(_("-S <num>  Set the saturation level"));
    puts(_("-n <num>  Set threshold for wavelet denoising"));
    puts(_("-H [0-9]  Highlight mode (0=clip, 1=unclip, 2=blend, 3+=rebuild)"));
    puts(_("-t [0-7]  Flip image (0=none, 3=180, 5=90CCW, 6=90CW)"));
    puts(_("-o [0-5]  Output colorspace (raw,sRGB,Adobe,Wide,ProPhoto,XYZ)"));
#ifndef NO_LCMS
    puts(_("-o <file> Apply output ICC profile from file"));
    puts(_("-p <file> Apply camera ICC profile from file or \"embed\""));
#endif
    puts(_("-d        Document mode (no color, no interpolation)"));
    puts(_("-D        Document mode without scaling (totally raw)"));
    puts(_("-W        Don't automatically brighten the image"));
    puts(_("-b <num>  Adjust brightness (default = 1.0)"));
    puts(_("-g <p ts> Set custom gamma curve (default = 2.222 4.5)"));
    puts(_("-q [0-3]  Set the interpolation quality"));
    puts(_("-h        Half-size color image (twice as fast as \"-q 0\")"));
    puts(_("-f        Interpolate RGGB as four colors"));
    puts(_("-m <num>  Apply a 3x3 median filter to R-G and B-G"));
    puts(_("-s [0..N-1] Select one raw image or \"all\" from each file"));
    puts(_("-6        Write 16-bit instead of 8-bit"));
    puts(_("-4        Linear 16-bit, same as \"-6 -W -g 1 1\""));
    puts(_("-T        Write TIFF instead of PPM"));
    puts("");
    return 1;
  }
  argv[argc] = "";
  for (arg=1; (((opm = argv[arg][0]) - 2) | 2) == '+'; ) {
    opt = argv[arg++][1];
    if ((cp = (char *) strchr (sp="nbrkStqmHACg", opt)))
      for (i=0; i < "114111111422"[cp-sp]-'0'; i++)
	if (!isdigit(argv[arg+i][0])) {
	  fprintf (stderr,_("Non-numeric argument to \"-%c\"\n"), opt);
	  return 1;
	}
    switch (opt) {
      case 'n':  threshold   = atof(argv[arg++]);  break;
      case 'b':  bright      = atof(argv[arg++]);  break;
      case 'r':
	   FORC4 user_mul[c] = atof(argv[arg++]);  break;
      case 'C':  aber[0] = 1 / atof(argv[arg++]);
		 aber[2] = 1 / atof(argv[arg++]);  break;
      case 'g':  gamm[0] =     atof(argv[arg++]);
		 gamm[1] =     atof(argv[arg++]);
		 if (gamm[0]) gamm[0] = 1/gamm[0]; break;
      case 'k':  user_black  = atoi(argv[arg++]);  break;
      case 'S':  user_sat    = atoi(argv[arg++]);  break;
      case 't':  user_flip   = atoi(argv[arg++]);  break;
      case 'q':  user_qual   = atoi(argv[arg++]);  break;
      case 'm':  med_passes  = atoi(argv[arg++]);  break;
      case 'H':  highlight   = atoi(argv[arg++]);  break;
      case 's':
	shot_select = abs(atoi(argv[arg]));
	multi_out = !strcmp(argv[arg++],"all");
	break;
      case 'o':
	if (isdigit(argv[arg][0]) && !argv[arg][1])
	  output_color = atoi(argv[arg++]);
#ifndef NO_LCMS
	else     out_profile = argv[arg++];
	break;
      case 'p':  cam_profile = argv[arg++];
#endif
	break;
      case 'P':  bpfile     = argv[arg++];  break;
      case 'K':  dark_frame = argv[arg++];  break;
      case 'z':  timestamp_only    = 1;  break;
      case 'e':  thumbnail_only    = 1;  break;
      case 'i':  identify_only     = 1;  break;
      case 'c':  write_to_stdout   = 1;  break;
      case 'v':  verbose           = 1;  break;
      case 'h':  half_size         = 1;  break;
      case 'f':  four_color_rgb    = 1;  break;
      case 'A':  FORC4 greybox[c]  = atoi(argv[arg++]);
      case 'a':  use_auto_wb       = 1;  break;
      case 'w':  use_camera_wb     = 1;  break;
      case 'M':  use_camera_matrix = 3 * (opm == '+');  break;
      case 'I':  read_from_stdin   = 1;  break;
      case 'E':  document_mode++;
      case 'D':  document_mode++;
      case 'd':  document_mode++;
      case 'W':  no_auto_bright    = 1;  break;
      case 'T':  output_tiff       = 1;  break;
      case '4':  gamm[0] = gamm[1] =
		 no_auto_bright    = 1;
      case '6':  output_bps       = 16;  break;
      default:
	fprintf (stderr,_("Unknown option \"-%c\".\n"), opt);
	return 1;
    }
  }
  if (arg == argc) {
    fprintf (stderr,_("No files to process.\n"));
    return 1;
  }
  if (write_to_stdout) {
    if (isatty(1)) {
      fprintf (stderr,_("Will not write an image to the terminal!\n"));
      return 1;
    }
  }
  for ( ; arg < argc; arg++) {
    status = 1;
    raw_image = 0;
    image = 0;
    oprof = 0;
    meta_data = ofname = 0;
    ofp = stdout;
    if (setjmp (failure)) {
      if (fileno(ifp) > 2) fclose(ifp);
      if (fileno(ofp) > 2) fclose(ofp);
      status = 1;
      goto cleanup;
    }
    ifname = argv[arg];
    if (!(ifp = fopen (ifname, "rb"))) {
      perror (ifname);
      continue;
    }
    status = (identify(),!is_raw);
    if (user_flip >= 0)
      flip = user_flip;
    switch ((flip+3600) % 360) {
      case 270:  flip = 5;  break;
      case 180:  flip = 3;  break;
      case  90:  flip = 6;
    }
    if (timestamp_only) {
      if ((status = !timestamp))
	fprintf (stderr,_("%s has no timestamp.\n"), ifname);
      else if (identify_only)
	printf ("%10ld%10d %s\n", (long) timestamp, shot_order, ifname);
      else {
	if (verbose)
	  fprintf (stderr,_("%s time set to %d.\n"), ifname, (int) timestamp);
	ut.actime = ut.modtime = timestamp;
	utime (ifname, &ut);
      }
      goto next;
    }
    write_fun = &CLASS write_ppm_tiff;
    if (thumbnail_only) {
      if ((status = !thumb_offset)) {
	fprintf (stderr,_("%s has no thumbnail.\n"), ifname);
	goto next;
      } else if (thumb_load_raw) {
	load_raw = thumb_load_raw;
	data_offset = thumb_offset;
	height = thumb_height;
	width  = thumb_width;
	filters = 0;
	colors = 3;
      } else {
	fseek (ifp, thumb_offset, SEEK_SET);
	write_fun = write_thumb;
	goto thumbnail;
      }
    }
    if (identify_only && verbose && make[0]) {
      printf (_("\nFilename: %s\n"), ifname);
      printf (_("Timestamp: %s"), ctime(&timestamp));
      printf (_("Camera: %s %s\n"), make, model);
      if (artist[0])
	printf (_("Owner: %s\n"), artist);
      if (dng_version) {
	printf (_("DNG Version: "));
	for (i=24; i >= 0; i -= 8)
	  printf ("%d%c", dng_version >> i & 255, i ? '.':'\n');
      }
      printf (_("ISO speed: %d\n"), (int) iso_speed);
      printf (_("Shutter: "));
      if (shutter > 0 && shutter < 1)
	shutter = (printf ("1/"), 1 / shutter);
      printf (_("%0.1f sec\n"), shutter);
      printf (_("Aperture: f/%0.1f\n"), aperture);
      printf (_("Focal length: %0.1f mm\n"), focal_len);
      printf (_("Embedded ICC profile: %s\n"), profile_length ? _("yes"):_("no"));
      printf (_("Number of raw images: %d\n"), is_raw);
      if (pixel_aspect != 1)
	printf (_("Pixel Aspect Ratio: %0.6f\n"), pixel_aspect);
      if (thumb_offset)
	printf (_("Thumb size:  %4d x %d\n"), thumb_width, thumb_height);
      printf (_("Full size:   %4d x %d\n"), raw_width, raw_height);
    } else if (!is_raw)
      fprintf (stderr,_("Cannot decode file %s\n"), ifname);
    if (!is_raw) goto next;
    shrink = filters && (half_size || (!identify_only &&
	(threshold || aber[0] != 1 || aber[2] != 1)));
    iheight = (height + shrink) >> shrink;
    iwidth  = (width  + shrink) >> shrink;
    if (identify_only) {
      if (verbose) {
	if (document_mode == 3) {
	  top_margin = left_margin = fuji_width = 0;
	  height = raw_height;
	  width  = raw_width;
	}
	iheight = (height + shrink) >> shrink;
	iwidth  = (width  + shrink) >> shrink;
	    if (pixel_aspect < 1) iheight = iheight / pixel_aspect + 0.5;
	    if (pixel_aspect > 1) iwidth  = iwidth  * pixel_aspect + 0.5;
	if (flip & 4)
	  SWAP(iheight,iwidth);
	printf (_("Image size:  %4d x %d\n"), width, height);
	printf (_("Output size: %4d x %d\n"), iwidth, iheight);
	printf (_("Raw colors: %d"), colors);
/* 	if (filters) { */
/* 	  printf (_("\nFilter pattern: ")); */
/* 	  for (i=0; i < 16; i++) */
/* 	    putchar (cdesc[fcol(i >> 1,i & 1)]); */
/* 	} */
	printf (_("\nDaylight multipliers:"));
	FORCC printf (" %f", pre_mul[c]);
	if (cam_mul[0] > 0) {
	  printf (_("\nCamera multipliers:"));
	  FORC4 printf (" %f", cam_mul[c]);
	}
	putchar ('\n');
      } else
	printf (_("%s is a %s %s image.\n"), ifname, make, model);
next:
      fclose(ifp);
      continue;
    }
    if (meta_length) {
      meta_data = (char *) malloc (meta_length);
      merror (meta_data, "main()");
    }
    if (filters || colors == 1) {
      raw_image = (ushort *) calloc ((raw_height+7), raw_width*2);
      merror (raw_image, "main()");
    } else {
      image = (ushort (*)[4]) calloc (iheight, iwidth*sizeof *image);
      merror (image, "main()");
    }
    if (verbose)
      fprintf (stderr,_("Loading %s %s image from %s ...\n"),
	make, model, ifname);

    fseeko (ifp, data_offset, SEEK_SET);
	(*load_raw)();

    if (is_foveon) {
	  if (load_raw ==&CLASS foveon_dp_load_raw)
		foveon_f20_interpolate();
	  else foveon_interpolate();
    }
    fclose(ifp);
     convert_to_rgb();
thumbnail:
    if (write_fun == &CLASS jpeg_thumb)
      write_ext = ".jpg";
    else if (output_tiff && write_fun == &CLASS write_ppm_tiff)
      write_ext = ".tiff";
    else
      write_ext = ".pgm\0.ppm\0.ppm\0.pam" + colors*5-5;
    ofname = (char *) malloc (strlen(ifname) + 64);
    merror (ofname, "main()");
    if (write_to_stdout)
      strcpy (ofname,_("standard output"));
    else {
      strcpy (ofname, ifname);
      if ((cp = strrchr (ofname, '.'))) *cp = 0;
      if (multi_out)
	sprintf (ofname+strlen(ofname), "_%0*d",
		snprintf(0,0,"%d",is_raw-1), shot_select);
      if (thumbnail_only)
	strcat (ofname, ".thumb");
      strcat (ofname, write_ext);
      ofp = fopen (ofname, "wb");
      if (!ofp) {
	status = 1;
	perror (ofname);
	goto cleanup;
      }
    }
    if (verbose)
      fprintf (stderr,_("Writing data to %s ...\n"), ofname);
    (*write_fun)();
   if (ofp != stdout) fclose(ofp);
cleanup:
    if (meta_data) free (meta_data);
    if (ofname) free (ofname);
    if (oprof) free (oprof);
    if (image) free (image);
    if (multi_out) {
      if (++shot_select < is_raw) arg--;
      else shot_select = 0;
    }
  }
  return status;
}
