#ifndef __LIB_X3F__
#define __LIB_X3F__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <sys/types.h>
#include <inttypes.h>
#include <math.h>
#include <wchar.h>
#include <locale.h>
#include <jpeglib.h>
/* #include <glib.h> */
/* #include "x3f_util.h" */

#define LIB_X3F_VERSION "0.2"

#ifndef bool
#define bool boolean
#endif

#ifdef LOCALEDIR
#include <libintl.h>
#define _(String) gettext(String)
#else
#define _(String) (String)
#endif

/* Main file identifier */
#define X3F_FOVb (uint32_t)(0x62564f46)
/* Directory identifier */
#define X3F_SECd (uint32_t)(0x64434553)
/* Property section identifiers */
#define X3F_PROP (uint32_t)(0x504f5250)
#define X3F_SECp (uint32_t)(0x70434553)
/* Image section identifiers */
#define X3F_IMAG (uint32_t)(0x47414d49)
#define X3F_IMA2 (uint32_t)(0x32414d49)
#define X3F_IMA (uint32_t)(0x414d4900)
#define X3F_SECi (uint32_t)(0x69434553)
/* CAMF identifiers */
#define X3F_CAMF (uint32_t)(0x464d4143)
#define X3F_SECc (uint32_t)(0x63434553)
/* CAMF entry identifiers */
#define X3F_CMbP (uint32_t)(0x50624d43)
#define X3F_CMbT (uint32_t)(0x54624d43)
#define X3F_CMbM (uint32_t)(0x4d624d43)
#define X3F_CMb  (uint32_t)(0x00624d43)

#define X3F_FILE_VERSION (((x3f->major)<<16) + x3f->minor)
#define X3F_HEADER_2_0_SIZE 40
#define X3F_HEADER_2_1_EXT_SIZE 192
#define UNIQUE_ID_SIZE 16
#define WB_STRING_SIZE 32
#define PVALUES_SIZE 1024
#define HUFFMAN_TABLE_SIZE 1024
#define TRUE_SEED_TABLE_SIZE 4
#define TRUE_PLANE_SIZE_COUNT 3

#define MAX_THUMB_WIDTH 300 /* verify SD1 thumbnail size */
#define PROPERTY_LIST_HEADER_SIZE 24
#define PROPERTY_ENTRY_SIZE 8
#define CAMF_HEADER_SIZE 28
#define CMb_HEADER_SIZE 20
#define IMA_HEADER_SIZE 28
#define BYTE_BOUNDARY 32
#define TRUE_HUF_MAX_TABLE_SIZE 28
#define MAX_CAMF_NAME_LENGTH 64

#define DECODED_IMAGE 1
#define FLOAT_IMAGE 2
#define WBCC_MAX_LENGTH 24
#define WBGAINS_MAX_LENGTH 24

#define FORC(cnt) for (c=0; c < cnt; c++)
#define FORC3 FORC(3)



#define X3F_MEM_ERROR(_where, _who)					\
  do{									\
    fprintf(stderr, "Memory allocation error in function " _where " for array " _who "\n"); \
  } while(0)

#define X3F_READ_ERROR(_where)					\
  do{								\
    fprintf(stderr, "Reading error in function " _where "\n");	\
  } while(0)

#define X3F_ERROR(_msg)				\
  do{						\
    fprintf(stderr, _msg "\n");			\
  } while(0)

#define X3F_MSG(_msg)				\
  do{						\
    printf(_msg "\n");				\
  } while(0)

typedef uint16_t utf16_t;

/* Some enums */

/* this are the PostProcessing values used by SPP */
/* used in X3F file main header */
/* we should store 3 different sets of PP values */
typedef enum extended_data_type{
  UNUSED =0,
  EXPOSURE = 1,
  CONTRAST = 2,
  SHADOW = 3,
  HIGHLIGHT = 4,
  SATURATION = 5,
  SHARPNESS = 6,
  COLOR_RED = 7,
  COLOR_GREEN = 8,
  COLOR_BLUE = 9,
  FILL_LIGHT = 10,
} EXTENDED_DATA_TYPES;

/* these is the image format */
/* Used in IMA section header */
typedef enum data_format {
	X3F_DATA_FORMAT_PNM_THUMBNAIL = 3,
	X3F_DATA_FORMAT_RAW = 6,
	X3F_DATA_FORMAT_HUFFMAN_PREVIEW = 11,
	X3F_DATA_FORMAT_JPEG = 18,
	X3F_DATA_FORMAT_TRUE_RAW = 30,
} DATA_FORMAT;

/* indicates wether data is raw or ready to display (pnm, jpg) */
/* used in IMA section header */
typedef enum data_type {
  X3F_DATA_TYPE_RAW_SD1 = 1,
  X3F_DATA_TYPE_PROCESSED = 2,
  X3F_DATA_TYPE_RAW = 3,
} DATA_TYPE;

/* used for huffman tree */
struct  decode_s{
  struct decode_s *branch[2];
  int leaf;
};

typedef struct decode_s decode;

/* HEADER struct should be read in one time from file: fread(HEADER_p, 1, sizeof(HEADER_p *), fp) */
typedef struct {
  uint32_t filetypeID; /* X3F_FOVb */
  uint32_t version;
  uint8_t uniqueID[UNIQUE_ID_SIZE];
  uint8_t flags[4]; /* we dispose here of 32 bitsflags: some are used to mark images as protected, flagged*/
  /* argh: Sigma does not seem to use this flags as a bitfield, but as a 4 char array => only 4 flags? */
  /* Sure 32 bits are way to much, but using 8 bits to store an on/off flag... */
  /* OK. pointed images have an integer value of 01000000 */
  /* strangely protected images is not a flag? */
  /* we have to look for values marked by SPP: if SPP do not use the 2, 3 and 4 bytes, then we should use them for our own personal use */
  /* otherwise, we will need to use an 8bit flag in the x3f struct */
  /* 8 bits will be enough to store:
     DECODED_THUMB,
     DECODED_preview,
     DECODED_RAW,
     DECODED_CAMF,
     APPLIED_CAMF,
     MARKED_FOR DELETION,
     MARKED_FOR_EXPORT,
     1 bit free
  */
  uint32_t columns;
  uint32_t rows;
  uint32_t rotation;

  /* The following was added in version 2.1 */
  uint8_t whiteBalanceString[WB_STRING_SIZE]; /* white balance string is null terminated */
  uint8_t extendedDataType[32];
  uint32_t extendedData[32];
} __attribute__ ((packed)) HEADER;

/* DIR_ENTRYs 3 first values should be read in one go from the file: fread(DIR_ENTRY_p, 1, sizeof(uint32_t)*3, fp) */
/* datas must be manually set */
typedef struct {
  uint32_t offset;
  uint32_t dataLength;
  uint32_t type;
  void *datas; /* can be IMA *, PROP_LIST *, or CAMF * */ 
} __attribute__ ((packed)) DIR_ENTRY;

/* DIR_SECTION should be read in one go from the file: fread(DIR_SECTION_p, 1, sizeof(DIR_SECTION_p *), fp) */
typedef struct {
  uint32_t sectionID; /* X3F_SECd */
  uint32_t section_version;
  int32_t dirEntryCount;
} __attribute__ ((packed)) DIR_SECTION;


/* Typedef for the sections type */

/* IMAx: structure to hold the images sections */
typedef struct {
  uint32_t sectionID; /* X3F_SECi */
  uint32_t section_version;
  uint32_t imageDataType;
  uint32_t dataFormat;
  uint32_t columns;
  uint32_t rows;
  uint32_t rowSize;
  void *imageData; /* This can be a small thumb, a big thumb or the raw data */
  uint32_t *rowOffsets;
  /* in the case of thumb or preview, image is of type char *
      in case of raw data, image is of type ushort *[3] */
  uint8_t flags; /* Used to set the state of the image: DECODED_IMAGE, MODIFIED_IMAGE, ... */
  uint16_t max[3];
  uint16_t min[3]; /* have to check if this correspond to DarkLevel values */
} IMA;

typedef struct X3F_PROPERTY{
  utf16_t *name; /* should point to PROP->datas+nameOffset */
  utf16_t *value; /* should point to PROP->datas+valueOffset */
  struct X3F_PROPERTY *next;
} PROPERTY;

/* Do we really need to keep the PROP struct? Probably the chained list of properties is enough */
/* but then we need to strdup the names and values */
typedef  struct X3F_PAIRS_OFFSETS{
  uint32_t nameOffset; 
  uint32_t valueOffset;
} PAIRS_OFFSETS;

typedef struct {
  uint32_t sectionID; /* X3F_SECp */
  uint32_t section_version;
  uint32_t propEntryCount;
  uint32_t characterFormat; /* 0= 16-bit unicode */
  uint32_t reserved;
  uint32_t dataLength; /* in 16-bits characters: in fact we should test the characterFormat to be sure */
  PAIRS_OFFSETS *offsetsTable;
  utf16_t *datas;
 /* contains all properties names and values. MUST BE INITIALIZED with an malloc call */
  /* We have a total of DIR_ENTRY->dataLength-PROP_HEADER_SIZE bytes, which are propEntryCount*2*uint32_t pairs of offsets and the rest for pairs of prop names/values */
} PROP;

/* Do we really need to keep the CAMF struct? Probably the chained list of camf_list_entry is enough */
/* if so, we need to allocate memory for */
/* CAMF: structure to hold all required infos to convert raw image data */
typedef struct camf_typeN_s {
  uint32_t val0;
  uint32_t val1;
  uint32_t val2;
  uint32_t val3;
} camf_typeN;

typedef struct camf_type2_s {
  uint32_t reserved;
  uint32_t infotype; /* CAMF_FCEb */
  uint32_t infotype_version;
  uint32_t crypt_key;
} camf_type2;

typedef struct camf_type4_s {
  uint32_t reserved;
  uint32_t decode_bias;
  uint32_t block_size;
  uint32_t block_count;
} camf_type4;

typedef struct {
  uint32_t sectionID; /* X3F_SECc */
  uint32_t camf_version;
  uint32_t infoDataType;
  union {
    camf_typeN tN;
    camf_type2 t2;
    camf_type4 t4;
  };
  char *camf_data; /*thes is the decoded
		      data */
  uint dataSize;
} CAMF;


/* CAMF entries */
/* known types of entry are CMbT, CMbP, CMbM */
/* all entries should be stored in a chained list */
/* SD1 x3fs do have lots of entries, but still less than 200 */
/* maybe we should have a particular attention for the CMbP IncludeBlocks entry */
/* as it tells us all block's names contained in this camf section */
/* PAY ATTENTION THAT IF WE HAVE MULTIPLE CAMF SECTION IN THE FILE */
/* INCLUDEBLOCKS ENTRY SHOULD BE UPDATED ACCORDINGLY */
typedef struct {
  uint32_t CAMFsubsectionID; /* should be CMbT, CMbP or CMbM*/
  uint32_t version; /* will we need to test for version? actually not*/
  uint32_t length; /* length of entry in byte*/
  uint32_t nameOffset; /* offset from start of entry to name of the entry */
  uint32_t dataOffset; /* offset from start of entry to value of the entry */
/*   uint8_t *data; */
} CMb_HEADER;

/* CMbT = Technical? */
typedef struct {
/*   CMb_HEADER header; */
/*   uint8_t CAMFentryName[16];  */
  uint32_t CAMFentryDataLength;
  uint8_t *data;
} CMbT;

/* CMbP = Parameters ?*/
/* CMbP consists of pairs of strings(?) storable as a HashTable or as a chained list */
typedef struct {
/*   CMb_HEADER header; */
/*   uint8_t *CAMFentryName;  *//* seems to be null terminated, but lenght vary */
  uint32_t numberOfParameters; /* How many pairs are stored in the entry */
  uint32_t dataOffset; /* start of data from the start of the entry */
/*   uint32_t [numberOfEntryData*2]; */ /* pairs of offset from start of data to values */
/* CMbP parameters names and values are null terminated strings. We will need atoi() to convert the numeric values */
} CMbP;

typedef struct {
  char *name;
  char *value;
} PARAMETERS;


/* CMbM = Matrix */
/* Matrix should be seen as a multidimensionnal array */
/* it's representation is matrix[x][y][z][etc] */

/* MATRIX_INFOS struct contains :
   the number of elements,
   the offset to the name of the dimension,
   the indice of this array in the matrix array (multidimensionnal array) */
/* MATRIX_INFOS should be considered as an array of count element */
typedef struct {
  uint32_t count;
  uint32_t nameOffset;
  uint32_t indice;
/*   uchar *name; */
} MATRIX_INFOS;

typedef union {
  uint32_t ui_32;
  int32_t i_32;
  uint16_t ui_16;
  int16_t i_16;
  uint8_t ui_8;
  int8_t i_8;
  float f;
} MATRIX;

/* typedef struct { */
/* /\*   CMb_HEADER header; *\/ */
/*  /\*  uint8_t *CAMFentryName[24]; *\/		/\* seems to be a char[24] NULL terminated string *\/ */
/*   uint32_t typeOfData;			/\* 3=float, 6=uint8_t?, 2=int8_t?, 1=uint32_t?, 0=uint16_t?, 5?*\/ */
/*   uint32_t matrixDimension;		/\* dimension of the matrix *\/ */
/*   uint32_t dataStartOffset; */
/*   MATRIX_INFOS *matrixInfos;	/\* This is an array of size matrixDimension *\/ */
/*   MATRIX *matrix;		/\* matrix is an array of size [matrixInfos[0]->count][...][matrixInfos[n]->count] *\/ */
/* } CMbM; */
typedef enum {
  UI_16 = 0,
  UI_32 =1,
  I_8 = 2,
  FL_32 =3,
  UI_8 = 6,
} CAMF_DATA_TYPE;

typedef struct {
  uint32_t dataType;
  uint32_t *planeElements;
  MATRIX *matrix;
} CMbM;

typedef struct CAMF_ENTRY{
  uint32_t CAMFtype;
  uint32_t count; /* length of string if CMbT, number of parameters if CMbP, matrix dimension if CMbM */
  char *name; /* should point to value+value->header->nameOffset */
  void *value; /* The whole datablock. Could be of type CMbT, CMbP or CMbM */
  /* in the case of CMbT, value will only be a string */
  /* in the case of CMbP, value will be a bidimensionnal array. We should store the number of parameters in the struct */
  struct CAMF_ENTRY *next;
} CAMF_LIST_ENTRY;

typedef struct {
  uint32_t height;
  uint32_t width;
  ushort curve[0x10000];
  int histogram[3][0x2000];
  double gamm[6];
  float rgb_cam[3][3];
  int bps;
  int colors;
  void *img;
} INTERPOLATED_IMG;

/* Main X3F struct */
typedef struct {
  HEADER *header;
  DIR_SECTION *dir_section;
  DIR_ENTRY **dir_entries; /* we do not know the number of entry in the file, this is a table of DIR_ENTRY *
			    dir_entries table contains all entries found in the file, ie all thumbnails, all previews,
			   all raw images, all prop_lists, all camf */
  DIR_ENTRY *thumbnail; /* a pointer to the last recorded thumbnail (processed image with rowSize!=0) */
  DIR_ENTRY *preview; /* a pointer to the last recorded preview (maybe not full size preview) */
  DIR_ENTRY *raw; /* a pointer to the last recorded raw */
  PROPERTY *property; /* only a chained list of pointer to prop
						 section->data */
  char **camf_list;
/*   CAMF_LIST_ENTRY *camf_list;  */
  uint32_t dir_offset; /* offset of the directory section (4 last bytes of the X3F file) */
}  X3F;

void X3F_foveon_camf_decoder(decode *first_decode, uint size, uint16_t code, uint16_t *huff);
void X3F_foveon_image_decoder(decode *first_decode, uint size, uint code, uint *huff);
void X3F_decode_camf2(CAMF *camf, uint dataSize);
void X3F_decode_camf4(CAMF *camf, uint dataSize);
int X3F_decode_preview(IMA *preview, uint32_t dataLength);
int X3F_decode_raw(IMA *raw);
int X3F_decode_thumbnail(IMA *image, uint32_t dataLength);
void X3F_free(X3F *x3f);
X3F *X3F_init(void);
unsigned char *X3F_jpeg_decompress_image(unsigned char *data, uint32_t datasize);
HEADER *X3F_parse_header(FILE *fp);
/* CAMF_LIST_ENTRY *X3F_fill_camf_list (uint dataSize, uint8_t *camf_data, CAMF_LIST_ENTRY *camf_list); */
CAMF *X3F_read_camf(FILE *fp, uint32_t dataLength);
DIR_ENTRY *X3F_read_dir_entry(FILE *fp);
DIR_SECTION *X3F_read_dir_section(FILE *fp);
IMA *X3F_read_ima(FILE *fp, uint dataLength);
PROP *X3F_read_prop(FILE *fp, X3F *x3f, uint32_t dataLength);
X3F *X3F_load_full_x3f(char *filename);
DIR_ENTRY *X3F_get_section(X3F *x3f, uint32_t sectionType);
/* void f20_color_correction(X3F *x3f); */
/* void color_correction(X3F *x3f); */
/* MATRIX *X3F_get_matrix(CAMF_LIST_ENTRY *camf_list, char *entry_name); */
int X3F_foveon_interpolate(X3F *x3f);
int X3F_foveon_TRUE_interpolate(X3F *x3f);
int X3F_foveon_F20_interpolate(X3F *x3f);
int X3F_raw_interpolate(X3F *x3f);
char *X3F_foveon_get_property(PROPERTY *prop_list, char *name);
void * X3F_foveon_camf_matrix (CAMF *camf, unsigned dim[3], const char *name);
const char * X3F_foveon_camf_param (CAMF *camf, const char *block, const char *param);

/* char *foveon_get_param(CAMF_LIST_ENTRY *camf_list, char *blockName, const char *name); */
/* CAMF_LIST_ENTRY *X3F_get_camf_entry(CAMF_LIST_ENTRY *camf_list, char *entry_name); */

#endif /* __LIB_X3F__ */
