/* This libx3f will only compile on little-endian processors, which are the vast majority of personnal cpu */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <sys/types.h>
#include <math.h>
#include <wchar.h>
#include <locale.h>
#include <jpeglib.h>
#include "raw_x3f.h"

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


X3F *X3F_init(void){
  X3F *x3f=NULL;

  if (!(x3f=malloc(sizeof(*x3f)))){
    X3F_MEM_ERROR("x3f_init","x3f");
    return NULL;
  }

  x3f->header=NULL;
  x3f->dir_section=NULL;
  x3f->dir_entries=NULL;
  x3f->thumbnail=NULL;
  x3f->preview=NULL;
  x3f->raw=NULL;
  x3f->property=NULL;
  x3f->camf_list=NULL;
  x3f->dir_offset=0;

  return x3f;
}

void X3F_free(X3F *x3f){
  int i;
  PROP *prop;
  IMA *image;
  CAMF *camf;
  /* no need to free thumbnail, preview or raw */
  x3f->thumbnail=x3f->preview=x3f->raw=NULL;
  /* But we should go down to header, dir_section, and dir_entries[]  to free all dynamically allocated variables */
  for (i=x3f->dir_section->dirEntryCount-1; i>=0; i--){
    if (x3f->dir_entries[i]->datas) {
      switch (x3f->dir_entries[i]->type) {
      case X3F_PROP:
	printf("Cleaning PROP");
	prop=(PROP *)x3f->dir_entries[i]->datas;
	free(prop->offsetsTable);
	free(prop->datas);
	printf("\tDone\n");
	break;
      case X3F_IMAG:
      case X3F_IMA2:
	printf("Cleaning IMA");
	image=(IMA *)x3f->dir_entries[i]->datas;
	if (image->imageData)
	  free(image->imageData);
	if (image->rowOffsets)
	  free(image->rowOffsets);
	printf("\tDone\n");
	break;
      case X3F_CAMF:
	printf("Cleaning CAMF");
	camf=(CAMF *)x3f->dir_entries[i]->datas;
/* 	if (camf->camf_data) */
/* 	  free(camf->camf_data); */
/* 	free(camf); */
	printf("\tDone\n");
	break;
      default: /* unknown entry_type: something bad happened cause entry->datas should not have been allocated */
	/* ok, let this in place */
	X3F_MEM_ERROR("X3F_free", "Unknown entry type");
	break;
      }
      if (x3f->dir_entries[i]->datas)
	free(x3f->dir_entries[i]->datas);
   }
  }
  free(x3f->dir_entries);
  free(x3f->dir_section);
  free(x3f->header);
  /* then we free x3f */
  i=0;
  if(x3f->property) {
    PROPERTY *property, *index;
    for (property=x3f->property, index=property->next;index!=NULL;index=property->next){
      free(property);
      property=index;
    }
    free(property); /* don't forget the last node */
  }
  if(x3f->camf_list) {
    CAMF_LIST_ENTRY *entry, *index;
    CMbM *cmbm;
    PARAMETERS *params;
    for (entry=x3f->camf_list, index=entry->next;index!=NULL;index=entry->next){
      if (entry->CAMFtype==X3F_CMbM){
	cmbm=(CMbM *)entry->value;
	free(cmbm->planeElements);
	free(cmbm->matrix);
      } else if (entry->CAMFtype==X3F_CMbP){
	params=(PARAMETERS *)entry->value;
	for (i=0; i<(int)entry->count; i++){
	  free(params[i].name);
	  free(params[i].value);
	}
      }
      free(entry->name);
      free(entry->value);
      free(entry);
      entry=index;
    }
    if (entry->CAMFtype==X3F_CMbM){
      cmbm=(CMbM *)entry->value;
      free(cmbm->planeElements);
      free(cmbm->matrix);
    } else if (entry->CAMFtype==X3F_CMbP){
      params=(PARAMETERS *)entry->value;
      for (i=0; i<(int)entry->count; i++){
	free(params[i].name);
	free(params[i].value);
      }
    }

    free(entry->name);
    free(entry->value);
    free(entry); /* don't forget the last node */
  }
  free(x3f);
}

HEADER *X3F_parse_header(FILE *fp){
  HEADER *header=NULL;

  if (!(header=malloc(sizeof(*header)))) {
    X3F_MEM_ERROR("X3F_parse_header","header");
    return NULL;
  }

  /* function calling us is responsible to set the offset correctly */
  fread(header, 1, X3F_HEADER_2_0_SIZE, fp);
  if (ferror(fp)){
    X3F_READ_ERROR("X3F_parse_header");
    return NULL;
  }
  if (header->filetypeID != X3F_FOVb){
    X3F_ERROR("Not a FOVb file");
    return NULL;
  }
  printf("File Version: %x\n", header->version);
  if (header->version > (2<<16)) {  // be sure file version is > 2.0
  printf("File Version: %x\n", header->version);
    fread(header->whiteBalanceString, 1, X3F_HEADER_2_1_EXT_SIZE, fp); 
  printf("White Balance String %s\n", (char *)header->whiteBalanceString);
    if (ferror(fp)){
      X3F_READ_ERROR("X3F_parse_header");
      return NULL;
    }
  } else { // initialize all with 0
    memset(header+X3F_HEADER_2_0_SIZE, 0, X3F_HEADER_2_1_EXT_SIZE);
  }
  X3F_MSG("Header has been parsed\n");
  return header;
}

DIR_SECTION *X3F_read_dir_section(FILE *fp){
  /* We read and store the content of directory section. This is
     needed if we want to be able to write an x3f file */
  DIR_SECTION *dir_section=NULL;
  int c=0;

  if (!(dir_section=malloc(sizeof(*dir_section))))
    X3F_MEM_ERROR("X3F_read_dir_section","dir_section");

  c=fread(dir_section, 1, sizeof(*dir_section), fp);
  if (ferror(fp)){
    X3F_READ_ERROR("X3F_read_dir_section");
    return NULL;
  }

  X3F_MSG("Directory section has been parsed\n");
  return dir_section;
}

DIR_ENTRY *X3F_read_dir_entry(FILE *fp){
  DIR_ENTRY *subsection=NULL;

  if (!(subsection=malloc(sizeof(DIR_ENTRY))))
    X3F_MEM_ERROR("X3F_read_dir_entry","subsection");

  fread(subsection, 1, sizeof(uint32_t)*3, fp);
  if (ferror(fp)){
    X3F_READ_ERROR("X3F_read_dir_entry");
    return NULL;
  }
  subsection->datas=NULL;

  X3F_MSG("Directory entry has been parsed\n");
  return subsection;
}

PROP *X3F_read_prop(FILE *fp, X3F *x3f, uint32_t dataLength) {
  uint i=0;
  PROP *prop=NULL;
  utf16_t *datas;

  if (!(prop=malloc(sizeof(*prop)))){
    X3F_MEM_ERROR("X3F_read_prop", "prop");
    return NULL;
  }
  fread(prop, 1, PROPERTY_LIST_HEADER_SIZE, fp);
  if (ferror(fp)){
    X3F_READ_ERROR("X3F_read_prop");
    free(prop);
    return NULL;
  }

  if (!(prop->offsetsTable=malloc(prop->propEntryCount * sizeof(*prop->offsetsTable)))){
    X3F_MEM_ERROR("X3F_read_prop", "prop->offsetsTable");

    /* FIXME: try to get through and guess the pairs based on null terminated strings? */

    free(prop);
    return NULL;
  }
  fread(prop->offsetsTable, sizeof(*prop->offsetsTable), prop->propEntryCount,fp);

  if (!(datas=malloc(dataLength-PROPERTY_LIST_HEADER_SIZE-prop->propEntryCount*sizeof(*prop->offsetsTable)))){
    X3F_MEM_ERROR("X3F_read_prop", "prop->datas");
    free(prop->offsetsTable);
    free(prop);
    return NULL;
  }
 
  fread(datas, 1, dataLength-PROPERTY_LIST_HEADER_SIZE-prop->propEntryCount*sizeof(*prop->offsetsTable), fp);
  if (ferror(fp)){
    X3F_READ_ERROR("X3F_read_prop");
    free(prop->offsetsTable);
    free(datas);
    free(prop);
    return NULL;
  }
  prop->datas=datas;
  if (x3f) {
    PROPERTY *new_prop=NULL;
    utf16_t *name;
    do {
      name=prop->datas+(prop->offsetsTable[i].nameOffset);
      for (new_prop=x3f->property; new_prop!=NULL;new_prop=new_prop->next){
	if (!wcscmp((wchar_t *)name, (wchar_t *)new_prop->name)) {
	  break;
	}
      }
      if (new_prop==NULL){
	if (!(new_prop=malloc(sizeof(*new_prop)))){
	  X3F_MEM_ERROR("X3F_read_prop", "new_prop");
	  free(prop->offsetsTable);
	  free(datas);
	  free(prop);
	  return NULL;
	}
	new_prop->name=name;
 	new_prop->next=x3f->property;
	new_prop->value=prop->datas+(prop->offsetsTable[i].valueOffset);
	x3f->property=new_prop;
      }
      /* We just found an existing property with the same name */
      /* Just ignore the new value, as we are parsing the file in reverse order (most up-to-date entries first) */
      /* 	new_prop->value=value; */ 
      i++;
    } 
    while (i<prop->propEntryCount);

  }
  return prop;
}

CAMF_LIST_ENTRY *X3F_fill_camf_list(uint dataSize, uint8_t *camf_data, CAMF_LIST_ENTRY *camf_list){
  CMb_HEADER *CMb_entry=NULL;
/*   CAMF_LIST_ENTRY *camf_list=NULL; */
  uint n=0;
  uint8_t *dataPtr=camf_data;

  while(dataPtr<camf_data+dataSize){
    char *name;
    CAMF_LIST_ENTRY *new_list_entry=NULL;

/*     dataPtr+=pos; */
    CMb_entry=(CMb_HEADER *)dataPtr;
    if ((CMb_entry->CAMFsubsectionID & 0xffffff)!= X3F_CMb) {
      printf("Oups: %#x\n",CMb_entry->CAMFsubsectionID);
      break;
    }

      /*If we do a memcpy, we will be able to free the camf_data array and only keep valid camf entries */
    /* Whatever, do we really need to store whole camf data blocks? Maybe easier if we want to write x3f files? */
/*     memcpy(CMb_entry, elem, entryLength); */

    if (!(name=malloc(sizeof(*name)*strlen((char *)dataPtr+CMb_entry->nameOffset)))){
	X3F_MEM_ERROR("X3F_fill_camf_list", "new_list_entry->name");
	free(new_list_entry);
	return camf_list;
    }
    name=strdup((char *)dataPtr+CMb_entry->nameOffset);
/*     printf("CAMF entry name is %s\n", name); */
    if (!(strcmp(name, "IncludeBlocks"))){
      free(name);
      dataPtr+=CMb_entry->length;
      continue;
    }
    /* we should first check for duplicate entries in the list */
    /* normaly, only one CAMF section is to be found in the X3F file. But the specs give the ability to have more than one */
    /* We handle this case here */
    for (new_list_entry=camf_list; new_list_entry!=NULL; new_list_entry=new_list_entry->next){
      if (!(strcmp(new_list_entry->name, name))){
	free(name);
	break;
      }
    }
    if (new_list_entry){ 
      printf("We found a camf block with duplicate name: %s,\nWe will ignore it\n", new_list_entry->name);
      dataPtr+=CMb_entry->length;
      continue;
    }

    /* Create a new list_entry */
    if (!(new_list_entry=malloc(sizeof(*new_list_entry)))) {
      X3F_MEM_ERROR("X3F_fill_camf_list", "CMb_entry");
      free(name);
      return camf_list;
    }
    new_list_entry->next=NULL;
    new_list_entry->CAMFtype=CMb_entry->CAMFsubsectionID;
    new_list_entry->name=name;
    new_list_entry->count=0;

    if ((CMb_entry->CAMFsubsectionID&0xffffffff)== X3F_CMbT){
      new_list_entry->count=*((uint32_t *)(dataPtr+CMb_entry->dataOffset));
      new_list_entry->value=(void *)strndup((char *)(dataPtr+CMb_entry->dataOffset+4), new_list_entry->count);
/*       printf("Technical: %s\n", (char *)new_list_entry->value); */
    } else if ((CMb_entry->CAMFsubsectionID&0xffffffff)== X3F_CMbP){
      /* we should handle the "IncludeBlocks" as a special case */
      uint i, paramCount;
      PAIRS_OFFSETS *offsets;
      PARAMETERS *params;
      uint32_t *cmbp_head, dataLength;

      cmbp_head=(uint32_t *)(dataPtr+CMb_entry->dataOffset);
      paramCount=cmbp_head[0];
      dataLength=CMb_entry->length-cmbp_head[1]; /* cmbp_head[1] is the start of parameters offset */
/*       printf("Founded %d cmbp->numberOfParameters for a size of %d bytes\n", cmbp_head[0], dataLength ); */
/*       if (!(offsets=malloc(sizeof(*offsets)*cmbp->numberOfParameters))){ */
      if (!(params=malloc(sizeof(*params)*paramCount))){
	X3F_MEM_ERROR("X3F_fill_camf_list", "CMbP offsets table");
	free(new_list_entry->name);
	free(new_list_entry);
	return camf_list;
      }
      
      offsets=(PAIRS_OFFSETS *)(dataPtr+CMb_entry->dataOffset+sizeof(cmbp_head)); 
      for (i=0;i<paramCount;i++){
	params[i].name=strdup((char *)dataPtr+cmbp_head[1]+offsets[i].nameOffset);
	params[i].value=strdup((char *)dataPtr+cmbp_head[1]+offsets[i].valueOffset);
      }
      new_list_entry->count=paramCount;
      new_list_entry->value=(void *)params;
    } else if ((CMb_entry->CAMFsubsectionID&0xffffffff)== X3F_CMbM){
      /* header and name are OK */
      CMbM *cmbm;
      uint i=0/* , c=0,v */;
      uint8_t*cmbmPtr=dataPtr+CMb_entry->dataOffset;
      uint32_t matrixDimension, dataStartOffset, dataType, valuesCount=1;
      int32_t *planeElements;
      MATRIX_INFOS *matrixInfos;
      MATRIX *matrix;

      dataType=*(uint32_t *)(cmbmPtr);
      matrixDimension=*(uint32_t *)(cmbmPtr+4);
      dataStartOffset=*(uint32_t *)(cmbmPtr+8);

       new_list_entry->count=matrixDimension;
      /* Do we need to keep infos such as name of arrays? */
      matrixInfos=(MATRIX_INFOS *)(cmbmPtr+12);

      planeElements=malloc(sizeof(*planeElements)*matrixDimension);
      for (i=0; i<matrixDimension;i++){
	planeElements[i]=matrixInfos[i].count;
	valuesCount*=planeElements[i];
      }

      matrix=malloc(sizeof(*matrix)*valuesCount);
      if (!dataType ||dataType==6)
 	for (i=0; i<valuesCount;i++)
	  matrix[i].ui_16=*((uint16_t *)(dataPtr+dataStartOffset+i*sizeof(uint16_t)))&0xffff;
     else
	memcpy(matrix, dataPtr+dataStartOffset, sizeof(*matrix)*valuesCount);
       
      if (!(cmbm=malloc(sizeof(*cmbm)))) {
	X3F_MEM_ERROR("X3F_fill_camf_list", "CMb_entry");
	free(name);
	return camf_list;
      }
      cmbm->dataType=dataType;
      cmbm->planeElements=planeElements;
      cmbm->matrix=matrix;

      new_list_entry->value=cmbm;

    }

    /* add the new entry to the chained list */
    new_list_entry->next=camf_list;
    camf_list=new_list_entry;
    dataPtr+=CMb_entry->length;

    /* just for stats */
    n++;
  }
  printf("%d camf entries were found\n", n);
  free(camf_data);
  return camf_list;
}


CAMF *X3F_read_camf(FILE *fp, uint32_t dataLength) {
  CAMF *camf;
  uint dataSize;

  if (!(camf=malloc(sizeof(*camf)))){
    X3F_MEM_ERROR("X3F_read_camf()", "camf");
    return NULL;
  }
  fread(camf,1,CAMF_HEADER_SIZE,fp); 
  if (ferror(fp)){
    X3F_READ_ERROR("X3F_read_camf");
    free(camf);
    return NULL;
  }
  dataSize=dataLength-CAMF_HEADER_SIZE;
  if (!(camf->camf_data=malloc(dataSize))){
    X3F_MEM_ERROR("X3F_read_camf()", "camf_data");
    free(camf);
    return NULL;
  }
  fread(camf->camf_data, 1, dataSize, fp);
  if (ferror(fp)){
    X3F_READ_ERROR("X3F_read_camf");
    free(camf->camf_data);
    free(camf);
    return NULL;
  }
  switch (camf->infoDataType) {
  case 2:
    X3F_decode_camf2(camf, dataSize);
    break;
  case 4:
    dataSize= camf->t4.block_count*camf->t4.block_size*3/2; /* dataSize of TRUE decoded camf is less than what's computed here */
    X3F_decode_camf4(camf, dataSize);
    break;
  default:
    X3F_MSG("Unknown camf type");
  }
  return camf;
}

/* to decode camf section from earlier SD9/10/14 */
void X3F_decode_camf2(CAMF *camf, uint dataSize){
  uint32_t key, val;
  uint i;

  key = camf->t2.crypt_key;
  for (i=0; i < dataSize; i++) {
    key = (key * 1597 + 51749) % 244944;
    val = key * (int64_t) 301593171 >> 24;
    camf->camf_data[i] ^= ((((key << 8) - val) >> 1) + val) >> 17;
  }
  /* camf->camf_data contains decoded datas */
}

/* Should be renamed X3F_foveon_TRUE_huffman_decoder */
void X3F_foveon_camf_decoder(decode *first_decode, uint size, uint16_t code, uint16_t *huff){
  decode *cur;
  decode *free_decode;
  uint i;
  uint8_t  n, len;

  if (!code){
    memset (first_decode, 0, sizeof(first_decode));
    free_decode = first_decode;
    first_decode->branch[0]=first_decode->branch[1]=NULL;
    first_decode->leaf=-1;
  }
  // Here we fill the tree
  for (i=0;i<size;i++){
    len=huff[i]&0x3f;
    code=huff[i]>>8&0xff;
    cur=first_decode;
    for (n=0; n<len;n++){
      if (!cur->branch[(code>>(7-n))&1]){
	free_decode++;
	if (free_decode > first_decode+2*size) {
	  fprintf (stderr,_("decoder table overflow\n"));
	  exit(1);
	}
	cur=cur->branch[(code>>(7-n))&1]=free_decode;
	free_decode->branch[0]=free_decode->branch[1]=NULL;
	free_decode->leaf=-1;
     } else{
	cur=cur->branch[(code>>(7-n))&1];
      }
    }
    cur->leaf=i;
  }
}

/* to decode TRUE encoded camf section */
void X3F_decode_camf4(CAMF *camf, uint dataSize){
  /* we first have a huffman table whose size is not known yet, but should be 32 bytes max*/
  /* rolkar decoder give a size of 13 -> 12 leaves +1 */
  /* strlen is 26 */
  /* the table is in the form 2*8bits integer (2 bytes): lower byte is the code size, higher byte is the code  */
  /* followed by the encrypted data */
  /* Probably the encrypted datas start is padded to a 32 bytes boundary */
/*   X3F_MSG("camf type 4 are not supported yet ;\("); */
  uint16_t *huff=(uint16_t *)camf->camf_data;
  uint8_t *tmp=camf->camf_data;
  /* TRUE huffman table is 0 terminated */
  int TRUE_tableSize=strlen((char *)camf->camf_data)/2;
  decode *first_decode;
  decode *dindex;
  int32_t row_start_acc[2][2]={{camf->t4.decode_bias,camf->t4.decode_bias},{camf->t4.decode_bias,camf->t4.decode_bias}};

  uint8_t *decoded_data=NULL, *ptr;
  uint32_t col,row, bitbuf=0;
  int32_t diff, bit, i;
  uint8_t bits, b, sign;

  first_decode=malloc(sizeof(*first_decode)*TRUE_tableSize*4);

  uint8_t *raw_data=(uint8_t *)camf->camf_data+BYTE_BOUNDARY;
  X3F_foveon_camf_decoder(first_decode, TRUE_tableSize, 0, huff);

  if (!(decoded_data=malloc(dataSize*sizeof(*decoded_data)))) {
    X3F_MEM_ERROR("X3F_decode_camf4","decoded_data");
    return;
  }
  ptr=decoded_data;
  bit=0;
  for (row = 0; row < camf->t4.block_count; row++) {
    int32_t acc[2]={row_start_acc[row&1][0],row_start_acc[row&1][1]};
    for (col = 0; col < camf->t4.block_size; col++) {
      for (dindex=first_decode; dindex->leaf<0; ) {
	if ((bit = (bit-1) & 31) == 31){
	  for (i=0; i < 4; i++,raw_data++){
	    bitbuf = (bitbuf << 8) + raw_data[0];
	  }
	}
	dindex = dindex->branch[bitbuf >> bit & 1];
      }/* dindex now points to the leaf */
      bits = dindex->leaf;
      if (bits==0)
	diff=0;
      else {
	if ((bit=(bit-1)&31)==31)
	  for (b=0; b < 4; b++){
	    bitbuf = (bitbuf << 8) + raw_data[0];
	    raw_data++;
	  }
	sign=diff=bitbuf>>bit &1;
	/* 	  Attention, on ne veut garder que 16 bit max dans diff */
	for (i=1;i<bits;i++){
	  if ((bit=(bit-1)&31)==31)
	    for (b=0; b < 4; b++){
	      bitbuf = (bitbuf << 8) + raw_data[0];
	      raw_data++;
	    }
	  diff=(diff<<1) + ((bitbuf>>bit)&1);
	}
	if (sign == 0)
	  diff -= (1<<bits) - 1;
      }

/* These comes from Roland Karlsson 's X3F tools */
      acc[col&1]+=diff;
      if (col<2) row_start_acc[row&1][col&1]=acc[col&1];

      switch(col&1) {
      case 0:
	*ptr++  = (uint8_t)((acc[col&1]>>4)&0xff);
	*ptr    = (uint8_t)((acc[col&1]<<4)&0xf0);
	break;
      case 1:
	*ptr++ |= (uint8_t)((acc[col&1]>>8)&0x0f);
	*ptr++  = (uint8_t)((acc[col&1]<<0)&0xff);
	break;
      }

    } /* end col */
  } /* end row */
  camf->camf_data=decoded_data;
  free(tmp);
}

IMA *X3F_read_ima(FILE *fp, uint dataLength){
  IMA *ima=NULL;
  uint32_t imageDataSize;

  if (!(ima=malloc(sizeof(IMA)))){
    X3F_MEM_ERROR("X3F_read_ima","ima");
    return NULL;
  }
  fread(ima,1,IMA_HEADER_SIZE, fp);
  if (ferror(fp)){
    X3F_READ_ERROR("X3F_read_ima");
    free(ima);
    return NULL;
  }

  imageDataSize=dataLength-IMA_HEADER_SIZE;
  if (ima->dataFormat == X3F_DATA_FORMAT_RAW || ima->dataFormat==X3F_DATA_FORMAT_HUFFMAN_PREVIEW)
    imageDataSize -= ima->rows*sizeof(*ima->rowOffsets);

  if (!(ima->imageData=(void *)malloc(imageDataSize))){
    X3F_MEM_ERROR("X3F_read_ima", "image");
    free(ima);
    return NULL;
  }
  fread(ima->imageData, 1, imageDataSize, fp);
  if (ferror(fp)){
    X3F_READ_ERROR("X3F_read_ima");
    free(ima->imageData);
    free(ima);
    return NULL;
  }
  ima->rowOffsets=NULL;
  if (ima->dataFormat == X3F_DATA_FORMAT_RAW || ima->dataFormat==X3F_DATA_FORMAT_HUFFMAN_PREVIEW){
    if (!(ima->rowOffsets=(uint32_t *)malloc(ima->rows*sizeof(*ima->rowOffsets)))){
      X3F_MEM_ERROR("X3F_read_ima", "rowOffsets");
      free(ima);
      return NULL;
    }
    fread(ima->rowOffsets, 1, ima->rows*sizeof(*ima->rowOffsets), fp);
    if (ferror(fp)){
      X3F_READ_ERROR("X3F_read_ima");
      free(ima->rowOffsets);
      free(ima->imageData);
      free(ima);
      return NULL;
    }
  }

  ima->flags=0;
  return ima;
}

/* should be renamed X3F_foveon_huffman_decoder */
/* This is stolen from DCraw */
void X3F_foveon_image_decoder(decode *first_decode, uint size, uint code, uint *huff){
  decode *cur;
  static decode *free_decode;
  unsigned i, len;

  if (!code){
    memset (first_decode, 0, sizeof first_decode);
    free_decode = first_decode;
  }
  cur=free_decode++;
  cur->leaf=-1;
  if (free_decode > first_decode+2*size) {
    fprintf (stderr,_("decoder table overflow\n"));
    exit(1);//    longjmp (failure, 2);
  }
  if (code)
    // Here we fill the tree
    for (i=0; i < size; i++)
      if (huff[i] == code) {
	cur->leaf = i;
	return;
      }
  /* this is just to isolate the bottom 27 bits */
  if ((len = code >> 27) > 26) return;
  code = (len+1) << 27 | (code & 0x3ffffff) << 1;

  cur->branch[0] = free_decode;
  X3F_foveon_image_decoder (first_decode, size, code, huff);
  cur->branch[1] = free_decode;
  X3F_foveon_image_decoder (first_decode, size, code+1, huff);
}
/* end dcraw */

unsigned char *X3F_jpeg_decompress_image(unsigned char *data, uint32_t datasize){
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  unsigned char *decoded_data=NULL;
  unsigned char *outbuf[1];
 
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
	
  jpeg_mem_src(&cinfo, data, datasize);
  jpeg_read_header(&cinfo, TRUE);
	
  jpeg_start_decompress(&cinfo);
  if (!(decoded_data=malloc(cinfo.output_width*cinfo.output_height*cinfo.output_components)))
    X3F_MEM_ERROR("X3F_jpeg_decompress_image","decoded_data");
  outbuf[0]=decoded_data;
  while (cinfo.output_scanline < cinfo.output_height){
    jpeg_read_scanlines(&cinfo, outbuf, 1);
    outbuf[0]+=cinfo.num_components * cinfo.output_width ;
  }
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  return decoded_data;
}

/* The X3F_decode functions returns 1 in case of success, 0 otherwise */
int X3F_decode_thumbnail(IMA *image, uint32_t dataLength){
  /* thumbnail (the last one found in the x3f file) is pointed by x3f->thumbnail */
  unsigned char *raw_data=(unsigned char*)image->imageData;
  unsigned char *decoded_data=NULL;

  if ((image->flags & DECODED_IMAGE) == DECODED_IMAGE)
    return 1;
  if (image->dataFormat==X3F_DATA_FORMAT_PNM_THUMBNAIL){ //pnm thumbnail
    /* hey, nothing to be done here */
    printf("PNM thumbnail\n");
    image->flags |= DECODED_IMAGE;
    return 1;
  } else if (image->dataFormat==X3F_DATA_FORMAT_JPEG) { // jpeg tumbnail
    /* decode the jpeg */
    decoded_data=X3F_jpeg_decompress_image(raw_data, dataLength-IMA_HEADER_SIZE);
    printf("JPEG thumbnail\n");
  } else {
    X3F_ERROR("Unknown thumbnail data format");
    free(decoded_data);
  }
  if (decoded_data){
    void *tmp=image->imageData;
    image->imageData=decoded_data;
/*     if (!image->rowSize) */
      image->rowSize=image->columns*3;
    image->flags |= DECODED_IMAGE;
    free(tmp);
    return 1;
  }
  return 0;
}

int X3F_decode_preview(IMA *preview,uint32_t dataLength){
  uint8_t *raw_data=(uint8_t *)preview->imageData;
  uint8_t *decoded_data=NULL;
  uint8_t *imageDataStart;

  if ((preview->flags & DECODED_IMAGE)==DECODED_IMAGE)
    return 1;
  if (preview->dataFormat==X3F_DATA_FORMAT_HUFFMAN_PREVIEW){ /* huffman compressed pnm preview */
    /* 256 first uint32 are the codeword table */
    /* following are the encoded data */
    /* we need a function to create the huffman tree */
    uint32_t *huff=(uint32_t *)raw_data;
    decode *first_decode;
    decode *dindex;
    uint8_t pixel[3];
    uint row, col,  c, i, n=0;
    int bitbuf=0, bit=-1;

    first_decode=malloc(sizeof(*first_decode)*512);
    X3F_foveon_image_decoder(first_decode, 256, 0, huff);
    raw_data+=256*sizeof(*huff);
    imageDataStart=raw_data;

    if (!(decoded_data=malloc(preview->columns*preview->rows*3*sizeof(*decoded_data)))) {
      X3F_MEM_ERROR("X3F_decode_preview","decoded_data");
      return 0;
    }
    /* This comes mainly from dcraw */
    for (row=0; row < preview->rows; row++) {
      memset (pixel, 0, sizeof(uint8_t)*3);
      raw_data=imageDataStart+preview->rowOffsets[row];

      for (bit=col=0; col < preview->columns; col++){
	for (c=0; c < 3; c++) {
	  for (dindex=first_decode; dindex->leaf<0; ) {
	    if ((bit = (bit-1) & 31) == 31){
	      for (i=0; i < 4; i++){
		bitbuf = (bitbuf << 8) + raw_data[0];
		raw_data++;
	      }
	    }
	    dindex = dindex->branch[bitbuf >> bit & 1];
	  }
	  pixel[c] += dindex->leaf;
	}
	for (c=0;c<3;c++)
	  decoded_data[n++]=pixel[c];
      }
    }
    /* end dcraw */
    free(first_decode);
  } else if (preview->dataFormat==X3F_DATA_FORMAT_JPEG) { /* jpeg preview */
    /* decode the jpeg */
    decoded_data=X3F_jpeg_decompress_image((unsigned char *)raw_data, dataLength-IMA_HEADER_SIZE);
     printf("JPEG preview\n");
  } else {
    X3F_ERROR("Unknown preview data format");
  }
  if (decoded_data){
    void *tmp=preview->imageData;
    preview->imageData=decoded_data;
/*     if (!preview->rowSize) */
      preview->rowSize=preview->columns*3;
    preview->flags|=DECODED_IMAGE;
    free(tmp);
    return 1;
  }
  return 0;

}

int X3F_decode_raw(IMA *raw){
  uint8_t *raw_data=(uint8_t *)raw->imageData;
  uint8_t *imageDataStart;
  unsigned row,col,c,i,b;
  int bitbuf=0, bit=-1;
  decode *dindex;
  uint16_t pix[3], *pvalues;
  uint16_t *decoded_raw;

  if ((raw->flags & DECODED_IMAGE)== DECODED_IMAGE)
    return 1;

  if (raw->dataFormat==X3F_DATA_FORMAT_RAW){
    /* Huffman compressed raw (SD9/10/14) */
    uint32_t *huff;
    decode *first_decode;
    uint n=0;

    /* first 1024 uint16 are the pixels values table */
    pvalues=(uint16_t *)raw_data;
    raw_data+= PVALUES_SIZE*sizeof(*pvalues);
    /* following 1024 uint32 are the huffman table */
    huff=(uint32_t *)raw_data;
    raw_data+=HUFFMAN_TABLE_SIZE*sizeof(*huff);
    imageDataStart=raw_data;

    first_decode=malloc(sizeof(*first_decode)*2048);
    /* next are the compressed image data */
    X3F_foveon_image_decoder(first_decode, 1024, 0, huff);
    if (!(decoded_raw=(uint16_t *) calloc (raw->rows*raw->columns*3, sizeof (*decoded_raw)))){
      X3F_MEM_ERROR("X3F_decode_raw", "decoded_raw");
      return 0;
    }
    /* This comes from dcraw */
    for (row=0; row < raw->rows; row++) {
      memset (pix, 0, sizeof(pix));
      raw_data=imageDataStart+raw->rowOffsets[row];
      for (bit=col=0; col < raw->columns; col++){
	if (raw->rowSize) {
	  for (c=0; c<3; c++) pix[2-c] += pvalues[raw_data[c] >> c*10 & 0x3ff];
	  raw_data+=3;
	} else {
	  for (c=0; c < 3; c++) {
	    for (dindex=first_decode; dindex->leaf<0; ) {
	      if ((bit = (bit-1) & 31) == 31){
		for (i=0; i < 4; i++){
		  bitbuf = (bitbuf << 8) + raw_data[0];
		  raw_data++;
		}
	      }
	      dindex = dindex->branch[bitbuf >> bit & 1];
	    }
	    pix[c] += pvalues[dindex->leaf];
	    pix[c]=(int16_t)pix[c]>0?pix[c]:0;
	    if (pix[c] >> 16 && ~pix[c] >> 16){
	      printf("Oups\n");
	      free(decoded_raw);
	      return 0; /* Corrupted data? */
	    }
	  }
	}
	for (c=0; c <3; c++) decoded_raw[n++]= pix[c];
      }
    }
    /* end dcraw */
    free(first_decode);
  } else if (raw->dataFormat == X3F_DATA_FORMAT_TRUE_RAW){
    /* TRUE compressed raw (SD15/1, DP1/2) */
    /* We have a TRUE header */
    /* uint16_t seed[4]; */
    /* huffman table */
    /* uint32_t planeSize [3] */
    uint16_t *seed, *huff=NULL;
    uint32_t *PlaneSize;
    decode *first_decode, *dindex;
    uint32_t col,row;
    uint8_t bits, sign;
    int16_t diff;
    int TRUE_tableSize;
    uint8_t *dataStart;

    seed=(uint16_t *)raw_data;
    raw_data+=TRUE_SEED_TABLE_SIZE*sizeof(uint16_t);
    huff=(uint16_t *)raw_data;
    TRUE_tableSize=strlen((char *)huff)/2; /* should be 13 */
    raw_data+=TRUE_HUF_MAX_TABLE_SIZE;
    PlaneSize=(uint32_t *)raw_data;
    raw_data+=TRUE_PLANE_SIZE_COUNT*sizeof(*PlaneSize);
    dataStart=raw_data;

    first_decode=malloc(sizeof(*first_decode)*TRUE_tableSize*4); /* FIXME: Why *4? */
    X3F_foveon_camf_decoder(first_decode, TRUE_tableSize, 0, huff);
    /* pixels are stored R0->Rn, G0->Gn, B0->Bn */
    /* this means we first decode the red channel, then the green one and finally the blue one */
    if (!(decoded_raw=(uint16_t *) calloc (raw->rows*raw->columns*3, sizeof (*decoded_raw)))){
      X3F_MEM_ERROR("X3F_decode_raw", "decoded_raw");
      free(first_decode);
      return 0;
    }
    /* Borrowed to Roland Karlsson X3FTools */
    for (c=0; c <3; c++) {
      int n=c;
      int16_t row_start_acc[2][2]={{seed[c],seed[c]},{seed[c],seed[c]}};
      if (c>0) {
	dataStart+=(((PlaneSize[c-1]+15)/16)*16);
	raw_data=dataStart;
      }
      bit=0;
      for (row=0; row < raw->rows; row++) {
	int16_t acc[2]={row_start_acc[row&1][0],row_start_acc[row&1][1]};
	for (col=0; col < raw->columns; col++){
	  for (dindex=first_decode; dindex->branch[0]||dindex->branch[1]; ) {
	    if ((bit = (bit-1) & 31) == 31){
	      for (i=0; i < 4; i++){
		bitbuf = (bitbuf << 8) + raw_data[0];
		raw_data++;
	      }
	    }
	    if (!(dindex = dindex->branch[bitbuf>>bit&1])){
	      printf("Something wrong happened when traversing the huffman tree :( \n");
	      break;
	    }
	  }
	  /* dindex->leaf is the number of bits to be read */
	  /* this number of bits are the actual diff with pixel color value n-2 */
	  bits = dindex->leaf;
	  if (bits==0)
	    diff=0;
	  else {
	    if ((bit=(bit-1)&31)==31)
	      for (b=0; b < 4; b++){
		bitbuf = (bitbuf << 8) + raw_data[0];
		raw_data++;
	      }
	    sign=diff=bitbuf>>bit &1;
	    /* 	  Attention, on ne veut garder que 16 bit max dans diff */
	    for (i=1;i<bits;i++){
	      if ((bit=(bit-1)&31)==31)
		for (b=0; b < 4; b++){
		  bitbuf = (bitbuf << 8) + raw_data[0];
		  raw_data++;
		}
	      diff=(diff<<1) + ((bitbuf>>bit)&1);
	    }
	    if (sign == 0)
	      diff -= (1<<bits) - 1;
	  }
	  acc[col&1]+=diff;
	  if (acc[col&1] >> 16 && ~acc[col&1] >>16) {/* Corrupted data? */
	    printf("Something wrong happened :( \n");
	    free(decoded_raw);
	    free(first_decode);
	    return 0;
	  }
	  if (col<2) row_start_acc[row&1][col&1]=acc[col&1];
	  decoded_raw[n+=3]=acc[col&1];
	}
      }
    }
    /* end X3FTools */
    free(first_decode);
  } else {
    X3F_ERROR("Unknown raw data format");
  }
  if (decoded_raw){
    void *tmp=raw->imageData;

    raw->imageData=decoded_raw;
    raw->rowSize=raw->columns*3;
    raw->flags|=DECODED_IMAGE;
    free(tmp);
    return 1;
  }
  return 0;
}

CAMF_LIST_ENTRY *find_camf_block(const char *name, CAMF_LIST_ENTRY *camf_list){
  CAMF_LIST_ENTRY *entry;

  for (entry=camf_list;entry!=NULL;entry=entry->next){
    if (!strcmp(entry->name, name))
      break;
  }
  return entry;

}

#define FORC(cnt) for (c=0; c < cnt; c++)
#define FORC3 FORC(3)

void X3F_compute_darkdrift(float ddft[3][3][2], X3F *x3f) {
  /* We first check if a darkdrift matrix is to be found in the camf datas */
  CAMF_LIST_ENTRY *darkdrift;
  CMbM *cmbm;
  uint16_t (*image)[3]=(uint16_t (*)[3])(((IMA*)(x3f->raw->datas))->imageData);
  uint c,d, i, j=0, row, col, width;

  width=x3f->header->columns;

  if ((darkdrift=find_camf_block("DarkDrift", x3f->camf_list))){
	cmbm=(CMbM *)darkdrift->value;

	for (d=1; d<3; d++)
	  FORC3
		for (i=0; i<2; i++)
		 ddft[d][c][i]=cmbm->matrix[j++].f;
/*   we just have to take the matrix from this entry */
  } else {
    /* compute the darkdrift from darkshieldtop and darkshieldbottom */
    for (i=0; i < 2; i++) {
	  if ((darkdrift=find_camf_block(i ? "DarkShieldBottom":"DarkShieldTop", x3f->camf_list))){
		cmbm=(CMbM *)darkdrift->value;
		for (row = cmbm->matrix[1].ui_32; row <= cmbm->matrix[3].ui_32; row++)
		  for (col = cmbm->matrix[0].ui_32; col <= cmbm->matrix[2].ui_32; col++)
			FORC3 ddft[i+1][c][1] +=  image[row*width+col][c];
		FORC3 ddft[i+1][c][1] /= (cmbm->matrix[3].ui_32-cmbm->matrix[1].ui_32+1) * (cmbm->matrix[2].ui_32-cmbm->matrix[0].ui_32+1);
	  }
    }
  }
  return;
}

char *foveon_get_param(const char *blockName, const char *name, CAMF_LIST_ENTRY *camf_list){
  CAMF_LIST_ENTRY *entry;
  PARAMETERS *params;
  uint i;

  entry=find_camf_block(blockName, camf_list);
  if (entry !=NULL){
    params=(PARAMETERS *)entry->value;
    for (i=0;i<entry->count;i++){
      if (!strcmp(params[i].name,name))
	return params[i].value;
    }
  }
  return NULL;
}

CMbM *foveon_get_matrix(const char *blockName, CAMF_LIST_ENTRY *camf_list){
  CAMF_LIST_ENTRY *entry;
  CMbM *cmbm;

  entry=find_camf_block(blockName, camf_list);
  if (entry !=NULL){
    cmbm=(CMbM *)entry->value;
    return cmbm;
  }
  return NULL;
}

int16_t *foveon_make_curve (double max, double mul, double filt)
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

void foveon_make_curves
	(int16_t **curvep, float dq[3], float div[3], float filt)
{
  double mul[3], max=0;
  int c;

  FORC3 mul[c] = dq[c]/div[c];
  FORC3 if (max < mul[c]) max = mul[c];
  FORC3 curvep[c] = foveon_make_curve (max, mul[c], filt);
}

int foveon_apply_curve (short *curve, int i)
{
  if (abs(i) >= curve[0]) return 0;
  return i < 0 ? -curve[1-i] : curve[1+i];
}

float foveon_avg (uint16_t *pix, uint range[2], float cfilt)
{
  uint i;
  float val, min=FLT_MAX, max=-FLT_MAX, sum=0;

  for (i=range[0]; i <= range[1]; i++) {
    sum += val = pix[i*4] + (pix[i*4]-pix[(i-1)*4]) * cfilt;
    if (min > val) min = val;
    if (max < val) max = val;
  }
  if (range[1] - range[0] == 1) return sum/2;
  return (sum - min - max) / (range[1] - range[0] - 1);
}

void X3F_compute_wbdiv(float div[3], const char *wbstring, CAMF_LIST_ENTRY*camf_list){
  CMbM *cmbm;
  char str[64];
  const char *cp;
  float cam_xyz[3][3], correct[3][3];
  float last[3][3], diag[3][3], num;
  int i, j, k, c;

  /* should we go now to XYZ colorspace? */
  /* try this */
  /* the transformation should be */
  /* x=cam_xyz[0][0]*R+cam_xyz[0][1]*G+cam_xyz[0][2]*B */
  /* y=cam_xyz[1][0]*R+cam_xyz[1][1]*G+cam_xyz[1][2]*B  */
  /* z=cam_xyz[2][0]*R+cam_xyz[2][1]*G+cam_xyz[2][2]*B  */

/*   X3F_cam_to_xyz(); */

  sprintf (str, "%sRGBNeutral", wbstring);
  if (!(cmbm=foveon_get_matrix(str, camf_list))) {
    /* compute the div only if %sRGBNeutral was not found */
    /* Find the camera WhiteBalance */
    if (!(cp = foveon_get_param ("WhiteBalanceIlluminants", wbstring, camf_list)))
	  printf ("Invalid white balance \"%s\"\n", wbstring);
    cmbm=foveon_get_matrix (cp, camf_list);
	for (i=0,k=0; i<3; i++)
	  for (j=0; j<3; j++)
		cam_xyz[i][j]=cmbm->matrix[k++].f;
	
	for (row=0; row<width; row++){
	  pix = image[row*width];
	  for (col=0; col<width; col++)
		FORC3 interpolation->pix[c]=cam_xyz[c][0]*pix[0]+cam_xyz[c][1]*pix[1]+cam_xyz[c][2]*pix[2];
	}

    cmbm=foveon_get_matrix (foveon_get_param ("WhiteBalanceCorrections",
											  wbstring, camf_list), camf_list);
	for (i=0,k=0; i<3; i++)
	  for (j=0; j<3; j++)
		correct[i][j]=cmbm->matrix[k++].f;
    
	memset (last, 0, sizeof last);
    for (i=0; i < 3; i++)
      for (j=0; j < 3; j++)
		FORC3 last[i][j] += correct[i][c] * cam_xyz[c][j];

#define LAST(x,y) last[(i+x)%3][(c+y)%3]
    for (i=0; i < 3; i++)
      FORC3 diag[c][i] = LAST(1,1)*LAST(2,2) - LAST(1,2)*LAST(2,1);
#undef LAST
    /* where do those constants come from? */
/*     FORC3 div[c] = diag[c][0]*0.3127 + diag[c][1]*0.329 + diag[c][2]*0.3583; */
    FORC3 div[c] = diag[c][0]*0.304 + diag[c][1]*0.346 + diag[c][2]*0.35;
    /* div[] has been computed from WBIlluminants and WBCorrections */
  } else {
	FORC3 div[c]=cmbm->matrix[c].f;
  }

  /* this is needed only if max div >1.0 */
  num = 0;
  FORC3 if (num < div[c]) num = div[c]; /* this place the max value of div[] in num */
  FORC3 div[c] /= num; /* this make the max value of div[] to 1, others may vary from 0 to 1 */
  return;
}

void X3F_compute_black_point(X3F_MATRIX_TRANSFORM *interpolation, uint width, uint height, uint16_t (*image)[3]){
  uint row, col, c, i;
  float last[3][3], fsum[3], val;
  int total[4];

#define ddft interpolation->ddft
#define black interpolation->black
#define dscr interpolation->dscr
#define cfilt interpolation->cfilt

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
  /* reset of last */
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
      FORC3 total[c] +=  image[row*width+col][c];
      total[3]++;
    }
  for (row=0; row < height; row++)
    FORC3 black[row][c] += fsum[c]/2 + total[c]/(total[3]*100.0);
/*   black is computed for each row */
/*   printf("Ready to interpolate?\n"); */

#undef ddft
#undef black
#undef dscr
#undef cfilt
}

void X3F_dcraw_interpolate_raw(X3F *x3f){
  /* what if 2 adjacent pixels are weird? */
  static const short hood[] = { -1,-1, -1,0, -1,1, 0,-1, 0,1, 1,-1, 1,0, 1,1 };
  static const int rgb_cam[3][3]={{1,0,0}, {0,1,0},{0,0,1}};

  short prev[3], *curve[8], (*shrink)[3];
  uint16_t *pix;
  float cfilt=0;
  float last[3][3], trans[3][3], diag[3][3], ddft[3][3][2], div[3];
  float cam_xyz[3][3], correct[3][3], sgain[33*49][3], chroma_dq[3], color_dq[3], ppm[3][3][3];
  float (*sgrow)[3];
  float fsum[3], val, frow, num;
  int row, width, height, col, c, i, j, k, diff, sgx, irow, sum, min, max, limit;
  int dscr[2][2], dstb[4];
  int (*smrow[7])[3], total[4], ipix[3];
  int work[3][3], smlast, smred, smred_p=0, dev[3];
  uint16_t satlev[3], keep[4], active[4];
  int32_t *dim, *badpix;
  double dsum=0, trsum[3];
  char str[128];
  const char* cp;
  IMA *ima=(IMA *)x3f->raw->datas;
  CAMF_LIST_ENTRY *camf_list=x3f->camf_list;
  uint16_t (*image)[3]=(uint16_t (*)[3])ima->imageData;
  CMbM *cmbm;
  X3F_MATRIX_TRANSFORM *interpolation=x3f->interpolation;

  fprintf (stderr,_("Foveon interpolation...\n"));

  width=ima->columns;
  height=ima->rows;

  /* foveon_fixed sets a copy of the matrix in the first arg ptr */
  /* foveon_camf_param returns TRUE if the CMbP is found, FALSE otherwise */
  /* we need to replace foveon_fixed with a func returning the camf matrix */
  /* we need to replace foveon_camf_param with find_camf_matrix */

/*   foveon_get_matrix(dscr, "DarkShieldColRange", camf_list); */
/*   dscr=(MATRIX[2][2])cmbm->matrix; */
/*   foveon_get_matrix (ppm[0][0], "PostPolyMatrix", camf_list); /\* not present in TRUE camf *\/ */
/*   foveon_get_matrix (satlev, "SaturationLevel", camf_list); */
/*   foveon_get_matrix (keep, "KeepImageArea", camf_list); */
/*   foveon_get_matrix (active, "ActiveImageArea", camf_list); */

/*   if (foveon_camf_param ("IncludeBlocks", "ColumnFilter")) */ /* not present in TRUE camf; have ColumnFilterRGB and ColumnFilterSquareRGB instead */
/*   if (foveon_get_matrix (&cfilt, "ColumnFilter", camf_list)){ */
    /* we should compute cfilt from ColumnFilterRGB which is float[3]
	   matrix (values so far are 0) */
/*     X3F_MSG("No ColumnFilter entry\n"); */
/*   } */
  if (!(interpolation=malloc(sizeof(*interpolation)))){
	X3F_MEM_ERROR("X3F_dcraw_interpolate_raw", "interpolation");
	return;
  }
  interpolation->cfilt=0.8;


  memset (interpolation->ddft, 0, sizeof interpolation->ddft);
  printf("Computing DarkDrift\n");
  X3F_compute_darkdrift(interpolation->ddft, x3f);
  for (j=0; j<3; j++)
	for (c=0; c<3; c++)
	  for (i=0; i<2; i++)
		printf("ddft[%d][%d][%d]: %f\n", j, c, i, interpolation->ddft[j][c][i]);

  cmbm=(CMbM*)((find_camf_block("DarkShieldColRange", camf_list))->value);
  for (c=0,j=0; c<2; c++)
	for (i=0; i<2; i++){
	  interpolation->dscr[c][i]=cmbm->matrix[j++].ui_32;
	  printf("dscr[%d][%d]: %d\n", c, i, interpolation->dscr[c][i]);
/* 	  printf("matrix[%d][%d]: %d\n", c, i, cmbm->matrix[c*2+i].ui_32); */
	}

  printf("Computing white balance divider\n");
  X3F_compute_wbdiv(interpolation->wbdiv, (char*)x3f->header->whiteBalanceString, camf_list);
  FORC3 printf("wbdiv[%d]: %f\n", c, interpolation->wbdiv[c]);

  /* rgb_cam is identity matrix */
  /* what is trans? */
  /* it's used to modify last then last is used to recompute trans values */
  memset (trans, 0, sizeof trans);
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      FORC3 trans[i][j] += rgb_cam[i][c] * last[c][j] * interpolation->wbdiv[j];
  FORC3 trsum[c] = trans[c][0] + trans[c][1] + trans[c][2];
  dsum = (6*trsum[0] + 11*trsum[1] + 3*trsum[2]) / 20;
  for (i=0; i < 3; i++)
    FORC3 last[i][c] = trans[i][c] * dsum / trsum[i];
  memset (trans, 0, sizeof trans);
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      FORC3 trans[i][j] += (i==c ? 32 : -1) * last[c][j] / 30;
  /* trans will be used for colorspace transformation */

  /* Here curves are computed */
  printf("Making curves\n");
  if (!(cmbm=foveon_get_matrix ("ColorDQ", camf_list)))
	cmbm=foveon_get_matrix("ColorDQCamRGB", camf_list);
  FORC3 color_dq[c]=cmbm->matrix[c].f;
  foveon_make_curves (interpolation->curve, color_dq, interpolation->wbdiv, interpolation->cfilt);

  cmbm=foveon_get_matrix  ("ChromaDQ", camf_list);
  FORC3 chroma_dq[c]=cmbm->matrix[c].f;
  FORC3 chroma_dq[c] /= 3;
  foveon_make_curves (interpolation->curve+3, chroma_dq, interpolation->wbdiv, interpolation->cfilt);
  FORC3 dsum += chroma_dq[c] / interpolation->wbdiv[c];
  curve[6] = foveon_make_curve (dsum, dsum, cfilt);
  curve[7] = foveon_make_curve (dsum*2, dsum*2, cfilt);

  /* spatial gain */
/*   cp = foveon_get_param ("SpatialGainTable",  (char*)x3f->header->whiteBalanceString, camf_list); */
/*   printf("%s\n", cp); */
/*   cmbm=foveon_get_matrix(cp?cp:"SpatialGain", camf_list); */
/*   for (i=0, k=0; i<33*49;i++) */
/* /\* 	for (j=0; j<49; j++) *\/ */
/* 	  FORC3 sgain[i][c]=cmbm->matrix[k++].f; */
/*   dim=cmbm->planeElements; */
/*   sgrow = (float (*)[3]) calloc (dim[1], sizeof *sgrow); */
/*   sgx = (width + dim[1]-2) / (dim[1]-1); */

  /* Black */
/*   X3F_compute_black_point(interpolation,width, height, image); */
/*   printf("Black\n"); */

  /* START OF INTERPOLATION */
/*   for (row=0; row < height; row++) { */
/*     for (i=0; i < 6; i++) */
/*       interpolation->ddft[0][0][i] = interpolation->ddft[1][0][i] + */
/* 	row / (height-1.0) * (interpolation->ddft[2][0][i] - interpolation->ddft[1][0][i]); */
/*     /\* pix point to start of row *\/ */
/*     pix = image[row*width]; */
/*     memcpy (prev, pix, sizeof prev); */
/* /\*     frow = row / (height-1.0) * (dim[0]-1); *\/ */
/* /\*     if ((irow = frow) == dim[0]-1) irow--; *\/ */
/* /\*     frow -= irow; *\/ */
/* /\*     for (i=0; i < dim[1]; i++) *\/ */
/* /\*       FORC3 sgrow[i][c] = sgain[ irow*dim[1]+i][c] * (1-frow) + *\/ */
/* /\* 			  sgain[(irow+1)*dim[1]+i][c] * frow; *\/ */
/*     for (col=0; col < width; col++) { */
/*       /\* FIRST ADJUST BLACK *\/ */
/*       FORC3 { */
/* 		diff = pix[c] - prev[c]; */
/* 		prev[c] = pix[c]; */
/* 		ipix[c] = pix[c] + floor ((diff + (diff*diff >> 14)) * interpolation->cfilt */
/* 								  - interpolation->ddft[0][c][1] - interpolation->ddft[0][c][0] * ((float) col/width - 0.5) */
/* 								  - interpolation->black[row][c] ); */
/*       } */
/* /\*       THEN ? *\/ */
/*       FORC3 { */
/* /\* 		work[0][c] = ipix[c] * ipix[c] >> 14; *\/ */
/* /\* 		work[2][c] = ipix[c] * work[0][c] >> 14; *\/ */
/* /\* 		work[1][2-c] = ipix[(c+1) % 3] * ipix[(c+2) % 3] >> 14; *\/ */
/* /\*       } *\/ */
/* /\*       /\\* APPLY SPATIAL GAIN (use PostPolyMatrix) *\\/ *\/ */
/* /\* 	  cmbm=foveon_get_matrix("PostPolyMatrix", camf_list); *\/ */
/* /\* 	  for (i=0,k=0; i<3; i++) *\/ */
/* /\* 		for (j=0; j<3; j++) *\/ */
/* /\* 		  FORC3 ppm[i][j][c]=cmbm->matrix[k++].f; *\/ */
/* /\*       FORC3 { *\/ */
/* /\* 	for (val=i=0; i < 3; i++) *\/ */
/* /\* 	  for (  j=0; j < 3; j++) *\/ */
/* /\* 	    val += ppm[j][i][c] * work[i][j]; *\/ */
/* /\* 	ipix[c] = floor ((ipix[c] + floor(val)) * *\/ */
/* /\* 		( sgrow[col/sgx  ][c] * (sgx - col%sgx) + *\/ */
/* /\* 		  sgrow[col/sgx+1][c] * (col%sgx) ) / sgx / div[c]); *\/ */
/* /\* 	if (ipix[c] > 32000) ipix[c] = 32000; *\/ */
/* 		pix[c] = ipix[c]; */
/* 	  } */
/*       pix += 3; */
/* 	} */
/*   } */
  /* All pixels have been interpolated */
/*   free (black); */
/*   free (sgrow); */
/* /\*   free (sgain); *\/ */

/* /\*   /\\* correction of badPixels (averaging neighboor pixels) *\\/ *\/ */
/* /\*   if ((badpix = (unsigned int *) foveon_get_matrix (dim, "BadPixels", camf_list))) { *\/ */
/* /\*     for (i=0; i < dim[0]; i++) { *\/ */
/* /\*       col = (badpix[i] >> 8 & 0xfff) - keep[0]; *\/ */
/* /\*       row = (badpix[i] >> 20       ) - keep[1]; *\/ */
/* /\*       if ((unsigned)(row-1) > height-3 || (unsigned)(col-1) > width-3) *\/ */
/* /\* 	continue; *\/ */
/* /\*       memset (fsum, 0, sizeof fsum); *\/ */
/* /\*       for (sum=j=0; j < 8; j++) *\/ */
/* /\* 	if (badpix[i] & (1 << j)) { *\/ */
/* /\* 	  FORC3 fsum[c] += (short) *\/ */
/* /\* 		image[(row+hood[j*2])*width+col+hood[j*2+1]+c]; *\/ */
/* /\* 	  sum++; *\/ */
/* /\* 	} *\/ */
/* /\*       if (sum) FORC3 image[row*width+col+c] = fsum[c]/sum; *\/ */
/* /\*     } *\/ */
/* /\*     free (badpix); *\/ */
/* /\*   } *\/ */


/*   /\* Smoothing and sharpening (only red channel) *\/ */

/*   /\* Array for 5x5 Gaussian averaging of red values *\/ */
/* /\*   smrow[6] = (int (*)[3]) calloc (width*5, sizeof **smrow); *\/ */
/* /\* /\\*   merror (smrow[6], "foveon_interpolate()"); *\\/ *\/ */
/* /\*   for (i=0; i < 5; i++) *\/ */
/* /\*     smrow[i] = smrow[6] + i*width; *\/ */

/* /\*   /\\* Sharpen the reds against these Gaussian averages *\\/ *\/ */
/* /\*   for (smlast=-1, row=2; row < height-2; row++) { *\/ */
/* /\*     while (smlast < row+2) { *\/ */
/* /\*       for (i=0; i < 6; i++) *\/ */
/* /\* 	smrow[(i+5) % 6] = smrow[i]; *\/ */
/* /\*       pix = image[++smlast*width+2]; *\/ */
/* /\*       for (col=2; col < width-2; col++) { *\/ */
/* /\* 	smrow[4][col][0] = *\/ */
/* /\* 	  (pix[0]*6 + (pix[-4]+pix[4])*4 + pix[-8]+pix[8] + 8) >> 4; *\/ */
/* /\* 	pix += 3; *\/ */
/* /\*       } *\/ */
/* /\*     } *\/ */
/* /\*     pix = image[row*width+2]; *\/ */
/* /\*     for (col=2; col < width-2; col++) { *\/ */
/* /\*       smred = ( 6 *  smrow[2][col][0] *\/ */
/* /\* 	      + 4 * (smrow[1][col][0] + smrow[3][col][0]) *\/ */
/* /\* 	      +      smrow[0][col][0] + smrow[4][col][0] + 8 ) >> 4; *\/ */
/* /\*       if (col == 2) *\/ */
/* /\* 	smred_p = smred; *\/ */
/* /\*       i = pix[0] + ((pix[0] - ((smred*7 + smred_p) >> 3)) >> 3); *\/ */
/* /\*       if (i > 32000) i = 32000; *\/ */
/* /\*       pix[0] = i; *\/ */
/* /\*       smred_p = smred; *\/ */
/* /\*       pix += 3; *\/ */
/* /\*     } *\/ */
/* /\*   } *\/ */

  /* Adjust the brighter pixels for better linearity */
  /* USE: SatLevel */
  min = 0xffff;
  cmbm=foveon_get_matrix("SaturationLevel", camf_list);
  FORC3 satlev[c]=cmbm->matrix[c].ui_16;
  FORC3 {
    i = satlev[c];
	printf("SaturationLevel[%d]: %d\n", c, satlev[c]);
    if (min > i) min = i;
  }
  limit = min * 9 >> 4;
  for (pix=image[0]; pix < image[height*width]; pix+=3) {
    if (pix[0] <= limit || pix[1] <=  limit || pix[2] <= limit)
      continue;
    min = max = pix[0];
    for (c=1; c < 3; c++) {
      if (min > pix[c]) min = pix[c];
      if (max < pix[c]) max = pix[c];
    }
    if (min >= limit*2) {
      pix[0] = pix[1] = pix[2] = max;
    } else {
/* 	  printf("Here\n"); */
      i = 0x4000 - ((min - limit) << 14) / limit;
      i = 0x4000 - (i*i >> 14);
      i = i*i >> 14;
      FORC3 pix[c] += (max - pix[c]) * i >> 14;
    }
  }

/*   /\* Denoising *\/ */
/*   /\* Use: curve[7] and curve[6] (relative to chromaDQ)*\/ */
/* /\* */
/*    Because photons that miss one detector often hit another, */
/*    the sum R+G+B is much less noisy than the individual colors. */
/*    So smooth the hues without smoothing the total. */
/*  *\/ */
/* /\*   for (smlast=-1, row=2; row < height-2; row++) { *\/ */
/* /\*     while (smlast < row+2) { *\/ */
/* /\*       for (i=0; i < 6; i++) *\/ */
/* /\* 	smrow[(i+5) % 6] = smrow[i]; *\/ */
/* /\*       pix = image[++smlast*width+2]; *\/ */
/* /\*       for (col=2; col < width-2; col++) { *\/ */
/* /\* 	FORC3 smrow[4][col][c] = (pix[c-4]+2*pix[c]+pix[c+4]+2) >> 2; *\/ */
/* /\* 	pix += 3; *\/ */
/* /\*       } *\/ */
/* /\*     } *\/ */
/* /\*     pix = image[row*width+2]; *\/ */
/* /\*     for (col=2; col < width-2; col++) { *\/ */
/* /\*       FORC3 dev[c] = -foveon_apply_curve (curve[7], pix[c] - *\/ */
/* /\* 	((smrow[1][col][c] + 2*smrow[2][col][c] + smrow[3][col][c]) >> 2)); *\/ */
/* /\*       sum = (dev[0] + dev[1] + dev[2]) >> 3; *\/ */
/* /\*       FORC3 pix[c] += dev[c] - sum; *\/ */
/* /\*       pix += 3; *\/ */
/* /\*     } *\/ */
/* /\*   } *\/ */
/* /\*   for (smlast=-1, row=2; row < height-2; row++) { *\/ */
/* /\*     while (smlast < row+2) { *\/ */
/* /\*       for (i=0; i < 6; i++) *\/ */
/* /\* 	smrow[(i+5) % 6] = smrow[i]; *\/ */
/* /\*       pix = image[++smlast*width+2]; *\/ */
/* /\*       for (col=2; col < width-2; col++) { *\/ */
/* /\* 	FORC3 smrow[4][col][c] = *\/ */
/* /\* 		(pix[c-8]+pix[c-4]+pix[c]+pix[c+4]+pix[c+8]+2) >> 2; *\/ */
/* /\* 	pix += 3; *\/ */
/* /\*       } *\/ */
/* /\*     } *\/ */
/* /\*     pix = image[row*width+2]; *\/ */
/* /\*     for (col=2; col < width-2; col++) { *\/ */
/* /\*       for (total[3]=375, sum=60, c=0; c < 3; c++) { *\/ */
/* /\* 	for (total[c]=i=0; i < 5; i++) *\/ */
/* /\* 	  total[c] += smrow[i][col][c]; *\/ */
/* /\* 	total[3] += total[c]; *\/ */
/* /\* 	sum += pix[c]; *\/ */
/* /\*       } *\/ */
/* /\*       if (sum < 0) sum = 0; *\/ */
/* /\*       j = total[3] > 375 ? (sum << 16) / total[3] : sum * 174; *\/ */
/* /\*       FORC3 pix[c] += foveon_apply_curve (curve[6], *\/ */
/* /\* 		((j*total[c] + 0x8000) >> 16) - pix[c]); *\/ */
/* /\*       pix += 3; *\/ */
/* /\*     } *\/ */
/* /\*   } *\/ */

/*   /\* Use: curve[0]->curve[2] (colorDQ) *\/ */
/*   /\* Transform the image to a different colorspace *\/ */
/*   for (pix=image[0]; pix < image[height*width]; pix+=3) { */
/*     FORC3 pix[c] -= foveon_apply_curve (interpolation->curve[c], pix[c]); */
/*     sum = (pix[0]+pix[1]+pix[1]+pix[2]) >> 2; */
/*     FORC3 pix[c] -= foveon_apply_curve (interpolation->curve[c], pix[c]-sum); */
/*     FORC3 { */
/*       for (dsum=i=0; i < 3; i++) */
/* 	dsum += trans[c][i] * pix[i]; */
/*       if (dsum < 0)  dsum = 0; */
/*       if (dsum > 24000) dsum = 24000; */
/*       ipix[c] = dsum + 0.5; */
/*     } */
/*     FORC3 pix[c] = ipix[c]; */
/*   } */

/*   /\* Smooth the image bottom-to-top and save at 1/4 scale *\/ */
/* /\*   shrink = (short (*)[3]) calloc ((width/4) * (height/4), sizeof *shrink); *\/ */
/* /\* /\\*   merror (shrink, "foveon_interpolate()"); *\\/ *\/ */
/* /\*   for (row = height/4; row--; ) *\/ */
/* /\*     for (col=0; col < width/4; col++) { *\/ */
/* /\*       ipix[0] = ipix[1] = ipix[2] = 0; *\/ */
/* /\*       for (i=0; i < 4; i++) *\/ */
/* /\* 	for (j=0; j < 4; j++) *\/ */
/* /\* 	  FORC3 ipix[c] += image[(row*4+i)*width+col*4+j][c]; *\/ */
/* /\*       FORC3 *\/ */
/* /\* 	if (row+2 > height/4) *\/ */
/* /\* 	  shrink[row*(width/4)+col][c] = ipix[c] >> 4; *\/ */
/* /\* 	else *\/ */
/* /\* 	  shrink[row*(width/4)+col][c] = *\/ */
/* /\* 	    (shrink[(row+1)*(width/4)+col][c]*1840 + ipix[c]*141 + 2048) >> 12; *\/ */
/* /\*     } *\/ */
/* /\*   /\\* From the 1/4-scale image, smooth right-to-left *\\/ *\/ */
/* /\*   for (row=0; row < (height & ~3); row++) { *\/ */
/* /\*     ipix[0] = ipix[1] = ipix[2] = 0; *\/ */
/* /\*     if ((row & 3) == 0) *\/ */
/* /\*       for (col = width & ~3 ; col--; ) *\/ */
/* /\* 	FORC3 smrow[0][col][c] = ipix[c] = *\/ */
/* /\* 	  (shrink[(row/4)*(width/4)+col/4][c]*1485 + ipix[c]*6707 + 4096) >> 13; *\/ */

/* /\*   /\\* Then smooth left-to-right *\\/ *\/ */
/* /\*     ipix[0] = ipix[1] = ipix[2] = 0; *\/ */
/* /\*     for (col=0; col < (width & ~3); col++) *\/ */
/* /\*       FORC3 smrow[1][col][c] = ipix[c] = *\/ */
/* /\* 	(smrow[0][col][c]*1485 + ipix[c]*6707 + 4096) >> 13; *\/ */

/* /\*   /\\* Smooth top-to-bottom *\\/ *\/ */
/* /\*     if (row == 0) *\/ */
/* /\*       memcpy (smrow[2], smrow[1], sizeof **smrow * width); *\/ */
/* /\*     else *\/ */
/* /\*       for (col=0; col < (width & ~3); col++) *\/ */
/* /\* 	FORC3 smrow[2][col][c] = *\/ */
/* /\* 	  (smrow[2][col][c]*6707 + smrow[1][col][c]*1485 + 4096) >> 13; *\/ */

/* /\*   /\\* Adjust the chroma toward the smooth values *\\/ *\/ */
/* /\*     /\\* SATURATION ? *\\/ *\/ */
/* /\*     /\\* Use: curve[3]->curve[5] (ChromaDQ) *\\/ *\/ */
/* /\*     for (col=0; col < (width & ~3); col++) { *\/ */
/* /\*       for (i=j=30, c=0; c < 3; c++) { *\/ */
/* /\* 	i += smrow[2][col][c]; *\/ */
/* /\* 	j += image[row*width+col][c]; *\/ */
/* /\*       } *\/ */
/* /\*       if (i!=0) *\/ */
/* /\* 	j = (j << 16) / i; *\/ */
/* /\*       else *\/ */
/* /\* 	j = (j<<16); *\/ */
/* /\*       for (sum=c=0; c < 3; c++) { *\/ */
/* /\* 	ipix[c] = foveon_apply_curve (curve[c+3], *\/ */
/* /\* 	  ((smrow[2][col][c] * j + 0x8000) >> 16) - image[row*width+col][c]); *\/ */
/* /\* 	sum += ipix[c]; *\/ */
/* /\*       } *\/ */
/* /\*       sum >>= 3; *\/ */
/* /\*       FORC3 { *\/ */
/* /\* 	i = image[row*width+col][c] + ipix[c] - sum; *\/ */
/* /\* 	if (i < 0) i = 0; *\/ */
/* /\* 	image[row*width+col][c] = i; *\/ */
/* /\*       } *\/ */
/* /\*     } *\/ */
/* /\*   } *\/ */
/* /\*   free (shrink); *\/ */
/* /\*   free (smrow[6]); *\/ */
/* /\*   for (i=0; i < 8; i++) *\/ */
/* /\*     free (curve[i]); *\/ */

/* /\*   /\\* Trim off the black border *\\/ *\/ */
/* /\*   active[1] -= keep[1]; *\/ */
/* /\*   active[3] -= 2; *\/ */
/* /\*   i = active[2] - active[0]; *\/ */
/* /\*   for (row=0; row < active[3]-active[1]; row++) *\/ */
/* /\*     memcpy (image[row*i], image[(row+active[1])*width+active[0]], *\/ */
/* /\* 	 i * sizeof *image); *\/ */
/* /\*   width = i; *\/ */
/* /\*   height = row;     *\/ */
}

X3F *X3F_load_full_x3f(char *filename){
  FILE *fp;
  uint i;
  DIR_ENTRY *subsection=NULL;
  X3F *x3f=NULL;

  X3F_MSG("Loading \n");
  if((fp = fopen(filename, "rb")) == NULL) {
    X3F_ERROR("Cannot open file\n");
    return NULL;
  }

  if (!(x3f=X3F_init())){
    X3F_ERROR("Error initializing X3F structure");
    return NULL;
  }

  if (!(x3f->header=X3F_parse_header( fp))){
    X3F_ERROR("Error parsing header");
    X3F_free(x3f);
    return NULL;
  }

  fseek(fp, -4, SEEK_END);
  fread(&x3f->dir_offset, 1, sizeof(x3f->dir_offset), fp);
  if (ferror(fp)){
    X3F_READ_ERROR("X3F_load_full_x3f");
    X3F_free(x3f);
    return NULL;
  }

  fseek(fp, x3f->dir_offset, SEEK_SET);
  /* Read the directory section */
  if (!(x3f->dir_section=X3F_read_dir_section(fp))){
    X3F_ERROR("Error reading X3F directory");
    X3F_free(x3f);
    return NULL;
  }
  
  /* Read each directory entry */
#define SECTION x3f->dir_section
/*   if (!(stored_entries=malloc(sizeof(*stored_entries)*SECTION->dirEntryCount))){ */
/*     X3F_MEM_ERROR("X3F_load_full_x3f","dir_entries"); */
/*     X3F_free(x3f); */
/*     return NULL; */
/*   } */
  if (!(x3f->dir_entries=malloc(sizeof(*x3f->dir_entries)*SECTION->dirEntryCount))){
    X3F_MEM_ERROR("X3F_load_full_x3f","dir_entries");
    X3F_free(x3f);
    return NULL;
  }
  /* parse each entry */
  for (i=SECTION->dirEntryCount;i>0;i--){
/*     int unneeded=0; */

    fseek(fp,x3f->dir_offset+12*i, SEEK_SET);
    if (!(subsection=X3F_read_dir_entry(fp))){
      X3F_MSG("Error reading directory entry");
 /*      free(stored_entries); */
      X3F_free(x3f);
      return NULL;
    }
    subsection->datas=NULL;
    fseek(fp, subsection->offset, SEEK_SET);
    if (subsection->type==X3F_PROP){
      PROP *prop=NULL;
      X3F_MSG("Directory entry of type PROP");
      prop=X3F_read_prop(fp,x3f,subsection->dataLength);
      subsection->datas=(void *)prop;
      /* prop does not need (actually) to be decoded */
      /* Theoretically, it's not necessary to store the PROP sections */
      /* We just create a hashtable in the main X3F struct, this way all property pairs will be available in one single table*/
      /* so just free the subsection */
/*       free(subsection); */
      /* and decrease the stored_sections_count */
/*       unneeded=1; */
    } else if (subsection->type==X3F_IMAG||subsection->type==X3F_IMA2){
      X3F_MSG("Directory entry of type IMA");
      IMA *image=NULL;
      image=X3F_read_ima(fp,subsection->dataLength);
      subsection->datas=(void *)image;
      if (image->imageDataType == X3F_DATA_TYPE_PROCESSED){
	if (!x3f->thumbnail && (image->dataFormat == X3F_DATA_FORMAT_PNM_THUMBNAIL ||
				image->rowSize!=0)) {
	  /* seems that jpeg thumbnail have a rowSize!=0 */
	  /* thumbnails should be available as soon as possible */
	  /* so decode them right now */
	  X3F_decode_thumbnail(image,subsection->dataLength);
	  x3f->thumbnail=subsection;
	} else if (!x3f->preview && (image->dataFormat == X3F_DATA_FORMAT_HUFFMAN_PREVIEW ||
				     image->rowSize==0)) {
	  /* seems that all jpeg preview have a rowSize==0 */
	  /* preview should be available quickly */
	  /* decode them now? */
	  X3F_decode_preview(image,subsection->dataLength);
	  printf("Preview size is: %dx%d\n", image->rows, image->columns);
	  x3f->preview=subsection;
	}
      } else if (!x3f->raw){
	/* raw should be available on demand */
	/* decode them later */

	x3f->raw=subsection;
      }
    } else if (subsection->type==X3F_CAMF){
      X3F_MSG("Directory entry of type CAMF");
      CAMF *camf;
      camf=X3F_read_camf(fp, subsection->dataLength);
      /* OK, camf datas will be decoded when parsing */
      x3f->camf_list=X3F_fill_camf_list(subsection->dataLength-CAMF_HEADER_SIZE,
					camf->camf_data, 
					x3f->camf_list);
      subsection->datas=(void *)camf;
    } else {
      X3F_MSG("Unknown entry type: -> skipping");
      printf("Entry type was: %x\n", subsection->type);
    }
/*     if(!unneeded){ */
/*       stored_entries[stored_sections_count]=subsection; */
/*       stored_sections_count++; */
    x3f->dir_entries[i-1]=subsection;
/*     } */
  }
/*   SECTION->dirEntryCount=stored_sections_count; */
/*   memcpy(x3f->dir_entries, stored_entries, sizeof(*x3f->dir_entries)*SECTION->dirEntryCount);  */
/*   free(stored_entries); */
#undef SECTION
  
  if( fclose( fp )) 
    X3F_ERROR("File: close error.\n");
  
  X3F_MSG("X3F has been loaded in memory\n");
  return x3f;
}

DIR_ENTRY *X3F_get_section(X3F *x3f, uint32_t sectionType){
  int i;
  for (i=0; i<x3f->dir_section->dirEntryCount;i++) {
    /* Hum, image entry should have type IMAG or IMA2; test for IMA */
    if (x3f->dir_entries[i]->type==sectionType||
	(sectionType==X3F_IMA && ((x3f->dir_entries[i]->type<<2)==sectionType)))
      return x3f->dir_entries [i];
  }
  return NULL;
}

