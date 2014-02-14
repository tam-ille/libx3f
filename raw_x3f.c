/* This libx3f will only compile on little-endian processors, which are the vast majority of personnal cpu */

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
#include <jpeglib.h>
#include "raw_x3f.h"

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
	if (camf->camf_data)
	  free(camf->camf_data);
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
      uint32_t *planeElements;
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
	  matrix[i].ui_32=*((uint32_t *)(dataPtr+dataStartOffset+i*sizeof(uint16_t)))&0xffff;
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
/*   free(camf_data); */
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
  camf->dataSize=dataSize;
  printf("dataSize: %d\n", dataSize);
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

  char *decoded_data=NULL, *ptr;
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
  ima->max[0]=ima->max[1]=ima->max[2]=0;
  ima->min[0]=ima->min[1]=ima->min[1]=65535;
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
  int16_t pix[3], *pvalues;
  int16_t (*decoded_raw)[4];

  if ((raw->flags & DECODED_IMAGE)== DECODED_IMAGE)
    return 1;
  memset(raw->max, 0, sizeof raw->max);
  memset(raw->min, 4096, sizeof raw->min);
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
    if (!(decoded_raw=(int16_t (*)[4]) calloc (raw->rows*raw->columns*4, sizeof (*decoded_raw)))){
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
/* 				    pix[c]=(int16_t)pix[c]>0?pix[c]:0; */
			if (pix[c] >> 16 && ~pix[c] >> 16){
			  printf("Oups\n");
			  free(decoded_raw);
			  return 0; /* Corrupted data? */
			}
		  }
		}
		for (c=0; c <3; c++) {
		  decoded_raw[row*raw->columns+col][c] = pix[c];
		  if (pix[c]>raw->max[c]&&pix[c]>0) raw->max[c]=pix[c];
		  if (pix[c]<raw->min[c]) raw->min[c]=pix[c]>0?pix[c]:0;
		}
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
    if (!(decoded_raw=(int16_t (*)[4]) calloc (raw->rows*raw->columns*3, sizeof (*decoded_raw)))){
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
		  decoded_raw[row*raw->columns+col][c] = acc[col&1];
		  if (acc[col&1]>raw->max[c]) raw->max[c]=acc[col&1];
		  if (acc[col&1]<raw->min[c]) raw->min[c]=acc[col&1];
		}
      }
    }
    /* end X3FTools */
    free(first_decode);
  } else {
    X3F_ERROR("Unknown raw data format");
  }
  if (decoded_raw){
	free(raw->imageData);
    raw->imageData=(void *)decoded_raw;
    raw->rowSize=raw->columns*4;
    raw->flags|=DECODED_IMAGE;
    printf ("Max R: %d\nMax G: %d\nMax B: %d\n", raw->max[0], raw->max[1], raw->max[2]);
    printf ("Min R: %d\nMin G: %d\nMin B: %d\n", raw->min[0], raw->min[1], raw->min[2]);
    return 1;
  }
  return 0;
}

void create_floatimage(IMA *ima)
{ 
	unsigned int size, i, c;
	uint16_t (*img)[4]=(uint16_t (*)[4])ima->imageData;
	double val, saturation, dark, (*image)[4], image_max=0, image_min=DBL_MAX;
	size = ima->rows*ima->columns;

	image=(double (*)[4])calloc(size*4, sizeof(*image));
	for (i=0; i < size*4; i++) 
		{ 
//type cast "integerimage" from integer to double and set val equal to integerimage
		val = (double)img[0][i]/4096;
		if (!val) continue;  

//keep track of maximum and minimum values
		if (val> image_max) image_max = val;
		if (val< image_min) image_min = val;

//set the double floating point array "image" equal to val
		image[0][i]=val;
		} 
	free(ima->imageData);
	ima->imageData=(void*)image;
	ima->flags |= FLOAT_IMAGE;
} 


CAMF_LIST_ENTRY *X3F_get_camf_entry(CAMF_LIST_ENTRY *camf_list, char *entry_name){
  /* go through the camf_entries to find the entry with name entry_name */
  CAMF_LIST_ENTRY *entry;

  if (!entry_name) return NULL;
  for (entry=camf_list; entry!=NULL; entry=entry->next)
    if (!(strcmp(entry->name, entry_name))){
      /* we found what we were looking for ! */
      printf ("Found %s\n", entry->name);
      break;
    }
  return entry;	
}

char *foveon_get_param(CAMF_LIST_ENTRY *camf_list, char *blockName, const char *name){
  CAMF_LIST_ENTRY *entry;
  PARAMETERS *params;
  uint i;

  entry=X3F_get_camf_entry(camf_list, blockName);
  if (entry !=NULL){
    params=(PARAMETERS *)entry->value;
    for (i=0;i<entry->count;i++){
      if (!strcmp(params[i].name,name)){
	printf("Found parameter %s: %s\n", name, params[i].value);
	return params[i].value;
      }
    }
  } else {
    printf("Oups\n");
  }
  return NULL;
}

char *foveon_get_property(PROPERTY *prop_list, char *name){
  char *value=NULL, *prop_name=NULL;
  PROPERTY *prop;

  for (prop=prop_list;prop!=NULL;prop=prop->next){
    prop_name=g_utf16_to_utf8(prop->name,-1,NULL,NULL,NULL);
    if (!strcmp(prop_name, name)){
      value=g_utf16_to_utf8(prop->value, -1, NULL,NULL,NULL);
      printf("Found property %s = %s\n", prop_name, value);
      break;
    }
  }
  g_free(prop_name);
  return value;
}

void get_matrix(CMbM *cmbm, void *ptr,  int valueCount){
  int i;
  int32_t *mat;

  mat=malloc(valueCount*4);
  for (i=0; i<valueCount;i++)
	switch (cmbm->dataType) {
	case 0:
	  mat[i]=cmbm->matrix[i].ui_16;
	  break;
	case 3:
	  mat[i]=cmbm->matrix[i].ui_32;
	  break;
	case 6:
	  mat[i]=cmbm->matrix[i].ui_8 & 0xffff;
	  break;
	default:
	  mat[i]=cmbm->matrix[i].ui_32;
	}
  memcpy(ptr, mat, valueCount*4);
  free(mat);
  return;
}


int X3F_raw_interpolate(X3F *x3f){
  IMA *ima=(IMA *)x3f->raw->datas;

/*   printf ("DATA Type:%s\n DATA Format: %s", ima->imageDataType, ima->dataFormat); */
  if (ima->imageDataType == X3F_DATA_TYPE_RAW) {
    if (ima->dataFormat == X3F_DATA_FORMAT_RAW)
      interpolate(x3f);
    else
      foveon_f20_interpolate(x3f);
  } else if (ima->imageDataType == X3F_DATA_TYPE_RAW_SD1)
    foveon_f20_interpolate(x3f);
  else
    X3F_ERROR("Unknown raw data type\n");
  return 1;
}

/* this is the interpolation process for pre-TRUEII x3f */
int X3F_foveon_interpolate(X3F *x3f){

/*   bad_pixels_correction(x3f,"BadPixels"); */
/*   compute_black_point(); */
/*   compute_white_point(); */
/*   apply_spatial_gain(); */
/*   gamma_correction(); */
/*   apply_colors_curves(); */
/*   sharpen_reds(); */
/*   apply_chroma_curve(); */
/*   X3F_cam2xyz(ima, (char *)x3f->header->whiteBalanceString, x3f->camf_list); */
/*   X3F_RGBNeutral(ima, (char *)x3f->header->whiteBalanceString, x3f->camf_list); */
  color_correction(x3f);
/*   colorspace_transformation(); */
/*   black_border_suppression(); */

  return 0;
}

/* this is the interpolation process for TRUEII x3f */
int X3F_foveon_TRUE_interpolate(X3F *x3f){

/*   bad_pixels_correction(x3f,"BadPixels"); */
/*   compute_black_point(); */
/*   compute_white_point(); */
/*   apply_spatial_gain(); */
/*   gamma_correction(); */
/*   apply_colors_curves(); */
/*   sharpen_reds(); */
/*   apply_chroma_curve(); */
/*   X3F_cam2xyz(ima, (char *)x3f->header->whiteBalanceString, x3f->camf_list); */
  f20_color_correction(x3f);
/*   colorspace_transformation(); */
/*   black_border_suppression(); */

  return 0;
}

/* this is the interpolation process for F20 captor x3f */
int X3F_foveon_F20_interpolate(X3F *x3f){

  /*   bad_pixels_SD1_correction(x3f,"BadPixelsF20"); */
  /*   compute_black_point(); */
  /*   compute_white_point(); */
  /*   apply_spatial_gain(); */
  /*   gamma_correction(); */
  /*   apply_colors_curves(); */
  /*   sharpen_reds(); */
  /*   apply_chroma_curve(); */
  f20_color_correction(x3f);
  /*   colorspace_transformation(); */
  /*   black_border_suppression(); */

  return 0;
}

/* color correction for f20 captors */
void f20_color_correction(X3F *x3f){
  IMA *ima=(IMA *)x3f->raw->datas;
  uint16_t (*image)[3]=(uint16_t (*)[3])ima->imageData;
  int row, col, c, height, width, i, d;
  uint16_t *pix;
  char wbcc[WBCC_MAX_LENGTH], wbgains[WBGAINS_MAX_LENGTH], cmcc[64], cmcm[64];
  char *prop_value;
  CAMF_LIST_ENTRY *camf_entry;
  CMbM *cmbm;
  float cc[3][3], gain[3];

  width=ima->columns;
  height=ima->rows;
 
  /* Pour SD1/1M & DP2/DPM
     Il faut regarder dans WhiteBalanceColorCorrections (camf) la matrice à récupérer
     et dans WhiteBalanceGains (camf) la table à récupérer
  */
  /* DP1 should be processed just like TRUEII raws, but the camf name
     differs
     WARNING!: SD15 should be processed just as DP1!!
     SD15 also have ColorMode!!
  */

    bool f20=TRUE;
    /* Argh!!! DP1 use TRUE RAW, but camf_name are not the same */
    /* We should use DP1_%sCCMatrix, DP1_%sWBGain, DP1_CP2_Matrix */
    /* Also there are no CMCM, CMCC, CorrectColorGain, FNumberGainFact SensorAdjustmentGainFact */
    /* but we have BaseGain and ISOGainFact */

    printf ("entered X3F_do_color_correction\n");
    if (ima->imageDataType==X3F_DATA_TYPE_RAW) f20=FALSE;

    camf_entry=X3F_get_camf_entry(x3f->camf_list, foveon_get_param(x3f->camf_list, "WhiteBalanceColorCorrections", (char *)x3f->header->whiteBalanceString));
    if (!camf_entry)
    camf_entry=X3F_get_camf_entry(x3f->camf_list, foveon_get_param(x3f->camf_list, "DP1_WhiteBalanceColorCorrections", (char *)x3f->header->whiteBalanceString));

    printf("Got %s camf entry\n", camf_entry->name);

    cmbm=camf_entry->value;
    /* compute color correction matrix */
    /* We know ccmatrix is [3][3] */
    for (c=0,i=0;c<3;c++){
      for(d=0;d<3;d++) {
	cc[c][d]=cmbm->matrix[i++].f;
      }
    }

    /* CP2_Matrix? */
/*     if (!f20) camf_entry=X3F_get_camf_entry(x3f->camf_list, "DP1_CP2_Matrix"); */
/*     else camf_entry=X3F_get_camf_entry(x3f->camf_list, "CP2_Matrix"); */
/*     printf("Got %s camf entry\n", camf_entry->name); */

/*     cmbm=camf_entry->value; */
/*     for (c=0,i=0;c<3;c++){ */
/*       for(d=0;d<3;d++) { */
/* 	cc[c][d]*=cmbm->matrix[i++].f; */
/*       } */
/*     } */

    /* color mode */
    prop_value=foveon_get_property (x3f->property, "CM_DESC");
    if (camf_entry=X3F_get_camf_entry(x3f->camf_list, foveon_get_param(x3f->camf_list, "ColorModeCompensations", prop_value))){
      printf("Got %s camf entry\n", camf_entry->name);

      cmbm=camf_entry->value;
      /* compute color correction matrix */
      /* We know ccmatrix is [3][3] */
      for (c=0,i=0;c<3;c++){
	for(d=0;d<3;d++) {
	  cc[c][d]*=cmbm->matrix[i++].f*.45; /* .45 = TCGamma ? */
	}
      }
    }
    g_free(prop_value);
      /*CorrectColorGain RR*/
    if (camf_entry=X3F_get_camf_entry(x3f->camf_list, "CorrectColorGain_RR")){
      printf("Got %s camf entry\n", camf_entry->name);

      cmbm=camf_entry->value;
      for (c=0;c<3;c++){
	cc[0][c]*=cmbm->matrix[c].f;
      }
      camf_entry=X3F_get_camf_entry(x3f->camf_list, "CorrectColorGain_GR");
      printf("Got %s camf entry\n", camf_entry->name);

      cmbm=camf_entry->value;
      for (c=0;c<3;c++){
	cc[1][c]*=cmbm->matrix[c].f;
      }
      camf_entry=X3F_get_camf_entry(x3f->camf_list, "CorrectColorGain_BR");
      printf("Got %s camf entry\n", camf_entry->name);

      cmbm=camf_entry->value;
      for (c=0;c<3;c++){
	cc[2][c]*=cmbm->matrix[c].f;
      }
    }

    printf("values of cc: %f %f %f %f %f %f %f %f %f\n", cc[0][0], cc[1][0], cc[2][0], cc[1][0], cc[1][1], cc[1][2], cc[2][0], cc[2][1], cc[2][2]);
		
    /* WBGains */
    camf_entry=X3F_get_camf_entry(x3f->camf_list, foveon_get_param(x3f->camf_list, "WhiteBalanceGains", (char *)x3f->header->whiteBalanceString));
    if (!camf_entry) camf_entry=X3F_get_camf_entry(x3f->camf_list, foveon_get_param(x3f->camf_list, "DP1_WhiteBalanceGains", (char *)x3f->header->whiteBalanceString));
    printf("Got %s camf entry\n", camf_entry->name);
    cmbm=camf_entry->value;
    /* compute gain matrix */
    /* We know gainMatrix is [3] */
    for (i=0; i<3; i++)
      gain[i]=cmbm->matrix[i].f;
    camf_entry=X3F_get_camf_entry(x3f->camf_list, "TempGainFact");
    printf("Got %s camf entry\n", camf_entry->name);
    cmbm=camf_entry->value;
    /* compute gain matrix */
    /* We know gainMatrix is [3] */
    for (i=0; i<3; i++)
      gain[i]*=cmbm->matrix[i].f;

    if (camf_entry=X3F_get_camf_entry(x3f->camf_list, "FNumberGainFact")){
      printf("Got %s camf entry\n", camf_entry->name);
      cmbm=camf_entry->value;
      /* compute gain matrix */
      /* We know gainMatrix is [3] */
      for (i=0; i<3; i++)
	gain[i]*=cmbm->matrix[i].f;
    }
    if (camf_entry=X3F_get_camf_entry(x3f->camf_list, "SensorAdjustmentGainFact")){
      printf("Got %s camf entry\n", camf_entry->name);
      cmbm=camf_entry->value;
      /* compute gain matrix */
      /* We know gainMatrix is [3] */
      for (i=0; i<3; i++)
	gain[i]*=cmbm->matrix[i].f;
    }
/*     } else { */
      //	camf_entry=X3F_get_camf_entry(x3f->camf_list, "BaseGain");
      //printf("Got %s camf entry\n", camf_entry->name);
      //  cmbm=camf_entry->value;
      /* compute gain matrix */
      /* We know gainMatrix is [3] */
      //  for (i=0; i<3; i++)
      //			gain[i]*=cmbm->matrix[i].f;
      //		camf_entry=X3F_get_camf_entry(x3f->camf_list, "ISOGainFact");
      //printf("Got %s camf entry\n", camf_entry->name);
      //  cmbm=camf_entry->value;
      /* compute gain matrix */
      /* We know gainMatrix is [3] */
      //  for (i=0; i<3; i++)
      //			gain[i]/=cmbm->matrix[i].f;

/*     } */
    /* ColorMode ColorCorrection*/
    /* We need to get CM_DESC property */

    //camf_entry=X3F_get_camf_entry(x3f->camf_list, "CMCC_landscape");
    //printf("Got %s camf entry\n", camf_entry->name);
    //  cmbm=camf_entry->value;
    /* compute gain matrix */
    /* We know gainMatrix is [3] */
    //  for (i=0; i<3; i++)
    //				gain[i]=cmbm->matrix[i].f;


    for (i=0;i<3;i++)
      for (c=0;c<3;c++) cc[i][c]*=gain[c];
    
    printf("values of cc: %f %f %f %f %f %f %f %f %f\n", cc[0][0], cc[1][0], cc[2][0], cc[1][0], cc[1][1], cc[1][2], cc[2][0], cc[2][1], cc[2][2]);



  /* On veut modifier les valeurs RVB de chaque pixel */

  for (row=0; row<height; row++){
    pix = image[row*width];
    for (col=0; col<width; col++){
      uint16_t val[3];
      for (c=0; c < 3; c++) val[c]=pix[c];
      for (c=0; c < 3; c++) {
	float tmp;
	tmp=cc[c][0]*val[0]+cc[c][1]*val[1]+cc[c][2]*val[2];
	if (tmp<0) pix[c]=0; else pix[c]=floor(tmp/* *4095/ima->max[c] */);
      }
      pix+=3;
    }
  }
}

void color_correction(X3F *x3f){
  static const float table[][12] =
    {{1.4032,-0.2231,-0.1016,-0.5263,1.4816,0.017,-0.0112,0.0183,0.9113 },
    {3.4032,-0.2231,-0.1016,-0.5263,3.4816,0.017, -0.0112,0.0183,0.4113 }};

  IMA *ima=(IMA *)x3f->raw->datas;
  uint16_t (*image)[3]=(uint16_t (*)[3])ima->imageData;
  int row, col, c, height, width, i, j, d;
  uint16_t *pix;
  char wbcc[WBCC_MAX_LENGTH], wbgains[WBGAINS_MAX_LENGTH], cmcc[64], cmcm[64];
  char *param_value;
  CAMF_LIST_ENTRY *camf_entry;
  CMbM *cmbm;
  float cc[3][3], wbc[3][3], last[3][3], div[3], trans[3][3];
  float  rgb_cam[3][3];
  double  dsum=0, trsum[3];
  uint  keep[4], active[4];


  for (i=0; i < 3; i++)
    for (c=0;c<3;c++) rgb_cam[i][c] = table[0][i*3+c];


  /* Trim off the black border */
/*   camf_entry=X3F_get_camf_entry(x3f->camf_list, "KeepImageArea"); */
/*   cmbm=camf_entry->value; */
/*   get_matrix(cmbm, keep,4); */
/*   printf("keep: %d %d %d %d\n", keep[0], keep[1], keep[2], keep[3]); */

/*   camf_entry=X3F_get_camf_entry(x3f->camf_list, "ActiveImageArea"); */
/*   cmbm=camf_entry->value; */
/*   get_matrix(cmbm, active, 4); */
/*   printf("active: %d %d %d %d\n", active[0], active[1], active[2], active[3]); */

/*   active[1] -= keep[1]; */
/*   active[3] -= 2; */
/*   i = active[2] - active[0]; */
/*   for (row=0; row < active[3]-active[1]; row++) */
/* 	memcpy (image[row*i], image[(row+active[1])*ima->columns+active[0]], */
/* 			i * sizeof *image); */
/*   ima->columns = width = i; */
/*   ima->rows = height = row; */
  width = ima->columns;
  height = ima->rows;

  /* /\* DP1 : WhiteBalanceIlluminants * WhiteBalanceCorrections *\/ */
  /* /\*   float mat[3][3]={{1.325910,0.084080,0.570090}, *\/ */
  /* /\* 		   {-2.328980,5.665190,-1.810120,}, *\/ */
  /* /\* 		   {2.327940,-8.260040,8.626610}}; *\/ */
  /* SD14:  WhiteBalanceIlluminants * WhiteBalanceCorrections */
  /*   float mat[3][3]={{1.281590,0.105849,0.523760}, */
  /* 		   {-2.273677,5.651217,-2.115273}, */
  /* 		   {2.101025,-7.713531,8.711900}}; */
  camf_entry=X3F_get_camf_entry(x3f->camf_list, foveon_get_param(x3f->camf_list, "WhiteBalanceIlluminants", (char *)x3f->header->whiteBalanceString));
  printf("Got %s camf entry\n", camf_entry->name);
  cmbm=camf_entry->value;
  get_matrix(cmbm, cc, 9);

  camf_entry=X3F_get_camf_entry(x3f->camf_list, foveon_get_param (x3f->camf_list, "WhiteBalanceCorrections", (char *)x3f->header->whiteBalanceString));
  if (camf_entry==NULL) {
	/* SD15 hack */
	camf_entry=X3F_get_camf_entry(x3f->camf_list, "CamToXYZ_Flash");
  }
  printf("Got %s camf entry\n", camf_entry->name);
  cmbm=camf_entry->value;
  get_matrix(cmbm, wbc, 9);

  memset(last, 0, sizeof last);
  for(i=0;i<3;i++) {
	for(d=0;d<3;d++) {
	  for (c=0;c<3;c++) last[i][c]+=cc[i][d]*wbc[d][c];
	}
  }

  sprintf(wbcc, "%sRGBNeutral", (char *)x3f->header->whiteBalanceString);
  camf_entry=X3F_get_camf_entry(x3f->camf_list,  wbcc);
  cmbm=camf_entry->value;
  get_matrix(cmbm, div, 3);

  memset (trans, 0, sizeof trans);
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      FORC3 { trans[i][j] += rgb_cam[i][c] * last[c][j]/*  * div[j] */;
	}
/*   FORC3 trsum[c] = trans[c][0] + trans[c][1] + trans[c][2]; */
/*   dsum = (6*trsum[0] + 11*trsum[1] + 3*trsum[2]) / 20; */
/*   for (i=0; i < 3; i++) */
/*     FORC3 last[i][c] = trans[i][c] * dsum / trsum[i]; */
/*   memset (trans, 0, sizeof trans); */
/*   for (i=0; i < 3; i++) */
/*     for (j=0; j < 3; j++){ */
/*       FORC3 trans[i][j] += (i==c ? 32 : -1) * last[c][j] / 30; */
/* 	  printf("trans[%d][%d]: %f\n", i, j, trans[i][j]); */
/* 	} */

  /* On veut modifier les valeurs RVB de chaque pixel */

  for (row=0; row<height; row++){
    pix = image[row*width];
    for (col=0; col<width; col++){
      uint16_t val[3];
      for (c=0; c < 3; c++) val[c]=pix[c];
      for (c=0; c < 3; c++) {
	float tmp;
	tmp=trans[c][0]*val[0]+trans[c][1]*val[1]+trans[c][2]*val[2];
	if (tmp<0) pix[c]=0; else pix[c]=floor(tmp*4095/(ima->max[c]));
      }
      pix+=3;
    }
  }
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
/* 	  x3f->camf=camf; */
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

