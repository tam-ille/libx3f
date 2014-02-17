/* This simple program is just a way to test the libx3f */
/* We should open a raw x3f file, display the thumbnail, the preview and the raw image */
/* We should also list all properties, and all camf entries */
/* */

#include <gtk/gtk.h>
#include "raw_x3f.h"
#include "interpolation.h"

/* #include "dcraw_func.h" */
#include <wchar.h>
#include <locale.h>
#include <math.h>
#include <unistd.h>
#include <arpa/inet.h>

#define FORC(cnt) for (c=0; c < cnt; c++)
#define FORC3 FORC(3)

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#define SQR(x) ((x)*(x))
#define LIM(x,min,max) MAX(min,MIN(x,max))
#define ULIM(x,y,z) ((y) < (z) ? LIM(x,y,z) : LIM(x,z,y))
#define CLIP(x) LIM(x,0,65535)
#define SWAP(a,b) { a=a+b; b=a-b; a=a-b; }


typedef   struct {
  char *filename;
  X3F *x3f_struct;
  INTERPOLATED_IMG *interpolated_img;
} X3F_FILE;

GtkBuilder  *  p_builder   = NULL;
GError      *  p_err       = NULL;


void add_property(utf16_t *key, utf16_t *val, GtkWidget *tree){
  GtkListStore *Metadata=  GTK_LIST_STORE(gtk_tree_view_get_model
										  (GTK_TREE_VIEW(tree)));
  GtkTreeIter Iter;

  gtk_list_store_append(Metadata, &Iter);
  gtk_list_store_set(Metadata, &Iter,
					 0, g_utf16_to_utf8(key, -1,NULL,NULL,NULL),
					 1, g_utf16_to_utf8(val, -1,NULL,NULL,NULL),
					 -1);

  return;
}

CMb_HEADER *get_camf_header (CAMF *camf, const char *name) {
  unsigned idx;
  char *pos;

  for (idx=0; idx < camf->dataSize; idx += sget4(pos+8)) {
    pos = camf->camf_data + idx;
    if (strncmp (pos, "CMb", 3)) break;
    if (strcmp (name, pos+sget4(pos+12))) continue;
	return (CMb_HEADER *)pos;
  }
/*   printf("camf entry %s not found!!\n", name); */
  return NULL;

}

void add_camf(X3F *x3f, GtkWidget *tree){
  GtkListStore *CAMFdata=  GTK_LIST_STORE(gtk_tree_view_get_model
										  (GTK_TREE_VIEW(tree)));
  GtkTreeIter Iter;
  union {
    char atype[5];
    uint32_t itype;
  } conv;


  /*  FOR TRUEII camera, lots of camf blocks are missing, so IncludeBlocks list is not accurate */
  memset(&conv, 0, sizeof(conv));
  unsigned idx, num;
  int i=0,j;
  char *pos, *cp, *dp, *entry_name;
  char *value; /* should be part of x3f */
  
  CAMF *camf;
  DIR_ENTRY *section;
  CMb_HEADER *header;

  section=X3F_get_section(x3f, X3F_CAMF);
  camf=section->datas;

  for (idx=0; idx < camf->dataSize; idx += sget4(pos+8)) {
    pos = camf->camf_data + idx;
    if (strncmp (pos, "CMb", 3)) break;
    if (pos[3] != 'P') continue;
    if (strcmp ("IncludeBlocks", pos+sget4(pos+12))) continue;
	/* IncludeBlocks found */
    cp = pos + sget4(pos+16);
    j = num = sget4(cp); /* number of entries */
    dp = pos + sget4(cp+4);

    while (j--) {
      cp += 8;
	  /* find the type of the entry */
	  if ((header=get_camf_header(camf, dp+sget4(cp))))
		conv.itype=header->CAMFsubsectionID;
	  else 
		conv.itype=X3F_CMb;

	  gtk_list_store_append(CAMFdata, &Iter);
	  gtk_list_store_set(CAMFdata, &Iter,
						 0, dp+sget4(cp),
						 1, conv.atype,
						 -1);

    }
  }

  return;
}

void update_camf_view(GtkWidget *widget, X3F *x3f){
  GtkWidget *camfTextview=(GtkWidget *)gtk_builder_get_object(p_builder, "CAMFdataTextview" );
  GtkTreeIter iter;
  GtkTreeModel *model;
  GtkTextBuffer *buffer;
  char *value, *pos, *dp, *cp;
  CAMF_LIST_ENTRY *camf_entry;
  uint32_t i, c, v, valueCount, count, type;
  GString *string;
/*   CMbM *cmbm; */
/*   PARAMETERS *params; */
  int32_t *matrix;
  float *fmatrix;
  CMb_HEADER *header;
  CAMF *camf;
  DIR_ENTRY *section;
  uint dim[3];

  section=X3F_get_section(x3f, X3F_CAMF);
  camf=section->datas;

  if (gtk_tree_selection_get_selected(GTK_TREE_SELECTION(widget), &model, &iter)) {

    gtk_tree_model_get(model, &iter, 0, &value,  -1);

    /* OK, we need to find the corresponding camf entry */
	/* obtain the Header to get the CAMFsubsectionID (CMbT, CMbP, CMbM) */
	/* place the pos to datastart */
	/* act accordingly depending of CAMFsubsectionID */
	if (!(header= get_camf_header(camf, value))) return;
	pos=(char *)header+header->dataOffset;

    switch (header->CAMFsubsectionID) {
    case X3F_CMbT:
      string=g_string_new("Technical Infos\n");
	  g_string_append_printf(string, "%s\n",
							 pos+4);
	  break;
    case X3F_CMbP:
	  /* CMbP: pos=> number of params to be found */
 	  cp=pos;
	  dp=(char*)header+sget4(cp+4);
	  count=sget4(pos);
	  string=g_string_new("Parameters ");
      g_string_append_printf(string, "%s\nParameters count: %d\n",
							 value,
							 count);
	  while (count--) {
		cp+=8;
		g_string_append_printf(string, "%s\t=>\t%s\n", dp+sget4(cp),dp+sget4(cp+4));
	  }
      break;
    case X3F_CMbM:
     string=g_string_new("Matrix ");
	  count=sget4(pos+4);
	  type=sget4(pos);
      g_string_append_printf(string, "%s\nMatrix dimension: %d\tMatrix data type: %d\n", value, count, type);
	  if (type && type != 3)
		matrix=X3F_foveon_camf_matrix(camf, dim, value);
	  else
		fmatrix=X3F_foveon_camf_matrix(camf, dim, value);

      g_string_append_printf(string, "Individual plane size:");
      for (i=0; i<count;i++)
		g_string_append_printf(string, " %d", dim[i]);
      g_string_append_printf(string, "\n\n");
	  /*       g_string_append_printf(string, "The following values are not properly formatted\n"); */
	  /*       g_string_append_printf(string, "I made the following assumption:\nfirst matrix dimension is the number of rows\n"); */
	  /*       g_string_append_printf(string, "second matrix dimension is the number of columns\n"); */
	  /*       g_string_append_printf(string, "third matrix dimension is the number of planes\n"); */

	  /*       /\* FIXME: we here assume the matrix dimensions will never exceed 3, which is right so far *\/ */
	  /*       /\* We will make the assumption first dimension is the number of row, second dim is the number of columns, third is the deep of the array *\/ */
	  /*       /\* this is really ugly code *\/ */

	  /*       /\* compute the number of values to display *\/ */
      valueCount=1;
      for (i=0; i<count; i++)
		valueCount*=dim[i];

      for (i=0; i<valueCount; i++){
		switch (count){
		case 1:
		  g_string_append_printf(string, "\n");
		  break;
		case 2:
		  if (i%dim[1]==0) g_string_append_printf(string, "\n");
		  break;
		case 3:
		  if (i%dim[1]==0) g_string_append_printf(string, "\n");
		  if (i%(dim[0]*dim[1])==0) g_string_append_printf(string, "\n");
		}
		if (type && type!=3)
		  g_string_append_printf(string, "%d ", matrix[i]);
		else
		  g_string_append_printf(string, "%f ", fmatrix[i]);
      }
	  if (type && type!=3){
		if (matrix)
		  free(matrix);
	  } else {
		if (fmatrix) free(fmatrix);
	  }
	  break;
	default:
	  printf("Unknnown CAMF subsection!²\n");
    }
    value=g_string_free(string, FALSE);
	gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(camfTextview), GTK_WRAP_WORD);
	buffer = gtk_text_view_get_buffer (GTK_TEXT_VIEW (camfTextview));
    gtk_text_buffer_set_text (buffer, value, -1);
  }

  g_free(value);

}


int sort_metadata(GtkTreeModel *model, GtkTreeIter *a, GtkTreeIter *b, gpointer userdata){
  gchar *str_a, *str_b;
  gint i;

  gtk_tree_model_get(model, a, 0, &str_a, -1);
  gtk_tree_model_get(model, b, 0, &str_b, -1);
  i=strcmp(str_a, str_b);
  g_free(str_a);
  g_free(str_b);
  return i;
}

void main_quit() {
  gtk_main_quit();
}

void free_img24(guchar *pixels, gpointer userdata){
  free(pixels);
}

gboolean load_file(gpointer userdata) {
  GtkWidget *splash = (GtkWidget *) gtk_builder_get_object (p_builder, "Splash" );
  GtkWidget * p_win = (GtkWidget *) gtk_builder_get_object (p_builder, "mainWindow" );
  GtkWidget *holder = (GtkWidget *) gtk_builder_get_object (p_builder, "leftbox" );
  GtkWidget *headerInfos = (GtkWidget *) gtk_builder_get_object (p_builder, "headerInfos" );
  GtkWidget *metadataTreeView = (GtkWidget *) gtk_builder_get_object (p_builder, "MetaDataTreeView" );
  GtkWidget *CAMFdataTreeView = (GtkWidget *) gtk_builder_get_object (p_builder, "CAMFdataTreeview" );
  GtkListStore *Metadata;
  GtkListStore *CAMFdata;
  GtkCellRenderer *renderer;
  GtkTreeViewColumn *column;
  GtkTreeSelection *selection= (GtkTreeSelection *) gtk_builder_get_object (p_builder, "CAMFtreeviewSelection" );

  GtkWidget *thumbnail, *preview, *raw;
  GdkPixbuf *thumbImg, *tmp, *rawImg_small;
  GdkPixbuf *previewImg, *rawImg;
  IMA *ima;
  PROPERTY *prop;
  /*   CAMF_LIST_ENTRY *camf_entry; */
  X3F_FILE *x3f_file=(X3F_FILE *)userdata;
  INTERPOLATED_IMG *interpolated=x3f_file->interpolated_img;
  uint count, row, col;
  uint8_t *img24=NULL;
  uint i,c,value, num;
  GString *string;
  uint16_t (*raw_img)[4];
  uint16_t (*i_img)[4]=NULL;

#define x3f x3f_file->x3f_struct

  x3f=X3F_load_full_x3f(x3f_file->filename);
  string=g_string_new("Header Infos : \n");
  g_string_append_printf(string, "Size: %d x %d\nWhite Balance: %s \n", x3f->header->columns, x3f->header->rows, (char *)x3f->header->whiteBalanceString);
  for (i=0;i<4;i++)
    g_string_append_printf(string, " flag[%d]: %d\n", i, x3f->header->flags[i]);
  gtk_label_set_text(GTK_LABEL(headerInfos), g_string_free(string, FALSE));

  X3F_decode_raw(x3f->raw->datas);

  ima=(IMA *)x3f->raw->datas;
  i_img= (uint16_t (*)[4]) calloc (ima->rows*ima->columns*4, sizeof (uint16_t));
  memcpy(i_img, ima->imageData, ima->rows*ima->columns*4*sizeof (uint16_t));
	  
  interpolated->width=ima->columns;
  interpolated->height=ima->rows;
  interpolated->bps=16;
  interpolated->colors=4;
  memset(interpolated->gamm, 0, sizeof interpolated->gamm);
  interpolated->gamm[0]=0.45;
  interpolated->gamm[1]=4.5;
  interpolated->img=i_img;
  for (i=0; i < 0x10000; i++) interpolated->curve[i] = i;
  if (ima->dataFormat == X3F_DATA_FORMAT_TRUE_RAW) {
	x3f_simple_coeff(interpolated, 1);
	x3f_trueII_interpolate(interpolated, x3f);
  } else {
	x3f_simple_coeff(interpolated, 0);
	x3f_interpolate(interpolated, x3f);
  }
  x3f_convert_to_rgb(interpolated, 1);
  x3f_apply_gamma(interpolated, 1.0);
  x3f_output_ppm(interpolated, 16, 0, x3f_file->filename);

  renderer = gtk_cell_renderer_text_new();
  column = gtk_tree_view_column_new_with_attributes("Property",
													renderer, "text", 0, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(metadataTreeView), column);
  column = gtk_tree_view_column_new_with_attributes("Value",
													renderer, "text", 1, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(metadataTreeView), column);

  Metadata=gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_STRING);
  gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(Metadata), 0, GTK_SORT_ASCENDING);
  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(Metadata), 0, sort_metadata, NULL, NULL);
  gtk_tree_view_set_model(GTK_TREE_VIEW(metadataTreeView),
						  GTK_TREE_MODEL(Metadata));
  g_object_unref(Metadata);

  for (prop=x3f->property;prop!=NULL;prop=prop->next){
    add_property(prop->name, prop->value, metadataTreeView);
  }

  column=gtk_tree_view_column_new_with_attributes("CAMF Block name", renderer, "text", 0, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(CAMFdataTreeView), column);
  column=gtk_tree_view_column_new_with_attributes("Type", renderer, "text", 1, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(CAMFdataTreeView), column);
  CAMFdata=gtk_list_store_new(2, G_TYPE_STRING,G_TYPE_STRING);
  gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(CAMFdata), 1, GTK_SORT_DESCENDING);
  gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(CAMFdata), 0, sort_metadata, NULL, NULL);
  gtk_tree_view_set_model(GTK_TREE_VIEW(CAMFdataTreeView),
						  GTK_TREE_MODEL(CAMFdata));
  g_object_unref(CAMFdata);
  /*   for (i=0; i<num; i++){ */
  add_camf(x3f, CAMFdataTreeView);
  /*   } */
  g_signal_connect(selection, "changed",
				   G_CALLBACK(update_camf_view), x3f);

 
  /* Display THUMBNAIL */
  ima=(IMA *)x3f->thumbnail->datas;
  /* create a pixbuf only if imageData has been decoded */
  if ((ima->flags & DECODED_IMAGE)==DECODED_IMAGE) {
    tmp=gdk_pixbuf_new_from_data(ima->imageData, GDK_COLORSPACE_RGB, FALSE, 8,
								 ima->columns, ima->rows,
								 ima->rowSize, NULL, NULL);
    if (x3f->header->rotation != 0) {
      thumbImg = gdk_pixbuf_rotate_simple(tmp,360-x3f->header->rotation);
    } else {
      thumbImg=tmp;
    }
    thumbnail=gtk_image_new_from_pixbuf(thumbImg);
    gtk_box_pack_start(GTK_BOX(holder), thumbnail, FALSE, TRUE, 0);
    gtk_box_reorder_child(GTK_BOX(holder), thumbnail, 1);
    gtk_widget_show(thumbnail);
    g_object_unref(tmp);
  }

  /* Display PREVIEW */
  ima=(IMA *)x3f->preview->datas;
  if ((ima->flags & DECODED_IMAGE)==DECODED_IMAGE) {
    previewImg=gdk_pixbuf_new_from_data(ima->imageData, GDK_COLORSPACE_RGB, FALSE, 8,
										ima->columns, ima->rows,
										ima->rowSize, NULL, NULL);

    if (x3f->header->rotation != 0) {
      tmp= gdk_pixbuf_rotate_simple(previewImg,360-x3f->header->rotation);
    } else {
      tmp=previewImg;
    }

    preview=gtk_image_new_from_pixbuf( tmp);
    holder= (GtkWidget *) gtk_builder_get_object (p_builder, "previewWindow" );
    gtk_container_add(GTK_CONTAINER(holder), preview);
    gtk_widget_show(preview);
  }

  /* Display RAW */
  /* Attention à la taille des pixels: sur DP1M, en résolution moyenne, le capteur garde 4927x1631 pixels
	 soit la largeur haute def et la hauteur basse def or l'image finale devrait faire 3264x2176*/
  ima=(IMA *)x3f->raw->datas;

  count=interpolated->height*interpolated->width;
  uint16_t (*temp)[4];
  temp=(uint16_t  (*)[4])interpolated->img;

  if ((ima->flags & DECODED_IMAGE)==DECODED_IMAGE) {

    img24=malloc(sizeof(*img24)*count*3);

    for (i=0;i<count;i++){
      for (c=0;c<3;c++){
		if (((temp[i][c] >> 8)>255))
		  img24[i*3+c]=255;
		else
		  img24[i*3+c]=temp[i][c] >>8;
	  }
    }
    tmp=gdk_pixbuf_new_from_data(img24, GDK_COLORSPACE_RGB, FALSE, 8,
								 interpolated->width, interpolated->height,
								 interpolated->width*3, free_img24, NULL);



	if (x3f->header->rotation != 0) {
	  rawImg= gdk_pixbuf_rotate_simple(tmp,360-x3f->header->rotation);
	  g_object_unref(tmp);
	} else {
	  rawImg=tmp;
	}
    
	rawImg_small=gdk_pixbuf_scale_simple(rawImg, gdk_pixbuf_get_width(rawImg)/4, gdk_pixbuf_get_height(rawImg)/4, GDK_INTERP_HYPER);

    raw=gtk_image_new_from_pixbuf(rawImg_small);
    holder= (GtkWidget *) gtk_builder_get_object (p_builder, "rawWindow" );
    gtk_container_add(GTK_CONTAINER(holder), raw);
    gtk_widget_show(raw);

  }

  /* Display RAW red channel */

  if ((ima->flags & DECODED_IMAGE)==DECODED_IMAGE) {
    img24=malloc(sizeof(*img24)*count*3);
    for (i=0;i<count;i++){
	  if ((temp[i][0]>>8)>255)
		value=255;
	  else
		value=temp[i][0]>>8;
	  img24[i*3]=img24[i*3+1]=img24[i*3+2]=value;
	}
    tmp=gdk_pixbuf_new_from_data(img24, GDK_COLORSPACE_RGB, FALSE, 8,
								 interpolated->width, interpolated->height,
								 interpolated->width*3, free_img24, NULL);

    if (x3f->header->rotation != 0) {
      rawImg= gdk_pixbuf_rotate_simple(tmp,360-x3f->header->rotation);
      g_object_unref(tmp);
    } else {
      rawImg=tmp;
    }
    raw=gtk_image_new_from_pixbuf(rawImg);
    holder= (GtkWidget *) gtk_builder_get_object (p_builder, "redWindow" );
    gtk_container_add(GTK_CONTAINER(holder), raw);
    gtk_widget_show(raw);
  }

  /* Display RAW green channel */
  if ((ima->flags & DECODED_IMAGE)==DECODED_IMAGE) {
    img24=malloc(sizeof(*img24)*count*3);
    for (i=0;i<count;i++){
	  if ((temp[i][1]>>8)>255)
		value=255;
	  else
		value=temp[i][1]>>8;
	  img24[i*3+c]=img24[i*3+c+1]=img24[i*3+c+2]=value;
    }
    tmp=gdk_pixbuf_new_from_data(img24, GDK_COLORSPACE_RGB, FALSE, 8,
								 interpolated->width, interpolated->height,
								 interpolated->width*3, free_img24, NULL);

    if (x3f->header->rotation != 0) {
      rawImg= gdk_pixbuf_rotate_simple(tmp,360-x3f->header->rotation);
      g_object_unref(tmp);
    } else {
      rawImg=tmp;
    }

    raw=gtk_image_new_from_pixbuf(rawImg);
    holder= (GtkWidget *) gtk_builder_get_object (p_builder, "greenWindow" );
    gtk_container_add(GTK_CONTAINER(holder), raw);
    gtk_widget_show(raw);
  }
  /* Display RAW blue channel */
  if ((ima->flags & DECODED_IMAGE)==DECODED_IMAGE) {
	img24=malloc(sizeof(*img24)*count*3);
	for (i=0;i<count;i++){
	  if ((temp[i][2]>>8)>255)
		value=255;
	  else
		value=temp[i][2]>>8;
	  img24[i*3+c]=img24[i*3+c+1]=img24[i*3+c+2]=value;
	}
    tmp=gdk_pixbuf_new_from_data(img24, GDK_COLORSPACE_RGB, FALSE, 8,
								 interpolated->width, interpolated->height,
								 interpolated->width*3, free_img24, NULL);

    if (x3f->header->rotation != 0) {
      rawImg= gdk_pixbuf_rotate_simple(tmp,360-x3f->header->rotation);
      g_object_unref(tmp);
    } else {
      rawImg=tmp;
    }
    raw=gtk_image_new_from_pixbuf(rawImg);
    holder= (GtkWidget *) gtk_builder_get_object (p_builder, "blueWindow" );
    gtk_container_add(GTK_CONTAINER(holder), raw);
    gtk_widget_show(raw);
  }


  /*   Signals */
  g_signal_connect (gtk_builder_get_object (p_builder, "button1"),
					"clicked", G_CALLBACK (main_quit),
					NULL);
  g_signal_connect (p_win, "delete_event", G_CALLBACK (main_quit),
					NULL);
  gtk_builder_connect_signals(p_builder, NULL);
#undef x3f
  gtk_widget_show_all (p_win);
  gtk_widget_destroy(splash);
  return FALSE;
}

int main(int argc, char *argv[]){
  X3F_FILE *x3f_file=malloc(sizeof(*x3f_file));
  INTERPOLATED_IMG *interpolated;
  if (argc!=2){
    fprintf(stderr, "usage: %s <X3F File>\n", argv[0]);
    return 1;
  }
  setlocale(LC_CTYPE, "");

  x3f_file->filename=argv[1];
  interpolated=malloc(sizeof(*interpolated));
  x3f_file->interpolated_img=interpolated;

  gtk_init (& argc, & argv);


  p_builder = gtk_builder_new ();

  if (p_builder != NULL)
    {
      gtk_builder_add_from_file (p_builder, "x3f_show.xml", & p_err);
       
      if (p_err == NULL)
		{
		  GtkWidget *splash=(GtkWidget *) gtk_builder_get_object (p_builder, "Splash");

		  g_idle_add(load_file, x3f_file);
		  gtk_widget_show_all(splash);
	  
		  gtk_main ();
		}
      else
		{
		  g_error ("%s", p_err->message);
		  g_error_free (p_err);
		}
    }
  X3F_free(x3f_file->x3f_struct);

  free(interpolated->img);
  free(interpolated);
  free(x3f_file);
  return 0;
   
}
