/* This simple program is just a way to test the libx3f */
/* We should open a raw x3f file, display the thumbnail, the preview and the raw image */
/* We should also list all properties, and all camf entries */
/* */

#include <gtk/gtk.h>
#include "raw_x3f.h"
/* #include "dcraw_func.h" */
#include <wchar.h>
#include <locale.h>

typedef   struct {
    char *filename;
    X3F *x3f_struct;
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

void add_camf(char *key, uint32_t type, GtkWidget *tree){
  GtkListStore *CAMFdata=  GTK_LIST_STORE(gtk_tree_view_get_model
					  (GTK_TREE_VIEW(tree)));
  GtkTreeIter Iter;
  union {
    char atype[5];
    uint32_t itype;
  } conv;

  memset(&conv, 0, sizeof(conv));
  conv.itype=type;

  gtk_list_store_append(CAMFdata, &Iter);
  gtk_list_store_set(CAMFdata, &Iter,
		     0, key,
		     1, conv.atype,
		     -1);

  return;
}

void update_camf_view(GtkWidget *widget, X3F_FILE *x3f_file){
  GtkWidget *camfTextview=(GtkWidget *)gtk_builder_get_object(p_builder, "CAMFdataTextview" );
  GtkTreeIter iter;
  GtkTreeModel *model;
  GtkTextBuffer *buffer;
  char *value;
  CAMF_LIST_ENTRY *camf_entry;
  uint32_t i, c, v, valueCount;
  GString *string;
  CMbM *cmbm;
  PARAMETERS *params;

  if (gtk_tree_selection_get_selected(GTK_TREE_SELECTION(widget), &model, &iter)) {

    gtk_tree_model_get(model, &iter, 0, &value,  -1);

    /* OK, we need to find the corresponding camf entry */
    for (camf_entry=x3f_file->x3f_struct->camf_list;camf_entry!=NULL;camf_entry=camf_entry->next)
      if (!(strcmp(camf_entry->name, value)))
	break;

    g_free(value);

    switch (camf_entry->CAMFtype) {
    case X3F_CMbT:
      string=g_string_new("Technical Infos\n");
       g_string_append_printf(string, "%s\n",
			      (char*)camf_entry->value);
       break;
    case X3F_CMbP:
      params=camf_entry->value;
      string=g_string_new("Parameters ");
      g_string_append_printf(string, "%s\nParameters count: %d\n",
			     camf_entry->name,
			     camf_entry->count);
      for (i=0; i<camf_entry->count; i++)
	g_string_append_printf(string, "%s\t=>\t%s\n", params[i].name,params[i].value); 
      break;
    case X3F_CMbM:
      cmbm=camf_entry->value;
      string=g_string_new("Matrix ");
      g_string_append_printf(string, "%s\nMatrix dimension: %d\tMatrix data type: %d\n",
			    camf_entry->name,
			    camf_entry->count,
			    cmbm->dataType);
      g_string_append_printf(string, "Individual plane size:");
      for (i=0; i<camf_entry->count;i++)
	g_string_append_printf(string, " %d", cmbm->planeElements[i]);
      g_string_append_printf(string, "\n\n");
      g_string_append_printf(string, "The following values are not properly formatted\n");
      g_string_append_printf(string, "I made the following assumption:\nfirst matrix dimension is the number of rows\n");
      g_string_append_printf(string, "second matrix dimension is the number of columns\n");
      g_string_append_printf(string, "third matrix dimension is the number of planes\n");

      /* FIXME: we here assume the matrix dimensions will never exceed 3, which is right so far */
      /* We will make the assumption first dimension is the number of row, second dim is the number of columns, third is the deep of the array */
      /* this is really ugly code */

      /* compute the number of values to display */
      valueCount=1;
      for (i=0; i<camf_entry->count; i++)
	valueCount*=cmbm->planeElements[i];

      for (i=0; i<valueCount; i++){
	switch (camf_entry->count){
	case 1:
	  g_string_append_printf(string, "\n");
	  break;
	case 2:
	  if (i%cmbm->planeElements[1]==0) g_string_append_printf(string, "\n");
	  break;
	case 3:
	  if (i%cmbm->planeElements[1]==0) g_string_append_printf(string, "\n");
	  if (i%(cmbm->planeElements[0]*cmbm->planeElements[1])==0) g_string_append_printf(string, "\n");
	}
	switch (cmbm->dataType) {
	case 0:
	  g_string_append_printf(string, "%d ", cmbm->matrix[i].ui_16);
	  break;
	case 3:
	  g_string_append_printf(string, "%f ", cmbm->matrix[i].f);
	  break;
	case 6:
	  g_string_append_printf(string, "%d ", cmbm->matrix[i].ui_8);
	  break;
	default:
	  g_string_append_printf(string, "%d ", cmbm->matrix[i].ui_32);
	}
      }
     break;
    }
    value=g_string_free(string, FALSE);
 
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
  GdkPixbuf *thumbImg, *tmp;
  GdkPixbuf *previewImg, *rawImg;
  IMA *ima;
  PROPERTY *prop;
  CAMF_LIST_ENTRY *camf_entry;
  X3F_FILE *x3f_file=(X3F_FILE *)userdata;
  int count;
  uint8_t *img24=NULL;
  int i,c,value;
  GString *string;

	  
#define x3f x3f_file->x3f_struct

  x3f=X3F_load_full_x3f(x3f_file->filename);
  string=g_string_new("Header Infos : \n");
  g_string_append_printf(string, "Size: %d x %d\nWhite Balance: %s \n", x3f->header->columns, x3f->header->rows, (char *)x3f->header->whiteBalanceString);
  for (i=0;i<4;i++)
    g_string_append_printf(string, " flag[%d]: %d\n", i, x3f->header->flags[i]);
  gtk_label_set_text(GTK_LABEL(headerInfos), g_string_free(string, FALSE));

  X3F_decode_raw(x3f->raw->datas);

 
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
  for (camf_entry=x3f->camf_list;camf_entry!=NULL;camf_entry=camf_entry->next){
    add_camf(camf_entry->name, camf_entry->CAMFtype, CAMFdataTreeView);
  }
  g_signal_connect(selection, "changed", 
      G_CALLBACK(update_camf_view), x3f_file);

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

    /* Scale down pixbuf */
/*     gdouble ratio=(gdouble) gdk_pixbuf_get_height(tmp)/700; */
/*     previewImg=gdk_pixbuf_scale_simple(tmp, gdk_pixbuf_get_width(tmp)/ratio, 700, GDK_INTERP_HYPER); */

/*     g_object_unref(tmp); */
    preview=gtk_image_new_from_pixbuf(/* previewImg */ tmp);
    holder= (GtkWidget *) gtk_builder_get_object (p_builder, "previewWindow" );
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(holder), preview);
    gtk_widget_show(preview);
  }

  /* Display RAW */
  ima=(IMA *)x3f->raw->datas;

  if ((ima->flags & DECODED_IMAGE)==DECODED_IMAGE) {
    uint16_t *rawData=(uint16_t *)ima->imageData;

    count=ima->rows*ima->columns;
    img24=malloc(sizeof(*img24)*count*3);
    for (i=0;i<count;i++){
      for (c=0;c<3;c++){
	if ((rawData[i*3+c]>>3)>255)
	  img24[i*3+c]=255;
	else
	  img24[i*3+c]=rawData[i*3+c]>>3;
      }
    }
    tmp=gdk_pixbuf_new_from_data(img24, GDK_COLORSPACE_RGB, FALSE, 8,
				 ima->columns, ima->rows,
				 ima->rowSize, free_img24, NULL);

    if (x3f->header->rotation != 0) {
      rawImg= gdk_pixbuf_rotate_simple(tmp,360-x3f->header->rotation);
      g_object_unref(tmp);
    } else {
      rawImg=tmp;
    }
    raw=gtk_image_new_from_pixbuf(rawImg);
    holder= (GtkWidget *) gtk_builder_get_object (p_builder, "rawWindow" );
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(holder), raw);
    gtk_widget_show(raw);

  }

  /* Display RAW red channel */

  if ((ima->flags & DECODED_IMAGE)==DECODED_IMAGE) {
    uint16_t *rawData=(uint16_t *)ima->imageData;

    count=ima->rows*ima->columns*3;
    img24=malloc(sizeof(*img24)*count*3);
    for (i=0;i<count;i+=3){
      if ((rawData[i]>>3)>255)
	value=255;
      else
	value=rawData[i]>>3;
      img24[i]=img24[i+1]=img24[i+2]=value;
    }
    tmp=gdk_pixbuf_new_from_data(img24, GDK_COLORSPACE_RGB, FALSE, 8,
				 ima->columns, ima->rows,
				 ima->rowSize, free_img24, NULL);

    if (x3f->header->rotation != 0) {
      rawImg= gdk_pixbuf_rotate_simple(tmp,360-x3f->header->rotation);
      g_object_unref(tmp);
    } else {
      rawImg=tmp;
    }
    raw=gtk_image_new_from_pixbuf(rawImg);
    holder= (GtkWidget *) gtk_builder_get_object (p_builder, "redWindow" );
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(holder), raw);
    gtk_widget_show(raw);
  }

  /* Display RAW green channel */
  if ((ima->flags & DECODED_IMAGE)==DECODED_IMAGE) {
   uint16_t *rawData=(uint16_t *)ima->imageData;

    count=ima->rows*ima->columns*3;
    img24=malloc(sizeof(*img24)*count*3);
    for (i=0;i<count;i+=3){
      if ((rawData[i+1]>>3)>255)
	value=255;
      else
	value=rawData[i+1]>>3;
      img24[i]=img24[i+1]=img24[i+2]=value;
    }
    tmp=gdk_pixbuf_new_from_data(img24, GDK_COLORSPACE_RGB, FALSE, 8,
				 ima->columns, ima->rows,
				 ima->rowSize, free_img24, NULL);

    if (x3f->header->rotation != 0) {
      rawImg= gdk_pixbuf_rotate_simple(tmp,360-x3f->header->rotation);
      g_object_unref(tmp);
    } else {
      rawImg=tmp;
    }
    raw=gtk_image_new_from_pixbuf(rawImg);
    holder= (GtkWidget *) gtk_builder_get_object (p_builder, "greenWindow" );
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(holder), raw);
    gtk_widget_show(raw);
  }
  /* Display RAW blue channel */
  if ((ima->flags & DECODED_IMAGE)==DECODED_IMAGE) {
    uint16_t *rawData=(uint16_t *)ima->imageData;

    count=ima->rows*ima->columns*3;
    img24=malloc(sizeof(*img24)*count*3);
    for (i=0;i<count;i+=3){
      if ((rawData[i+2]>>3)>255)
	value=255;
      else
	value=rawData[i+2]>>3;
      img24[i]=img24[i+1]=img24[i+2]=value;
    }
    tmp=gdk_pixbuf_new_from_data(img24, GDK_COLORSPACE_RGB, FALSE, 8,
				 ima->columns, ima->rows,
				 ima->rowSize, free_img24, NULL);

    if (x3f->header->rotation != 0) {
      rawImg= gdk_pixbuf_rotate_simple(tmp,360-x3f->header->rotation);
      g_object_unref(tmp);
    } else {
      rawImg=tmp;
    }
    raw=gtk_image_new_from_pixbuf(rawImg);
    holder= (GtkWidget *) gtk_builder_get_object (p_builder, "blueWindow" );
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(holder), raw);
    gtk_widget_show(raw);
  }


  /* Signals */
/*   g_signal_connect (gtk_builder_get_object (p_builder, "button1"), */
/* 		    "clicked", G_CALLBACK (main_quit), */
/* 		    NULL); */
/*   g_signal_connect (p_win, "delete_event", G_CALLBACK (main_quit), */
/* 		    NULL);	    */
  gtk_builder_connect_signals(p_builder, NULL);
#undef x3f	   
  gtk_widget_show_all (p_win);
  gtk_widget_destroy(splash);
  return FALSE;
}

int main(int argc, char *argv[]){
  X3F_FILE *x3f_file=malloc(sizeof(*x3f_file));

  if (argc!=2){
    fprintf(stderr, "usage: %s <X3F File>\n", argv[0]);
    return 1;
  }
  setlocale(LC_CTYPE, "");

  x3f_file->filename=argv[1];


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
  free(x3f_file);
  return 0;
   
}
