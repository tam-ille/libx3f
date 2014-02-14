
CC=gcc

CFLAGS=-W  -Wl,--export-dynamic\
	`pkg-config gtk+-3.0 gmodule-2.0 --cflags` \
	-I/usr/lib/glib-2.0/include -I/usr/include/glib-2.0

LDFLAGS=-L./ -lx3f -lm -ljpeg `pkg-config gtk+-3.0 gmodule-2.0 --libs`

libx3f: raw_x3f.o interpolation.o
	$(CC) -shared -Wl,-soname,libx3f.so.0 -o libx3f.so.0.0.1 -lc -lm -lglib-2.0 -fPIC raw_x3f.o interpolation.o
# 	ln -s libx3f.so.0.0.1 libx3f.so.0
# 	ln -s libx3f.so.0 libx3f.so

raw_x3f.o: raw_x3f.c raw_x3f.h
	$(CC) -o raw_x3f.o -c raw_x3f.c -fPIC $(CFLAGS)
interpolation.o: interpolation.c interpolation.h
	$(CC) -o interpolation.o -c interpolation.c -fPIC $(CFLAGS)

x3f_show: x3f_show.o 
	$(CC) -o x3f_show x3f_show.o  $(LDFLAGS)

x3f_show.o: x3f_show.c
	$(CC) -o x3f_show.o -c x3f_show.c $(CFLAGS)

#dcraw_func.o: dcraw_func.c dcraw_func.h
#	$(CC) -o dcraw_func.o -c dcraw_func.c $(CFLAGS)

clean:
	rm -rf *.o

mrproper: clean
	rm -rf libx3f.so*
	rm -rf x3f_show
