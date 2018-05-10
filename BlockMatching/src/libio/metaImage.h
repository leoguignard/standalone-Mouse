
#ifndef METAIMAGE_H
#define METAIMAGE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <ImageIO.h>


extern int writeMetaImageHeader( char *name, _image *im );
extern int writeMetaImage( char *name, _image *im );

extern PTRIMAGE_FORMAT createMetaImageFormat();

#ifdef __cplusplus
}
#endif

#endif
