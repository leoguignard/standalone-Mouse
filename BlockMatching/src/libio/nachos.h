/* ad-hoc reader/writer for Nachos raw data
 */

#ifndef NACHOS_H
#define NACHOS_H


#include <ImageIO.h>


#ifdef __cplusplus
extern "C" {
#endif



int readNachosImage(const char *basename,_image *im);

/** Writes the given image body in an already opened file.
    @param im image descriptor */
int _writeNachosData(const _image *im);

/** test if an image file has the raw format
    @return -1 */
int testNachosHeader(char *magic,const char *filename);

/** calls _writeInrimageHeader and _writeInrimageData
    @param basename the file name without any extension 
    @param im structure
    @return same as  _writeInrimageHeader*/
int writeNachosImage(char *basename,_image *im);

/** creates an return the file format structure associated with the Inrimage file format */
PTRIMAGE_FORMAT createNachosFormat();

#ifdef __cplusplus
}
#endif

#endif
