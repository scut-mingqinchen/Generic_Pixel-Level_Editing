#ifndef IMAGE_HDR
#define IMAGE_HDR

//
//  imageHDR.cpp
//  app_tonemapping
//
//  Created by Ana Cambra on 28/10/16.
//
//
#include <iostream>
//Opencv
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <chrono>
#include <strings.h>
#include <sstream>
#include <iomanip>

/* flags indicating which fields in an rgbe_header_info are valid */
#define RGBE_VALID_PROGRAMTYPE 0x01
#define RGBE_VALID_GAMMA       0x02
#define RGBE_VALID_EXPOSURE    0x04

/* return codes for rgbe routines */
#define RGBE_RETURN_SUCCESS 0
#define RGBE_RETURN_FAILURE -1

/* offsets to red, green, and blue components in a data (float) pixel */
#define RGBE_DATA_RED    0
#define RGBE_DATA_GREEN  1
#define RGBE_DATA_BLUE   2
/* number of floats per pixel */
#define RGBE_DATA_SIZE   3


enum rgbe_error_codes {
    rgbe_read_error,
    rgbe_write_error,
    rgbe_format_error,
    rgbe_memory_error,
};
/* default error routine.  change this to change error handling */
int rgbe_error(int rgbe_error_code, const char *msg)
;


typedef struct {
    int valid;            /* indicate which fields are valid */
    char programtype[16]; /* listed at beginning of file to identify it
                           * after "#?".  defaults to "RGBE" */
    float gamma;          /* image has already been gamma corrected with
                           * given gamma.  defaults to 1.0 (no correction) */
    float exposure;       /* a value of 1.0 in an image corresponds to
                           * <exposure> watts/steradian/m^2.
                           * defaults to 1.0 */
} rgbe_header_info;


/* default error routine.  change this to change error handling */
//int rgbe_error(int rgbe_error_code, const char *msg);

/* standard conversion from float pixels to rgbe pixels */
/* note: you can remove the "inline"s if your compiler complains about it */
inline void
float2rgbe(unsigned char rgbe[4], float red, float green, float blue);

/* standard conversion from rgbe to float pixels */
/* note: Ward uses ldexp(col+0.5,exp-(128+8)).  However we wanted pixels */
/*       in the range [0,1] to map back into the range [0,1].            */
inline void
rgbe2float(float *red, float *green, float *blue, unsigned char rgbe[4]);
/* default minimal header. modify if you want more information in header */
int RGBE_WriteHeader(FILE *fp, int width, int height, rgbe_header_info *info);

/* minimal header reading.  modify if you want to parse more information */
int RGBE_ReadHeader(FILE *fp, int *width, int *height, rgbe_header_info *info);

/* simple write routine that does not use run length encoding */
/* These routines can be made faster by allocating a larger buffer and
 fread-ing and fwrite-ing the data in larger chunks */
int RGBE_WritePixels(FILE *fp, float *data, int numpixels);

/* simple read routine.  will not correctly handle run length encoding */
int RGBE_ReadPixels(FILE *fp, float *data, int numpixels);

/* The code below is only needed for the run-length encoded files. */
/* Run length encoding adds considerable complexity but does */
/* save some space.  For each scanline, each channel (r,g,b,e) is */
/* encoded separately for better compression. */

int RGBE_WriteBytes_RLE(FILE *fp, unsigned char *data, int numbytes);

int RGBE_WritePixels_RLE(FILE *fp, float *data, int scanline_width,
                         int num_scanlines);

int RGBE_ReadPixels_RLE(FILE *fp, float *data, int scanline_width,
                        int num_scanlines);


std::string type2str(int type) ;

std::string cvMat_description(const char* name, const cv::Mat& m)
;

cv::Mat load_hdr(const char* filename)
;
#endif // SUPERPIXELS_H

