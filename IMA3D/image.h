#pragma once
/*
*   $Source: /net/users/roe/onderzoek/morphology/skeletons/code/democode/RCS/image.h,v $
*   $Revision: 1.1 $
*   $Date: 2009/06/16 14:19:09 $
*   $State: Exp $
*   $Author: roe $
*
* ------------------------------------------------------------*/
#ifndef _IMAGE_H
#define _IMAGE_H
#endif

#define MAXSTR 200
#define SQR(a) ((a)*(a))
#define PIXEL short
#define BYTE unsigned char

#ifndef TRUE
#define	TRUE 1
#endif
#ifndef FALSE
#define	FALSE 0
#endif

typedef struct {
	int xdim, ydim;
	unsigned char **data;
} ByteImage;

ByteImage *ByteImage_New(int xdim, int ydim);
void ByteImage_GetMinMax(ByteImage *im, BYTE *im_min, BYTE *im_max);
ByteImage *ByteImage_Crop(ByteImage *im, int factor);
void ByteImage_Delete(ByteImage *im);
void ByteImage_CropWrite(ByteImage *im, char name[], int cropfactor);
void ByteImage_WritePGM(ByteImage *im, char name[]);

