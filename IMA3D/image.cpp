/*
*   $Source: /net/users/roe/onderzoek/morphology/skeletons/code/democode/RCS/image.c,v $
*   $Revision: 1.1 $
*   $Date: 2009/06/16 14:18:59 $
*   $State: Exp $
*   $Author: roe $
*
* ------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "image.h"
#define F_MAX 255
#define F_MIN 0

ByteImage *ByteImage_New(int xdim, int ydim)
{
	ByteImage *im;
	int i;

	im = (ByteImage *)malloc(sizeof(ByteImage));
	if (im == NULL) {
		printf("Out of memory!\n");
		exit(1);
	}

	im->xdim = xdim;
	im->ydim = ydim;

	im->data = (BYTE **)calloc(ydim, sizeof(BYTE*));

	if (!im->data) {
		printf("Out of memory!\n");
		exit(1);
	}

	im->data[0] = (BYTE *)calloc(xdim*ydim, sizeof(BYTE));

	if (!im->data[0]) {
		printf("Out of memory!\n");
		exit(1);
	}

	for (i = 1; i<ydim; i++)
		im->data[i] = im->data[i - 1] + xdim;

	return im;
}


void ByteImage_Delete(ByteImage *im)
{
	free(im->data[0]);
	free(im->data);
}


void ByteImage_GetMinMax(ByteImage *im, BYTE *im_min, BYTE *im_max)
{
	int i, size;
	BYTE min, max, val;
	BYTE *imptr;

	size = (im->xdim) * (im->ydim);

	max = min = im->data[0][0];

	imptr = im->data[0];

	for (i = 0; i<size; i++) {
		val = imptr[i];
		if (val > max) max = val;
		if (val < min) min = val;
	}

	*im_min = min;
	*im_max = max;

}

void ByteImage_WritePGM(ByteImage *im, char name[])
/* Convert a double field to a PGM gray value image. */
/* ByteImage_WritePGM rescales values between 0 and 255 */
{
	int i, xdim, ydim;
	BYTE min, max, *data;
	double scalefac;
	FILE *f;

	xdim = im->xdim;
	ydim = im->ydim;

	ByteImage_GetMinMax(im, &min, &max);
	/* If all values are 1, map them to 255 */
	if (min == max) {
		if (min == 1) min = 0;
		max = min + 1;
	}
	scalefac = 255.0 / (max - min);

	data = im->data[0];
	f = fopen(name, "w");
	fprintf(f, "P5 \n%d %d\n255\n", ydim, xdim);
	for (i = 0; i<xdim*ydim; i++)
		fprintf(f, "%c", (char)(floor((data[i] - min)*scalefac)));
	fclose(f);
}

ByteImage *ByteImage_Crop(ByteImage *im, int factor)
{
	int i, j;
	int xsize, ysize, xnew, ynew, xoffset, yoffset;
	ByteImage *imnew;

	xsize = im->xdim;
	ysize = im->ydim;

	xnew = xsize / factor;
	ynew = ysize / factor;
	xoffset = (xsize - xnew) / 2;
	yoffset = (ysize - ynew) / 2;

	imnew = ByteImage_New(xnew, ynew);

	for (i = 0; i < imnew->xdim; i++)
		for (j = 0; j < imnew->ydim; j++)
			imnew->data[i][j] = im->data[xoffset + i][yoffset + j];

	return imnew;
}

void ByteImage_CropWrite(ByteImage *im, char name[], int cropfactor)
{
	ByteImage *imnew;

	imnew = ByteImage_Crop(im, cropfactor);
	ByteImage_WritePGM(imnew, name);

	ByteImage_Delete(imnew);
}


