#pragma once
/*
*   $Source: /net/users/roe/onderzoek/morphology/skeletons/code/democode/RCS/volume.h,v $
*   $Revision: 1.1 $
*   $Date: 2009/06/16 14:19:09 $
*   $State: Exp $
*   $Author: roe $
*
* ------------------------------------------------------------*/
#ifndef _VOLUME_H
#define _VOLUME_H

#include <stdlib.h>
#include <stdio.h>
#include "image.h"
#define VOXEL unsigned short
#define BYTE unsigned char
#define PI 3.1415926535897932384626433
#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))


#define PR2( f, p)	(void) fprintf ( stdout, f, p)
#define MAXSTR 200
#define MAX(a,b) ((a)>(b) ? (a) : (b))

typedef struct {
	int xdim, ydim, zdim;
	VOXEL ***data;
} Volume;

typedef struct {
	int xdim, ydim, zdim;
	BYTE ***data;
} ByteVolume;

typedef struct {
	int xdim, ydim, zdim;
	int ***data;
} IntVolume;

typedef struct {
	double u, v;
} PlanarPoint;

Volume *Volume_New(int xdim, int ydim, int zdim);
Volume *Volume_ReadSFF(char filename[]);
ByteVolume *ByteVolume_ReadSFF(char* filename);
ByteVolume *ByteVolume_ReadVTK(char* filename);
void Volume_MinMax(Volume *v, VOXEL *min, VOXEL *max);
Volume *Volume_ReadAVS(char filename[]);
void Volume_WriteAVS(Volume *v, const char filename[]);
void Volume_WriteSFF(Volume *v, const char* filename);
void ByteVolume_WriteVTK(ByteVolume *v, const char* filename);
void ByteVolume_WriteSFF(ByteVolume *v, const char* filename);
ByteVolume *ByteVolume_New(int x, int y, int z);
void Volume_Delete(Volume *v);
void ByteVolume_Delete(ByteVolume *v);
//void ByteVolume_Print(char *volume_name, ByteVolume *v);
void ByteVolume_MinMax(ByteVolume *v, BYTE *min, BYTE *max);
IntVolume *IntVolume_New(int x, int y, int z);
void IntVolume_Delete(IntVolume * iv);
char *basename_no_ext(const char* filename, const char *desired_ext);
float *vector(int nl, int nh);
int *ivector(int nl, int nh);
void free_vector(float *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
float mean(float *v, int nl, int nh);
float stdev(float *v, int nl, int nh);
void ByteVolume_Clear(ByteVolume *v);
void ByteVolume_Clamp(ByteVolume *v, BYTE clampval);
ByteImage *ByteVolume_MIP(ByteVolume *v, double phi, double theta, double alpha, int size);
PlanarPoint ProjectPoint(double x, double y, double z, double phi, double theta, double alpha);
void makemip(char *file_name, ByteVolume *v, double phi, double theta, double alpha, int cropfactor);
#endif
