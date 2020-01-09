#include <math.h>
#include "volume.h"

#include <iostream>
#include <string>

#include <stdio.h>
#include <string.h>
#include "avs_io.h"



inline unsigned int to_int(unsigned char bytes[4])
{
	return bytes[3] + 256 * (bytes[2] + 256 * (bytes[1] + 256 * bytes[0]));
}

inline unsigned short to_short(unsigned char bytes[2])
{
	return bytes[1] + 256 * bytes[0];
}

inline void int_to_bytes(int nr, unsigned char bytes[4])
{

	bytes[3] = nr % 256;
	nr = nr / 256;
	bytes[2] = nr % 256;
	nr = nr / 256;
	bytes[1] = nr % 256;
	bytes[0] = nr / 256;

}

inline void from_short(unsigned short nr, unsigned char bytes[2])
{
	bytes[1] = nr % 256;
	bytes[0] = nr / 256;
}

Volume *Volume_New(int xdim, int ydim, int zdim)
{
	long i, j;
	Volume *v;

	v = (Volume *)malloc(sizeof(Volume));
	if (v == NULL) {
		printf("Out of memory!\n");
		exit(0);
	}

	v->xdim = xdim;
	v->ydim = ydim;
	v->zdim = zdim;

	v->data = (VOXEL ***)calloc(xdim, sizeof(VOXEL **));
	if (!v->data) {
		printf("Out of memory!\n");
		exit(-1);
	}

	v->data[0] = (VOXEL **)calloc(xdim * ydim, sizeof(VOXEL *));
	if (!v->data[0]) {
		printf("Out of memory!\n");
		exit(-1);
	}

	v->data[0][0] = (VOXEL *)calloc(xdim * ydim * zdim, sizeof(VOXEL));
	if (!v->data[0][0]) {
		printf("Out of memory!\n");
		exit(-1);
	}

	for (j = 1; j<ydim; j++)
		v->data[0][j] = v->data[0][j - 1] + zdim;

	for (i = 1; i<xdim; i++) {
		v->data[i] = v->data[i - 1] + ydim;
		v->data[i][0] = v->data[i - 1][0] + ydim * zdim;
		for (j = 1; j<ydim; j++)
			v->data[i][j] = v->data[i][j - 1] + zdim;
	}

	return v;
}

Volume *Volume_ReadSFF(char filename[])
{
	int i, j, k;
	Volume *temp;
	FILE *f;
	unsigned int txdim, tydim, tzdim;
	unsigned char bytes[4], bytes2[2];
	int xdim, ydim, zdim;
	int xoff, yoff, zoff;
	unsigned char ftype, dims, delim, byteval;
	unsigned short shortval;
	VOXEL val, min, max;

	f = fopen(filename, "rb");

	if (f == NULL) {
		printf("Error opening the file %s\n", filename);
		exit(0);
		/* return temp;*/
	}

	fread(&ftype, sizeof(unsigned char), 1, f);


	fread(&dims, sizeof(unsigned char), 1, f);
	if (dims != 3) {
		printf("Wrong dimensions\n");
		exit(0);
		/* return temp;*/
	}

	fread(bytes, sizeof(unsigned char), 4, f);
	txdim = to_int(bytes);
	fread(bytes, sizeof(unsigned char), 4, f);
	tydim = to_int(bytes);
	fread(bytes, sizeof(unsigned char), 4, f);
	tzdim = to_int(bytes);

	xdim = txdim;
	ydim = tydim;
	zdim = tzdim;

	printf("In Volume_ReadSFF: xdim=%d,ydim=%d,zdim=%d\n", xdim, ydim, zdim);

	temp = Volume_New(xdim, ydim, zdim);

	xoff = (xdim - txdim) / 2;
	yoff = (ydim - tydim) / 2;
	zoff = (zdim - tzdim) / 2;

	/* skip comments */
	delim = fgetc(f);
	while (delim != '\f')
		delim = fgetc(f);

	if (ftype == 1) {
		printf("Reading SFF file: %s... DataType: Bytes\n", filename);
		for (k = 0; k<tzdim; k++)
			for (j = 0; j<tydim; j++)
				for (i = 0; i<txdim; i++) {
					fread(&byteval, sizeof(unsigned char), 1, f);
					val = (VOXEL)byteval;
					temp->data[i + xoff][j + yoff][k + zoff] = val;
				}
	}
	else if (ftype == 2) {
		printf("Reading SFF file: %s... DataType: Shorts\n", filename);
		for (k = 0; k<tzdim; k++)
			for (j = 0; j<tydim; j++)
				for (i = 0; i<txdim; i++) {
					fread(bytes2, sizeof(unsigned char), 2, f);
					shortval = to_short(bytes2);
					temp->data[i + xoff][j + yoff][k + zoff] = (VOXEL)shortval;
				}
	}
	else {
		printf("Unknown datatype\n");
		exit(0);
		/* return temp;*/
	}

	fclose(f);

	Volume_MinMax(temp, &min, &max);
	printf("In Volume_ReadSFF of %s: min=%d, max=%d\n", filename, min, max);

	return temp;
}


void ByteVolume_WriteVTK(ByteVolume* v, const char* filename)
{
	FILE* f = fopen(filename, "wb");

	char buf[1024];

	strcpy(buf, "# vtk DataFile Version 3.0\n");
	fwrite(buf, sizeof(char), strlen(buf), f);

	strcpy(buf, "#generated by IMA3D\n");
	fwrite(buf, sizeof(char), strlen(buf), f);

	strcpy(buf, "BINARY\n");
	fwrite(buf, sizeof(char), strlen(buf), f);

	strcpy(buf, "DATASET STRUCTURED_POINTS\n");
	fwrite(buf, sizeof(char), strlen(buf), f);

	sprintf(buf, "DIMENSIONS %d %d %d\n", v->xdim + 1, v->ydim + 1, v->zdim + 1);	//Needed to add +1 if we want cell-data
	fwrite(buf, sizeof(char), strlen(buf), f);

	strcpy(buf, "ORIGIN 0 0 0\n");
	fwrite(buf, sizeof(char), strlen(buf), f);

	strcpy(buf, "SPACING 1 1 1\n");
	fwrite(buf, sizeof(char), strlen(buf), f);

	int size = v->xdim*v->ydim*v->zdim;
	sprintf(buf, "CELL_DATA %d\n", size);							//Replace this by POINT_DATA if you want point data
	fwrite(buf, sizeof(char), strlen(buf), f);

	strcpy(buf, "SCALARS voxel_data unsigned_char\n");
	fwrite(buf, sizeof(char), strlen(buf), f);

	strcpy(buf, "LOOKUP_TABLE default\n");
	fwrite(buf, sizeof(char), strlen(buf), f);


	for (int k = 0; k<v->zdim; k++)
		for (int j = 0; j<v->ydim; j++)
			for (int i = 0; i<v->xdim; i++)
			{
				unsigned char byteval = v->data[i][j][k];
				fwrite(&byteval, sizeof(unsigned char), 1, f);
			}

	fclose(f);
}

void Volume_WriteSFF(Volume *v, const char* filename)
{
	unsigned char bytes[4];
	unsigned char ftype, dims, delim;
	int x, y, z;
	int i, j, k;
	VOXEL min, max;

	Volume_MinMax(v, &min, &max);
	printf("In  VolumeData_Write of %s: min=%d, max=%d\n", filename, min, max);

	FILE* f = fopen(filename, "wb");

	dims = 3;

	if (max <= 255)						//write volume as bytes
	{
		printf("Writing SFF file: %s... DataType: Bytes\n", filename);

		ftype = 1;
		fwrite(&ftype, sizeof(unsigned char), 1, f);
		fwrite(&dims, sizeof(unsigned char), 1, f);

		int_to_bytes(v->xdim, bytes);
		fwrite(bytes, sizeof(unsigned char), 4, f);
		int_to_bytes(v->ydim, bytes);
		fwrite(bytes, sizeof(unsigned char), 4, f);
		int_to_bytes(v->zdim, bytes);
		fwrite(bytes, sizeof(unsigned char), 4, f);

		delim = '\f';
		fwrite(&delim, sizeof(unsigned char), 1, f);

		for (k = 0; k<v->zdim; k++)
			for (j = 0; j<v->ydim; j++)
				for (i = 0; i<v->xdim; i++)
				{
					unsigned char byteval = v->data[i][j][k];
					fwrite(&byteval, sizeof(unsigned char), 1, f);
				}
	}
	else									//write volume as shorts
	{
		printf("Writing SFF file: %s... DataType: Shorts\n", filename);

		ftype = 2;
		fwrite(&ftype, sizeof(unsigned char), 1, f);
		fwrite(&dims, sizeof(unsigned char), 1, f);

		int_to_bytes(v->xdim, bytes);
		fwrite(bytes, sizeof(unsigned char), 4, f);
		int_to_bytes(v->ydim, bytes);
		fwrite(bytes, sizeof(unsigned char), 4, f);
		int_to_bytes(v->zdim, bytes);
		fwrite(bytes, sizeof(unsigned char), 4, f);

		delim = '\f';
		fwrite(&delim, sizeof(unsigned char), 1, f);

		for (z = 0; z<v->zdim; z++)
			for (y = 0; y<v->ydim; y++)
				for (x = 0; x<v->xdim; x++)
				{
					unsigned short shortval = v->data[x][y][z];/* VOL(z,y,x);*/
					from_short(shortval, bytes);
					fwrite(bytes, sizeof(unsigned char), 2, f);
				}
	}

	fclose(f);
}

void Volume_Delete(Volume *v)
{
	free(v->data[0][0]);
	free(v->data[0]);
	free(v->data);
}

ByteVolume *ByteVolume_New(int xdim, int ydim, int zdim)
{
	long i, j;
	ByteVolume *v;

	v = (ByteVolume *)malloc(sizeof(ByteVolume));
	if (v == NULL) {
		printf("Out of memory!\n");
		exit(0);
	}

	v->xdim = xdim;
	v->ydim = ydim;
	v->zdim = zdim;

	v->data = (BYTE ***)calloc(xdim, sizeof(BYTE **));
	if (!v->data) {
		printf("Out of memory!\n");
		exit(-1);
	}

	v->data[0] = (BYTE **)calloc(xdim * ydim, sizeof(BYTE *));
	if (!v->data[0]) {
		printf("Out of memory!\n");
		exit(-1);
	}

	v->data[0][0] = (BYTE *)calloc(xdim * ydim * zdim, sizeof(BYTE));
	if (!v->data[0][0]) {
		printf("Out of memory!\n");
		exit(-1);
	}

	for (j = 1; j<ydim; j++)
		v->data[0][j] = v->data[0][j - 1] + zdim;

	for (i = 1; i<xdim; i++) {
		v->data[i] = v->data[i - 1] + ydim;
		v->data[i][0] = v->data[i - 1][0] + ydim * zdim;
		for (j = 1; j<ydim; j++)
			v->data[i][j] = v->data[i][j - 1] + zdim;
	}

	return v;
}

void ByteVolume_Delete(ByteVolume *v)
{
	free(v->data[0][0]);
	free(v->data[0]);
	free(v->data);
}

ByteVolume *ByteVolume_ReadSFF(char* filename)
{
	int i, j, k;
	ByteVolume *temp;
	FILE *f;
	unsigned char bytes[4];
	unsigned int txdim, tydim, tzdim;
	int xdim, ydim, zdim;
	int xoff, yoff, zoff;
	unsigned char ftype, dims, delim, byteval;

	f = fopen(filename, "rb");

	if (f == NULL) {
		printf("Error opening the file %s\n", filename);
		exit(0);
	}

	fread(&ftype, sizeof(unsigned char), 1, f);

	printf("ftype %d\n", ftype);


	fread(&dims, sizeof(unsigned char), 1, f);
	if (dims != 3) {
		printf("Wrong dimensions\n");
		exit(0);
	}

	fread(bytes, sizeof(unsigned char), 4, f);
	txdim = to_int(bytes);
	fread(bytes, sizeof(unsigned char), 4, f);
	tydim = to_int(bytes);
	fread(bytes, sizeof(unsigned char), 4, f);
	tzdim = to_int(bytes);

	xdim = txdim;
	ydim = tydim;
	zdim = tzdim;

	temp = ByteVolume_New(xdim, ydim, zdim);

	xoff = (xdim - txdim) / 2;
	yoff = (ydim - tydim) / 2;
	zoff = (zdim - tzdim) / 2;

	/* skip comments */
	delim = fgetc(f);
	while (delim != '\f')
		delim = fgetc(f);

	if (ftype == 1) {
		for (k = 0; k<tzdim; k++)
			for (j = 0; j<tydim; j++)
				for (i = 0; i<txdim; i++) {
					fread(&byteval, sizeof(unsigned char), 1, f);
					temp->data[i + xoff][j + yoff][k + zoff] = byteval;
				}
	}
	else {
		printf("Input is not of type byte\n");
		exit(0);
	}

	fclose(f);

	return temp;
}



ByteVolume *ByteVolume_ReadVTK(char* filename)
{
	int i, j, k;
	ByteVolume *temp;
	int xdim, ydim, zdim;
	int xoff, yoff, zoff;
	int dx, dy, dz, size;
	unsigned char byteval;

	FILE* f = fopen(filename, "rb");

	if (!f)
	{
		printf("Error opening the file %s\n", filename);
		exit(0);
	}

	char buf[1024], dset_type[100], stype[100];
	bool bin;

	for (;;)
	{
		fgets(buf, 1024, f);
		if (buf[0] != '#') break;
	}

	if (!strcmp(buf, "BINARY\n")) bin = true;
	else if (!strcmp(buf, "ASCII\n")) bin = false;
	else
	{
		printf("Error: expected BINARY/ASCII, got: %s\n", buf);
		exit(1);
	}

	fscanf(f, "%s %s", buf, dset_type);
	if (!strcmp(buf, "DATASET"))
	{
		if (strcmp(dset_type, "STRUCTURED_POINTS"))
		{
			printf("Error: expected structured points, got: %s\n", dset_type);
			exit(1);
		}
	}
	else
	{
		printf("Error: expected DATASET, got: %s\n", buf);
		exit(1);
	}

	fscanf(f, "%s %u %u %u", buf, &xdim, &ydim, &zdim);
	if (strcmp(buf, "DIMENSIONS"))
	{
		printf("Error: expected DIMENSIONS, got: %s\n", buf);
		exit(1);
	}

	fscanf(f, "%s %d %d %d", buf, &xoff, &yoff, &zoff);
	if (strcmp(buf, "ORIGIN"))
	{
		printf("Error: expected ORIGIN, got: %s\n", buf);
		exit(1);
	}

	fscanf(f, "%s %d %d %d", buf, &dx, &dy, &dz);
	if (strcmp(buf, "SPACING"))
	{
		printf("Error: expected SPACING, got: %s\n", buf);
		exit(1);
	}

	fscanf(f, "%s %d", buf, &size);
	if (strcmp(buf, "POINT_DATA"))
	{
		printf("Error: expected POINT_DATA, got: %s\n", buf);
		exit(1);
	}

	if (size != xdim*ydim*zdim)
	{
		printf("Error: expected to have %d point scalars, file has only %d scalars", xdim*ydim*zdim, size);
		exit(1);
	}

	fscanf(f, "%s %*s %s", buf, stype);
	if (strcmp(buf, "SCALARS"))
	{
		printf("Error: expected SCALARS, got: %s\n", buf);
		exit(1);
	}
	if (strcmp(stype, "unsigned_char"))
	{
		printf("Error: expected scalars of unsigned char type, got: %s\n", stype);
		exit(1);
	}

	fscanf(f, "%s %*s", buf);
	if (strcmp(buf, "LOOKUP_TABLE"))
	{
		printf("Error: expected LOOKUP_TABLE, got: %s\n", buf);
		exit(1);
	}

	xoff = 0;
	yoff = 0;
	zoff = 0;

	printf("Volume: %d %d %d\n", xdim, ydim, zdim);

	temp = ByteVolume_New(xdim, ydim, zdim);


	for (k = 0; k<zdim; k++)
		for (j = 0; j<ydim; j++)
			for (i = 0; i<xdim; i++)
			{
				fread(&byteval, sizeof(unsigned char), 1, f);
				temp->data[i + xoff][j + yoff][k + zoff] = byteval;
			}

	fclose(f);

	return temp;
}



void ByteVolume_MinMax(ByteVolume *v, BYTE *min, BYTE *max)
{
	int x, y, z;
	BYTE val, minval, maxval;

	minval = maxval = v->data[0][0][0];

	for (x = 0; x<v->xdim; x++)
		for (y = 0; y<v->ydim; y++)
			for (z = 0; z<v->zdim; z++) {
				val = v->data[x][y][z];
				if (val > maxval) maxval = val;
				if (val < minval) minval = val;
			}

	*min = minval;
	*max = maxval;

}

void ByteVolume_Clamp(ByteVolume *v, BYTE clampval)
{
	int x, y, z;
	BYTE val;

	for (x = 0; x<v->xdim; x++)
		for (y = 0; y<v->ydim; y++)
			for (z = 0; z<v->zdim; z++) {
				val = v->data[x][y][z];
				if (val > clampval)
					v->data[x][y][z] = clampval;
			}
}


void ByteVolume_WriteSFF(ByteVolume *v, const char filename[])
{
	unsigned char bytes[4], byteval;
	FILE *f;
	unsigned char ftype, dims, delim;
	int i, j, k;
	BYTE min, max;


	ByteVolume_MinMax(v, &min, &max);
	f = fopen(filename, "wb");


	printf("Writing SFF file: %s... DataType: Bytes\n", filename);

	ftype = 1;
	fwrite(&ftype, sizeof(unsigned char), 1, f);

	dims = 3;

	fwrite(&dims, sizeof(unsigned char), 1, f);
	int_to_bytes(v->xdim, bytes);
	fwrite(bytes, sizeof(unsigned char), 4, f);
	int_to_bytes(v->ydim, bytes);
	fwrite(bytes, sizeof(unsigned char), 4, f);
	int_to_bytes(v->zdim, bytes);
	fwrite(bytes, sizeof(unsigned char), 4, f);

	delim = '\f';
	fwrite(&delim, sizeof(unsigned char), 1, f);


	for (k = 0; k<v->zdim; k++)
		for (j = 0; j<v->ydim; j++)
			for (i = 0; i<v->xdim; i++) {
				byteval = v->data[i][j][k];
				fwrite(&byteval, sizeof(unsigned char), 1, f);
			}

	fclose(f);
}


void nrerror(char error_text[])
{
	fprintf(stderr, "run-time error...\n");
	fprintf(stderr, "%s\n", error_text);
	exit(1);
}


float *vector(int nl, int nh)
{
	float *v;

	v = (float *)malloc((unsigned)(nh - nl + 1) * sizeof(float));
	if (!v) nrerror("allocation failure in vector()");
	return v - nl;
}

int *ivector(int nl, int nh)
{
	int *v;

	v = (int *)malloc((unsigned)(nh - nl + 1) * sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v - nl;
}

void free_vector(float* v, int nl, int nh)
{
	free((char*)(v + nl));
}

void free_ivector(int* v, int nl, int nh)
{
	free((char*)(v + nl));
}

float mean(float *v, int nl, int nh)
/* Return the mean of vector v */
{
	int k;
	float m = 0.0;
	float sum = 0.0;

	for (k = nl; k <= nh; k++) sum += v[k];
	if (nh >= nl) m = sum / (nh - nl + 1);
	else nrerror("zero length of vector in function: mean");

	return m;
}

float stdev(float *v, int nl, int nh)
/* Return the standard deviation of vector v */
{
	int k;
	float m = 0.0, std = 0.0;
	float sum = 0.0;

	m = mean(v, nl, nh);

	for (k = nl; k <= nh; k++) sum += (v[k] - m)*(v[k] - m);
	if (nh >= nl) 	std = sqrt(sum / (nh - nl + 1));
	else nrerror("zero length of vector in function: stdev");

	return std;
}


void ByteVolume_Clear(ByteVolume *v) {
	int x, y, z;

	for (x = 0; x<v->xdim; x++)
		for (y = 0; y<v->ydim; y++)
			for (z = 0; z<v->zdim; z++) {
				v->data[x][y][z] = 0;
			}
}

ByteImage *ByteVolume_MIP(ByteVolume *v, double phi, double theta, double alpha, int size)
/* ByteVolume_MIP assumes square images */
{
	ByteImage *im;
	int x, y, z, u1f, u2f;
	int xdim, ydim, zdim;
	int xhalf, yhalf, zhalf;
	double offset1, offset2;
	double cphi, sphi, ctheta, stheta, u1, u2;
	PlanarPoint P;
	BYTE val;


	ctheta = cos(theta);
	stheta = sin(theta);
	cphi = cos(phi);
	sphi = sin(phi);
	xdim = v->xdim;
	ydim = v->ydim;
	zdim = v->zdim;
	xhalf = xdim / 2;
	yhalf = ydim / 2;
	zhalf = zdim / 2;

	offset1 = size / 2 + 0.5;
	offset2 = size / 2 + 0.5;

	im = ByteImage_New(size, size);	/* im is initialized to zero */

	for (x = 0; x<v->xdim; x++)
		for (y = 0; y<v->ydim; y++) {
			for (z = 0; z<v->zdim; z++) {
				P = ProjectPoint((double)(x - xhalf), (double)(y - yhalf), (double)(z - zhalf), phi, theta, alpha);
				u1 = P.u;
				u2 = P.v;
				u1f = floor(u1 + offset1);
				u2f = floor(u2 + offset2);
				val = v->data[x][y][z];
				im->data[u2f][u1f] = MAX(im->data[u2f][u1f], val);
			}
		}

	return im;
}



PlanarPoint ProjectPoint(double x, double y, double z, double phi, double theta, double alpha)
{
	PlanarPoint P;
	double cphi, sphi, ctheta, stheta, calpha, salpha;

	ctheta = cos(theta);
	stheta = sin(theta);
	cphi = cos(phi);
	sphi = sin(phi);
	calpha = cos(alpha);
	salpha = sin(alpha);


	P.u = x*ctheta*(cphi*calpha + sphi*salpha) + y*(sphi*ctheta*calpha + cphi*salpha)\
		- z*stheta*calpha;
	P.v = x*(-cphi*ctheta*salpha - sphi*calpha) + y*(-sphi*ctheta*salpha + cphi*calpha)\
		+ z*stheta*salpha;

	return P;
}

void makemip(char *input_file_name, ByteVolume *v, double phi, double theta, double alpha, int cropfactor)
{
	ByteImage *projv;
	char mip_file_name[MAXSTR];
	int mipsize;

	printf("\nMIP of %s\n", input_file_name);
	mipsize = MAX(MAX(v->xdim, v->ydim), v->zdim);
	sprintf(mip_file_name, "%s_%s", input_file_name, "mip.pgm");
	projv = ByteVolume_MIP(v, phi, theta, alpha, mipsize);
	ByteImage_CropWrite(projv, mip_file_name, cropfactor);
	printf("Wrote %s\n", mip_file_name);
	ByteImage_Delete(projv);
}


IntVolume *IntVolume_New(int xdim, int ydim, int zdim)
{
	long i, j;
	IntVolume *v;

	v = (IntVolume *)malloc(sizeof(IntVolume));
	if (v == NULL) {
		printf("Out of memory!\n");
		exit(0);
	}

	v->xdim = xdim;
	v->ydim = ydim;
	v->zdim = zdim;

	v->data = (int ***)calloc(xdim, sizeof(int **));
	if (!v->data) {
		printf("Out of memory!\n");
		exit(-1);
	}

	v->data[0] = (int **)calloc(xdim * ydim, sizeof(int *));
	if (!v->data[0]) {
		printf("Out of memory!\n");
		exit(-1);
	}

	v->data[0][0] = (int *)calloc(xdim * ydim * zdim, sizeof(int));
	if (!v->data[0][0]) {
		printf("Out of memory!\n");
		exit(-1);
	}

	for (j = 1; j<ydim; j++)
		v->data[0][j] = v->data[0][j - 1] + zdim;

	for (i = 1; i<xdim; i++) {
		v->data[i] = v->data[i - 1] + ydim;
		v->data[i][0] = v->data[i - 1][0] + ydim * zdim;
		for (j = 1; j<ydim; j++)
			v->data[i][j] = v->data[i][j - 1] + zdim;
	}

	return v;
}

void IntVolume_Delete(IntVolume *v)
{
	free(v->data[0][0]);
	free(v->data[0]);
	free(v->data);
}

Volume *Volume_ReadAVS(char filename[])
{

	avs_header header;
	FILE *fp;
	unsigned char byteValue;
	short value;
	int i, j, k;
	Volume *temp;

	fp = fopen(filename, "r");
	if (!fp)
	{
		printf("Cannot open the input AVS field file\n");
	}

	/* Read AVS Header */
	_avs_read_header(fp, &header);

	printf("In Volume_ReadAVS: xdim=%d,ydim=%d,zdim=%d\n", header.dim1, header.dim2, header.dim3);

	/* Create a volume */
	temp = Volume_New(header.dim1, header.dim2, header.dim3);

	switch (header.datatype) {
	case 1: /* bytes */
		for (k = 0; k<header.dim3; k++) {
			for (j = 0; j<header.dim2; j++) {
				for (i = 0; i<header.dim1; i++) {
					fread(&byteValue, sizeof(unsigned char), 1, fp);
					value = (VOXEL)byteValue;
					temp->data[i][j][k] = value;
				}
			}
		}
		printf("Reading AVS file: %s... DataType: Bytes\n", filename);
		break;
	case 2: /* shorts */
	case 3: /* shorts */
		for (k = 0; k<header.dim3; k++) {
			for (j = 0; j<header.dim2; j++) {
				for (i = 0; i<header.dim1; i++) {
					fread(&value, sizeof(short), 1, fp);
					temp->data[i][j][k] = (VOXEL)value;
				}
			}
		}
		printf("Reading AVS file: %s... DataType: Shorts\n", filename);
		break;
	default:
		printf("read AVS: Wrong datatype\n");
		fclose(fp);
	}
	fclose(fp);


	return temp;

}

void Volume_WriteAVS(Volume *v, const char filename[])
{

	avs_header avs_head;
	VOXEL minVal, maxVal;
	FILE *fp;
	int i, j, k;
	unsigned char byteVal;
	VOXEL voxelVal;

	avs_head.ndim = 3;
	avs_head.dim1 = v->xdim;
	avs_head.dim2 = v->ydim;
	avs_head.dim3 = v->zdim;
	avs_head.min_x = 0;
	avs_head.min_y = 0;
	avs_head.min_z = 0;
	avs_head.max_x = v->xdim - 1;
	avs_head.max_y = v->ydim - 1;
	avs_head.max_z = v->zdim - 1;
	avs_head.filetype = 0;
	avs_head.skip = 0;
	avs_head.nspace = 3;
	avs_head.veclen = 1;
	avs_head.dataname[0] = '\0';


	Volume_MinMax(v, &minVal, &maxVal);
	/*	    if ((minVal >= 0) && (maxVal <= 255)) */
	if (maxVal <= 255)
		avs_head.datatype = 1;
	else
		avs_head.datatype = 2;



	fp = fopen(filename, "w");
	if (!fp)
	{
		printf("Cannot open file: %s for writing the AVS field\n", filename);
		return;
	}

	switch (avs_head.datatype) {
	case 1:
		avs_head.datatype = 1;
		avs_write_header(fp, &avs_head);
		for (k = 0; k<avs_head.dim3; k++) {
			for (j = 0; j<avs_head.dim2; j++) {
				for (i = 0; i<avs_head.dim1; i++) {
					byteVal = (unsigned char)v->data[i][j][k];
					fwrite(&byteVal, sizeof(unsigned char), 1, fp);
				}
			}
		}
		printf("Writing AVS file: %s... DataType: Bytes\n", filename);
		break;
	case 2:
	case 3:
		avs_head.datatype = 3;
		avs_write_header(fp, &avs_head);
		for (k = 0; k<avs_head.dim3; k++) {
			for (j = 0; j<avs_head.dim2; j++) {
				for (i = 0; i<avs_head.dim1; i++) {
					voxelVal = v->data[i][j][k];
					fwrite(&voxelVal, sizeof(short), 1, fp);
				}
			}
		}
		printf("Writing AVS file: %s... DataType: Shorts\n", filename);
		break;
	default:
		printf("writeAVS: Datatype is not supported\n");
		fclose(fp);
		return;
	}
	fclose(fp);
}

void Volume_MinMax(Volume *v, VOXEL *min, VOXEL *max)
{
	int x, y, z;
	VOXEL val, minval, maxval;

	minval = maxval = v->data[0][0][0];

	for (x = 0; x<v->xdim; x++)
		for (y = 0; y<v->ydim; y++)
			for (z = 0; z<v->zdim; z++) {
				val = v->data[x][y][z];
				if (val > maxval) maxval = val;
				if (val < minval) minval = val;
			}

	*min = minval;
	*max = maxval;

}


char *basename_no_ext(const char* filename, const char *desired_ext)
{
	
	//Simple C++ basename extraction
	std::string fullname(filename);
	size_t lastindex = fullname.find_last_of(".");
	std::string fname_base = fullname.substr(0, lastindex);

	//Cast back to C
	char *basefilename = new char[fname_base.size() + 1];
	std::copy(fname_base.begin(), fname_base.end(), basefilename);
	basefilename[fname_base.size()] = '\0'; // don't forget the terminating 0

	/*
	char *dot;
	char filename_copy[1024];
	int i, match;
	char *ext;

	strcpy(filename_copy, filename);


	basefilename = basename(filename_copy);

	dot = strrchr(basefilename, '.');
	ext = dot + 1;

	if (dot && ((dot - basefilename) < strlen(basefilename)))
	{
		match = 1;

		for (i = 0; i < strlen(desired_ext); i++)
		{
			if ((ext[i] == '\0') || (tolower(ext[i]) != tolower(desired_ext[i])))
			{
				match = 0;
				break;
			}
		}
		if (match) {
			*dot = '\0';
		}
		else {
			printf("This is not a file of type: %s\n", desired_ext);
			exit(0);
		}

	}
	*/

	return basefilename;
}