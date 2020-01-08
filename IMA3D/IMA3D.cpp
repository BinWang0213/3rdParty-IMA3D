/*
SYNOPSIS

This program computes the IMA skeleton of a 3D data set as described in
the paper: W. H. Hesselink, J. B. T. M. Roerdink: Euclidean skeletons of
digital image and volume data in linear time by integer medial axis
transform. IEEE Trans. Pattern Anal. Machine Intell. 30 (2008), pp. 2204--2217.

The skeleton is computed using a feature transform. The feature transform
algorithm is an adaptation of the distance transform algorithm described
in "A general algorithm for computing distance transforms in linear time",
by Arnold Meijster, Jos B.T.M. Roerdink and Wim H. Hesselink, in:
Proceedings Mathematical Morphology and its Applications to Image and
Signal Processing, 2000, pp. 331-340.

The time complexity of the several algorithms is linear in the number of
data elements.

AUTHORS
W. H. Hesselink, J. B. T. M. Roerdink

DESCRIPTION
The construction of the 3D feature transform, applied to the construction
of skeletons.  The boundary is extended with a plane at x = xdim - 1 +
INFTY. The "empty boundary" is therefore allowed.  The pruning parameter gamma
for the skeleton is chosen in such a way that the smaller the parameter,
the more points are admitted.  No filtering (exact IMA) corresponds to the
choice: gamma = 1.

NOTES

The size of the volume is xdimy*dimension*zdim. The feature transform is
implicitly encoded in a single integer. In the skelton routine, a test is
inserted to check that from the two voxels to be compared, at least one is
a FOREGROUND voxel.  To save memory, the skeleton values are written back
into the input array

INPUT:	data set SFF file format.

OUTPUT: data set SFF file format, with the nonzero data elements being part
of the skeleton.

USAGE:  ima3d [INPUTFILE] [-g gamma]

*/

#include <iostream>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "volume.h"
//#include <sys/time.h>
//#include <sys/resource.h>

const int GAMMA = 1;
const int FOREGROUND = 1;
const int BACKGROUND = 0;
const int SKEL = 255;

#define sqr(x) ((x)*(x))
#define sqrlength(i, j, k) ((i)*(i)+(j)*(j)+(k)*(k))

static int			gamma_val = GAMMA;
static IntVolume*	ft;
static ByteVolume*	indata;
static char			input_file[MAXSTR];
static char*		basefilename;
static char			basename[MAXSTR];
static char			skel_file[MAXSTR];

/*************** SKELETONIZATION *****************/
void doubleScan(int ft[], int lim, int ftline[],
	int dtline[], int ss[], int tt[])
{
	int q = 0, j, w;
	ss[0] = tt[0] = 0;
	for (j = 1; j < lim; j++) {
		while (q >= 0 &&
			(j - ss[q])*(j + ss[q] - 2 * tt[q]) < dtline[ss[q]] - dtline[j]) {
			q--;
		}
		if (q < 0) {
			q = 0;
			ss[0] = j;
		}
		else {
			w = 1 +
				((j + ss[q])*(j - ss[q]) + dtline[j] - dtline[ss[q]]) / (2 * (j - ss[q]));
			if (w < lim) {
				q++;
				ss[q] = j;
				tt[q] = w;
			}
		}
	}
	for (j = lim - 1; j >= 0; j--) {
		ft[j] = ftline[ss[q]] * lim + ss[q]; /* encoding */
		if (j == tt[q]) q--;
	}
}

IntVolume *featTrans(ByteVolume *boundary)
{
	/* interpretation (x,y,z) in boundary IFF boundary[x][y][z] == 0
	first phase: construction of feature transform in the x-direction
	Vectors (x, y, z) are encoded as integers by
	encode(x, y, z) = z + zdim * (y + ydim * x) */
	int x, y, z, xy, right, left;
	int xdim = boundary->xdim;
	int ydim = boundary->ydim;
	int zdim = boundary->zdim;
	IntVolume *ftx, *ftxy, *ftxyz;
	int INFTY = 1 + (int)sqrt(sqr(xdim) + sqr(ydim) + sqr(zdim));
	int LIM = MAX(MAX(xdim, ydim), zdim);
	int *ftline, *dtline, *ss, *tt;

	ftline = ivector(0, LIM - 1);
	dtline = ivector(0, LIM - 1);
	ss = ivector(0, LIM - 1);
	tt = ivector(0, LIM - 1);

	/* The pure algorithm require a nonempty boundary; the encoding
	* requires all coordinates nonnegative. We therefore extend the
	* boundary with points with the plane x = xdim - 1 + INFTY.
	* Conditions: (xdim-1)^2 + ... + (zdim-1)^2 < INFTY^2
	* and (xdim-1+INFTY) * ydim * zdim <= Integer.MAX_VALUE */
	ftx = IntVolume_New(ydim, zdim, xdim);
	for (y = 0; y < ydim; y++) {
		for (z = 0; z < zdim; z++) {
			dtline[xdim - 1] = right =
				(boundary->data[xdim - 1][y][z] == 0 ? 0 : INFTY);
			/* INFTY simulates a boundary outside the image */
			for (x = xdim - 2; x >= 0; x--) {
				dtline[x] = right =
					(boundary->data[x][y][z] == 0 ? 0 : right + 1);
			}
			ftx->data[y][z][0] = left = dtline[0];
			for (x = 1; x < xdim; x++) {
				right = dtline[x];
				ftx->data[y][z][x] = left =
					(x - left <= right ? left : x + right);
			}
		}
	}

	/* second phase: construction of feature transform in the xy-direction
	* based on feature transform in the x-direction */
	ftxy = IntVolume_New(zdim, xdim, ydim);
	for (z = 0; z < zdim; z++) {
		for (x = 0; x < xdim; x++) {
			for (y = 0; y < ydim; y++) {
				ftline[y] = xy = ftx->data[y][z][x];
				dtline[y] = sqr(xy - x);
			}
			doubleScan(ftxy->data[z][x], ydim, ftline,
				dtline, ss, tt);
		}
	}
	IntVolume_Delete(ftx);

	/* third phase: construction of feature transform in the xyz-direction
	* based on feature transform in the xy-direction */
	ftxyz = IntVolume_New(xdim, ydim, zdim);
	for (x = 0; x < xdim; x++) {
		for (y = 0; y < ydim; y++) {
			for (z = 0; z < zdim; z++) {
				ftline[z] = xy = ftxy->data[z][x][y];
				dtline[z] = sqr(xy / ydim - x) + sqr(xy % ydim - y);
			}
			doubleScan(ftxyz->data[x][y], zdim, ftline,
				dtline, ss, tt);
		}
	}
	IntVolume_Delete(ftxy);
	free_ivector(ftline, 0, LIM - 1);
	free_ivector(dtline, 0, LIM - 1);
	free_ivector(ss, 0, LIM - 1);
	free_ivector(tt, 0, LIM - 1);

	return ftxyz;
}

inline double euclidnorm(int x, int y, int z)
{
	return sqrt(x * x + y * y + z*z);
}


void compare(BYTE *xskel, BYTE *pskel,
	int x, int y, int z, int p, int q, int r,
	int xf, int pf, int ydim, int zdim)
{
	/* compute feature transform vectors vf = (xf, yf, zf), uf = (pf, qf, rf) */
	int yf, zf, qf, rf, dif, crit;
	zf = xf % zdim; xf = xf / zdim;
	yf = xf % ydim; xf = xf / ydim;
	rf = pf % zdim; pf = pf / zdim;
	qf = pf % ydim; pf = pf / ydim;
	dif = sqr(xf - pf) + sqr(yf - qf) + sqr(zf - rf);

	if (dif > 1 && dif >
		(gamma_val > 0 ? gamma_val * gamma_val /* Constant pruning */
			: gamma_val < 0 ?         /* Linear pruning */
			(sqr(x - xf + p - pf) + sqr(y - yf + q - qf) + sqr(z - zf + r - rf)) * gamma_val * gamma_val
			: /* Square root pruning */
			euclidnorm(x - xf + p - pf, y - yf + q - qf, z - zf + r - rf)
			+ 2 * ((x - p)*(xf - pf) + (y - q)*(yf - qf) + (z - r)*(zf - rf)) + 1.5)
		)
		/*	 C_d=1.5 for d=3; C_d=1 for d=2;  */
	{
		/* the point xskel is in the skeleton
		iff ||m-uf|| <= ||m-vf|| iff crit >= 0, where */
		crit = (xf - pf)*(xf + pf - x - p) + (yf - qf)*(yf + qf - y - q) + (zf - rf)*(zf + rf - z - r);
		if (crit >= 0) *xskel = SKEL;
		if (crit <= 0) *pskel = SKEL;
	}
}


void skeleton(IntVolume *ft, ByteVolume *skel)
{
	int xdim = ft->xdim;
	int ydim = ft->ydim;
	int zdim = ft->zdim;
	int x, y, z;

	for (x = 1; x < xdim; x++)
		for (y = 0; y < ydim; y++)
			for (z = 0; z < zdim; z++)
				if (skel->data[x][y][z] == FOREGROUND || skel->data[x - 1][y][z] == FOREGROUND)
					compare(&skel->data[x][y][z], &skel->data[x - 1][y][z],
						x, y, z, x - 1, y, z,
						ft->data[x][y][z], ft->data[x - 1][y][z], ydim, zdim);

	for (x = 0; x < xdim; x++) {
		for (y = 1; y < ydim; y++)
			for (z = 0; z < zdim; z++)
				if (skel->data[x][y][z] == FOREGROUND || skel->data[x][y - 1][z] == FOREGROUND)
					compare(&skel->data[x][y][z], &skel->data[x][y - 1][z],
						x, y, z, x, y - 1, z,
						ft->data[x][y][z], ft->data[x][y - 1][z], ydim, zdim);
		for (y = 0; y < ydim; y++)
			for (z = 1; z < zdim; z++)
				if (skel->data[x][y][z] == FOREGROUND || skel->data[x][y][z - 1] == FOREGROUND)
					compare(&skel->data[x][y][z], &skel->data[x][y][z - 1],
						x, y, z, x, y, z - 1,
						ft->data[x][y][z], ft->data[x][y][z - 1], ydim, zdim);
	}
}


/*************** MAIN PROGRAM *****************/


int main(int argc, const char **argv)
{
	int i;
	int x, y, z;
	int infile = -1;
	int gammavalue = -1;
	int xdim, ydim, zdim;
	BYTE min, max;


	/* Parse command line params */
	for (i = 1; i<argc; i++)
	{
		if (strcmp(argv[i], "--help") == 0) {
			printf("\nUsage: %s INFILE [-g gamma]\n", argv[0]);
			printf("INFILE is the VTK file (unsigned char ()) to use as input.\n");
			printf("gamma is a value for the pruning parameter (default=1)\n");
			printf("gamma>1: constant pruning; gamma<1: linear pruning; gamma=0: square-root pruning.\n");
			return 0;
		}
		else if (strcmp(argv[i], "-g") == 0) {
			if (i + 1<argc) {
				gammavalue = i + 1;
				i++;
			}
			else printf("Missing value for gamma.\n");
		}
		else infile = i;
	}

	if (infile == -1) {
		printf("Missing input file_name. Use '%s --help' for full help.\n", argv[0]);
		return 0;
	}

	if (gammavalue != -1)
		gamma_val = atof(argv[gammavalue]);

	strcpy(input_file, argv[infile]);

	basefilename = basename_no_ext(input_file, "vtk");

	sprintf(skel_file, "%s_%s%d_%s%s", basefilename, "g=", gamma_val, "skel", ".vtk");
	fprintf(stdout, "gamma = %d\n", gamma_val);
	fflush(stdout);

	// Create datastructures
	indata = ByteVolume_ReadVTK(input_file);
	ByteVolume_MinMax(indata, &min, &max);
	xdim = indata->xdim;
	ydim = indata->ydim;
	zdim = indata->zdim;
	if (MAX(MAX(xdim, ydim), zdim)>1024)
	{
		printf("Maximum input dimension = 1024!\n");
		exit(0);
	}

	if (max > FOREGROUND)
	{
		//  Clamp all input values to FOREGROUND
		ByteVolume_Clamp(indata, FOREGROUND);
	}

	// Construct the feature transform 

	ft = featTrans(indata);

	// Construct the skeleton from the feature transform into 'indata'

	skeleton(ft, indata);

	// Put skeleton points to FOREGROUND, and non-skeleton points to BACKGROUND
	for (x = 0; x<xdim; x++)
		for (y = 0; y<ydim; y++)
			for (z = 0; z<zdim; z++)
			{
				if (indata->data[x][y][z] == FOREGROUND) indata->data[x][y][z] = BACKGROUND;
				if (indata->data[x][y][z] == SKEL) indata->data[x][y][z] = FOREGROUND;
			}

	// Write skeleton data
	std::string output_fname = std::string(basefilename) + "_skel.vtk";
	ByteVolume_WriteVTK(indata, output_fname.c_str());
	printf("Skeleton Output = %s", output_fname.c_str());

	// Cleanup
	ByteVolume_Delete(indata);
	IntVolume_Delete(ft);

	return 0;
}
