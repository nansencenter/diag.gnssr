/*
 * This program creates anomalies in sigo relative to a
 * moving spatial average (following correspondence between
 * sigo anomalies and optical MSS in Kudryavtsev et al. 2012
 * on Imaging mesoscale upper ocean dynamics) - September 2016.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "/home/ricani/prog/include.netcdf/include/cdcnetcdf.h"

#define LEN        100
#define LOTS       2000
#define MISSING    -9e9                      /* generic missing value */
#define DELMISS    -8e9                      /* generic missing value comparison */

float **sigo, **smoo, *array;

double **alloc_dblmatrix(int nrows, int ncols);
float  **alloc_flmatrix(int nrows, int ncols);
void free_dblmatrix(double **matrix, int nrows, int ncols);
void free_flmatrix(float **matrix, int nrows, int ncols);

main(int argc, char *argv[])
{
    FILE *fpa, *fpb, *fpc, *fopen();
    char date[LEN], line[LOTS], command[LOTS];
    char infila[LEN], infilb[LEN], infilc[LEN], outfile[LEN], stema[LEN];
    int a, b, c, d, latind, lonind, lat_end, lon_end, midline, midpixel, gridlen;
    float minlat, minlon, midlat, midlon, maxlat, maxlon, sum, smoosum, subscenesize;
    float maskresol, minang, midang, maxang, refminang, refmidang, refmaxang, smoolen;
    float pixel_spacing, tmpval;

    if (argc != 3) {
      printf("Usage: %s 2004-03-30_43859_01.00040.hdr smoothing_length_km\n",argv[0]);
      printf("e.g.,  %s 2004-03-30_43859_01.00040.hdr 20.0\n",argv[0]);
      exit(1);
    }
    sscanf(argv[2],"%f",&smoolen);

    strcpy(infila,argv[1]);                                                   /* read the header info */
    if ((fpa = fopen(infila,"r")) == NULL) {
      printf("couldn't open %s\n",argv[1]);
      exit(1);
    }
    for (a = 0; a < 27; a++) {
      fgets(line,LEN,fpa);
      if (a == 9)  sscanf(line,"%*s %*s %*s %s",date);
      if (a == 24) sscanf(line,"%*s %d",&lat_end);
      if (a == 25) sscanf(line,"%*s %d",&lon_end);
    }
    fclose(fpa);
    strcpy(stema,argv[1]);  stema[25] = '\0';                                 /* but get resolution from */
    pixel_spacing = atoi(stema+20);                                           /* the file name and compute */
    gridlen = smoolen * 500.0 / pixel_spacing;                                /* a pixel smoothing distance */
    printf("smoothing will be done over %d points NSEW\n", gridlen);

    sigo  = alloc_flmatrix(lat_end,lon_end);                                  /* allocate some space */
    smoo  = alloc_flmatrix(lat_end,lon_end);
    array = (float *)calloc(lat_end*lon_end,sizeof(float));

    strcpy(infila,argv[1]);                                                   /* read the input data */
    infila[25] = '\0';
    strcat(infila,".sar.nc");
    sarread(infila, "sigo", date, array, lat_end*lon_end);
    d = 0;
    for (a = 0; a < lat_end; a++)
      for (b = 0; b < lon_end; b++) {
        sigo[a][b] = array[d];
        d++;
      }

    for (a = 0; a < lat_end; a++)                                             /* initialize the smoothed values */
      for (b = 0; b < lon_end; b++)
        smoo[a][b] = 0.0;

    for (a = 0; a < lat_end; a++)                                             /* then sum over moving domains */
      for (b = 0; b < lon_end; b++) {
        sum = smoosum = 0;
        if (sigo[a][b] > DELMISS)
          for (c = -gridlen; c <= gridlen; c++)
            if (a + c >= 0 && a + c < lat_end)
              for (d = -gridlen; d <= gridlen; d++)
                if (b + d >= 0 && b + d < lon_end) {
                  tmpval = sigo[a+c][b+d];
                  if (tmpval > DELMISS) {
                    sum += tmpval;
                    smoosum += 1.0;
                  }
                }
        if (smoosum == 0)
          smoo[a][b] = MISSING;
        else
          smoo[a][b] = sum / smoosum;
      }

    strcpy(outfile,argv[1]);                                                  /* write the backscatter */
    outfile[25] = '\0';                                                       /* and the land/sea mask */
    strcat(outfile,".sar.nc");
    d = 0;
    for (a = 0; a < lat_end; a++)
      for (b = 0; b < lon_end; b++)
        array[d++] = sigo[a][b] - smoo[a][b];
    sarwrite(outfile, "sice", date, array, lat_end*lon_end);

    free_flmatrix(sigo,lat_end,lon_end);                                      /* and free the storage */
    free_flmatrix(smoo,lat_end,lon_end);
    free(array);
    exit(0);
}


double **alloc_dblmatrix(int nrows, int ncols)
{
    int i, j;
    double **matrix;

    if ((matrix = (double **)calloc((size_t)(nrows),sizeof(double *))) == NULL) {
      printf("*** Error 1 using calloc in alloc_dblmatrix\n");
      exit(1);
    }

    for (i = 0; i < nrows; i++) {
      if ((matrix[i] = (double *)calloc((size_t)(ncols),sizeof(double))) == NULL) {
        printf("*** Error 2 using calloc in alloc_dblmatrix\n");
        for (j = 0; j < i; j++)
          free(matrix[j]);
        free(matrix);
        exit(1);
      }
    }
    return matrix;
}

float **alloc_flmatrix(int nrows, int ncols)
{
    int i, j;
    float **matrix;

    if ((matrix = (float **)calloc((size_t)(nrows),sizeof(float *))) == NULL) {
      printf("*** Error 1 using calloc in alloc_flmatrix\n");
      exit(1);
    }

    for (i = 0; i < nrows; i++) {
      if ((matrix[i] = (float *)calloc((size_t)(ncols),sizeof(float))) == NULL) {
        printf("*** Error 2 using calloc in alloc_flmatrix\n");
        for (j = 0; j < i; j++)
          free(matrix[j]);
        free(matrix);
        exit(1);
      }
    }
    return matrix;
}

void free_dblmatrix(double **matrix, int nrows, int ncols)
{
    int i;

    for (i = nrows-1; i >= 0; i--)
      free(matrix[i]);
    free(matrix);
}

void free_flmatrix(float **matrix, int nrows, int ncols)
{
    int i;

    for (i = 0; i < nrows; i++ )
      free(matrix[i]);
    free(matrix);
}
