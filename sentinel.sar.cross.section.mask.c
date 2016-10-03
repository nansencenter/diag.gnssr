/*
 * This program creates a land mask for a SAR acquisition
 * using the GMT command gmtselect.  Included is the option
 * to mask local maxima in sigma-naught, including connected
 * regions.  Conversion of sigma-naught from raw units to dB
 * is also done here - RD August 2005, September 2016.
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

float **sigo, **angl, **mask, **land, **lats, **lons, *array;

void landswit(int a, int b);
int  sigoswit(int a, int b, float cutoff);
void     swit(int a, int b);
double **alloc_dblmatrix(int nrows, int ncols);
float  **alloc_flmatrix(int nrows, int ncols);
void free_dblmatrix(double **matrix, int nrows, int ncols);
void free_flmatrix(float **matrix, int nrows, int ncols);

main(int argc, char *argv[])
{
    FILE *fpa, *fpb, *fpc, *fopen();
    char date[LEN], line[LOTS], command[LOTS];
    char infila[LEN], infilb[LEN], infilc[LEN], outfile[LEN];
    int a, b, c, d, latind, lonind, lat_end, lon_end, flag;
    float minlat, minlon, midlat, midlon, maxlat, maxlon, crita, critb;
    float maskresol, minang, midang, maxang, refminang, refmidang, refmaxang;

    if (argc != 4) {
      printf("\nUsage: %s 2004-03-30_43859_01.00400.hdr hotspot_cutoff connected_cutoff\n",argv[0]);
      printf("e.g.,  %s 2004-03-30_43859_01.00400.hdr -15.0 -22.0\n\n",argv[0]);
      exit(1);
    }
    sscanf(argv[2],"%f",&crita);
    sscanf(argv[3],"%f",&critb);

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

    sigo  = alloc_flmatrix(lat_end,lon_end);                                  /* allocate some space */
    angl  = alloc_flmatrix(lat_end,lon_end);
    mask  = alloc_flmatrix(lat_end,lon_end);
    land  = alloc_flmatrix(lat_end,lon_end);
    lats  = alloc_flmatrix(lat_end,lon_end);
    lons  = alloc_flmatrix(lat_end,lon_end);
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
    sarread(infila, "angl", date, array, lat_end*lon_end);
    d = 0;
    for (a = 0; a < lat_end; a++)
      for (b = 0; b < lon_end; b++) {
        angl[a][b] = array[d];
        d++;
      }
    sarread(infila, "lats", date, array, lat_end*lon_end);
    d = 0;
    for (a = 0; a < lat_end; a++)
      for (b = 0; b < lon_end; b++) {
        lats[a][b] = array[d];
        d++;
      }
    sarread(infila, "lons", date, array, lat_end*lon_end);
    d = 0;
    for (a = 0; a < lat_end; a++)
      for (b = 0; b < lon_end; b++) {
        lons[a][b] = array[d];
        d++;
      }

    for (a = 0; a < lat_end; a++)                                             /* convert sigo from raw to dB */
      for (b = 0; b < lon_end; b++) {
        if (sigo[a][b] > 0)
          sigo[a][b] = 10.0 * log10f(sigo[a][b]);
        else
          sigo[a][b] = MISSING;
      }

/*  for (a = 0; a < lat_end; a++)                                              * initialize the mask *
      for (b = 0; b < lon_end; b++) {
        land[a][b] = 1.0;
        if (sigo[a][b] < DELMISS)
          mask[a][b] = 1;
        else
          mask[a][b] = 0;
      }

    strcpy(infilb,argv[1]);                                                    * write the lats and lons *
    strcpy(infilc,argv[1]);
    strcat(infilb,".gmtdata");
    strcat(infilc,".gmtland");
    if ((fpb = fopen(infilb,"w")) == NULL) {
      printf("couldn't open %s\n",infilb);
      exit(1);
    }
    printf("writing lat/lon pairs to %s\n",infilb);
    for (a = 0; a < lat_end; a++)
      for (b = 0; b < lon_end; b++)
        fprintf(fpb,"%f %f %d %d\n",lons[a][b],lats[a][b],a,b);
    fclose(fpb);

    if (lats[0][0]         > lats[0][lon_end-1])         maxlat = lats[0][0]                 + 1.0;
    else                                                 maxlat = lats[0][lon_end-1]         + 1.0;
    if (lats[lat_end-1][0] < lats[lat_end-1][lon_end-1]) minlat = lats[lat_end-1][0]         - 1.0;
    else                                                 minlat = lats[lat_end-1][lon_end-1] - 1.0;
    if (lons[0][lon_end-1] > lons[lat_end-1][lon_end-1]) maxlon = lons[0][lon_end-1]         + 1.0;
    else                                                 maxlon = lons[lat_end-1][lon_end-1] + 1.0;
    if (lons[0][0]         < lons[lat_end-1][0])         minlon = lons[0][0]                 - 1.0;
    else                                                 minlon = lons[lat_end-1][0]         - 1.0;
    midlat = (maxlat + minlat) / 2.0;
    midlon = (maxlon + minlon) / 2.0;
    sprintf(command,"/usr/lib/gmt/bin/gmtselect %s -R%.2f/%.2f/%.2f/%.2f -Jq%.2f/1 -Dh -Ns/k/k/k/k > %s\n",
      infilb,minlon,maxlon,minlat,maxlat,midlon,infilc);
    printf("%s", command);
    system(command);

    if ((fpc = fopen(infilc,"r")) == NULL) {                                   * and read the land points *
      printf("couldn't open %s\n",infilc);
      exit(1);
    }
    printf("reading land points from %s\n",infilc);
    while (fscanf(fpc,"%*s %*s %d %d",&a,&b) != EOF) {
                                                landswit(a  ,b);
        * landswit(a-2,b-2); landswit(a-2,b-1); landswit(a-2,b); landswit(a-2,b+1); landswit(a-2,b+2);
          landswit(a-1,b-2); landswit(a-1,b-1); landswit(a-1,b); landswit(a-1,b+1); landswit(a-1,b+2);
          landswit(a  ,b-2); landswit(a  ,b-1); landswit(a  ,b); landswit(a  ,b+1); landswit(a  ,b+2);
          landswit(a+1,b-2); landswit(a+1,b-1); landswit(a+1,b); landswit(a+1,b+1); landswit(a+1,b+2);
          landswit(a+2,b-2); landswit(a+2,b-1); landswit(a+2,b); landswit(a+2,b+1); landswit(a+2,b+2);
                             landswit(a-1,b-1); landswit(a-1,b); landswit(a-1,b+1);
                             landswit(a  ,b-1); landswit(a  ,b); landswit(a  ,b+1);
                             landswit(a+1,b-1); landswit(a+1,b); landswit(a+1,b+1); *
    }
    fclose(fpc);
    sprintf(command,"rm .gmtcommands4 %s %s",infilb,infilc);
    system(command);  */

    for (a = 0; a < lat_end; a++)                                             /* then seed the hotspots and mask */
      for (b = 0; b < lon_end; b++)                                           /* the contiguous regions by */
        sigoswit(a, b, crita);

    flag = 1;                                                                 /* expanding the subdomain outward */
    while (flag != 0) {
      flag = 0;
      for (a = 1; a < lat_end-1; a++)
        for (b = 1; b < lon_end-1; b++)
          if (sigo[a][b] < DELMISS) {
            flag += sigoswit(a+1, b,   critb);
            flag += sigoswit(a,   b+1, critb);
            flag += sigoswit(a,   b-1, critb);
            flag += sigoswit(a-1, b,   critb);
            flag += sigoswit(a+1, b+1, critb);
            flag += sigoswit(a+1, b-1, critb);
            flag += sigoswit(a-1, b+1, critb);
            flag += sigoswit(a-1, b-1, critb);
          }
    }

/*  for (a = 0; a < lat_end; a++)                                              * also trim the edges *
      for (b = 0; b < lon_end; b++) {                                          * of the acquistion and *
        if (mask[a][b] == 1) {
          if (a-2 >= 0 && a-2 < lat_end && b-2 >= 0 && b-2 < lon_end) swit(a-2,b-2);
          if (a-2 >= 0 && a-2 < lat_end && b-1 >= 0 && b-1 < lon_end) swit(a-2,b-1);
          if (a-2 >= 0 && a-2 < lat_end && b   >= 0 && b   < lon_end) swit(a-2,b);
          if (a-2 >= 0 && a-2 < lat_end && b+1 >= 0 && b+1 < lon_end) swit(a-2,b+1);
          if (a-2 >= 0 && a-2 < lat_end && b+2 >= 0 && b+2 < lon_end) swit(a-2,b+2);
          if (a-1 >= 0 && a-1 < lat_end && b-2 >= 0 && b-2 < lon_end) swit(a-1,b-2);
          if (a-1 >= 0 && a-1 < lat_end && b-1 >= 0 && b-1 < lon_end) swit(a-1,b-1);
          if (a-1 >= 0 && a-1 < lat_end && b   >= 0 && b   < lon_end) swit(a-1,b);
          if (a-1 >= 0 && a-1 < lat_end && b+1 >= 0 && b+1 < lon_end) swit(a-1,b+1);
          if (a-1 >= 0 && a-1 < lat_end && b+2 >= 0 && b+2 < lon_end) swit(a-1,b+2);
          if (a   >= 0 && a   < lat_end && b-2 >= 0 && b-2 < lon_end) swit(a  ,b-2);
          if (a   >= 0 && a   < lat_end && b-1 >= 0 && b-1 < lon_end) swit(a  ,b-1);
          if (a   >= 0 && a   < lat_end && b   >= 0 && b   < lon_end) swit(a  ,b);
          if (a   >= 0 && a   < lat_end && b+1 >= 0 && b+1 < lon_end) swit(a  ,b+1);
          if (a   >= 0 && a   < lat_end && b+2 >= 0 && b+2 < lon_end) swit(a  ,b+2);
          if (a+1 >= 0 && a+1 < lat_end && b-2 >= 0 && b-2 < lon_end) swit(a+1,b-2);
          if (a+1 >= 0 && a+1 < lat_end && b-1 >= 0 && b-1 < lon_end) swit(a+1,b-1);
          if (a+1 >= 0 && a+1 < lat_end && b   >= 0 && b   < lon_end) swit(a+1,b);
          if (a+1 >= 0 && a+1 < lat_end && b+1 >= 0 && b+1 < lon_end) swit(a+1,b+1);
          if (a+1 >= 0 && a+1 < lat_end && b+2 >= 0 && b+2 < lon_end) swit(a+1,b+2);
          if (a+2 >= 0 && a+2 < lat_end && b-2 >= 0 && b-2 < lon_end) swit(a+2,b-2);
          if (a+2 >= 0 && a+2 < lat_end && b-1 >= 0 && b-1 < lon_end) swit(a+2,b-1);
          if (a+2 >= 0 && a+2 < lat_end && b   >= 0 && b   < lon_end) swit(a+2,b);
          if (a+2 >= 0 && a+2 < lat_end && b+1 >= 0 && b+1 < lon_end) swit(a+2,b+1);
          if (a+2 >= 0 && a+2 < lat_end && b+2 >= 0 && b+2 < lon_end) swit(a+2,b+2);
        * if (a-1 >= 0 && a-1 < lat_end && b-1 >= 0 && b-1 < lon_end) swit(a-1,b-1);
          if (a-1 >= 0 && a-1 < lat_end && b   >= 0 && b   < lon_end) swit(a-1,b);
          if (a-1 >= 0 && a-1 < lat_end && b+1 >= 0 && b+1 < lon_end) swit(a-1,b+1);
          if (a   >= 0 && a   < lat_end && b-1 >= 0 && b-1 < lon_end) swit(a  ,b-1);
          if (a   >= 0 && a   < lat_end && b   >= 0 && b   < lon_end) swit(a  ,b);
          if (a   >= 0 && a   < lat_end && b+1 >= 0 && b+1 < lon_end) swit(a  ,b+1);
          if (a+1 >= 0 && a+1 < lat_end && b-1 >= 0 && b-1 < lon_end) swit(a+1,b-1);
          if (a+1 >= 0 && a+1 < lat_end && b   >= 0 && b   < lon_end) swit(a+1,b);
          if (a+1 >= 0 && a+1 < lat_end && b+1 >= 0 && b+1 < lon_end) swit(a+1,b+1); *
        }
      }

    strcpy(infilc,argv[1]);                                                    * then open the masking *
    infilc[19] = '\0';                                                         * file (if it exists) *
    strcat(infilc,".mask");
    if ((fpc = fopen(infilc,"r")) == NULL) {
      printf("\nno mask file was found : %s\n",infilc);
      printf("and no beam seams will be masked...\n\n");
    }
    else {                                                                     * and for each mask range *
      printf("\nreading mask file %s\n",infilc);                               * check (linearly) for *
      printf("and masking beam seams...\n\n");                                 * gridbox overlap using *
      while (fgets(line,LOTS,fpc) != NULL) {                                   * incidence angle *
        sscanf(line,"%f %f %f %f",&maskresol,&refminang,&refmidang,&refmaxang);
        for (a = 0; a < lon_end; a++) {
          if (angl[0][0] < angl[0][lon_end-1]) {
            if (a == 0)
              minang = MISSING;
            else
              minang = (angl[0][a-1] + angl[0][a]) / 2.0;
            if (a == lon_end-1)
              maxang = -MISSING;
            else
              maxang = (angl[0][a+1] + angl[0][a]) / 2.0;
          }
          else {
            if (a == 0)
              maxang = -MISSING;
            else
              maxang = (angl[0][a-1] + angl[0][a]) / 2.0;
            if (a == lon_end-1)
              minang = MISSING;
            else
              minang = (angl[0][a+1] + angl[0][a]) / 2.0;
          }
          if ((refminang <= minang && refmaxang >  minang) ||                  * where there is overlap *
              (refminang >= minang && refmaxang <= maxang) ||                  * mask the beam seam(s) *
              (refminang <  maxang && refmaxang >= maxang))
            for (b = 0; b < lat_end; b++)
              swit(b,a);
        }
      }
      fclose(fpc);
    }  */

    strcpy(outfile,argv[1]);                                                  /* write the backscatter */
    outfile[25] = '\0';                                                       /* and the land/sea mask */
    strcat(outfile,".sar.nc");
    d = 0;
    for (a = 0; a < lat_end; a++)
      for (b = 0; b < lon_end; b++)
        array[d++] = sigo[a][b];
    sarwrite(outfile, "sigo", date, array, lat_end*lon_end);
/*  d = 0;
    for (a = 0; a < lat_end; a++)
      for (b = 0; b < lon_end; b++)
        array[d++] = land[a][b];
    sarwrite(outfile, "land", date, array, lat_end*lon_end);  */

    free_flmatrix(sigo,lat_end,lon_end);                                      /* and free the storage */
    free_flmatrix(angl,lat_end,lon_end);
    free_flmatrix(mask,lat_end,lon_end);
    free_flmatrix(land,lat_end,lon_end);
    free_flmatrix(lats,lat_end,lon_end);
    free_flmatrix(lons,lat_end,lon_end);
    free(array);
    exit(0);
}

void landswit(int a, int b)
{
    if (sigo[a][b] > DELMISS)
      sigo[a][b] = MISSING;
    land[a][b] = -1.0;
    mask[a][b] =  1.0;
}

int sigoswit(int a, int b, float cutoff)
{
    int flag;

    flag = 0;                                                                 /* mask any hotspot (assuming */
    if (sigo[a][b] > DELMISS && sigo[a][b] > cutoff) {                        /* it is not already masked) */
      sigo[a][b] = MISSING;
      flag = 1;
    }

    return(flag);
}

void swit(int a, int b)
{
    if (sigo[a][b] > DELMISS)
      sigo[a][b] = MISSING;
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
