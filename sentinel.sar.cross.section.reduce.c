/*
 * This program is designed to smooth a radar cross section field using
 * a binomial smoother similar to that of Koch (TGRS 2004).  The Sentinel
 * scene is assumed to be in netCDF format (i.e., following a conversion
 * from SAFE format) - RD July 2005, September 2016.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "/home/ricani/prog/include.netcdf/include/cdcnetcdf.h"

#define LEN        200
#define LOTS       2000
#define MISSING    -9e9                      /* generic missing value */
#define DELMISS    -8e9                      /* generic missing value comparison */

double **alloc_dblmatrix(int nrows, int ncols);
float  **alloc_flmatrix(int nrows, int ncols);
void free_dblmatrix(double **matrix, int nrows, int ncols);
void free_flmatrix(float **matrix, int nrows, int ncols);

main(int argc, char *argv[])
{
    FILE *fpa, *fpb, *fpc, *fopen();
    char infila[LEN], infilb[LEN], infilc[LEN], outfile[LEN], line[LEN],
         stema[LEN], stemb[LEN], date[LEN], resolution[LEN], command[LEN];
    int a, b, c, d, nlines, npixels, subscenesize, reduce, gridint, aziblocks, rngblocks,
        startline_valid,  numlines_valid,  stopline_valid, midline,  linenum, subscene_centreline,
        startpixel_valid, numpixels_valid, stopixel_valid, midpixel, joffset, subscene_centrepixel,
        nlines2skipbeg, nlines2skipend, npixels2skipbeg, npixels2skipend;
    double pixel_spacing, sum, count, resol, ullat, ullon, urlat, urlon, lllat, lllon, lrlat, lrlon,
           firstline_latitude_firstpixel, firstline_latitude_lastpixel, lastline_latitude_firstpixel,
           lastline_latitude_lastpixel, firstline_longitude_firstpixel, firstline_longitude_lastpixel,
           lastline_longitude_firstpixel, lastline_longitude_lastpixel, **smoo, smoosum;
    float tmpval, **sigo, **angl, **land, **lats, **lons, **ssigo, **sangl, **sland, **slats, **slons, *array;

    if (argc != 3) {
      printf("\nUsage: %s 2005-05-16_49750_01.00800.hdr reduction_factor\n",argv[0]);
      printf("       where a reduction factor of 2 (3) yields 1/4 (1/8) the number of gridboxes\n\n");
      exit(1);
    }
    sscanf(argv[2],"%d",&reduce);
    if (reduce == 0) {
      printf("\nUsage: %s 2005-05-16_49750_01.00800.hdr reduction_factor\n",argv[0]);
      printf("       where the reduction factor should be greater than 0\n\n");
      exit(1);
    }

    strcpy(stema,argv[1]);  stema[25] = '\0';                                 /* create the file names */
    strcpy(stemb,argv[1]);  stemb[19] = '\0';
    strcpy(infila,argv[1]);
    strcpy(infilb,stema);
    strcat(infilb,".sar.nc");
    pixel_spacing = atoi(stema+20);
    if (pixel_spacing == 12 || pixel_spacing == 13)
      pixel_spacing = 12.5;

    strcpy(infila,argv[1]);                                                   /* read the header info */
    if ((fpa = fopen(infila,"r")) == NULL) {
      printf("couldn't open %s\n",argv[1]);
      exit(1);
    }
    for (a = 0; a < 27; a++) {
      fgets(line,LEN,fpa);
      if (a == 9)  sscanf(line,"%*s %*s %*s %s",date);
      if (a == 24) sscanf(line,"%*s %d",&nlines);
      if (a == 25) sscanf(line,"%*s %d",&npixels);
    }
    fclose(fpa);

    gridint = pow(2.0,reduce);                                                /* determine the new resol */
    subscenesize = 1;                                                         /* relative to the old, and */
    for (a = 1; a <= reduce; a++)                                             /* the size of the desired */
      subscenesize += 2.0 * pow(2.0,(double)a);                               /* smoothing matrix */
    aziblocks = (nlines - subscenesize) / gridint + 1;                        /* then get the desired */
    rngblocks = (npixels - subscenesize) / gridint + 1;                       /* number of subscenes */
    numlines_valid   = (aziblocks - 1) * gridint + subscenesize;              /* in azimuth and range */
    numpixels_valid  = (rngblocks - 1) * gridint + subscenesize;              /* and allocate some memory */
    startline_valid  = (nlines - numlines_valid) / 2 + 1;
    stopline_valid   = numlines_valid + startline_valid - 1;
    startpixel_valid = (npixels - numpixels_valid) / 2 + 1;
    stopixel_valid   = numpixels_valid + startpixel_valid - 1;
    nlines2skipbeg   = startline_valid - 1;
    nlines2skipend   = nlines - stopline_valid;
    npixels2skipbeg  = startpixel_valid - 1;
    npixels2skipend  = npixels - stopixel_valid;
    printf("\nreduce power      %d\n",reduce);
    printf("new grid int      %d\n",gridint);
    printf("subscenesize      %d\n",subscenesize);
    printf("aziblocks         %d\n",aziblocks);
    printf("rngblocks         %d\n",rngblocks);
    printf("nlines2skipbeg   %10d   npixels2skipbeg  %10d\n",nlines2skipbeg,npixels2skipbeg);
    printf("startline_valid  %10d   startpixel_valid %10d\n",startline_valid,startpixel_valid);
    printf("numlines_valid   %10d   numpixels_valid  %10d\n",numlines_valid,numpixels_valid);
    printf("stopline_valid   %10d   stopixel_valid   %10d\n",stopline_valid,stopixel_valid);
    printf("nlines2skipend   %10d   npixels2skipend  %10d\n",nlines2skipend,npixels2skipend);
    printf("nlines           %10d   npixels          %10d\n\n",nlines,npixels);

    smoo   = alloc_dblmatrix(subscenesize,subscenesize);
    sigo   = alloc_flmatrix(nlines,npixels);
    angl   = alloc_flmatrix(nlines,npixels);
    land   = alloc_flmatrix(nlines,npixels);
    lats   = alloc_flmatrix(nlines,npixels);
    lons   = alloc_flmatrix(nlines,npixels);
    ssigo  = alloc_flmatrix(aziblocks,rngblocks);
    sangl  = alloc_flmatrix(aziblocks,rngblocks);
    sland  = alloc_flmatrix(aziblocks,rngblocks);
    slats  = alloc_flmatrix(aziblocks,rngblocks);
    slons  = alloc_flmatrix(aziblocks,rngblocks);
    array  = (float *)calloc(nlines*npixels,sizeof(float));

    sprintf(infilc,                                                           /* next read the */
      "/home/ricani/prog/format.ccrs/include/conv.kernel.%d",reduce);         /* convolution kernel */
    if ((fpb = fopen(infilc,"r")) == NULL) {                                  /* (and normalize it) */
      printf("ERROR opening %s\n",infilc);
      exit(1);
    }
    printf("reading %s\n\n",infilc);
    count = 0;
    for (a = 0; a < subscenesize; a++)
      for (b = 0; b < subscenesize; b++) {
        fscanf(fpb,"%lf",&smoo[a][b]);
        count += smoo[a][b];
      }
    fclose(fpb);
/*  for (a = 0; a < subscenesize; a++)
      for (b = 0; b < subscenesize; b++)
        smoo[a][b] /= count;  */

    sarread(infilb, "sigo", date, array, nlines*npixels);
    d = 0;  for (a = 0; a < nlines; a++)  for (b = 0; b < npixels; b++) {sigo[a][b] = array[d];  d++;}
    sarread(infilb, "angl", date, array, nlines*npixels);
    d = 0;  for (a = 0; a < nlines; a++)  for (b = 0; b < npixels; b++) {angl[a][b] = array[d];  d++;}
    sarread(infilb, "land", date, array, nlines*npixels);
    d = 0;  for (a = 0; a < nlines; a++)  for (b = 0; b < npixels; b++) {land[a][b] = array[d];  d++;}
    sarread(infilb, "lats", date, array, nlines*npixels);
    d = 0;  for (a = 0; a < nlines; a++)  for (b = 0; b < npixels; b++) {lats[a][b] = array[d];  d++;}
    sarread(infilb, "lons", date, array, nlines*npixels);
    d = 0;  for (a = 0; a < nlines; a++)  for (b = 0; b < npixels; b++) {lons[a][b] = array[d];  d++;}
    printf("\n");

    for (a = 0; a < aziblocks; a++)                                           /* calculate smoothed values */
      for (b = 0; b < rngblocks; b++) {                                       /* for each range and azimuth */
        sum = smoosum = 0;                                                    /* set sigo to missing if */
        midline  = nlines2skipbeg  + (subscenesize - 1) / 2 + gridint * a;    /* central value is undefined */
        midpixel = npixels2skipbeg + (subscenesize - 1) / 2 + gridint * b;
        if (sigo[midline][midpixel] > DELMISS) {
          for (c = 0; c < subscenesize; c++)
            for (d = 0; d < subscenesize; d++) {
              tmpval = sigo[c+midline-(subscenesize-1)/2][d+midpixel-(subscenesize-1)/2];
              if (tmpval > DELMISS) {
                sum += smoo[c][d] * pow(10.0,tmpval/10.0);
                smoosum += smoo[c][d];
              }
            }
        }
        if (sum <= 0)
          ssigo[a][b] = MISSING;
        else
          ssigo[a][b] = 10.0 * log10(sum/smoosum);
      }

    for (a = 0; a < aziblocks; a++)                                           /* loop through the data */
      for (b = 0; b < rngblocks; b++) {                                       /* again to get locations */
        midline  = nlines2skipbeg  + (subscenesize - 1) / 2 + gridint * a;    /* of the subsampled data */
        midpixel = npixels2skipbeg + (subscenesize - 1) / 2 + gridint * b;
        sland[a][b] = land[midline][midpixel];
        slats[a][b] = lats[midline][midpixel];
        slons[a][b] = lons[midline][midpixel];
        sangl[a][b] = angl[midline][midpixel];
      }

    resol = pixel_spacing * pow(2.0,reduce);
    if (resol >= 10000)                                                       /* create an output file */
      sprintf(resolution,"%.0f",(float)resol);                                /* and store the fields */
    else if (resol >= 1000)                                                   /* correctly oriented */
      sprintf(resolution,"0%.0f",(float)resol);
    else if (resol >= 100)
      sprintf(resolution,"00%.0f",(float)resol);
    else if (resol >= 10)
      sprintf(resolution,"000%.0f",(float)resol);
    else if (resol >= 1)
      sprintf(resolution,"0000%.0f",(float)resol);
    else strcpy(resolution,"00000");
    sprintf(outfile,"%s.%s.sar.nc",stemb,resolution);
    sprintf(command,"/home/ricani/bin/nc.template.sar %s %s %d %d\n",
      outfile,date,aziblocks,rngblocks);
    printf("%s", command);
    system(command);
    printf("\n");

    d = 0;  for (a = 0; a < aziblocks; a++)  for (b = 0; b < rngblocks; b++) {array[d] = ssigo[a][b];  d++;}
    sarwrite(outfile, "sigo", date, array, aziblocks*rngblocks);
    d = 0;  for (a = 0; a < aziblocks; a++)  for (b = 0; b < rngblocks; b++) {array[d] = sangl[a][b];  d++;}
    sarwrite(outfile, "angl", date, array, aziblocks*rngblocks);
    d = 0;  for (a = 0; a < aziblocks; a++)  for (b = 0; b < rngblocks; b++) {array[d] = sland[a][b];  d++;}
    sarwrite(outfile, "land", date, array, aziblocks*rngblocks);
    d = 0;  for (a = 0; a < aziblocks; a++)  for (b = 0; b < rngblocks; b++) {array[d] = slats[a][b];  d++;}
    sarwrite(outfile, "lats", date, array, aziblocks*rngblocks);
    d = 0;  for (a = 0; a < aziblocks; a++)  for (b = 0; b < rngblocks; b++) {array[d] = slons[a][b];  d++;}
    sarwrite(outfile, "lons", date, array, aziblocks*rngblocks);
    ullat = slats[0][0];                      ullon = slons[0][0];
    urlat = slats[0][rngblocks-1];            urlon = slons[0][rngblocks-1];
    lllat = slats[aziblocks-1][0];            lllon = slons[aziblocks-1][0];
    lrlat = slats[aziblocks-1][rngblocks-1];  lrlon = slons[aziblocks-1][rngblocks-1];

    sarread(infilb, "lat", date, array,  nlines);                             /* also transfer arbitrary netCDF */
    d = 0;  for (a = 0; a <  nlines; a++) {lats[a][0] = array[d];  d++;}      /* lats and lons (subsampled) in */
    sarread(infilb, "lon", date, array, npixels);                             /* order to maintain hdr mappings */
    d = 0;  for (b = 0; b < npixels; b++) {lons[0][b] = array[d];  d++;}
    for (a = 0; a < aziblocks; a++) {
      midline  = nlines2skipbeg  + (subscenesize - 1) / 2 + gridint * a;
      slats[a][0] = lats[midline][0];
    }
    for (b = 0; b < rngblocks; b++) {
      midpixel = npixels2skipbeg + (subscenesize - 1) / 2 + gridint * b;
      slons[0][b] = lons[0][midpixel];
    }
    d = 0;  for (a = 0; a < aziblocks; a++) {array[d] = slats[a][0];  d++;}
    sarwrite(outfile, "lat", date, array, aziblocks);
    d = 0;  for (b = 0; b < rngblocks; b++) {array[d] = slons[0][b];  d++;}
    sarwrite(outfile, "lon", date, array, rngblocks);
    printf("\n");

    free_dblmatrix(smoo,subscenesize,subscenesize);                           /* free allocated memory */
    free_flmatrix(sigo,nlines,npixels);
    free_flmatrix(angl,nlines,npixels);
    free_flmatrix(land,nlines,npixels);
    free_flmatrix(lats,nlines,npixels);
    free_flmatrix(lons,nlines,npixels);
    free_flmatrix(ssigo,aziblocks,rngblocks);
    free_flmatrix(sangl,aziblocks,rngblocks);
    free_flmatrix(sland,aziblocks,rngblocks);
    free_flmatrix(slats,aziblocks,rngblocks);
    free_flmatrix(slons,aziblocks,rngblocks);
    free(array);

    sprintf(outfile,"%s.%s.hdr",stemb,resolution);                            /* and store hdr metadata */
    printf("reading %s\nwriting %s\n\n",infila,outfile);
    if ((fpa = fopen(infila,"r")) == NULL || (fpc = fopen(outfile,"w")) == NULL) {
      printf("\n\n***Error opening %s or %s\n\n",infila,outfile);
      exit(1);
    }
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56s\n","ProductID",infilb);
    for (a = 0; a < 9; a++) {
      fgets(line,LEN,fpa);  fprintf(fpc,"%s",line);
    }
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56f\n","UpperLeft_Latitude",ullat);
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56f\n","UpperLeft_Longitude",ullon);
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56f\n","UpperRight_Latitude",urlat);
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56f\n","UpperRight_Longitude",urlon);
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56f\n","LowerLeft_Latitude",lllat);
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56f\n","LowerLeft_Longitude",lllon);
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56f\n","LowerRight_Latitude",lrlat);
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56f\n","LowerRight_Longitude",lrlon);
    for (a = 0; a < 2; a++) {
      fgets(line,LEN,fpa);  fprintf(fpc,"%s",line);
    }
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56d\n","StartLineValid",startline_valid);
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56d\n","StopLineValid",stopline_valid);
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56d\n","StartPixelValid",startpixel_valid);
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56d\n","StopPixelValid",stopixel_valid);
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56d\n","NumAzimuthSubScenes",aziblocks);
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56d\n","NumRangeSubScenes",rngblocks);
    fgets(line,LEN,fpa);  fprintf(fpc,"%-23s%-56d\n","SubSceneSize",subscenesize);
    for (a = 0; a < 4; a++) {
      fgets(line,LEN,fpa);  fprintf(fpc,"%s",line);
    }
    fclose(fpa);
    fclose(fpc);
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
