#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <unistd.h>
#include <ctype.h>
#include <fcntl.h>
/*----------------------------------------------------------------------*/
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/stat.h>
#include <sys/shm.h>
/*----------------------------------------------------------------------*/
#define DEBUG printf("File: %s\tCompiled: %s %s\tLine: %d\n", __FILE__, __DATE__, __TIME__, __LINE__);
/*----------------------------------------------------------------------*/
#define OUTPUT_FILE "output.dat"
/*----------------------------------------------------------------------
 * %
 * % NTOOLS_DECIMATE(INPUT_FILENAME, OUTPUT_FILENAME, DECIMATION)
 * %
 * % NTOOLS_DECIMATE decimates a data file on disk and saves to another file
 * %
 * % Inputs:        DECIMATION      =  Scalar, Decimation factor 
 * %
 * % Outputs: None 
 * %
 *  ----------------------------------------------------------------------*/
int main(int argc, char *argv[])
{

  FILE *output_fp;
  short int *buf;
  long int n;
  long int data_len;
  int i,k;

  int CH = 256;
  int SS = 2;
  int DECIMATION_FACTOR = 15;
  long int NUM_SAMPLES = 10;

  /*------ Parse input options ---------------------*/
  while ((i=getopt(argc,argv,"c:d:h")) !=EOF)
    {
      switch(i)
        {
          case 'd': DECIMATION_FACTOR = atoi(optarg); break;
          case 'c': CH = atoi(optarg); break;
          case 'h':
          default:
	  exit(1);
        } /* end switch statement */
    }  /* end while statement */

printf("Factor is %d\n", DECIMATION_FACTOR);
printf("Num channels is %d\n", CH);

  output_fp = fopen(OUTPUT_FILE, "w");
  if(output_fp == NULL) { fprintf(stderr, "Open %s failure\n", OUTPUT_FILE);  exit(EXIT_FAILURE);}

  data_len = CH*DECIMATION_FACTOR*NUM_SAMPLES;

  buf=malloc(data_len);

  while ((k=read(0, buf, data_len)) != 0) 
  {
    for (n = 0; n < k/SS; n+=CH*DECIMATION_FACTOR)
    {
      fwrite(&buf[n],SS,CH,output_fp);
    }
  }

  fclose(output_fp);
  return EXIT_SUCCESS;

}
