/*---------------------------------------------------------------------------
 * broker_child.c - Child process that is forked from broker. This process
 * will read stdin in this case pipe from broker and calculate status and
 * time and write it into share memory circula buffer.
 * IPC communications used 1.) Semaphores 2.) Share Memory 3.) Circular Buffer
 *---------------------------------------------------------------------------*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/*---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
/*---------------------------------------------------------------------------*/
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include "myshare.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{

  int n;
  int i = 0;
  FILE *raw1_fd, *raw2_fd, *raw3_fd, *raw4_fd;

  /*------- Data Buffer ------------*/
  char *buf, *raw1, *raw2, *raw3, *raw4;
  long int  k;
  long int DD = 100000;
  long int   data_len = DD*256*2;
 
  /*------ set output filenames --------------------

  /*-----  Initialization -------------------------*/
  
  /*------ intialize eye position output file ------*/
  /*
  raw1_fd = creat(RAW1_FILE, 0666);
  if(raw1_fd == -1) { fprintf(stderr, "Open %s failure\n", RAW1_FILE);  exit(EXIT_FAILURE);}

  raw2_fd = creat(RAW2_FILE, 0666);
  if(raw2_fd == -1) { fprintf(stderr, "Open %s failure\n", RAW2_FILE);  exit(EXIT_FAILURE);}
      
  raw3_fd = creat(RAW3_FILE, 0666);
  if(raw3_fd == -1) { fprintf(stderr, "Open %s failure\n", RAW3_FILE);  exit(EXIT_FAILURE);}
      
  raw4_fd = creat(RAW4_FILE, 0666);
  if(raw4_fd == -1) { fprintf(stderr, "Open %s failure\n", RAW4_FILE);  exit(EXIT_FAILURE);}
  */
  raw1_fd = fopen(RAW1_FILE, "w");
  raw2_fd = fopen(RAW2_FILE, "w");
  raw3_fd = fopen(RAW3_FILE, "w");
  raw4_fd = fopen(RAW4_FILE, "w");
    
  buf=malloc(data_len);
  raw1 = malloc(60*DD*2); 
  raw2 = malloc(60*DD*2);
  raw3 = malloc(60*DD*2); 
  raw4 = malloc(60*DD*2);
  while ((k=read(0, buf, data_len)) != 0) 
  {
      //printf("%ld\n",k);
  /*memset(raw1,0,sizeof(raw1));
  memset(raw2,0,sizeof(raw1));
  memset(raw3,0,sizeof(raw1));
  memset(raw4,0,sizeof(raw1)); */
      for (n = 0; n < k; n += 256*SS)
      {
        memcpy(raw1 + (n/256)*60,buf + n,60*SS); // copy 1 sample from 60 channels to first raw buffer
        memcpy(raw2 + (n/256)*60,buf+(n+64*SS),60*SS); // 
        memcpy(raw3 + (n/256)*60,buf+(n+128*SS),60*SS); //
        memcpy(raw4 + (n/256)*60,buf+(n+192*SS),60*SS); //
      }
      fwrite(raw1,1,(k/256)*60,raw1_fd);  // write DD samples from first buffer to disk
      fwrite(raw2,1,(k/256)*60,raw2_fd);
      fwrite(raw3,1,(k/256)*60,raw3_fd);
      fwrite(raw4,1,(k/256)*60,raw4_fd);

  }  /* ======================== END OF WHILE LOOP ==============================*/
  fclose(raw1_fd);
  fclose(raw2_fd);
  fclose(raw3_fd);
  fclose(raw4_fd);
  printf("%d %d\n",i,sizeof(raw1));

  return EXIT_SUCCESS;
}
/*------------------------------------------------------------*/






