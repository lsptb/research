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
  long int ch1[8], ch2[8]; 
  FILE *raw1_fd, *raw2_fd;

  /*------- Data Buffer ------------*/
  char *buf, *raw1, *raw2;
  int chan_count = 128;
  long int  k, ich;
  long int DD = 100000;
  long int   data_len = DD*chan_count*2;
 
  /*------ set output filenames --------------------

  /*-----  Initialization -------------------------*/
  for (ich = 0; ich < 8; ich++)
    {
      ch1[ich] = 15-2*ich;
      ch2[ich] = 16+2*ich;
    }
  
  /*------ intialize eye position output file ------*/
  raw1_fd = fopen(RAW1_FILE, "w");
  raw2_fd = fopen(RAW2_FILE, "w");
    
  buf=malloc(data_len);
  raw1 = malloc(8*DD*2); 
  raw2 = malloc(8*DD*2);
  while ((k=read(0, buf, data_len)) != 0) 
  {
      //printf("%ld\n",k);
      for (n = 0; n < k; n += chan_count*SS)
      {
        for (ich = 0; ich < 8; ich++)
          {
          memcpy(raw1 + (n/chan_count)*8+ich*SS,buf+n+ch1[ich]*SS,SS); // copy sample from raw buffer
          memcpy(raw2 + (n/chan_count)*8+ich*SS,buf+n+ch2[ich]*SS,SS); // 
          }
      }
      fwrite(raw1,1,(k/chan_count)*8,raw1_fd);  // write DD samples from first buffer to disk
      fwrite(raw2,1,(k/chan_count)*8,raw2_fd);

  }  /* ======================== END OF WHILE LOOP ==============================*/
  fclose(raw1_fd);
  fclose(raw2_fd);
  printf("%d %d\n",i,sizeof(raw1));

  return EXIT_SUCCESS;
}
/*------------------------------------------------------------*/

