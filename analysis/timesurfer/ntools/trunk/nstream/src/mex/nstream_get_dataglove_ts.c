#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/shm.h>
/*----------------------------------------------------------------------*/
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
/*----------------------------------------------------------------------*/
#include "mex.h"
#include "nstream.h"
/*----------------------------------------------------------------------*/
#define DEBUG printf("File: %s\tCompiled: %s %s\tLine: %d\n", __FILE__, __DATE__, __TIME__, __LINE__);
/*----------------------------------------------------------------------
% Y = NSTREAM_GET_DATAGLOVE_TS(TS)
%
% NSTREAM_GET_DATAGLOVE_TS returns all dataglove samples > TS
%
% Inputs:	 N	=  Scalar, number of data points to return 
%
% Outputs:   Y  =  [M,N] Matrix
% M = timestamp
% N = [thumb index middle ring little gesture]
%
----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*------ Semaphore Init. Var ---------*/
  int i, j, read_pos;
  int shm_id = 0;
  int write_pos = 0;

  /*------ Share Memory Init. Var ------*/
  int running = 1;
  void *shared_memory = (void *) 0;
  struct shared_state_vars *share;
  /*--------- End of Init. Var ---------*/
  
  double N;
  double *data, *timestamp;

  shm_id = shmget((key_t)NSTREAM_STATE_VARS_SHM_KEY, sizeof(struct shared_state_vars), 0666);
	
  /* Check share memory */
  if(shm_id == -1) { perror("shmget failed: "); exit(EXIT_FAILURE); }
	
  /* Make shared memory accessible by program */
  shared_memory = (struct shared_state_vars*)shmat(shm_id,(void *)0, 0);
  if(shared_memory == (void *)-1)
  { 
    fprintf(stderr, "shmat failed in %s\n", __FILE__);
    return; 
  } 
	
  /*  Make pointer to the share memory */
  share = (struct shared_state_vars *)shared_memory;

  /*====================================================================*/

  /* Check for proper number of arguments. */
  if(nrhs!=1) { mexErrMsgTxt("One input required."); } 
  else if(nlhs > 2) { mexErrMsgTxt("Too many output arguments"); }
  
  int start_ts = mxGetScalar(prhs[0]);
  int number_of_samples = 0;

  /* first determine how many events we have to return and bail out if zero */
  write_pos = share->dataglove_buffer_index - 1;
  if(write_pos < 0) {
    write_pos = NSTREAM_DATAGLOVE_NUM_SAMPLES + write_pos;
  }

  for (number_of_samples = 0; number_of_samples < NSTREAM_DATAGLOVE_NUM_SAMPLES; number_of_samples++)
  {
    if (share->dataglove_buffer[write_pos].timestamp/30 <= start_ts)
      break;
    write_pos--;
	if(write_pos < 0) { 	 
	  write_pos = NSTREAM_DATAGLOVE_NUM_SAMPLES + write_pos; 
    }
  }
 
  /* Create matrix for the return argument. [mxREAL|mxCOMPLEX] */
  plhs[0] = mxCreateDoubleMatrix(7,number_of_samples, mxREAL);
  data = mxGetPr(plhs[0]);
  /*========================== Critical Section =======================*/

  write_pos = share->dataglove_buffer_index - 1;
  if(write_pos < 0) {
    write_pos = NSTREAM_DATAGLOVE_NUM_SAMPLES + write_pos;
  }
  
  /*writing backwards into the matrix*/
  for(i=0; i < number_of_samples; i++)
  { 
     data[7*number_of_samples - 1 - 7*i] = share->dataglove_buffer[write_pos].gesture;
     data[7*number_of_samples - 2 - 7*i] = share->dataglove_buffer[write_pos].little;
     data[7*number_of_samples - 3 - 7*i] = share->dataglove_buffer[write_pos].ring;
     data[7*number_of_samples - 4 - 7*i] = share->dataglove_buffer[write_pos].middle;
     data[7*number_of_samples - 5 - 7*i] = share->dataglove_buffer[write_pos].index;
     data[7*number_of_samples - 6 - 7*i] = share->dataglove_buffer[write_pos].thumb;
     data[7*number_of_samples - 7 - 7*i] = share->dataglove_buffer[write_pos].timestamp/30;
     write_pos--;
	 if(write_pos < 0) { 	 
	   write_pos = NSTREAM_DATAGLOVE_NUM_SAMPLES + write_pos; 
     }
  }  
  /*===================== End of Critical Section ======================*/

  if (shmdt(share) == -1) { fprintf(stderr, "shmdt function failed\n"); exit(EXIT_FAILURE); } 
}  

/*----------------------------------------------------------------------*/
