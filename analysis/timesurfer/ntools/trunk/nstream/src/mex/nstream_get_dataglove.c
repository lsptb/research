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
% Y = NSTREAM_GET_DATAGLOVE(N)
%
% NSTREAM_GET_DATAGLOVE returns N dataglove event points back from last written position 
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
  
  int number_of_samples = mxGetScalar(prhs[0]);
  
  /* Create matrix for the return argument. [mxREAL|mxCOMPLEX] */
  plhs[0] = mxCreateDoubleMatrix(1,number_of_samples, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(6,number_of_samples, mxREAL);
  timestamp = mxGetPr(plhs[0]); 
  data = mxGetPr(plhs[1]);
  /*========================== Critical Section =======================*/

  write_pos = share->dataglove_buffer_index - 1;
  if(write_pos < 0) {
    write_pos = NSTREAM_DATAGLOVE_NUM_SAMPLES + write_pos;
  }
  
  /*writing backwards into the matrix*/
  for(i=0; i < number_of_samples; i++)
  { 
     timestamp[number_of_samples - 1 - i] = share->dataglove_buffer[write_pos].timestamp/30;
     data[6*number_of_samples - 1 - 6*i] = share->dataglove_buffer[write_pos].gesture;
     data[6*number_of_samples - 2 - 6*i] = share->dataglove_buffer[write_pos].little;
     data[6*number_of_samples - 3 - 6*i] = share->dataglove_buffer[write_pos].ring;
     data[6*number_of_samples - 4 - 6*i] = share->dataglove_buffer[write_pos].middle;
     data[6*number_of_samples - 5 - 6*i] = share->dataglove_buffer[write_pos].index;
     data[6*number_of_samples - 6 - 6*i] = share->dataglove_buffer[write_pos].thumb;
     write_pos--;
	 if(write_pos < 0) { 	 
	   write_pos = NSTREAM_DATAGLOVE_NUM_SAMPLES + write_pos; 
     }
  }  
  /*===================== End of Critical Section ======================*/

  if (shmdt(share) == -1) { fprintf(stderr, "shmdt function failed\n"); exit(EXIT_FAILURE); } 
}  

/*----------------------------------------------------------------------*/
