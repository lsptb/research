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
%
% Y = NSTREAM_GETRAW_NSPIKE_DECIMATE(N,D)
%
% NSTREAM_GETRAW_NSPIKE_DECIMATE returns N data points back from last 
%                                written position in shm decimated by 
%                                a factor of D
%
% Inputs:	N  =  Scalar, number of data points to return 
%           D  =  Scalar, decimation factor
%
% Outputs:  Y  =  [M,N] Matrix, channel, data
%
*/
/*----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*------ Semaphore Init. Var ---------*/
  int i, j, read_pos;
  int shm_id = 0;

  /*------ Share Memory Init. Var ------*/
  int running = 1;
  void *shared_memory = (void *) 0;
  struct shared_state_vars *share;
  /*--------- End of Init. Var ---------*/
  
  double N;
  double *raw_matrix;

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
  if(nrhs!=2) { mexErrMsgTxt("Two inputs required."); } 
  else if(nlhs > 2) { mexErrMsgTxt("Too many output arguments"); }
  
  int number_of_samples = mxGetScalar(prhs[0]);
  int decimation_factor = mxGetScalar(prhs[1]);

  if ((NSTREAM_SR % decimation_factor) != 0)
  {
      mexErrMsgTxt("Decimation factor must be an integer factor of 30000");
  }
  int decimated_number_of_samples = number_of_samples/decimation_factor;
 
  /* Create matrix for the return argument. [mxREAL|mxCOMPLEX] */
  plhs[0] = mxCreateDoubleMatrix(NSTREAM_NUM_NSPIKE_CHANNELS,decimated_number_of_samples, mxREAL);
  raw_matrix = mxGetPr(plhs[0]); 

  /*========================== Critical Section =======================*/

  unsigned long time = share->nspike_time;

  /*writing backwards into the matrix*/
  int write_pos = (time*NSTREAM_NUM_NSPIKE_CHANNELS)%(NSTREAM_NSPIKE_NUM_SAMPLES); 
  if (write_pos == 0) { write_pos = NSTREAM_NSPIKE_NUM_SAMPLES; }
  write_pos--; 
 
  for(i=0; i < decimated_number_of_samples*NSTREAM_NUM_NSPIKE_CHANNELS; i++)
  { 
     raw_matrix[NSTREAM_NUM_NSPIKE_CHANNELS*decimated_number_of_samples - 1 - i] = share->nspike_buffer[write_pos];
     write_pos--;
 
     if (((i+1)%NSTREAM_NUM_NSPIKE_CHANNELS) == 0)
     {
        write_pos -= NSTREAM_NUM_NSPIKE_CHANNELS*(decimation_factor-1);
     }
	 
    if(write_pos < 0) {
        printf("WRAP!\n");
	    write_pos = NSTREAM_NSPIKE_NUM_SAMPLES - (write_pos*-1); 
     }
  }  
  /*===================== End of Critical Section ======================*/

  if (shmdt(share) == -1) { fprintf(stderr, "shmdt function failed\n"); exit(EXIT_FAILURE); } 
}  

/*----------------------------------------------------------------------*/
