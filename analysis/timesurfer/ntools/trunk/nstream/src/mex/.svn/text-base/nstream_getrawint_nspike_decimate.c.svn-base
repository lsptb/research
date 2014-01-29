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
% Y = NSTREAM_GETRAWINT_NSPIKE_DECIMATE(T1,T2,D)
%
% NSTREAM_GETRAWINT_NSPIKE_DECIMATE returns data for a given time interval 
% decimated by a factor of D from the shared memory buffer
%
% Inputs:	 T1	=  Time in ms to start from 
%       	 T2	=  Time in ms to end at 
%             D =  Decimation Factor (MUST BE INTEGER FACTOR OF 30)
%
% Outputs:  Y  =  [M,N] Matrix, channel, data
%
*/ 
/*----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*------ Semaphore Init. Var ---------*/
  int i;
  int shm_id = 0;

  /*------ Share Memory Init. Var ------*/
  int running = 1;
  void *shared_memory = (void *) 0;
  struct shared_state_vars *share;
  /*--------- End of Init. Var ---------*/
  
  double *return_matrix;
  /* int ROW=RAW_OUTER, COL, n, array_row; */

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
  if(nrhs!=3) { mexErrMsgTxt("Three inputs required."); } 
  else if(nlhs > 2) { mexErrMsgTxt("Too many output arguments"); }
  
  int from_ms = mxGetScalar(prhs[0]);
  int to_ms = mxGetScalar(prhs[1]);
  int decimation_factor = mxGetScalar(prhs[2]);

  if ((NSTREAM_SR % decimation_factor) != 0)
  {
      mexErrMsgTxt("Decimation factor must be an integer factor of 30000");
  }

  if (from_ms > to_ms)
      mexErrMsgTxt("T1 must be < T2");

  int dur_ms = to_ms - from_ms;
 
  /* convert ms (1k) to samples (30k) */ 
  unsigned long dur_samp = (dur_ms * 30)/decimation_factor; 
  unsigned long from_samp = from_ms * 30;
  unsigned long to_samp = to_ms * 30;
  unsigned long nspike_buf_time = share->nspike_time;

  /* verify that the range requested is available in the buffer */
  if (nspike_buf_time < to_samp)
      mexErrMsgTxt("T2 > nspike_time");

  if (from_samp < nspike_buf_time - (NSTREAM_DURATION * NSTREAM_SR))
    if (nspike_buf_time > nspike_buf_time - (NSTREAM_DURATION * NSTREAM_SR))
      mexErrMsgTxt("T1 < nspike_time - buffer length in samples");

  /* ok, we're lookin good.  create the return matrix and do the copy! */
  
  /* Create matrix for the return argument. [mxREAL|mxCOMPLEX] */
  /* plhs[0] = mxCreateDoubleMatrix(NSTREAM_NUM_NSPIKE_CHANNELS, num_samples_to_get, mxREAL); */
  plhs[0] = mxCreateDoubleMatrix(NSTREAM_NUM_NSPIKE_CHANNELS, dur_samp, mxREAL);
  return_matrix = mxGetPr(plhs[0]); 

  /*========================== Critical Section =======================*/

  /*writing backwards into the matrix*/
  int write_pos = ((to_samp*NSTREAM_NUM_NSPIKE_CHANNELS)%(NSTREAM_NSPIKE_NUM_SAMPLES));
  if (write_pos == 0) { write_pos = NSTREAM_NSPIKE_NUM_SAMPLES; }
  write_pos--;

  
  for(i=0; i < dur_samp*NSTREAM_NUM_NSPIKE_CHANNELS; i++)
  {  
     return_matrix[NSTREAM_NUM_NSPIKE_CHANNELS*dur_samp - 1 - i] = share->nspike_buffer[write_pos];
     write_pos--;

     if (((i+1)%NSTREAM_NUM_NSPIKE_CHANNELS) == 0)
     {
        write_pos -= NSTREAM_NUM_NSPIKE_CHANNELS*(decimation_factor-1);
     }
 	 
	 if(write_pos < 0) {
	    write_pos = NSTREAM_NSPIKE_NUM_SAMPLES - (write_pos*-1); 
     }
  }  
  /*===================== End of Critical Section ======================*/

  if (shmdt(share) == -1) { fprintf(stderr, "shmdt function failed\n"); exit(EXIT_FAILURE); } 
}  

/*----------------------------------------------------------------------*/
