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
% Y = NSTREAM_GETRAW_NSPIKE_ALL()
%
% NSTREAM_GETRAW_NSPIKE_ALL returns all data from the nspike shm buffer 
%
% Inputs:   None	
%
% Outputs:  Y  =  [M,N] Matrix, nspike data 
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
  if(nrhs!=0) { mexErrMsgTxt("No inputs required."); } 
  else if(nlhs > 2) { mexErrMsgTxt("Too many output arguments"); }
  
  /* Create matrix for the return argument. [mxREAL|mxCOMPLEX] */
  plhs[0] = mxCreateDoubleMatrix(NSTREAM_NUM_NSPIKE_CHANNELS,NSTREAM_DURATION * NSTREAM_SR, mxREAL);
  raw_matrix = mxGetPr(plhs[0]); 


  for(i=0; i < NSTREAM_NSPIKE_NUM_SAMPLES; i++)
  {  
     raw_matrix[i] = share->nspike_buffer[i];
  }  

  if (shmdt(share) == -1) { fprintf(stderr, "shmdt function failed\n"); exit(EXIT_FAILURE); } 
}  

/*----------------------------------------------------------------------*/
