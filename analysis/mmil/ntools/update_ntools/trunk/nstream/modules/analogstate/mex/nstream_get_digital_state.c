#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/*----------------------------------------------------------------------*/
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>
/*----------------------------------------------------------------------*/
#include "nstream.h"
#include "mex.h"
/*----------------------------------------------------------------------*/
#define DEBUG printf("File: %s\tCompiled: %s %s\tLine: %d\n", __FILE__, __DATE__, __TIME__, __LINE__);
/*----------------------------------------------------------------------
%
% Y = NSTREAM_GET_DIGITAL_STATE(N)
%
% NSTREAM_GET_DIGITAL_STATE returns state data from shared memory buffer
%
% Inputs:	 N	=  Scalar, Length of state data to return
%
% Outputs:  Y  =  [3,N] Matrix, state data returned
%
 ----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*------ Semaphore Init. Var ---------*/
  int i, j, read_pos;
  
  /*------ Share Memory Init. Var ------*/
  int running = 1;
  void *shared_memory = (void *) 0;
  struct shared_state_vars *share;
  int shm_id = 0;

  /*--------- End of Init. Var ---------*/
  
  double N;
  double *state_matrix;
  int ROW=3, COL, n, array_row;

  /*================ Get semaphore and share memory ======================*/
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
  
  COL = mxGetScalar(prhs[0]);

  
  /* Create matrix for the return argument. [mxREAL|mxCOMPLEX] */
  plhs[0] = mxCreateDoubleMatrix(ROW, COL, mxREAL); 
  state_matrix = mxGetPr(plhs[0]); 

  /*========================== Critical Section =======================*/

  read_pos = share->analog_state_index-1;
  if(read_pos == -1) { read_pos = NSTREAM_ANALOG_STATE_NUM_SAMPLES-1; } 

  /*writing backwards into the matrix*/
  for(i=ROW*COL-1; i > -1; i--)
  { 
    /* if share memory is at begining just reset it to end */
    state_matrix[i] = share->analog_state_buffer[i % ROW][read_pos];
    if (i%ROW==0)
      { --read_pos; if(read_pos == -1) { read_pos = NSTREAM_ANALOG_STATE_NUM_SAMPLES-1; }}
  }  
  /*===================== End of Critical Section ======================*/

  if (shmdt(share) == -1) { fprintf(stderr, "shmdt function failed\n"); exit(EXIT_FAILURE); } 
}  

/*----------------------------------------------------------------------*/
