#include <stdio.h>
#include <stdlib.h>
/*----------------------------------------------------------------------*/
/* #include <sys/types.h> */
#include <sys/shm.h>
#include <sys/ipc.h>
/* #include <sys/sem.h> */
/*----------------------------------------------------------------------*/
#include "mex.h"
#include "nstream.h"
/*----------------------------------------------------------------------*/
#define DEBUG1 printf("File: %s\tCompiled: %s %s\tLine: %d\n", __FILE__, __DATE__, __TIME__, __LINE__);
/*----------------------------------------------------------------------
%
% Y = NSTREAM_GETTIME_COMEDI
%
% NSTREAM_GETTIME_COMEDI returns time in _milliseconds_ from comedi shm buffer
%
% Inputs:	 No inputs
%
% Outputs:  Y  =  scalar, time data
*/
/*----------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *data;
  int shm_id = 0;

  /*------ Share Memory Init. Var ------*/
  void *shared_memory = (void *) 0;
  struct shared_state_vars *share;
  /*--------- End of Init. Var ---------*/

  /* Create matrix for the return argument. [mxREAL|mxCOMPLEX] */
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);  
  data = mxGetPr(plhs[0]);

  /*====================== Get share memory ====================*/
  shm_id = shmget((key_t)NSTREAM_STATE_VARS_SHM_KEY, sizeof(struct shared_state_vars), 0666);
	
  /* Check share memory */
  if(shm_id == -1) { perror("shmget failed: "); exit(EXIT_FAILURE); }
	
  /* Make shared memory accessible by program */
  shared_memory = (struct shared_state_vars*)shmat(shm_id,(void *)0, 0);
  if(shared_memory == (void *)-1)
  { 
    fprintf(stderr, "shmat failed in %s\n", __FILE__);
    data[0] = -1;
    return; 
  } 
	
  /*  Make pointer to the share memory */
  share = (struct shared_state_vars *)shared_memory;
        
  /*------------------- Start of Critical Section ---------------------*/
  data[0] = share->comedi_time / 30;
  /*------------------- End of Critical Section ---------------------*/

  if (shmdt(share) == -1) { fprintf(stderr, "shmdt function failed\n"); exit(EXIT_FAILURE); }
  return;
}
/*----------------------------------------------------------------------*/
