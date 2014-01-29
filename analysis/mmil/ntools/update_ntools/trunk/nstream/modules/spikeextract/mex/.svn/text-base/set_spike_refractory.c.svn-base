#include <stdio.h>
#include <stdlib.h>
/*----------------------------------------------------------------------*/
/* #include <sys/types.h> */
#include <sys/shm.h>
#include <sys/ipc.h>
/* #include <sys/sem.h> */
/*----------------------------------------------------------------------*/
#include "mex.h"
/* #include "mysem.h" */
#include "nstream.h"
/*----------------------------------------------------------------------*/
#define DEBUG1 printf("File: %s\tCompiled: %s %s\tLine: %d\n", __FILE__, __DATE__, __TIME__, __LINE__);
/*----------------------------------------------------------------------
%
% Y = SET_SPIKE_REFRACTORY()
%
% Sets refractory period for spike detection.
%
% SET_SPIKE_REFRACTORY(R) 
%
% Inputs:	 R = Refractory period in ms 
%
% Outputs:  None 
%
 ----------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *data;
  int shm_id = 0;

  /*------ Share Memory Init. Var ------*/
  void *shared_memory = (void *) 0;
  struct shared_state_vars *share;
  /*--------- End of Init. Var ---------*/

 
  /* Check for proper number of arguments. */
  if(nrhs!=1) { mexErrMsgTxt("One inputs required."); } 
  else if(nlhs != 0) { mexErrMsgTxt("No output arguments"); }
 
  int refractory_period_ms = mxGetScalar(prhs[0]); 
 
  /*====================== Get share memory ====================*/
  shm_id = shmget((key_t)NSTREAM_STATE_VARS_SHM_KEY, sizeof(struct shared_state_vars), 0666);
	
  /* Check share memory */
  if(shm_id == -1) { perror("shmget failed: "); exit(EXIT_FAILURE); }
	
  /* Make shared memory accessible by program */
  shared_memory = (struct shared_use_st*)shmat(shm_id,(void *)0, 0);
  if(shared_memory == (void *)-1)
  { 
    fprintf(stderr, "shmat failed in %s\n", __FILE__);
    data[0] = -1;
    return; 
  } 
	
  /*  Make pointer to the share memory */
  share = (struct shared_state_vars *)shared_memory;

  share->spike_detect_refractory_period = refractory_period_ms * 30;        
  if (shmdt(share) == -1) { fprintf(stderr, "shmdt function failed\n"); exit(EXIT_FAILURE); }
  return;
}
/*----------------------------------------------------------------------*/
