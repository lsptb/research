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
% Y = NSTREAM_GETMU_NSPIKE_16(N)
%
% NSTREAM_GETMU_NSPIKE_16 returns N data points back from last written position in shm 
%
% Inputs:	 N	=  Scalar, number of data points to return 
%
% Outputs:  Y  =  [M,N] Matrix, channel, data
%
*/
/*----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*------ Semaphore Init. Var ---------*/
  int i, j, read_pos, ich;
  int shm_id = 0;
  long int ch[16]; 
  int output_pos = 0;

  /*------ Share Memory Init. Var ------*/
  int running = 1;
  void *shared_memory = (void *) 0;
  struct shared_state_vars *share;

  for (ich = 0; ich < 8; ich++)
    {
      ch[ich] = 15-2*ich+128;
      ch[ich+8] = 16+2*ich+128;
    }
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
  if(nrhs!=1) { mexErrMsgTxt("One input required."); } 
  else if(nlhs > 2) { mexErrMsgTxt("Too many output arguments"); }
  
  int number_of_samples = mxGetScalar(prhs[0]);
  
  /* Create matrix for the return argument. [mxREAL|mxCOMPLEX] */
  plhs[0] = mxCreateDoubleMatrix(16,number_of_samples, mxREAL);
  raw_matrix = mxGetPr(plhs[0]); 

  /*========================== Critical Section =======================*/


  unsigned long time = share->nspike_time; 

  /*writing backwards into the matrix*/

  int write_pos = (time*NSTREAM_NUM_NSPIKE_CHANNELS)%(NSTREAM_NSPIKE_NUM_SAMPLES); 
  if (write_pos == 0) { write_pos = NSTREAM_NSPIKE_NUM_SAMPLES; }
  write_pos-=NSTREAM_NUM_NSPIKE_CHANNELS;
  if (write_pos < 0) { write_pos += NSTREAM_NSPIKE_NUM_SAMPLES; }

/*  write_pos--; */

 

  for(i=0; i < number_of_samples*NSTREAM_NUM_NSPIKE_CHANNELS; i+=NSTREAM_NUM_NSPIKE_CHANNELS)
  {  
    for(ich = 0; ich < 16; ich++)
      {
       raw_matrix[16*(number_of_samples-1) - output_pos + ich ] = share->nspike_buffer[write_pos+ch[ich]];
      }
    output_pos+=16;
    write_pos -= NSTREAM_NUM_NSPIKE_CHANNELS;
    if(write_pos < 0) {
       write_pos += NSTREAM_NSPIKE_NUM_SAMPLES;
    }
  }  

  /*===================== End of Critical Section ======================*/

  if (shmdt(share) == -1) { fprintf(stderr, "shmdt function failed\n"); exit(EXIT_FAILURE); } 
}  

/*----------------------------------------------------------------------*/
