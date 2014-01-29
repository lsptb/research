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
% Y = GET_SPIKE_DATA(CH, CELL_ID, NPOINTS)
%
% GET_SPIKE_DATA returns N spike buffers from channel CH which match CELL_ID 
%
% Inputs:	 CH	     =  Channel 
%            N	     =  Scalar, number of data points to return 
%            CELL_ID = Cell ID to match, or none for all 
%
% Outputs:  Y  =  [M,N] Matrix, channel,data
*/
/*----------------------------------------------------------------------*/
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
  double *data, *timestamp, *cell_id;

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
  else if(nlhs > 3) { mexErrMsgTxt("Too many output arguments"); }
 
  int channel           = mxGetScalar(prhs[0]) - 1; 
  int number_of_buffers = mxGetScalar(prhs[1]);
  int input_cell_id     = mxGetScalar(prhs[2]);
 
  printf("Channel: %d, N: %d, Cell_ID: %d\n",channel, number_of_buffers, input_cell_id); 
  /* Create matrix for the return argument. [mxREAL|mxCOMPLEX] */
  plhs[0] = mxCreateDoubleMatrix(1,number_of_buffers, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(NSTREAM_SPIKE_TOTAL_SAMPLES,number_of_buffers, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1,number_of_buffers, mxREAL);

  timestamp = mxGetPr(plhs[0]); 
  data = mxGetPr(plhs[1]);
  cell_id = mxGetPr(plhs[2]);

  /*========================== Critical Section =======================*/

  write_pos = share->spike_index[channel] - 1;
  if(write_pos < 0) {
    write_pos = NSTREAM_SPIKE_NUM_ELEMENTS + write_pos;
  }

  /*writing backwards into the matrix*/
  for(i=0; i < number_of_buffers; i++)
  { 
     timestamp[number_of_buffers - 1 - i] = share->spike_event_buffer[channel][write_pos].timestamp/30;
     printf("Timestamp: %d\n",  share->spike_event_buffer[channel][write_pos].timestamp/30);

     cell_id[number_of_buffers - 1 - i] = share->spike_event_buffer[channel][write_pos].cell_id;
     for (j=0; j < NSTREAM_SPIKE_TOTAL_SAMPLES; j++)
     {
        data[i*NSTREAM_SPIKE_TOTAL_SAMPLES + j] = share->spike_event_buffer[channel][write_pos].waveform[j];
        /* data[(number_of_buffers - 1 - i)*NSTREAM_SPIKE_TOTAL_SAMPLES + j] = share->spike_event_buffer[channel - 1][write_pos].waveform[j]; */
     }
     write_pos--;
	 if(write_pos < 0) { 	 
	   write_pos = NSTREAM_SPIKE_NUM_ELEMENTS + write_pos; 
     }
  }  
  /*===================== End of Critical Section ======================*/

  if (shmdt(share) == -1) { fprintf(stderr, "shmdt function failed\n"); exit(EXIT_FAILURE); } 
}  

/*----------------------------------------------------------------------*/
