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
% Y = NSTREAM_GETDIO_TS(TS)
%
% NSTREAM_GETDIO_TS returns all known digital io data points with a timestamp > ts 
%
% Inputs:	 N	=  Scalar, number of data points to return 
%
% Outputs:   Y  =  [M,N] Matrix, channel,data
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
  int number_of_events = 0;

  /* first determine how many events we have to return and bail out if zero */
  write_pos = share->dio_message_index - 1;
  if(write_pos < 0) {
    write_pos = NSTREAM_DIO_STATE_NUM_ELEMENTS + write_pos;
  }

  for (number_of_events = 0; number_of_events < NSTREAM_DIO_STATE_NUM_ELEMENTS; number_of_events++)
  {
    if (share->dio_message_buffer[write_pos].nspike_global_timestamp/30 <= start_ts)
      break;
    write_pos--;
	if(write_pos < 0) { 	 
	  write_pos = NSTREAM_DIO_STATE_NUM_ELEMENTS + write_pos; 
    }
  }

   
  /* Create matrix for the return argument. [mxREAL|mxCOMPLEX] */
  plhs[0] = mxCreateDoubleMatrix(5,number_of_events, mxREAL);
  data = mxGetPr(plhs[0]);
  /*========================== Critical Section =======================*/

  write_pos = share->dio_message_index - 1;
  if(write_pos < 0) {
    write_pos = NSTREAM_DIO_STATE_NUM_ELEMENTS + write_pos;
  }

  /*writing backwards into the matrix*/
  for(i=0; i < number_of_events; i++)
  { 
     data[5*number_of_events - 1 - 5*i] = share->dio_message_buffer[write_pos].state[1] & 0xF;
     data[5*number_of_events - 2 - 5*i] = share->dio_message_buffer[write_pos].state[1] >> 8;
     data[5*number_of_events - 3 - 5*i] = share->dio_message_buffer[write_pos].state[0] & 0xF;
     data[5*number_of_events - 4 - 5*i] = share->dio_message_buffer[write_pos].state[0] >> 8;
     data[5*number_of_events - 5 - 5*i] = share->dio_message_buffer[write_pos].nspike_global_timestamp/30;
     write_pos--;
	 if(write_pos < 0) { 	 
	   write_pos = NSTREAM_DIO_STATE_NUM_ELEMENTS + write_pos; 
     }
  }  
  /*===================== End of Critical Section ======================*/

  if (shmdt(share) == -1) { fprintf(stderr, "shmdt function failed\n"); exit(EXIT_FAILURE); } 
}  

/*----------------------------------------------------------------------*/
