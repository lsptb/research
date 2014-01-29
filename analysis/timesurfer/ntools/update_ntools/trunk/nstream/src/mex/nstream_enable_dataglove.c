#include <stdio.h>
#include <stdlib.h>
/*----------------------------------------------------------------------*/
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <string.h>
#include "nstreamclient.h"
#include "mex.h"

#define NAME_LENGTH 256

/*----------------------------------------------------------------------
%
% Y = nstream_enable_dataglove(serial_port, delay)
%
% enables data glove support on serial_port, waiting a minimum of
% delay ticks of the nspike clock between sample capture 
%
% Inputs:	 serial_port  = serial port device (ie: /dev/ttyS0) 
%		     delay        = nspike clock ticks to wait between sampling 
%
% Outputs:  Y  =  
*/
/*----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char error_text[NSTREAMCLIENT_ERROR_LENGTH];
  char serial_port[NAME_LENGTH];
  int delay = 0;

  /*----- Check for proper number of arguments. -----*/ 
  if(nrhs!=2) { mexErrMsgTxt("two inputs required"); } 
  
  mxGetString(prhs[0], serial_port, NAME_LENGTH);           
  delay = mxGetScalar(prhs[1]);

  NStreamClient *nsc = NStreamClient_create();
    
  if (!NStreamClient_enable_dataglove(nsc, &serial_port, delay))
  {
      NStreamClient_get_error(nsc, (char*)&error_text, NSTREAMCLIENT_ERROR_LENGTH);
      mexErrMsgTxt(error_text);
  }

  NStreamClient_destroy(nsc);
                
}

/*----------------------------------------------------------------------*/

