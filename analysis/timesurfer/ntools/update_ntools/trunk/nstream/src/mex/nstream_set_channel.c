#include <stdio.h>
#include <stdlib.h>
/*----------------------------------------------------------------------*/
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <string.h>
#include "nstreamclient.h"
#include "mex.h"

/*----------------------------------------------------------------------
%
% Y = nstream_set_channel(sw_channel, hw_channel)
%
% maps hw_channel to sw_channel 
%
% Inputs:	 sw_channel   = logical channel number 
%		     hw_channel   = hardware channel number
%
% Outputs:  Y  =  
*/
/*----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char error_text[NSTREAMCLIENT_ERROR_LENGTH];
  int sw_channel, hw_channel;

  sw_channel = mxGetScalar(prhs[0]);
  hw_channel = mxGetScalar(prhs[1]);

  NStreamClient *nsc = NStreamClient_create();
    
  if (!NStreamClient_set_channel(nsc, sw_channel, hw_channel))
  {
      NStreamClient_get_error(nsc, (char*)&error_text, NSTREAMCLIENT_ERROR_LENGTH);
      mexErrMsgTxt(error_text);
  }

  NStreamClient_destroy(nsc);
                
}

/*----------------------------------------------------------------------*/

