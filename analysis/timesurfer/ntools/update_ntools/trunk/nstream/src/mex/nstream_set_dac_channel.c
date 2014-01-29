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
% Y = nstream_set_dac_channel(dac_channel, data_channel)
%
% selects channels for nspike d/a output 
%
% Inputs: dac_channel  = d/a number (0-1 or 0-3) 
%         data_channel = neural data channel number
%
% NOTE: on a 256 channel system, you can only use D/A 0,1 for ch 1-128
%       and D/A 2,3 for ch 129-256.
%
% Outputs: none.  fails on error.
*/
/*----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char error_text[NSTREAMCLIENT_ERROR_LENGTH];
  int dac_channel, data_channel;

  dac_channel = mxGetScalar(prhs[0]);
  data_channel = mxGetScalar(prhs[1]);

  NStreamClient *nsc = NStreamClient_create();
    
  if (!NStreamClient_set_dac_channel(nsc, dac_channel, data_channel))
  {
      NStreamClient_get_error(nsc, (char*)&error_text, NSTREAMCLIENT_ERROR_LENGTH);
      mexErrMsgTxt(error_text);
  }

  NStreamClient_destroy(nsc);
                
}

/*----------------------------------------------------------------------*/

