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
% Y = nstream_set_dac_gain(dac_channel, gain)
%
% sets gain for nspike d/a output 
%
% Inputs: dac_channel  = d/a number (0-1 or 0-3) 
%         gain         = gain (1-100) 
%
% Outputs: none.  fails on error. 
*/
/*----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int dac_channel, gain;
  char error_text[NSTREAMCLIENT_ERROR_LENGTH];
  /*----- Check for proper number of arguments. -----*/ 
  if(nrhs!=2) { mexErrMsgTxt("two inputs required"); } 
 
  dac_channel = mxGetScalar(prhs[0]);
  gain = mxGetScalar(prhs[1]);

  NStreamClient *nsc = NStreamClient_create();
    
  if (!NStreamClient_set_dac_gain(nsc, dac_channel, gain))
  {
      NStreamClient_get_error(nsc, (char*)&error_text, NSTREAMCLIENT_ERROR_LENGTH);
      mexErrMsgTxt(error_text);
  }

  NStreamClient_destroy(nsc);
}

/*----------------------------------------------------------------------*/

