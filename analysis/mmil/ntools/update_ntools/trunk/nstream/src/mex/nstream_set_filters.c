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
% Y = nstream_set_filters(data_channel, high_pass, low_pass)
%
% sets filters for an nspike data channel 
%
% Inputs:	 data_channel = logical data channel number (ie: 1-256) 
%		     high_pass    = high pass filter 
%            low_pass     = low pass filter 
%
% Outputs:  Y  =  
*/
/*----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char error_text[NSTREAMCLIENT_ERROR_LENGTH];
  int data_channel, high_pass, low_pass;

  /*----- Check for proper number of arguments. -----*/ 
  if(nrhs!=3) { mexErrMsgTxt("three inputs required"); } 
 
  data_channel = mxGetScalar(prhs[0]);
  high_pass = mxGetScalar(prhs[1]);
  low_pass = mxGetScalar(prhs[2]);

  NStreamClient *nsc = NStreamClient_create();
    
  if (!NStreamClient_set_channel_filters(nsc, data_channel, high_pass, low_pass))
  {
      NStreamClient_get_error(nsc, (char*)&error_text, NSTREAMCLIENT_ERROR_LENGTH);
      mexErrMsgTxt(error_text);
  }

  NStreamClient_destroy(nsc);
                
}

/*----------------------------------------------------------------------*/

