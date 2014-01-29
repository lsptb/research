#include <stdio.h>
#include <stdlib.h>
/*----------------------------------------------------------------------*/
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <string.h>
#include "nstreamclient.h"
#include "mex.h"

#define ADDR_LEN 256

/*----------------------------------------------------------------------
%
% Y = nstream_add_auxdsp(ip_address)
%
% configures an auxdsp 
%
% Inputs:	 ip_address   = address of auxdsp
%
% Outputs:  Y  =  
%
% NOTE: Order matters when calling this function.  Adding them in
%       non-strictly-increasing order may result in network traffic
%       patterns that may result in packet loss.
%
*/
/*----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char error_text[NSTREAMCLIENT_ERROR_LENGTH];
  char ip_address[ADDR_LEN];
  
  mxGetString(prhs[0], ip_address, ADDR_LEN);           

  NStreamClient *nsc = NStreamClient_create();
    
  if (!NStreamClient_add_auxdsp(nsc, (char*)&ip_address))
  {
      NStreamClient_get_error(nsc, (char*)&error_text, NSTREAMCLIENT_ERROR_LENGTH);
      mexErrMsgTxt(error_text);
  }

  NStreamClient_destroy(nsc);
                
}

/*----------------------------------------------------------------------*/

