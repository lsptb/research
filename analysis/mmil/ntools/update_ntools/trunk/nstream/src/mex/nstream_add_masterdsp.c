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
% Y = nstream_add_masterdsp(ip_address)
%
% configures a masterdsp 
%
% Inputs: ip_address   = address of masterdsp
%
% Outputs: none.  fails on error. 
%
% NOTE: You must add the masterdsp from the master system first in a
%       two system setup
%
% NOTE2: You must add at least one masterdsp before calling
%        nstream_start_acquire 
%
*/
/*----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char error_text[NSTREAMCLIENT_ERROR_LENGTH];
  char ip_address[ADDR_LEN];
  
  mxGetString(prhs[0], ip_address, ADDR_LEN);           

  NStreamClient *nsc = NStreamClient_create();
    
  if (!NStreamClient_add_masterdsp(nsc, (char*)&ip_address))
  {
      NStreamClient_get_error(nsc, (char*)&error_text, NSTREAMCLIENT_ERROR_LENGTH);
      mexErrMsgTxt(error_text);
  }

  NStreamClient_destroy(nsc);
                
}

/*----------------------------------------------------------------------*/

