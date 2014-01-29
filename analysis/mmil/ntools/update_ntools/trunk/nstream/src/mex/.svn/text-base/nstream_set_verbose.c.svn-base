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
% Y = nstream_set_verbose(enable)
%
% enable verbosity for nstream engine messages 
%
% Inputs:	 enable  = int (0 or 1, off or on) 
%
% Outputs:  Y  =  
%
 ----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int enable;
  char error_text[NSTREAMCLIENT_ERROR_LENGTH];
  
  /*----- Check for proper number of arguments. -----*/ 
  if(nrhs!=1) { mexErrMsgTxt("one inputs required"); } 
 
  enable = mxGetScalar(prhs[0]);

  NStreamClient *nsc = NStreamClient_create();
    
  if (!NStreamClient_set_verbose(nsc, enable?1:0))
  {
      NStreamClient_get_error(nsc, (char*)&error_text, NSTREAMCLIENT_ERROR_LENGTH);
      mexErrMsgTxt(error_text);
  }

  NStreamClient_destroy(nsc);
                
}

/*----------------------------------------------------------------------*/

