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
% Y = nstream_exit()
%
% nstream_exit stops the nstream engine process 
%
% Inputs:	 (none) 
%
% Outputs:  Y  =  
*/
/*----------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char error_text[256];
    /*----- Check for proper number of arguments. -----*/ 
    if(nrhs!=0) { mexErrMsgTxt("no inputs required"); } 

    NStreamClient *nsc = NStreamClient_create();
    
    if (!NStreamClient_exit_nstream(nsc))
    {
        NStreamClient_get_error(nsc, (char*)&error_text, 256);
        mexErrMsgTxt(error_text);
    }

    NStreamClient_destroy(nsc);
                    
}

/*----------------------------------------------------------------------*/

