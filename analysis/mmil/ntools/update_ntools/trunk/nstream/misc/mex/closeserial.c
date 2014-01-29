#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <termios.h>
#include <stdio.h>
#include "mex.h"

/* 
 * closeserial.c - a mex function to close an open serial port.
 *
 *
% closeserial(serial_fid);
%
% closes the open serial port referenced by serial_fid
%
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int serial_fd = 0;

  /*----- Check for proper number of arguments. -----*/ 
  if(nrhs!=1) { mexErrMsgTxt("file descriptor number required"); } 
  
  serial_fd = mxGetScalar(prhs[0]);           
  
  close(serial_fd);
}

