#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/io.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <termios.h>
#include <stdio.h>
#include "mex.h"

/* 
% setparallel(port, byte)
% 
%  port: number of parallel port (1, 2, 3).  
%  NOTE: parallel port bases used are machine dependent
%  and may need to be changed when moving across machines. 
%
%  byte: one byte representing 8 bits of state for the 8 pins of the port
%
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int value = 1;  
  int port = 0;
  unsigned long base = 0;

    /*----- Check for proper number of arguments. -----*/ 
  if(nrhs!=2) { mexErrMsgTxt("port and value required"); } 
	        
  port  =  mxGetScalar(prhs[0]);          
  value =  mxGetScalar(prhs[1]);

  if (port == 0)
    base = 0x378;
  else if (port == 1)
    base = 0x6000;
  else if (port == 2)
    base = 0x6800;
  
 /* printf("value is %d\n",value); */

  if (ioperm(base,1,1))
    /*fprintf(stderr, "Couldn't get the port at %x\n", base), exit(1);*/
    perror("ioperm"); 
    outb(value, base);
		                  
}
