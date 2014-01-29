#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <termios.h>
#include <stdio.h>
#include "mex.h"

#define SERIAL_STRING_LENGTH 8192 /* max length string we'll send */ 

/* 
 * writeserial.c - a mex function to write a string to a serial port.
 *
%
% writeserial(serial_fd, "string");
%
% serial_fd: file descriptor returned by open serial
%    string: string to write to serial port
%
 */

bool write_data(int fd, char *data, int len);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char serial_string[SERIAL_STRING_LENGTH];
  int serial_fd;

  /*----- Check for proper number of arguments. -----*/ 
  if(nrhs!=2) { mexErrMsgTxt("serial_fd and string required"); } 
  
  mxGetString(prhs[1], serial_string, SERIAL_STRING_LENGTH);           
  serial_fd = mxGetScalar(prhs[0]);           

  if (write_data(serial_fd, (char*)&serial_string, strlen(serial_string)) == false)
  {
      mexErrMsgTxt("write serial failed");
  }
                
}

bool write_data(int fd, char *data, int len)
{
  
  if (write(fd, data, len) != len)
  {
     perror("write");
     return false;
  }

  return true; 
}
