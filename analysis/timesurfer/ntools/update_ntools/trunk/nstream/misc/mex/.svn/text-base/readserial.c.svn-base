#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <termios.h>
#include <stdio.h>
#include "mex.h"

#define SERIAL_STRING_LENGTH     16384 /* max length string we'll read */ 
#define TERMINATOR_STRING_LENGTH 8  

/* 
 * readserial.c - a mex function to read a string from a serial port.
 *
 *
%  x = readserial(serial_fd, "terminator string", timeout);
%
%  reads a string from a serial port terminated by "terminator string"
%  blocks until timeout reached (secs) or terminator read from port
%
 */


bool read_data(int fd, char *terminator, char *serial_data, int serial_data_len, int timeout);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char terminator_string[TERMINATOR_STRING_LENGTH];
    char *serial_data;
    int serial_fd;
    int timeout;
    int flags;

    /*----- Check for proper number of arguments. -----*/ 
    if(nrhs!=3) { mexErrMsgTxt("serial_fd, terminator string and timeout required"); } 

    serial_fd = mxGetScalar(prhs[0]);           
    mxGetString(prhs[1], terminator_string, TERMINATOR_STRING_LENGTH);           
    timeout = mxGetScalar(prhs[2]);           


    serial_data     = mxCalloc(SERIAL_STRING_LENGTH, sizeof(char));
   
    flags = fcntl(serial_fd, F_GETFL); 
    fcntl(serial_fd, F_SETFL, flags|O_NONBLOCK);
    if (read_data(serial_fd, terminator_string, serial_data, SERIAL_STRING_LENGTH, timeout) == false)
    {
        mexErrMsgTxt("read serial failed or timed out");
    }
    else
    {
        plhs[0] = mxCreateString(serial_data);
    }
    fcntl(serial_fd, F_SETFL, flags);
    
            
}

bool read_data(int fd, char *terminator, char *serial_data, int serial_data_len, int timeout)
{

    fd_set serial_fd;
    char *serial_data_ptr = serial_data;
    
    int start = gettimeinsec();
    int bytes_remaining = 0;
    int bytes_read = 0;
    int terminator_length = strlen(terminator);
    struct timeval spin_delay;

    fcntl(fd, F_SETFL, O_NONBLOCK);
    
    while (1)
    {
        FD_ZERO(&serial_fd);
        FD_SET(fd, &serial_fd);
        spin_delay.tv_sec  = 1;
        spin_delay.tv_usec = 0;
   
        if (gettimeinsec() - start > timeout)
        {
            printf("read timed out\n");
            return false; 
        }
       
        bytes_remaining = (SERIAL_STRING_LENGTH - (serial_data_ptr - serial_data)); 
        if (bytes_remaining <= 0)
        {
            printf("buffer overrun\n");
            return false;
        }

        select(fd + 1, &serial_fd, NULL, NULL, &spin_delay);
        if (FD_ISSET(fd, &serial_fd))
        {
            bytes_read = read(fd, serial_data_ptr, bytes_remaining); 
            serial_data_ptr += bytes_read;
        
            if ((serial_data_ptr - serial_data) >= terminator_length &&
                strncmp(terminator, serial_data_ptr - terminator_length, terminator_length) == 0)
            {
                *--serial_data_ptr = 0;
                return true;
            }
        }
            
    }

}

int gettimeinsec()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec;
}

