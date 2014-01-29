#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <termios.h>
#include <stdio.h>
#include "mex.h"

/* 
 * openserial.c - a mex function to open a serial port.
 *                port settings are hardwired in source for now.
 *
 *
% serial_fid = openserial();
%
% opens a serial port and returns a file descriptor for it.
% port settings are currently hardcoded for now
%
 */

/* 
 * to configure baudrate: change BAUDRATE define below
 * to configure port: change SERIALDEVICE define below
 * to configure flow control: configure newtio.c_cflag below, valid settings:
 *
 *

/* baudrate settings are defined in <asm/termbits.h>, which is
included by <termios.h> */
/*
  valid baudrates:
        B0
        B50
        B75
        B110
        B134
        B150
        B200
        B300
        B600
        B1200
        B1800
        B2400
        B4800
        B9600
        B19200
        B38400
        B57600
        B115200
        B230400
*/
#define BAUDRATE B115200           
/* change this definition for the correct port */
/* ttyS0 -> COM1
   ttyS1 -> COM2
   etc..
*/ 
#define SERIALDEVICE "/dev/ttyS0"
#define _POSIX_SOURCE 1 /* POSIX compliant source */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int fd = 0;
  double *ret;

  /*----- Check for proper number of arguments. -----*/ 
  if(nrhs!=0) { mexErrMsgTxt("no inputs required"); } 
  if(nlhs!=1) { mexErrMsgTxt("this function is useless if you ignore it's return value"); } 
  
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  ret = mxGetPr(plhs[0]);
  
  if ((fd = open_serial()) == -1)
  {
      mexErrMsgTxt("open serial failed");
      ret[0] = -1;
  }
  
  ret[0] = fd;  
}

int open_serial()
{
  int fd,c, res;
  struct termios oldtio,newtio;
  char buf[255];
/* 
  Open modem device for reading and writing and not as controlling tty
  because we don't want to get killed if linenoise sends CTRL-C.
*/
 fd = open(SERIALDEVICE, O_RDWR | O_NOCTTY ); 
 if (fd <0) {perror(SERIALDEVICE); close(fd); return -1; }

 tcgetattr(fd,&oldtio); /* save current serial port settings */
 bzero(&newtio, sizeof(newtio)); /* clear struct for new port settings */

/* 
  BAUDRATE: Set bps rate. You could also use cfsetispeed and cfsetospeed.
  CRTSCTS : output hardware flow control (only used if the cable has
            all necessary lines. See sect. 7 of Serial-HOWTO)
  CS8     : 8n1 (8bit,no parity,1 stopbit)
  CLOCAL  : local connection, no modem contol
  CREAD   : enable receiving characters
*/
 newtio.c_cflag = BAUDRATE | CRTSCTS | CS8 | CLOCAL | CREAD;
 
/*
  IGNPAR  : ignore bytes with parity errors
  ICRNL   : map CR to NL (otherwise a CR input on the other computer
            will not terminate input)
  otherwise make device raw (no other input processing)
*/
 newtio.c_iflag = IGNPAR | ICRNL;
 
/*
 Raw output.
*/
 newtio.c_oflag = 0;
 
/*
  raw input.
*/
newtio.c_lflag &= ~(ICANON | ECHO | ECHOE | ISIG); 

/* 
  now clean the modem line and activate the settings for the port
*/
 if (tcflush(fd, TCIFLUSH) == -1)
 {
    perror("tcflush");
    close(fd); 
    return -1;
 }


 if (tcsetattr(fd,TCSANOW,&newtio) == -1)
 {
    perror("tcsetattr");
    close(fd); 
    return -1;
 }

 return fd;
}

