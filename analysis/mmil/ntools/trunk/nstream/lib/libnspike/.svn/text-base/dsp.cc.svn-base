#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/select.h>
#include "cfulist.h"
#include "dsp.h"

DSP *DSP_create(int dspnum, UDPSocket *socket)
{
    DSP *dsp;
   
    if ((dsp = (DSP*)malloc(sizeof(DSP))) == NULL)
        return false;

    memset(dsp, 0, sizeof(DSP));

    dsp->ex = Exception_create("");
    dsp->dspnum = dspnum;
    dsp->control_socket = socket;

    return dsp;
}

bool DSP_get_address(DSP *dsp, char *ip, int len)
{
    if (!dsp->control_socket)
        Exception_set(dsp->ex, "DSP not connected");

    if(!UDPSocket_get_address(dsp->control_socket, ip, len))
    {
        Exception_propagate_with_text(dsp->control_socket->ex, dsp->ex, "DSP_get_address ");
        return false;
    }

    return true;
}


int DSP_write_data(DSP *dsp, unsigned short address, short n_words, unsigned short *data)
{
    int size;
    unsigned short dataout[1000]; 
    size = n_words * sizeof(unsigned short);

    int i;

    dataout[0] = NSPIKE_DEVICE_CONTROL_ADDR;
    dataout[1] = NSPIKE_SHORT_WRITE;
    dataout[2] = NSPIKE_DSP_SRAM;
    dataout[3] = address;
    dataout[4] = size;

    /* append data to the end of dataout */
    memcpy(dataout + 5, data, size);

    /* increment size to be the size of the full dataout array */
    size += 5 * sizeof(unsigned short);

    DSP__byte_swap(dataout, 5 + n_words);

    int count = 0;
    /* send out the command */
    for (count = 0; count <= 25; count++)
    {
        if (UDPSocket_write(dsp->control_socket, dataout, size) != size) {
            Exception_propagate_with_text(dsp->control_socket->ex, dsp->ex, "DSP_write_data ");
            return 0;
        }
        if (DSP_get_response(dsp, 0, NULL))
            return 1;
        Exception_printf("DSP_write_data: communications timeout, retrying, %d/25\n",count); 
        // usleep(100000);
    }

    return 0;
}



int DSP_read_data(DSP *dsp, unsigned short address, short n, 
                unsigned short *data)
    /* write data to the dsp. data must contain < 24 words */
{
    int 		size;
    int 		nread;
    unsigned short 	dataout[6], datain[26];
    (void)data; // prevent compiler warnings
    int i;

    if (!dsp->control_socket)
        Exception_set(dsp->ex, "DSP_read_data: not connected");
    
    size = (n+2)*sizeof(unsigned short);

    /* create the programming command */
    dataout[0] = NSPIKE_DEVICE_CONTROL_ADDR;
    dataout[1] = NSPIKE_SHORT_READ;
    dataout[2] = NSPIKE_DSP_SRAM;
    dataout[3] = address;
    dataout[4] = size;

    if (n > 24) {
	    Exception_set(dsp->ex, "DSP_read_data: can't read > 24 words from dsp");
	    return -1;
    }

    int 	dspacq, acqmessage, error;

	/* send the command to the DSP */
	DSP__byte_swap(dataout, 5);
    if (UDPSocket_write(dsp->control_socket, dataout, 5 * sizeof(unsigned short)) != 5 * sizeof(unsigned short)) {
        Exception_propagate_with_text(dsp->control_socket->ex, dsp->ex, "DSP_read_data ");
        return -1;
    }

	/* read from the DSP */
	usleep(10000);
    nread = UDPSocket_read(dsp->control_socket, datain, size);    
    
    /* swap bytes and check for a valid response */
	DSP__byte_swap(datain, n + 2);
	
    if (datain[0] != NSPIKE_DEVICE_CONTROL_ADDR) {
	    Exception_set(dsp->ex, "DSP_read_data: bad pipe_id: %u", datain[0]);
	    return -1;
	}

	if (datain[1] != 1) {
	    if (datain[1] == 3) {
		    Exception_set(dsp->ex, "DSP_read_data: bad address");
		    return -1;
	    }
	    else if (datain[1] == 5) {
		    Exception_set(dsp->ex, "DSP_read_data: unsupported command");
		    return -1;
	    }
	}
	/* copy the returned data to the output array */
	size -= 2 * sizeof(unsigned short);
	memcpy(data, datain + 2, size);

    return 1;
}

int DSP_get_response(DSP *dsp, int n, unsigned short *data)
{
    /* read a response to a programming command */
    unsigned short resp[NSPIKE_MAX_PACKET_SIZE];
    int size, dsize, readsize;
    int error = 0, done = 0;
    int ntries = 1;
    int i;
    struct timeval tv;
    tv.tv_sec  = 0;
    tv.tv_usec = 400000;
   // tv.tv_usec = 0;

    size = 2 * sizeof(unsigned short);
    dsize = n * sizeof(unsigned short);
    while (!done) {
        if (ntries > 0) {
            // printf("try %d\n", ntries);
            //readsize = UDPSocket_read(dsp->control_socket, resp, sizeof(resp));
            readsize = UDPSocket_blockedread(dsp->control_socket, resp, sizeof(resp), &tv);
            // printf("read %d bytes!\n",readsize);
        }
        else {
            break;
        }

        DSP__byte_swap(resp, NSPIKE_MAX_PACKET_SIZE);
        // printf("read %d bytes\n", readsize);
        /* check for status codes */
        if ((readsize > 0)) {
            done = 1;
            /* this is not a left over data packet, so parse it */
            if (resp[0] != NSPIKE_DEVICE_CONTROL_ADDR) {
                Exception_set(dsp->ex, "Error writing to dsp bad PortID %u\n", resp[0]);
                error = 1;
                break;
            }
            if (resp[1] != 1) {
                if (resp[1] == 3) {
                    Exception_set(dsp->ex, "Error writing to dsp: bad address\n");
                    error = 1;
                    break;
                }
                else if (resp[1] == 5) {
                    Exception_set(dsp->ex, "DSP_get_response: unsupported command\n");
                    error = 1;
                    break;
                }
            }
        }
        else if (readsize <= 0) {
            ntries--;
        }
    }
    
    if (readsize <= 0)
    {
        error = 1;
        // Exception_set(dsp->ex, "no response from DSP\n");
    }

    if (error)
    {
        return 0;
    }

    /* copy the data if it is available */
    if (n) {
	    memcpy(data, resp + 2, readsize - size);
        return 1;
    }
    
    return 1;
}


void DSP_destroy(DSP *dsp)
{
    if (dsp->control_socket != 0)
        UDPSocket_destroy(dsp->control_socket);

    Exception_destroy(dsp->ex);
    free(dsp);
}
    
void DSP__byte_swap(unsigned short *sptr, int nelem)
{
    int i;
    char *cptr, ctmp;
    /* swap the first and second bytes of each short to create a shorts for the
     * dsps */
    cptr = (char *) sptr;
    for (i = 0; i < nelem; i++, cptr+=2) {
        ctmp = *cptr;
        *cptr = *(cptr+1);
        *(cptr+1) = ctmp;
    }
    return;
}

