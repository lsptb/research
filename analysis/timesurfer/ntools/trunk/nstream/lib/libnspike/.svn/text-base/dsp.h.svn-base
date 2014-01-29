#ifndef __DSP_H__
#define __DSP_H__

#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdbool.h>
#include <pthread.h>
#include "cfulist.h"
#include "hw_defines.h"
#include "udpsocket.h"

struct DSP {
    int dspnum;
    UDPSocket *control_socket;
    Exception *ex;
};

DSP *DSP_create(int dspnum, UDPSocket *socket);

bool DSP_get_address(DSP *dsp, char *ip, int len);

int DSP_write_data(DSP *dsp, unsigned short address, short n_words, unsigned short *data);
int DSP_read_data(DSP *dsp, unsigned short address, short n, unsigned short *data);
int DSP_get_response(DSP *dsp, int n, unsigned short *data);
void DSP__byte_swap(unsigned short *sptr, int nelem);
void DSP__byte_swap_long(u32 *lptr, int nelem);

void DSP_destroy(DSP *dsp);

#endif

