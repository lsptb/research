#ifndef __NETKEEPALIVE_H__
#define __NETKEEPALIVE_H__

#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdbool.h>
#include "hw_defines.h"
#include "libnspike.h"

/* the switch will forget which port this machine is on unless
 * a packet is sent out from time to time.  this class creates
 * a thread and sends a junk packet to the given address and
 * port every interval seconds */

struct NetKeepAlive {
    pthread_t keepalive_thread;
    bool threads_started;
    bool shutdown;
    UDPSocket *keepalive_socket;
    int send_interval;
    Exception *ex;
};

NetKeepAlive *NetKeepAlive_create(char *address, int port, int interval);
bool NetKeepAlive_run(NetKeepAlive *nka, char *address, int port, int interval);
void NetKeepAlive_destroy(NetKeepAlive *nka);
void NetKeepAlive__keepalive_thread_main(void *void_nka_ptr);

#endif

