
#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include "libnspike.h"
#include "testmacros.h"

#define RECV     1 
#define SEND     2

typedef struct fakehw_msg {
    unsigned short msg[32];
    int len;
    int type;
} fakehw_msg;

typedef struct fakehw {
    int socket;
    int port;
    struct sockaddr_in address;
    struct sockaddr_in peeraddr;
    socklen_t peerlen;
    char testname[128];
    int ptr;
    fakehw_msg msgs[256];
    pthread_t fakehw_thread;
    bool failed;
    char failreason[256];
    char last_error[256];

    fd_set fdset;
} fakehw;

fakehw *fakehw_create(char *address, int port);
void fakehw_destroy(fakehw *fhw);
void fakehw_start_seq(fakehw *fhw, char *name);
void fakehw_add_msg(fakehw *fhw, unsigned short msg[32], int len, int type);
void fakehw_recv(fakehw *fhw, unsigned short msg[32], int len);
void fakehw_send(fakehw *fhw, unsigned short msg[32], int len);
void fakehw_run(fakehw *fhw);
bool fakehw_check(fakehw *fhw);
void fakehw_thread_main(void *fhw_ptr);


