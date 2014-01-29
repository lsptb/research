#ifndef __UDPSOCKET_H__
#define __UDPSOCKET_H__

#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdbool.h>
#include <pthread.h>
#include "cfulist.h"
#include "hw_defines.h"
#include "libnspike.h"

typedef int (*udpsocket_read_func_ptr_type)(UDPSocket *, void *, int);
typedef int (*udpsocket_write_func_ptr_type)(UDPSocket *, void *, int);
typedef int (*udpsocket_blockedread_func_ptr_type)(UDPSocket *, void *, int, struct timeval *);

struct UDPSocket {
    struct sockaddr_in address;
    struct timeval timeout;
    Exception *ex;
    int socket;
    fd_set fdset;
    udpsocket_read_func_ptr_type udpsocket_read_func_ptr;
    udpsocket_write_func_ptr_type udpsocket_write_func_ptr;
    udpsocket_blockedread_func_ptr_type udpsocket_blockedread_func_ptr;
    void *read_func_ptr_userdata;
    void *write_func_ptr_userdata;
    void *blockedread_func_ptr_userdata;
};

UDPSocket *UDPSocket_create();

/* for outbound data */
bool UDPSocket_connect(UDPSocket *udps, char *address, int port);

/* for listening to inbound data */
bool UDPSocket_bind(UDPSocket *udps, int port);
bool UDPSocket_get_address(UDPSocket *udps, char *ip, int len);

bool UDPSocket_set_recv_buf_size(UDPSocket *udps, int size);
void UDPSocket_set_read_func(UDPSocket *udps, udpsocket_read_func_ptr_type *func, void *userptr);
void UDPSocket_set_write_func(UDPSocket *udps, udpsocket_write_func_ptr_type *func, void *userptr);
void UDPSocket_set_blockedread_func(UDPSocket *udps, udpsocket_blockedread_func_ptr_type *func, void *userptr);

void *UDPSocket_get_read_userptr(UDPSocket *udps);
void *UDPSocket_get_write_userptr(UDPSocket *udps);
void *UDPSocket_get_blockedread_userptr(UDPSocket *udps);

int UDPSocket_read(UDPSocket *udps, void *buf, int count);
int UDPSocket_write(UDPSocket *udps, void *buf, int count);
int UDPSocket_blockedread(UDPSocket *udps, void *buf, int count, struct timeval *tv);

int UDPSocket__native_read(UDPSocket *udps, void *buf, int count);
int UDPSocket__native_write(UDPSocket *udps, void *buf, int count);
int UDPSocket__native_blockedread(UDPSocket *udps, void *buf, int count, struct timeval *tv);

void UDPSocket_destroy(UDPSocket *udps);

#endif

