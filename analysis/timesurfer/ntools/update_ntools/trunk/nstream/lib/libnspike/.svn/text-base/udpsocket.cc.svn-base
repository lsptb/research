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
#include <sys/select.h>
#include <unistd.h>
#include "udpsocket.h"

UDPSocket *UDPSocket_create()
{
    UDPSocket *udps;
   
    if ((udps = (UDPSocket*)malloc(sizeof(UDPSocket))) == NULL)
        return false;

    memset(udps, 0, sizeof(UDPSocket));

    udps->ex = Exception_create("");

    // initialize function pointers to native versions
    // which can be overrided later by with set_xxx_func
   
    udps->udpsocket_read_func_ptr = UDPSocket__native_read;
    udps->udpsocket_write_func_ptr = UDPSocket__native_write;
    udps->udpsocket_blockedread_func_ptr = UDPSocket__native_blockedread;

    return udps;
}


bool UDPSocket_connect(UDPSocket *udps, char *address, int port)
{
    if (!inet_pton(AF_INET, address, &udps->address.sin_addr))
    {
        Exception_set(udps->ex, "Unable to parse IP address");
        return false;
    }
    
    udps->address.sin_port   = htons(port);
    udps->address.sin_family = AF_INET;

    /* create the socket */
    if ((udps->socket = socket(PF_INET, SOCK_DGRAM, 0)) < 0)
    {
        Exception_set(udps->ex, "Unable to create UDP socket: %s", strerror(errno));
        return false;
    }

    /* connect the socket */
    if (connect(udps->socket, (struct sockaddr *) &udps->address, sizeof(udps->address))==-1) {
        Exception_set(udps->ex, "Cannot connect() UDP socket: %s", strerror(errno));
        return false;
    }

    /* set non-blocking */
    fcntl(udps->socket, F_SETFL, O_NONBLOCK);

    return true;
}

bool UDPSocket_bind(UDPSocket *udps, int port)
{
    udps->address.sin_port        = htons(port);
    udps->address.sin_family      = AF_INET;
    udps->address.sin_addr.s_addr = htonl(INADDR_ANY);

    /* create the socket */
    if ((udps->socket = socket(PF_INET, SOCK_DGRAM, 0)) < 0)
    {
        Exception_set(udps->ex, "Unable to create UDP socket: %s", strerror(errno));
        return false;
    }

    /* bind the socket */
    if (bind(udps->socket, (struct sockaddr *) &udps->address, sizeof(udps->address))==-1) {
        Exception_set(udps->ex, "Cannot bind() UDP socket: %s", strerror(errno));
        return false;
    }

    /* set non-blocking */
    fcntl(udps->socket, F_SETFL, O_NONBLOCK);

    return true;
}

bool UDPSocket_get_address(UDPSocket *udps, char *ip, int len)
{
    if (!inet_ntop(AF_INET, &udps->address.sin_addr, ip, len) ||
         udps->address.sin_port == 0)
    {
        Exception_set(udps->ex, "Socket has not yet been connected to an address");
        return false;
    }

    return true;
}

bool UDPSocket_set_recv_buf_size(UDPSocket *udps, int size)
{
    int check_size = 0;
    socklen_t check_size_len = sizeof(int);

    if (setsockopt(udps->socket, SOL_SOCKET, SO_RCVBUF, &size, sizeof(size)) == -1)
    {
        Exception_set(udps->ex, "Cannot set recv buf size for UDP socket: %s", strerror(errno));
        return false;
    }

    if (getsockopt(udps->socket, SOL_SOCKET, SO_RCVBUF, &check_size, &check_size_len) == -1)
    {
        Exception_set(udps->ex, "Cannot get recv buf size for UDP socket: %s", strerror(errno));
        return false;
    }

    /* linux is weird, it reports the socket buffer size to be 2x what we set it to */
    if (check_size != size * 2)
    {
        Exception_set(udps->ex, "Requested socket buffer size %d too large, truncated to %d", size, check_size);
        return false;
    }
    
    return true;
}

void UDPSocket_set_read_func(UDPSocket *udps, udpsocket_read_func_ptr_type *func, void *userptr)
{
    udps->udpsocket_read_func_ptr = (udpsocket_read_func_ptr_type)func;
    udps->read_func_ptr_userdata = userptr;
}

void UDPSocket_set_write_func(UDPSocket *udps, udpsocket_write_func_ptr_type *func, void *userptr)
{
    udps->udpsocket_write_func_ptr = (udpsocket_write_func_ptr_type)func;
    udps->write_func_ptr_userdata = userptr;
}

void UDPSocket_set_blockedread_func(UDPSocket *udps, udpsocket_blockedread_func_ptr_type *func, void *userptr)
{
    udps->udpsocket_blockedread_func_ptr = (udpsocket_blockedread_func_ptr_type)func;
    udps->blockedread_func_ptr_userdata = userptr;
}

int UDPSocket_read(UDPSocket *udps, void *buf, int count)
{
    return (udps->udpsocket_read_func_ptr)(udps, buf, count);
}

int UDPSocket_write(UDPSocket *udps, void *buf, int count)
{
    return (udps->udpsocket_write_func_ptr)(udps, buf, count);
}

int UDPSocket_blockedread(UDPSocket *udps, void *buf, int count, struct timeval *tv)
{
    return (udps->udpsocket_blockedread_func_ptr)(udps, buf, count, tv);
}

void *UDPSocket_get_read_userptr(UDPSocket *udps)
{
    return udps->read_func_ptr_userdata;
}

void *UDPSocket_get_write_userptr(UDPSocket *udps)
{
    return udps->write_func_ptr_userdata;
}

void *UDPSocket_get_blockedread_userptr(UDPSocket *udps)
{
    return udps->blockedread_func_ptr_userdata;
}

int UDPSocket__native_read(UDPSocket *udps, void *buf, int count)
{
    int bytesread;
    char tmp[1024];

    if (-1 == (bytesread = read(udps->socket, buf, count)))
    {
        // Exception_set(udps->ex, "read error: %s", strerror_r(errno, tmp, sizeof(tmp)));
        return -1;
    }

    return bytesread;
}

int UDPSocket__native_write(UDPSocket *udps, void *buf, int count)
{
    int byteswritten;

    if (count != (byteswritten = write(udps->socket, buf, count)))
    {
        char tmp[200];
        Exception_set(udps->ex, "write error: %s", strerror_r(errno, tmp, sizeof(tmp)));
        return -1;
    }

    return byteswritten;
}


int UDPSocket__native_blockedread(UDPSocket *udps, void *buf, int count, struct timeval *tv)
{
    int bytesread;

    FD_ZERO(&udps->fdset);
    FD_SET(udps->socket, &udps->fdset);

    select(udps->socket + 1, &udps->fdset, NULL, NULL, tv);
    if (!FD_ISSET(udps->socket, &udps->fdset))
    {
    // return UDPSocket__native_read(udps, buf, count);
        return 0;
    }

    return UDPSocket__native_read(udps, buf, count);
}

void UDPSocket_destroy(UDPSocket *udps)
{
    if (udps->socket != 0)
        close(udps->socket);

    Exception_destroy(udps->ex);
    free(udps);
}

