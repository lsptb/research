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
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include "libnspike.h"

NetKeepAlive *NetKeepAlive_create(char *address, int port, int interval)
{
    NetKeepAlive *nka;

    if ((nka = (NetKeepAlive*)malloc(sizeof(NetKeepAlive))) == NULL)
        return false;

    memset(nka, 0, sizeof(NetKeepAlive));

    nka->ex = Exception_create("NetKeepAlive: ");
    nka->threads_started = false;

    return nka;
}

bool NetKeepAlive_run(NetKeepAlive *nka, char *address, int port, int interval)
{
    nka->send_interval = interval;

    if (nka->keepalive_socket)
    {
        Exception_set(nka->ex, "NetKeepAlive does not support more than one call to run per instance");
        return false;
    }

    nka->keepalive_socket = UDPSocket_create();
    if (!UDPSocket_connect(nka->keepalive_socket, address, port))
    {
        Exception_propagate(nka->keepalive_socket->ex, nka->ex);
        UDPSocket_destroy(nka->keepalive_socket);
        return false;
    }

    pthread_create(&nka->keepalive_thread, NULL, (void* (*) (void*))&NetKeepAlive__keepalive_thread_main, (void*)nka);
    nka->threads_started = true;
    
    return true;
}


void NetKeepAlive_destroy(NetKeepAlive *nka)
{
    if (nka->threads_started)
    {
        nka->shutdown = true;
        
        pthread_join(nka->keepalive_thread, NULL);
    }
    
    Exception_destroy(nka->ex);
    free(nka);
}


void NetKeepAlive__keepalive_thread_main(void *void_nka_ptr)
{
    NetKeepAlive *nka;
    nka = (NetKeepAlive*) void_nka_ptr;

    struct timeval tv; 
    int last_send_time = 0;
    
    char bogus_data[64];
    memset(bogus_data, 0, 64);

    for(;;)
    {
        sleep(1);

        if (nka->shutdown)
        {
            pthread_exit(NULL);
        }
       
        gettimeofday(&tv, NULL);

        if (tv.tv_sec - last_send_time > nka->send_interval)
        {
            if (UDPSocket_write(nka->keepalive_socket, &bogus_data, sizeof(bogus_data)) != sizeof(bogus_data)) {
                Exception_printf("WARNING: unable to send bogus packet for netkeepalive\n");
            }

            last_send_time = tv.tv_sec;
        }  
        
    }
}

