
#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include <fcntl.h>
#include <errno.h>
#include "libnspike.h"
#include "udpsocket.h"
#include "fakehw.h"

fakehw *fakehw_create(char *address, int port)
{
    fakehw *fhw = (fakehw*)malloc(sizeof(fakehw));

    fail_if(!inet_pton(AF_INET, address, &fhw->address.sin_addr));

    fhw->address.sin_port   = htons(port);
    fhw->address.sin_family = AF_INET;
    fhw->address.sin_addr.s_addr = htonl(INADDR_ANY);
    fhw->peerlen = sizeof(fhw->peeraddr);

    /* create the socket */
    fail_if((fhw->socket = socket(PF_INET, SOCK_DGRAM, 0)) < 0);
    fail_if(bind(fhw->socket, (struct sockaddr *) &fhw->address, sizeof(fhw->address))==-1);
 
    /* set non-blocking */
    fcntl(fhw->socket, F_SETFL, O_NONBLOCK);

    fhw->failed = false;

    return fhw;
}

void fakehw_destroy(fakehw *fhw)
{
    close(fhw->socket);
    free(fhw);
}

int fakehw_write(fakehw *fhw, void *buf, int count)
{
    int byteswritten;

    if(count != (byteswritten = sendto(fhw->socket, buf, count, 0, &fhw->peeraddr, fhw->peerlen)))
    {
        char tmp[200];
        sprintf(fhw->last_error, "fakehw: write error: %s", strerror_r(errno, tmp, sizeof(tmp)));
        printf("%s\n", fhw->last_error);
        return -1;
    }

    return byteswritten;
}


int fakehw_blockedread(fakehw *fhw, void *buf, int count, struct timeval *tv)
{
    int bytesread;

    FD_ZERO(&fhw->fdset);
    FD_SET(fhw->socket, &fhw->fdset);

    select(fhw->socket + 1, &fhw->fdset, NULL, NULL, tv);
    if (!FD_ISSET(fhw->socket, &fhw->fdset))
    {
        return -1;
    }

    return fakehw_read(fhw, buf, count);
}

int fakehw_read(fakehw *fhw, void *buf, int count)
{
    int bytesread;
    char tmp[1024];

    if(-1 == (bytesread = recvfrom(fhw->socket, buf, count, 0, &fhw->peeraddr, &fhw->peerlen)))
    {
        sprintf(fhw->last_error, "fakehw: read error: %s", strerror_r(errno, tmp, sizeof(tmp)));
        printf("%s\n", fhw->last_error);
        return -1;
    }

    return bytesread;
}

void fakehw_start_seq(fakehw *fhw, char *name)
{
    fhw->ptr    = 0;
    strncpy(fhw->testname, name, sizeof(fhw->testname));
    memset(fhw->msgs, 0, sizeof(fhw->msgs));
}

void fakehw_add_msg(fakehw *fhw, unsigned short msg[32], int len, int type)
{
    memcpy(fhw->msgs[fhw->ptr].msg, msg, len);
    fhw->msgs[fhw->ptr].len = len;
    fhw->msgs[fhw->ptr].type = type;
    fhw->ptr++;
}

void fakehw_recv(fakehw *fhw, unsigned short msg[32], int len)
{
    fakehw_add_msg(fhw, msg, len, RECV);
}

void fakehw_send(fakehw *fhw, unsigned short msg[32], int len)
{
    fakehw_add_msg(fhw, msg, len, SEND);
}

void fakehw_run(fakehw *fhw)
{
    fhw->ptr = 0;
    pthread_create(&fhw->fakehw_thread, NULL, &fakehw_thread_main, (void*)fhw);
}

bool fakehw_check(fakehw *fhw)
{
    pthread_join(fhw->fakehw_thread, NULL);
    return fhw->failed ? false : true;
}
    

void fakehw_thread_main(void *fhw_ptr)
{
     fakehw *fhw;
     fhw = (fakehw*) fhw_ptr;
     int readsize, j, i;
     struct timeval tv;
     unsigned short readbuf[NSPIKE_MAX_PACKET_SIZE];
     char expected[256];
     char got[256];


     tv.tv_sec = 10;
     tv.tv_usec = 0;

     while(1)
     {
        if (fhw->msgs[fhw->ptr].type == RECV)
        {
            if (-1 == (readsize = fakehw_blockedread(fhw, readbuf, 
                       NSPIKE_MAX_PACKET_SIZE * sizeof(short), &tv)))
            {
                fhw->failed = true;
                sprintf(fhw->failreason, "timed out, ptr=%d", fhw->ptr);
                printf("fhw: %s\n", fhw->failreason);
                break;
            }

           // printf("fakehw got packet!\n");
            
            if (fhw->msgs[fhw->ptr].len != readsize)
            {
                fhw->failed = true;
                sprintf(fhw->failreason, "got wrong packet size: %s: ptr=%d", fhw->testname, fhw->ptr);
                printf("fhw: %s\n", fhw->failreason);
                break;
            }
 
            if (memcmp(&fhw->msgs[fhw->ptr].msg, readbuf, readsize) != 0)
            {
                j = 0; 
                for (i = 0; i < readsize; i+=2)
                {
                    j += sprintf(got + j, "0x%hx ", *(unsigned short*)(readbuf + (i/2)));
                }
    
                j = 0; 
                for (i = 0; i < fhw->msgs[fhw->ptr].len / sizeof(unsigned short); i++)
                {
                    j += sprintf(expected + j, "0x%hx ", fhw->msgs[fhw->ptr].msg[i]);
                }

                fhw->failed = true;
                sprintf(fhw->failreason, "packets don't match: %s: ptr=%d got: %s, expected: %s", 
                        fhw->testname, fhw->ptr, got, expected);
                printf("fhw: %s\n", fhw->failreason);
                break;                
            }

            // we're good!
            // printf("fhw: packet ok!\n");
        } 
        
        if (fhw->msgs[fhw->ptr].type == SEND)
        {
	        if (fakehw_write(fhw, fhw->msgs[fhw->ptr].msg, fhw->msgs[fhw->ptr].len) != fhw->msgs[fhw->ptr].len) {
                fhw->failed = true;
                sprintf(fhw->failreason, "unable to write to socket in fakehw: %s", fhw->last_error);
                printf("fhw: %s\n", fhw->failreason);
                break;
            }

            // printf("fhw: sent reply!\n");
        }

        fhw->ptr++;
        if (fhw->msgs[fhw->ptr].type == 0)
        {
            // we're at the end so break..
            break;
        }
     }
     pthread_exit(NULL);
} 


