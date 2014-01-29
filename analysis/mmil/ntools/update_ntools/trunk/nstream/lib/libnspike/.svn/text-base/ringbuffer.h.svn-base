#ifndef __RINGBUFFER_H__
#define __RINGBUFFER_H__

#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdbool.h>
#include "hw_defines.h"
#include "libnspike.h"

struct RingBuffer {
    int size;
    char *data_buffer;
    char *data_buffer_end_ptr;
    char *data_buffer_write_ptr;
    char *data_buffer_read_ptr;
    int data_in_buffer;
    Exception *ex;
};

RingBuffer *RingBuffer_create(int size);
void RingBuffer_destroy(RingBuffer *rb);
bool RingBuffer_push_data(RingBuffer *rb, char *data, int len);
int RingBuffer_data_available(RingBuffer *rb);
int RingBuffer_get_data(RingBuffer *rb, char *data, int len);
void RingBuffer_reset(RingBuffer *rb);

#endif

