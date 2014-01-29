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
#include "libnspike.h"

RingBuffer *RingBuffer_create(int size)
{
    RingBuffer *rb;

    if ((rb = (RingBuffer*)malloc(sizeof(RingBuffer))) == NULL)
        return false;

    memset(rb, 0, sizeof(RingBuffer));
    
    rb->data_buffer           = (char*)malloc(size);
    
    if (rb->data_buffer == NULL)
    {
        free(rb);
        return false;
    }
        
    memset(rb->data_buffer, 0, size);

    rb->size                  = size;
    rb->data_buffer_end_ptr   = rb->data_buffer + size;

    RingBuffer_reset(rb);

    return rb;
}

void RingBuffer_destroy(RingBuffer *rb)
{
    free(rb->data_buffer);
    free(rb);
}

void RingBuffer_reset(RingBuffer *rb)
{
    rb->data_buffer_write_ptr = rb->data_buffer;
    rb->data_buffer_read_ptr  = rb->data_buffer;
    rb->data_in_buffer        = 0;
}

bool RingBuffer_push_data(RingBuffer *rb, char *data, int len)
{
    int len_to_end = 0;
    
    if (len > rb->size - rb->data_in_buffer)
        return false;

    rb->data_in_buffer += len;

    /* if we can't fit this whole chunk in the buffer, calculate
     * how much space is remaining (len_to_end) fill it and reset
     * the ptr to the beginning.  */
    if (rb->data_buffer_write_ptr + len > rb->data_buffer_end_ptr)
    {
        // printf("writer wrapping\n");
        len_to_end = rb->data_buffer_end_ptr - (rb->data_buffer_write_ptr);
        memcpy(rb->data_buffer_write_ptr, data, len_to_end);
        rb->data_buffer_write_ptr = rb->data_buffer;
    }

    /* if we wrapped, only copy the remaining bytes by subtracting
     * len_to_end, if we didn't len_to_end is zero so we copy the
     * whole thing and advance the pointer */
    memcpy(rb->data_buffer_write_ptr, data + len_to_end, len - len_to_end);
    rb->data_buffer_write_ptr += len - len_to_end;

    return true;
}

int RingBuffer_data_available(RingBuffer *rb)
{
    return rb->data_in_buffer;
}

int RingBuffer_get_data(RingBuffer *rb, char *data, int len)
{
    int bytes_copied = 0;

    /* if the writer has wrapped */
    if (rb->data_buffer_write_ptr < rb->data_buffer_read_ptr)
    {
        /* if the number of bytes asked for takes us beyond the end
         * copy to the end and reset the ptr */
        if (len > rb->data_buffer_end_ptr - (rb->data_buffer_read_ptr))
        {
            memcpy(data, rb->data_buffer_read_ptr, rb->data_buffer_end_ptr - (rb->data_buffer_read_ptr));
            bytes_copied = rb->data_buffer_end_ptr - (rb->data_buffer_read_ptr);
            rb->data_buffer_read_ptr = rb->data_buffer;
            len -= bytes_copied;
        }
    }

    /* copy the number of bytes between where the reader is and where the
     * writer is up to the length requested */
    if (len > 0)
    {
        /* this will shorten the returned bytes, preventing the reader from
         * overstepping the writer */
        if (len > rb->data_buffer_write_ptr - rb->data_buffer_read_ptr &&
            rb->data_buffer_write_ptr > rb->data_buffer_read_ptr)
        {
            len = (rb->data_buffer_write_ptr - rb->data_buffer_read_ptr);
        }

        memcpy(data + bytes_copied, rb->data_buffer_read_ptr, len);
        rb->data_buffer_read_ptr += len;
        bytes_copied += len;
    }
    
    rb->data_in_buffer -= bytes_copied;
    return bytes_copied;
}


