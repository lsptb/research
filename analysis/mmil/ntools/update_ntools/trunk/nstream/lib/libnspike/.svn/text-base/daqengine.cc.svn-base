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
#include <unistd.h>
#include "libnspike.h"

FILE *fp[16];

DAQEngine *DAQEngine_create()
{
    DAQEngine *daq;

    if ((daq = (DAQEngine*)malloc(sizeof(DAQEngine))) == NULL)
        return false;

    memset(daq, 0, sizeof(DAQEngine));

    daq->ex = Exception_create("DAQEngine: ");
    daq->shutdown = false;
    daq->threads_started = false;
    
    return daq;
}

void DAQEngine_register_auxdsp(DAQEngine *daq, AuxDSP *adsp)
{
    int dsp_num = AuxDSP_get_number(adsp);

    daq->auxdsps[dsp_num] = adsp;
    daq->auxdsp_count++;
}

bool DAQEngine_initialize(DAQEngine *daq)
{
    int i;
    char fname[128];
    int socket_fd[2];

    daq->data_socket = UDPSocket_create();
    if (!UDPSocket_bind(daq->data_socket, NSPIKE_DATA_PORT))
    {
        Exception_propagate(daq->data_socket->ex, daq->ex);
        UDPSocket_destroy(daq->data_socket);
        return false;
    }
    
    if (!UDPSocket_set_recv_buf_size(daq->data_socket, NSPIKE_SOCKET_RCVBUF_SIZE))
    {
        Exception_propagate(daq->data_socket->ex, daq->ex);
        UDPSocket_destroy(daq->data_socket);
        return false;
    }

    /*
    for (i = 1; i < 17; i++)
    {
        sprintf(fname, "dataout.%d", i);
        fp[i] = fopen(fname, "w");
    }
*/
    daq->state = DAQENGINE_ACQ_STOPPED;

    /* create unix data socket for passing data back to the application */
    if (socketpair(AF_UNIX, SOCK_STREAM, 0, socket_fd) != 0)
    {
        Exception_set(daq->ex, "unable to create data unix socket");
        return false;
    }

    daq->daqengine_unix_socket = socket_fd[0];
    daq->clientapp_unix_socket = socket_fd[1];
    
    /* start reader thread */
    pthread_create(&daq->reader_thread, NULL, (void* (*)(void*))&DAQEngine__reader_thread_main, (void*)daq);

    return true;
}


void DAQEngine_destroy(DAQEngine *daq)
{
    Channel *ch;

    /*
    if (daq->threads_started)
    {
        daq->shutdown = true;
        
        pthread_join(daq->reader_thread, NULL);
        pthread_join(daq->dio_notify_thread, NULL);
    }
    
    cfulist_destroy(daq->inbound_message_queue);
    cfulist_destroy(daq->dio_message_queue);
    */

    /*
    pthread_mutex_destroy(&daq->inbound_queue_mutex);
    pthread_mutex_destroy(&daq->dio_queue_mutex);
    pthread_cond_destroy(&daq->dio_messageready);
    */

    Exception_destroy(daq->ex);
    free(daq);
}

int DAQEngine_start_acquire(DAQEngine *daq)
{
    int dsp_index;
    int timestamp = -1;
    int i;

    daq->first_timestamp_for_first_dsp = 0;
    daq->first_dsp_index = 0;
    daq->latest_timestamp = 0;

    for (dsp_index = 1; dsp_index <= daq->auxdsp_count; dsp_index++)
    {
        if (!AuxDSP_start_acquire(daq->auxdsps[dsp_index]))
        {
            printf("daqengine start acquire failed, auxdsp %d\n", dsp_index);
            Exception_propagate(daq->auxdsps[dsp_index]->ex, daq->ex);
            exit(1);
            for (dsp_index = 1; dsp_index <= daq->auxdsp_count; dsp_index++)
            {
                AuxDSP_stop_acquire(daq->auxdsps[dsp_index]);
            }
            return -1;
        }
    }
    
    daq->state = DAQENGINE_ACQ_STARTUP;
    
    for (i = 0; i < 20; i++)
    {
        if (daq->state == DAQENGINE_ACQ_RUN)
        {
            timestamp = daq->first_timestamp_for_first_dsp;
            break;
        }
        usleep(250000);
    }

    return timestamp;
}

bool DAQEngine_stop_acquire(DAQEngine *daq)
{
    int dsp_index;

    for (dsp_index = 1; dsp_index <= daq->auxdsp_count; dsp_index++)
    {
        AuxDSP_stop_acquire(daq->auxdsps[dsp_index]);
        Exception_printf("Aux DSP %d: lost %d packets\n", dsp_index, AuxDSP_get_lost_packet_count(daq->auxdsps[dsp_index]));
    }

    daq->state = DAQENGINE_ACQ_STOPPED;

    return true;
}

bool DAQEngine_is_acquiring(DAQEngine *daq)
{
    if (daq->state == DAQENGINE_ACQ_RUN)
        return true;
    else
        return false;
}

int DAQEngine_get_client_socket(DAQEngine *daq)
{
    return daq->clientapp_unix_socket;
}

long int DAQEngine_get_latest_timestamp(DAQEngine *daq)
{
    return daq->latest_timestamp;
}

int DAQEngine_get_first_timestamp(DAQEngine *daq)
{
    return daq->first_timestamp_for_first_dsp;
}

void DAQEngine__reader_thread_main(void *void_daq_ptr)
{
    DAQEngine *daq;
    daq = (DAQEngine*) void_daq_ptr;

    fd_set rfds;
    struct timeval timeout;
    unsigned short readbuf[NSPIKE_MAX_PACKET_SIZE];
    unsigned short *message;
    int readsize;
    char tmp[256];
    int pipeid;
    u32 timestamp;
    int data_length;
    unsigned short *read_data_ptr;

    for(;;)
    {
        // printf(">>>>>>>READER LOOP\n\n");
        if (daq->shutdown)
        {
            pthread_exit(NULL);
        }

        timeout.tv_usec = 0; // 50ms
        timeout.tv_sec  = 1;
        if (-1 == (readsize = UDPSocket_blockedread(daq->data_socket, readbuf, NSPIKE_MAX_PACKET_SIZE * sizeof(short), &timeout)))
        {
            Exception_get(daq->data_socket->ex, tmp, sizeof(tmp));
            Exception_printf("udp read error: %s\n", tmp);
            continue;
        }

        if (readsize == 0)
            continue;
       
        DSP__byte_swap(readbuf, readsize / sizeof(unsigned short));
      
        pipeid = readbuf[0];
        
        memcpy(&timestamp, readbuf + 1, 2 * sizeof(unsigned short));
        if (timestamp < 50)
            Exception_printf("WARNING: timestamp less than 50 %u\n",timestamp); 
       
        data_length = readsize - 3 * sizeof(unsigned short);
        read_data_ptr = readbuf + 3;
    
        if (daq->state == DAQENGINE_ACQ_STARTUP)
        {
            DAQEngine__process_packet_in_startup(daq, pipeid, timestamp, read_data_ptr, data_length); 
        }
        
        else if (daq->state == DAQENGINE_ACQ_RUN)
        {
            // printf(".");
            DAQEngine__process_packet(daq, pipeid, timestamp, read_data_ptr, data_length);
        }

        // fwrite(&timestamp, 1, sizeof(unsigned short) * 2, fp[pipeid]);
        // fwrite(readbuf + 3, 1, readsize - 3 * sizeof(unsigned short), fp[pipeid]);
    }
}

void DAQEngine__process_packet_in_startup(DAQEngine *daq, int pipeid, int timestamp, unsigned short *read_data_ptr, int data_length)
{
    // memcpy(daq->auxdsps[pipeid]->first_packet, read_data_ptr, data_length);

    daq->auxdsps[pipeid]->last_packet_timestamp = timestamp;
    AuxDSP_set_first_timestamp(daq->auxdsps[pipeid], timestamp);

    int dsp_index;
    int last_timestamp          = INT_MAX;
    int last_timestamp_dsp      = 0;
    int first_timestamp_for_dsp = 0;
    int timestamp_difference    = 0;
    short zeros[256000];
    memset(&zeros, 0, sizeof(zeros));
    
    // now check and see if all dsps have received their first packet
    // once they have, set the global first timestamp so we can start acq
    for (dsp_index = 1; dsp_index <= daq->auxdsp_count; dsp_index++)
    {
        first_timestamp_for_dsp = AuxDSP_get_first_timestamp(daq->auxdsps[dsp_index]);

        // haven't received packets for all dsps yet
        if (first_timestamp_for_dsp == 0)
        {
            return;
        }

        if (first_timestamp_for_dsp < last_timestamp)
        {
            last_timestamp     = first_timestamp_for_dsp;
            last_timestamp_dsp = dsp_index;
        }
    }
    
    daq->first_timestamp_for_first_dsp = last_timestamp;
    daq->first_dsp_index               = last_timestamp_dsp;
    
    for (dsp_index = 1; dsp_index <= daq->auxdsp_count; dsp_index++)
    {
        if (dsp_index != daq->first_dsp_index)
        {
            timestamp_difference = AuxDSP_get_first_timestamp(daq->auxdsps[dsp_index]) - daq->first_timestamp_for_first_dsp;
//          printf("dsp: %d timestamp diff: %d\n", dsp_index, timestamp_difference);
            AuxDSP_push_data(daq->auxdsps[dsp_index], (char*)zeros, timestamp_difference * sizeof(short) * 16);
            // AuxDSP_push_data(daq->auxdsps[dsp_index], daq->auxdsps[dsp_index]->first_packet, 928);
        } 
    } 
        
        
    daq->state                        = DAQENGINE_ACQ_RUN;


    
    // printf("all packets received.  first timestamp = %d on dsp %d\n",
    //        daq->first_timestamp_for_first_dsp, daq->first_dsp_index);
    return;
}


short buffer[17][29][16];

// switch these lines to toggle between 128 and 256 channels 
// short unrolled_buffer2[29][128];
short unrolled_buffer2[29][256];

void DAQEngine__process_packet(DAQEngine *daq, int pipeid, int timestamp, unsigned short *read_data_ptr, int data_length)
{
    bool data_on_all_dsps = true;
    int dsp_index;
    int bytes_read = 0;
    int ch_index = 0;
    int short_index = 0;
    int index = 0;
    char zeros[928];
    memset(zeros, 0, sizeof(zeros));

    if (daq->auxdsps[pipeid]->last_packet_timestamp + 29 != timestamp)
    {
        Exception_printf("PACKET LOSS: dsp %d, expected %d, got %d\n", pipeid, daq->auxdsps[pipeid]->last_packet_timestamp + 29, timestamp);
        int packets_lost = (timestamp - daq->auxdsps[pipeid]->last_packet_timestamp) / 29;
        
        AuxDSP_add_to_lost_packet_count(daq->auxdsps[pipeid], packets_lost);

        if (daq->auxdsps[pipeid]->last_packet_timestamp + 29 > timestamp)
            Exception_printf("PACKET LOSS: packets out of order!\n");

        for (index = 0; index < packets_lost - 1; index++)
        {
            if (!AuxDSP_push_data(daq->auxdsps[pipeid], zeros, 928))
            {
                for (;;)
                {
                    Exception_printf("FATAL ERROR: auxdsp data ringbuffer overflow in packetloss handling, restart!\n");
                    sleep(1);
                } 
            }
        }
        Exception_printf("WARNING: packet loss, dsp %d: %d packets zeroed out\n",pipeid, packets_lost);
    }

    daq->auxdsps[pipeid]->last_packet_timestamp = timestamp;
    
    // printf("%d %d\n",pipeid,timestamp);
    AuxDSP_push_data(daq->auxdsps[pipeid], (char*)read_data_ptr, data_length);
    
    for (dsp_index = 1; dsp_index <= daq->auxdsp_count; dsp_index++)
    {
        if (AuxDSP_data_available(daq->auxdsps[dsp_index]) <= 928)
        {
            data_on_all_dsps = false;
            continue;
        }
    }

    if (data_on_all_dsps == true)
    {
        for (dsp_index = 1; dsp_index <= daq->auxdsp_count; dsp_index++)
        {
            if ((bytes_read = AuxDSP_get_data(daq->auxdsps[dsp_index], (char*)&buffer[dsp_index], 928)) == 0)
                Exception_printf("no data!\n");
        }
        
        for (short_index = 0; short_index <= 28; short_index++)
        {
            // toggle these lines to switch between 128 and 256 channels 
            /* for (dsp_index = 1; dsp_index <= daq->auxdsp_count - 8; dsp_index++) */
            for (dsp_index = 1; dsp_index <= daq->auxdsp_count; dsp_index++)
            {
                for (ch_index = 0; ch_index < daq->auxdsp_count; ch_index++)
                {
                    index = (((dsp_index - 1) * 16) + ch_index);
                    unrolled_buffer2[short_index][index] = buffer[dsp_index][short_index][ch_index];
                }
            }
        }

        if (timestamp > daq->latest_timestamp)
            daq->latest_timestamp = timestamp;

        write(daq->daqengine_unix_socket, &unrolled_buffer2, sizeof(unrolled_buffer2));
    } 
}





                

