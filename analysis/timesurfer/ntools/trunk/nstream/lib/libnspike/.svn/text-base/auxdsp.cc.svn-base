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


AuxDSP *AuxDSP_create(int dsp_num)
{
    AuxDSP *adsp;

    if ((adsp = (AuxDSP*)malloc(sizeof(AuxDSP))) == NULL)
        return false;

    memset(adsp, 0, sizeof(AuxDSP));

    adsp->ex = Exception_create("Aux DSP %d: ", dsp_num);
    adsp->shutdown = false;
    adsp->threads_started = false;

    adsp->ring_buffer = RingBuffer_create(NSPIKE_AUXDSP_DATA_BUFFER_SIZE);

    /* set some reasonable defaults */
    adsp->sample_rate = 30000;
   
    adsp->dsp_num = dsp_num;
    
    return adsp;
}

bool AuxDSP_set_sample_rate(AuxDSP *adsp, int sample_rate)
{
    adsp->sample_rate = sample_rate;
    return true;
}

int AuxDSP_get_sample_rate(AuxDSP *adsp)
{
    return adsp->sample_rate;
}

int AuxDSP_get_number(AuxDSP *adsp)
{
    return adsp->dsp_num;
}

bool AuxDSP_set_channel(AuxDSP *adsp, Channel *ch)
{
    int slot_num = ((Channel_get_number(ch) -1) % NSPIKE_MAX_CHAN_PER_DSP);
    
    /* the hardware seems to get pissed if you have two copies of the
     * same channel on the same dsp */
    int i;
    Channel *ch2;
    for (i = 0; i < NSPIKE_MAX_CHAN_PER_DSP; i++)
    {
        ch2 = (Channel*)adsp->channels[i];
        if (ch2 == NULL)
            continue;
        if (Channel_get_hw_number(ch2) == Channel_get_hw_number(ch) && i != slot_num)
        {
            Exception_set(adsp->ex, "This DSP already handles a copy of hw_channel %d, a second copy must be configured on an alternate DSP (channel number < %d or > %d)",
                          Channel_get_hw_number(ch), ((adsp->dsp_num - 1) * NSPIKE_MAX_CHAN_PER_DSP), (adsp->dsp_num * NSPIKE_MAX_CHAN_PER_DSP));
            return false;
        }
    }
    Channel_set_dsp_slot_num(ch, slot_num);
    
    if (adsp->channels[slot_num] != NULL)
    {
        Channel *old_chan = adsp->channels[slot_num];
        Channel_destroy(old_chan);
    }
        
    adsp->channels[slot_num] = ch;

    Channel_set_auxdsp(ch, adsp);

    adsp->channel_count++;
    adsp->samples_per_packet = (short) (NSPIKE_MAX_PACKET_SIZE - 3 * sizeof(unsigned short)) / AuxDSP_get_channel_count(adsp) - 1;
    
    return true;
}

/* does this dsp use the right sample rate and does it have space for
 * another channel */
bool AuxDSP_can_add_channel(AuxDSP *adsp, Channel *ch)
{
    if (AuxDSP_get_channel_count(adsp) >= NSPIKE_MAX_CHAN_PER_DSP)
        return false;

    if (Channel_get_sample_rate(ch) != AuxDSP_get_sample_rate(adsp))
        return false;

    return true;
}

Channel *AuxDSP_get_channel_by_number(AuxDSP *adsp, int number)
{
    Channel *ch;
    int i;

    for (i = 0; i < NSPIKE_MAX_CHAN_PER_DSP; i++)
    {
        ch = (Channel*)adsp->channels[i];
        if (ch == NULL)
            continue;
        if (Channel_get_number(ch) == number)
            return ch;
    }
    return NULL;
}

int AuxDSP_get_channel_count(AuxDSP *adsp)
{
    return adsp->channel_count;
}

bool AuxDSP_all_channels_configured(AuxDSP *adsp)
{
    Channel *ch;
    int i;

    for (i = 0; i < NSPIKE_MAX_CHAN_PER_DSP; i++)
    {
        ch = (Channel*)adsp->channels[i];
        if (ch == NULL)
        {
            Exception_set(adsp->ex, "Channel %d not configured, unable to start.", (adsp->dsp_num - 1) * NSPIKE_MAX_CHAN_PER_DSP + i);
            return false;
        }
    }
    return true;
}

void AuxDSP_set_first_timestamp(AuxDSP *adsp, int timestamp)
{
    // printf("AuxDSP %d: first timestamp: %d\n", adsp->dsp_num, timestamp);
    adsp->first_packet_timestamp = timestamp;
}

int AuxDSP_get_first_timestamp(AuxDSP *adsp)
{
    return adsp->first_packet_timestamp;
}

bool AuxDSP_connect(AuxDSP *adsp, char *address)
{
    UDPSocket *socket;
    socket = UDPSocket_create();
    if (!UDPSocket_connect(socket, address, NSPIKE_MESSAGE_PORT))
    {
        Exception_propagate(socket->ex, adsp->ex);
        UDPSocket_destroy(socket);
        return false;
    }

    // UDPSocket_set_read_func(socket, &AuxDSP__read_message_from_queue, adsp);

    adsp->dsp = DSP_create(0, socket);
            
    return true;
}

bool AuxDSP_initialize(AuxDSP *adsp)
{
    unsigned short data[24];

    /* turn off acq in case it's on from an unclean shutdown */
    if (!AuxDSP_stop_acquire(adsp))
        return false;

    /* set the PipeID */
    data[0] = adsp->dsp_num;
    if(!DSP_write_data(adsp->dsp, NSPIKE_PIPE_ID_ADDR, 1, data))
        return false;
    // printf("libnspike: init DSP %d: ",adsp->dsp_num);

    /* set the number of channels per sample */
    data[0] = (unsigned short) AuxDSP_get_channel_count(adsp);
    if(!DSP_write_data(adsp->dsp, NSPIKE_NUM_CHAN_ADDR, 1, data))
        return false;
    // printf("num_chan=%d ",AuxDSP_get_channel_count(adsp));

    /* set the number of samples per block */
    data[0] = adsp->samples_per_packet;
    if(!DSP_write_data(adsp->dsp, NSPIKE_NUM_SAMPLES_ADDR, 1, data))
        return false;
    //printf("samples_packet=%d ",adsp->samples_per_packet);

    /* set the decimation factor */
    data[0] = NSPIKE_DSP_BASE_SAMP_RATE / adsp->sample_rate;
    if(!DSP_write_data(adsp->dsp, NSPIKE_DECIMATION_ADDR, 1, data))
        return false;
    // printf("decimate=%d\n", NSPIKE_DSP_BASE_SAMP_RATE / adsp->sample_rate);

    /* iterate through this dsp's channels and configure them */
    Channel *ch;
    int i;

    for (i = 0; i < AuxDSP_get_channel_count(adsp); i++)
    {
        ch = (Channel*)adsp->channels[i];
        if (!Channel_configure_on_home_auxdsp(ch))
            return false;
    }

    return true;
}

bool AuxDSP_start_acquire(AuxDSP *adsp)
{
    /* reset the state of this auxdsp */
    adsp->first_packet_timestamp = 0;
    adsp->last_packet_timestamp  = 0;
    adsp->lost_packets           = 0;
    RingBuffer_reset(adsp->ring_buffer);

    unsigned short data[1];

    /* set to a number > SHRT_MAX to tell the dsp to send packets forever */
    data[0] = USHRT_MAX;
    
    if (!DSP_write_data(adsp->dsp, NSPIKE_BLOCKS_TO_SEND_ADDR, 1, data)) {
        Exception_propagate_with_text(adsp->dsp->ex, adsp->ex, "Failed to start acq: ");
	    return false;
    }

    return true;
}

bool AuxDSP_stop_acquire(AuxDSP *adsp)
{
    unsigned short data[1];

    /* set to 0 to tell dsp to stop */
    data[0] = 0;
    
    if (!DSP_write_data(adsp->dsp, NSPIKE_BLOCKS_TO_SEND_ADDR, 1, data)) {
        Exception_propagate_with_text(adsp->dsp->ex, adsp->ex, "Failed to stop acq: ");
	    return false;
    }

    return true;
}

bool AuxDSP_get_address(AuxDSP *adsp, char *ip, int len)
{
    if (!DSP_get_address(adsp->dsp, ip, len))
    {
        Exception_propagate(adsp->dsp->ex, adsp->ex);
        return false;
    }

    return true;
}

void AuxDSP_destroy(AuxDSP *adsp)
{
    Channel *ch;
    int i;

    for (i = 0; i < NSPIKE_MAX_CHAN_PER_DSP; i++)
    {
        ch = (Channel*)adsp->channels[i];
        if (ch == NULL)
            continue;
        Channel_destroy(ch);
    }

    if (adsp->dsp)
        DSP_destroy(adsp->dsp);

    if (adsp->ring_buffer)
        RingBuffer_destroy(adsp->ring_buffer);

    free(adsp);
}

bool AuxDSP_push_data(AuxDSP *adsp, char *data, int len)
{
    if (RingBuffer_data_available(adsp->ring_buffer) + len > NSPIKE_AUXDSP_DATA_BUFFER_SIZE)
    {
        Exception_printf("***** WARNING: auxdsp data ringbuffer overwrite\n");
        return false;
    }
        
    return RingBuffer_push_data(adsp->ring_buffer, data, len);
}

int AuxDSP_data_available(AuxDSP *adsp)
{
    return RingBuffer_data_available(adsp->ring_buffer);
}

int AuxDSP_get_data(AuxDSP *adsp, char *data, int len)
{
    return RingBuffer_get_data(adsp->ring_buffer, data, len);
} 
    
void AuxDSP_add_to_lost_packet_count(AuxDSP *adsp, int number)
{
    adsp->lost_packets += number;
}

int AuxDSP_get_lost_packet_count(AuxDSP *adsp)
{
    return adsp->lost_packets;
}


        


