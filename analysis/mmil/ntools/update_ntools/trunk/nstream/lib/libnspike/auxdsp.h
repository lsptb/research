#ifndef __AUXDSP_H__
#define __AUXDSP_H__

#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdbool.h>
#include "hw_defines.h"
#include "libnspike.h"

struct AuxDSP {
    DSP *dsp;
    bool shutdown;   // flag for threads to terminate
    bool threads_started;
    int dsp_num;
    Channel *channels[NSPIKE_MAX_CHAN_PER_DSP];
    int samples_per_packet;
    int sample_rate;
    int first_packet_timestamp;
    int channel_count;
    int last_packet_timestamp;
    int lost_packets;
    RingBuffer *ring_buffer;
    Exception *ex;
};


AuxDSP *AuxDSP_create(int dsp_num);
bool AuxDSP_connect(AuxDSP *adsp, char *address);
bool AuxDSP_set_sample_rate(AuxDSP *adsp, int sample_rate);
int AuxDSP_get_number(AuxDSP *adsp);
int AuxDSP_get_sample_rate(AuxDSP *adsp);
bool AuxDSP_can_add_channel(AuxDSP *adsp, Channel *ch);
Channel *AuxDSP_get_channel_by_number(AuxDSP *adsp, int number);
bool AuxDSP_get_address(AuxDSP *adsp, char *ip, int len);
bool AuxDSP_add_channel(AuxDSP *adsp, Channel *ch);
int AuxDSP_get_channel_count(AuxDSP *adsp);
bool AuxDSP_all_channels_configured(AuxDSP *adsp);
void AuxDSP_set_first_timestamp(AuxDSP *adsp, int timestamp);
bool AuxDSP_set_channel(AuxDSP *adsp, Channel *ch);
int AuxDSP_get_first_timestamp(AuxDSP *adsp);
bool AuxDSP_initialize(AuxDSP *adsp);
bool AuxDSP_start_acquire(AuxDSP *adsp);
bool AuxDSP_stop_acquire(AuxDSP *adsp);
void AuxDSP_destroy(AuxDSP *adsp);
bool AuxDSP_push_data(AuxDSP *adsp, char *data, int len);
int AuxDSP_data_available(AuxDSP *adsp);
int AuxDSP_get_data(AuxDSP *adsp, char *data, int len);
void AuxDSP_add_to_lost_packet_count(AuxDSP *adsp, int number);
int AuxDSP_get_lost_packet_count(AuxDSP *adsp);




#endif

