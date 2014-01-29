#ifndef __NSTREAMCLIENT_H__
#define __NSTREAMCLIENT_H__

#include <stdbool.h>
#include "exception.h"


/* header file for nstream control */

/* message types
 * X -> exit
 * S -> start record (message_s)
 * T -> stop record 
 * G -> set dac gain
 * L -> set dac channel
 * C -> set channel mapping
 * V -> verbose
 * Q -> start acquisition
 * A -> add auxdsp
 * M -> add masterdsp
 * D -> enable dataglove
 */


#define NSTREAMCLIENT_ERROR_LENGTH 512 

/*
 * this sucks and will be changed to be generic, one day...
 */
typedef struct control_message
{
    char message_type;
    char recording_path[256];
    char recording_filename_root[256];
    char ip_address[256];
    char serial_port[256];
    int data_channel;
    int hw_channel;
    int sw_channel;
    int dac_channel;
    int dac_gain;
    int high_pass;
    int low_pass;
    int enable;
    int delay;
    int comedi_num_channels_to_write;
    int nspike_num_channels_to_write_low;
    int nspike_num_channels_to_write_high;
    int nspike_high_channels_offset;
    int decimation_factor;
} control_message;

typedef struct control_response
{
    bool success;
    char error_message[NSTREAMCLIENT_ERROR_LENGTH];
} control_response;


typedef struct NStreamClient NStreamClient;

struct NStreamClient {
    Exception *ex;
};

NStreamClient *NStreamClient_create();
void NStreamClient_destroy(NStreamClient *nsc);
bool NStreamClient_add_masterdsp(NStreamClient *nsc, char *masterdsp_address);
bool NStreamClient_add_auxdsp(NStreamClient *nsc, char *auxdsp_address);
bool NStreamClient_set_channel(NStreamClient *nsc, int sw_channel, int hw_channel);
bool NStreamClient_start_record(NStreamClient *nsc, char *recording_path, char *recording_filename_root, int decimation_factor, int comedi_num_channels_to_write, int nspike_num_channels_to_write_low, int nspike_num_channels_to_write_high, int nspike_high_channels_offset);
bool NStreamClient_stop_record(NStreamClient *nsc);
bool NStreamClient_start_acquisition(NStreamClient *nsc);
bool NStreamClient_exit_nstream(NStreamClient *nsc);
bool NStreamClient_set_verbose(NStreamClient *nsc, int enable);
bool NStreamClient_set_dac_channel(NStreamClient *nsc, int dac_channel, int data_channel);
bool NStreamClient_set_dac_gain(NStreamClient *nsc, int dac_channel, int dac_gain);
bool NStreamClient_set_channel_filters(NStreamClient *nsc, int data_channel, int high_pass, int low_pass);
bool NStreamClient_enable_dataglove(NStreamClient *nsc, char *serial_port, int delay);
void NStreamClient_get_error(NStreamClient *nsc, char *buf, int len);
bool NStreamClient__send_message(NStreamClient *nsc, control_message *cm);
bool NStreamClient__send_message_timeout(NStreamClient *nsc, control_message *cm, int timeout);
bool NStreamClient__wait_socket(NStreamClient *nsc, int socket, int timeout_seconds);

#endif

