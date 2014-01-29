#include <stdbool.h>
#include <sys/types.h>
#include <sys/un.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <stdio.h>
#include <fcntl.h>

#include "nstreamclient.h"
#include "exception.h"
#include "nstream.h"

NStreamClient *NStreamClient_create()
{
    NStreamClient *nsc;

    if ((nsc = (NStreamClient*)malloc(sizeof(NStreamClient))) == NULL)
        return false;

    memset(nsc, 0, sizeof(NStreamClient));

    nsc->ex                 = Exception_create("");
    Exception_enable_print(nsc->ex, false);
    return nsc;
}


void NStreamClient_destroy(NStreamClient *nsc)
{
    Exception_destroy(nsc->ex);
    free(nsc);
}


bool NStreamClient_start_record(NStreamClient *nsc, char *recording_path, char *recording_filename_root, int decimation_factor, int comedi_num_channels_to_write, int nspike_num_channels_to_write_low, int nspike_num_channels_to_write_high, int nspike_high_channels_offset)
{
    control_message cm;
    cm.message_type = 'S';
    strncpy((char*)&cm.recording_path, recording_path, sizeof(cm.recording_path));
    strncpy((char*)&cm.recording_filename_root, recording_filename_root, sizeof(cm.recording_filename_root));
    cm.decimation_factor = decimation_factor;
    cm.comedi_num_channels_to_write = comedi_num_channels_to_write;
    cm.nspike_num_channels_to_write_low = nspike_num_channels_to_write_low;
    cm.nspike_num_channels_to_write_high = nspike_num_channels_to_write_high;
    cm.nspike_high_channels_offset = nspike_high_channels_offset;

    return NStreamClient__send_message(nsc, &cm);
}


bool NStreamClient_stop_record(NStreamClient *nsc)
{
    control_message cm;
    memset(&cm, 0, sizeof(cm));
    
    cm.message_type = 'T';

    return NStreamClient__send_message(nsc, &cm);
}


bool NStreamClient_exit_nstream(NStreamClient *nsc)
{
    control_message cm;
    memset(&cm, 0, sizeof(cm));
    
    cm.message_type = 'X';
    return NStreamClient__send_message(nsc, &cm);
}

bool NStreamClient_start_acquire(NStreamClient *nsc)
{
    control_message cm;
    memset(&cm, 0, sizeof(cm));
    
    cm.message_type = 'Q';
    return NStreamClient__send_message_timeout(nsc, &cm, 15);
}

bool NStreamClient_set_verbose(NStreamClient *nsc, int enable)
{
    control_message cm;
    memset(&cm, 0, sizeof(cm));
    
    cm.message_type = 'V';
    cm.enable = enable;

    return NStreamClient__send_message(nsc, &cm); 
}


bool NStreamClient_set_dac_channel(NStreamClient *nsc, int dac_channel, int data_channel)
{
    control_message cm;
    memset(&cm, 0, sizeof(cm));
    
    cm.message_type = 'L';
    cm.data_channel = data_channel;
    cm.dac_channel  = dac_channel;

    return NStreamClient__send_message(nsc, &cm); 
}


bool NStreamClient_set_dac_gain(NStreamClient *nsc, int dac_channel, int dac_gain)
{
    control_message cm;
    memset(&cm, 0, sizeof(cm));
    
    cm.message_type = 'G';
    cm.dac_gain     = dac_gain;
    cm.dac_channel  = dac_channel;

    return NStreamClient__send_message(nsc, &cm); 
}


bool NStreamClient_set_channel_filters(NStreamClient *nsc, int data_channel, int high_pass, int low_pass)
{
    control_message cm;
    memset(&cm, 0, sizeof(cm));
    
    cm.message_type = 'F';
    cm.data_channel = data_channel;
    cm.high_pass    = high_pass;
    cm.low_pass     = low_pass;

    return NStreamClient__send_message(nsc, &cm); 
}

bool NStreamClient_enable_dataglove(NStreamClient *nsc, char *serial_port, int delay)
{
    control_message cm;
    memset(&cm, 0, sizeof(cm));
    
    cm.message_type = 'D';
    strncpy((char*)&cm.serial_port, serial_port, sizeof(cm.serial_port));
    cm.delay = delay;
 
    return NStreamClient__send_message(nsc, &cm); 
}

bool NStreamClient_set_channel(NStreamClient *nsc, int sw_channel, int hw_channel)
{
    control_message cm;
    memset(&cm, 0, sizeof(cm));
    
    cm.message_type = 'C';
    cm.sw_channel = sw_channel;
    cm.hw_channel = hw_channel;

    return NStreamClient__send_message(nsc, &cm); 
}

bool NStreamClient_add_auxdsp(NStreamClient *nsc, char *auxdsp_address)
{
    control_message cm;
    memset(&cm, 0, sizeof(cm));
    
    cm.message_type = 'A';
    strncpy((char*)&cm.ip_address, auxdsp_address, sizeof(cm.ip_address));

    return NStreamClient__send_message(nsc, &cm); 
}

bool NStreamClient_add_masterdsp(NStreamClient *nsc, char *masterdsp_address)
{
    control_message cm;
    memset(&cm, 0, sizeof(cm));
    
    cm.message_type = 'M';
    strncpy((char*)&cm.ip_address, masterdsp_address, sizeof(cm.ip_address));

    return NStreamClient__send_message(nsc, &cm); 
}

void NStreamClient_get_error(NStreamClient *nsc, char *buf, int len)
{
    Exception_get(nsc->ex, buf, len);
}

bool NStreamClient__send_message(NStreamClient *nsc, control_message *cm)
{
    return NStreamClient__send_message_timeout(nsc, cm, 5);
}

bool NStreamClient__send_message_timeout(NStreamClient *nsc, control_message *cm, int timeout)
{
    int sock;
    int len;
    struct sockaddr_in address;
    control_response cr;
    memset(&cr, 0, sizeof(cr));
    
    if ((sock = socket(AF_INET, SOCK_DGRAM, 0)) < 0)
    {
        perror("can't create control socket");
        return false;
    }

    fcntl(sock,F_SETFL,O_NONBLOCK);
    
    memset(&address, 0, sizeof(struct sockaddr_in));
    address.sin_family = AF_INET;
    address.sin_port   = htons(NSTREAM_CONTROL_SOCKET_NUMBER);
    len = sizeof(address);

    if(inet_pton(AF_INET, "127.0.0.1", &address.sin_addr) < 0)
    {
        Exception_set_strerror(nsc->ex, "can't resolve localhost");
        return false;
    }

    /* connect */
    if (connect(sock, (struct sockaddr_in *)&address, len) < 0)
    {
        Exception_set_strerror(nsc->ex, "can't connect control socket");
        return false;
    }

    /* send the message */
    if (send(sock, cm, sizeof(control_message), 0) < 0)
    {
        Exception_set_strerror(nsc->ex, "can't send control message");
        close(sock);
        return false;
    }
    
    /* wait 2 seconds for reply */
    if (!NStreamClient__wait_socket(nsc, sock, timeout))
    {
        Exception_set(nsc->ex, "no response from engine");
        close(sock);
        return false;
    } 
   
    /* read reply */
    if ((len = read(sock, &cr, sizeof(control_response))) < 0)
    {
        Exception_set_strerror(nsc->ex, "unable to read control message response");
        close(sock);
        return false;
    }

    /* check reply */
    if (!cr.success)
    {
        Exception_set(nsc->ex, &cr.error_message);
        close(sock);
        return false;
    }
    
    close(sock);
    return true;
}


bool NStreamClient__wait_socket(NStreamClient *nsc, int socket, int timeout_seconds)
{
    int bytesread;
    fd_set fdset;

    struct timeval timeout;
    timeout.tv_sec = timeout_seconds;
    timeout.tv_usec = 0;

    FD_ZERO(&fdset);
    FD_SET(socket, &fdset);

    select(socket + 1, &fdset, NULL, NULL, &timeout);
    if (FD_ISSET(socket, &fdset))
    {
        return true;
    }

    return false;
}






