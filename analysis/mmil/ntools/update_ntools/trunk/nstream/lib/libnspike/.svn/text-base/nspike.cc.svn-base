#include <stdio.h>
#include <stdlib.h>
#include <arpa/inet.h>
#include <fcntl.h>
#include <string.h>
#include "libnspike.h"
#include <sys/socket.h>
#include <sys/un.h>
#include <math.h>

void testfunc(int timestamp, unsigned short data[4], void *x)
{
    int pnum, i;

    Exception_printf("libnspike debug: dio state change: ts=%d ",timestamp);
    for (pnum = 0; pnum < 2; pnum++) {
        for (i = 0; i < 16; i++) {
            if (data[pnum] & (1 << 15-i)) {
                Exception_printf("1");
            }
            else {
                Exception_printf("0");
            }
        }
    }
    Exception_printf("\n");
}
        



NSpike *NSpike_create()
{
    NSpike *nsp;
   
    if ((nsp = (NSpike*)malloc(sizeof(NSpike))) == NULL)
        return false;

    memset(nsp, 0, sizeof(NSpike));
    
    nsp->auxdsps        = cfulist_new();
    nsp->daq_engine     = DAQEngine_create();
    // nsp->net_keep_alive = NetKeepAlive_create();
    nsp->ex             = Exception_create("");

    // EXCEPTION_OUTPUT_PIPE_fp = NULL;

    return nsp;
}


bool NSpike_set_master_address(NSpike *nsp, char *master_address)
{
    if (nsp->master_dsp)
        MasterDSP_destroy(nsp->master_dsp);
        
    nsp->master_dsp = MasterDSP_create();
    
    if (!MasterDSP_connect(nsp->master_dsp, master_address))
    {
        Exception_propagate(nsp->master_dsp->ex, nsp->ex);
        MasterDSP_destroy(nsp->master_dsp);
        nsp->master_dsp = NULL;
        return false;
    }

    return true;
}

bool NSpike_add_masterdsp(NSpike *nsp, char *master_address)
{
    MasterDSP **new_masterdsp;

    if (nsp->master_dsp == NULL)
        new_masterdsp = &nsp->master_dsp;
    else if (nsp->secondary_master_dsp == NULL)
        new_masterdsp = &nsp->secondary_master_dsp;
    else
    {
        Exception_set(nsp->ex, "Both primary and secondary MasterDSP units already configured");
        return false;
    }

    *new_masterdsp = MasterDSP_create();
    
    if (!MasterDSP_connect(*new_masterdsp, master_address))
    {
        MasterDSP *tmp;
        tmp = *new_masterdsp;
        Exception_propagate(tmp->ex, nsp->ex);
        MasterDSP_destroy(*new_masterdsp);
        *new_masterdsp = NULL;
        return false;
    }

    return true;
}


bool NSpike_set_dioport_count(NSpike *nsp, int dioports)
{
    nsp->dioports = dioports;
    return true;
}

bool NSpike_add_auxdsp(NSpike *nsp, char *ip_address, int sample_rate)
{
    /* argument to create is the "number" for this dsp */
    AuxDSP *adsp = AuxDSP_create(NSpike_get_num_auxdsps(nsp) + 1);

    if (!AuxDSP_set_sample_rate(adsp, sample_rate))
    {
        Exception_set(nsp->ex, "Unable to set sample rate for AuxDSP to %d", sample_rate);
        goto FAIL;
    }

    if (!AuxDSP_connect(adsp, ip_address))
    {
        Exception_set(nsp->ex, "Unable to attach AuxDSP to IP %s", ip_address);
        goto FAIL;
    }

    if (!cfulist_push(nsp->auxdsps, adsp))
    {
        Exception_set(nsp->ex, "Internal Error: Unable to add AuxDSP to internal list");
        goto FAIL;
    }
            
    return true; 

    /* on failure destroy the dsp and return false */
    FAIL:
    AuxDSP_destroy(adsp);
    return false;    
}

int NSpike_get_num_auxdsps(NSpike *nsp)
{
    return cfulist_num_entries(nsp->auxdsps);
}

bool NSpike_set_channel(NSpike *nsp, Channel *ch)
{
    int new_chan_number = Channel_get_number(ch);
    int dsp_count = NSpike_get_num_auxdsps(nsp);
    int dsp_num = floor((new_chan_number - 1) / NSPIKE_MAX_CHAN_PER_DSP);
    
    if (dsp_num >= dsp_count)
    {   
        Exception_set(nsp->ex, "Only %d AuxDSPs making for a max of %d channels configured.  Channel %d invalid.", 
                      dsp_count, dsp_count * NSPIKE_MAX_CHAN_PER_DSP, new_chan_number); 
        return false;
    }

    AuxDSP *adsp;
    size_t dummy;
    cfulist_nth_data(nsp->auxdsps, (void**)&adsp, &dummy, dsp_num);
 
    if (!AuxDSP_set_channel(adsp, ch))
    {
        Exception_propagate_with_text(adsp->ex, nsp->ex, "Failed setting channel: ");
        return false;
    }
   
    return true;
} 
    

Channel *NSpike_get_channel_by_number(NSpike *nsp, int number)
{
    size_t dummy;
    AuxDSP *adsp;
    Channel *ch;

    /* iterate through auxdsps looking for one that
     * owns this channel */
    
    cfulist_reset_each(nsp->auxdsps);

    while (cfulist_next_data(nsp->auxdsps, (void**)&adsp, &dummy))
    {
        if (ch = AuxDSP_get_channel_by_number(adsp, number))
            return ch;
    }

    Exception_set(nsp->ex, "Error: Unknown channel %d", number);
    return NULL;
}

/* apply changes to a given channel */
bool NSpike_apply_channel_changes(NSpike *nsp, Channel *ch)
{
    MasterDSP **master_dsp;
    int ch_num = Channel_get_number(ch);
    
    if (ch_num > 127)
        master_dsp = &nsp->secondary_master_dsp;
    else
        master_dsp = &nsp->master_dsp;    
    
    if (!MasterDSP_configure_dac_channel_if_selected(*master_dsp, ch))
    {
        Exception_set(nsp->ex, "Error: Unable to update channel %d settings for DAC", Channel_get_number(ch));
        return false;
    }

    if (!Channel_configure_on_home_auxdsp(ch))
    {
        Exception_set(nsp->ex, "Error: Unable to update channel %d settings", Channel_get_number(ch));
        return false;
    }

    return true;
}

/* set the dac to listen to a specific channel */
bool NSpike_set_dac_channel(NSpike *nsp, int channel_number, int dac_channel_number)
{
    Channel *ch;
    MasterDSP **master_dsp;

    if (!(ch = NSpike_get_channel_by_number(nsp, channel_number)))
        return false;

    if (dac_channel_number > 1)
    {
        dac_channel_number -= 2;
        master_dsp = &nsp->secondary_master_dsp;
    }
    else
        master_dsp = &nsp->master_dsp;

    if (!MasterDSP_set_dac_channel(*master_dsp, ch, dac_channel_number))
        return false;

    return true;
}

/* set the dac gain */
bool NSpike_set_dac_gain(NSpike *nsp, int channel_number, int gain)
{

    MasterDSP **master_dsp;

    if (channel_number > 1)
    {
        channel_number -= 2;
        master_dsp = &nsp->secondary_master_dsp;
    }
    else
        master_dsp = &nsp->master_dsp;

    if (!MasterDSP_set_dac_gain(*master_dsp, channel_number, gain))
        return false;

    return true;
}

bool NSpike_configure(NSpike *nsp)
{
    if (!nsp->master_dsp)
    {
        Exception_set(nsp->ex, "Cannot configure hardware until master_address is set");
        return false;
    }

    if (nsp->dioports == 0)
    {
        Exception_set(nsp->ex, "Cannot configure hardware until dioports is set");
        return false;
    }

    MasterDSP_set_dio_port_count(nsp->master_dsp, 2);
    MasterDSP_set_dio_mode(nsp->master_dsp, 0, true);
    MasterDSP_set_dio_mode(nsp->master_dsp, 1, false);
    MasterDSP_set_dio_mode(nsp->master_dsp, 2, true);
    MasterDSP_set_dio_mode(nsp->master_dsp, 3, true);

   
    // MasterDSP_register_dio_state_change_callback(nsp->master_dsp, (dio_status_funcptr*)&testfunc, NULL);
    
    if (!MasterDSP_initialize(nsp->master_dsp))
    {
        Exception_set(nsp->ex, "MasterDSP initialization failed.");
        return false;
    }

    /* initialize secondary master dsp if present */
    if (nsp->secondary_master_dsp)
    {
        MasterDSP_set_dio_port_count(nsp->master_dsp, 2);
        MasterDSP_set_dio_mode(nsp->master_dsp, 0, true);
        MasterDSP_set_dio_mode(nsp->master_dsp, 1, false);
        MasterDSP_set_dio_mode(nsp->master_dsp, 2, true);
        MasterDSP_set_dio_mode(nsp->master_dsp, 3, true);

        if (!MasterDSP_initialize(nsp->secondary_master_dsp))
        {
            Exception_set(nsp->ex, "Secondary MasterDSP initialization failed.");
            return false;
        }
    }
    
    size_t dummy;
    AuxDSP *adsp;

    cfulist_reset_each(nsp->auxdsps);

    while (cfulist_next_data(nsp->auxdsps, (void**)&adsp, &dummy))
    {
        if (!AuxDSP_all_channels_configured(adsp))
        {
            Exception_propagate(adsp->ex, nsp->ex);
            return false;
        }
    }
 
    cfulist_reset_each(nsp->auxdsps);

    while (cfulist_next_data(nsp->auxdsps, (void**)&adsp, &dummy))
    {
        AuxDSP_initialize(adsp);
        DAQEngine_register_auxdsp(nsp->daq_engine, adsp);
    }
    DAQEngine_initialize(nsp->daq_engine);
 
   /* 
    printf(">>>>>>>>>>DIO ON!\n\n");
    int i; 
    for (i=0; i<16; i++)
    { 
     MasterDSP_change_dio_out_state(nsp->master_dsp, i, 1);
     sleep(100);
    }
    */
   
    /*
    char address[128];
    MasterDSP_get_address(nsp->master_dsp, (char*)&address, sizeof(address));
    */ 

    // temporarily disabled due to nspike firmware bug - if we use it's echo facility
    // it will start sending digital io status messages to the wrong port!
    //
    //NetKeepAlive_run(nsp->net_keep_alive, address, NSPIKE_DSP_ECHO_PORT, NSPIKE_DSP_ECHO_INTERVAL); 
    return true;
}


int NSpike_get_data_socket(NSpike *nsp)
{
    return DAQEngine_get_client_socket(nsp->daq_engine);
}

int NSpike_get_messages_socket(NSpike *nsp)
{
    int socket_fd[2];
    /* create unix data socket for passing data back to the application */
    if (socketpair(AF_UNIX, SOCK_STREAM, 0, socket_fd) != 0)
    {
        Exception_set(nsp->ex, "unable to create data unix socket");
        return false;
    }

    // EXCEPTION_OUTPUT_PIPE_fp = fdopen(socket_fd[0], "w");
    return socket_fd[1];     
}

void NSpike_register_dio_state_change_callback(NSpike *nsp, dio_status_funcptr *callback, void *user_data_ptr)
{
    MasterDSP_register_dio_state_change_callback(nsp->master_dsp, callback, user_data_ptr);
}

long int NSpike_get_latest_timestamp(NSpike *nsp)
{
    return DAQEngine_get_latest_timestamp(nsp->daq_engine);
}

int NSpike_raise_dio(NSpike *nsp, int digio_num)
{
    int timestamp;

    if (!MasterDSP_change_dio_out_state(nsp->master_dsp, digio_num, 1))
    {
        Exception_propagate(nsp->daq_engine->ex, nsp->ex);
        return -1;
    }

    if ((timestamp = MasterDSP_wait_for_dio_change(nsp->master_dsp, 7, 1)) == -1)
    {
        Exception_set(nsp->ex, "Timed out waiting for dio state change");
    }
    
    if (!MasterDSP_change_dio_out_state(nsp->master_dsp, digio_num, 0))
    {
        Exception_propagate(nsp->daq_engine->ex, nsp->ex);
        return -1;
    }

    return timestamp;
}

/* start data acquisition */
int NSpike_start_acquire(NSpike *nsp)
{
    int timestamp = 0;
    if(!MasterDSP_reset_clock(nsp->master_dsp))
    {
       Exception_propagate(nsp->daq_engine->ex, nsp->ex);
       return -1;
    }

    if ((timestamp = DAQEngine_start_acquire(nsp->daq_engine)) < 0)
    {
       Exception_propagate(nsp->daq_engine->ex, nsp->ex);
       return -1;
    }

    return timestamp;
}

/* stop data acquisition */
bool NSpike_stop_acquire(NSpike *nsp)
{
    if (!DAQEngine_stop_acquire(nsp->daq_engine))
    {
       Exception_propagate(nsp->daq_engine->ex, nsp->ex);
       return false;
    }

    return true;
}

/* query acquire state */
bool NSpike_is_acquiring(NSpike *nsp)
{
    return DAQEngine_is_acquiring(nsp->daq_engine);
}


bool NSpike_get_master_address(NSpike *nsp, char *ip, int len)
{
    if (!nsp->master_dsp || !MasterDSP_get_address(nsp->master_dsp, ip, len))
    {
        Exception_set(nsp->ex, "Master IP address not set");
        return false;
    }

    return true;
}


void NSpike_get_error(NSpike *nsp, char *msg, int len)
{
    if (NSpike__is_valid_object(nsp))
        Exception_get(nsp->ex, msg, len);
    else
        snprintf(msg, len, "NSpike object invalid or not created");
}


int NSpike_get_num_dioports(NSpike *nsp)
{
    return nsp->dioports;
}


void NSpike_destroy(NSpike *nsp)
{
    if (nsp->master_dsp)
        MasterDSP_destroy(nsp->master_dsp);

    if (nsp->auxdsps)
        cfulist_destroy(nsp->auxdsps);
    
    if (nsp->net_keep_alive)
        NetKeepAlive_destroy(nsp->net_keep_alive);

    Exception_destroy(nsp->ex);
    free(nsp);
}


bool NSpike__is_valid_object(NSpike *nsp)
{
    if (nsp == NULL)
        return false;

    return true;
} 




