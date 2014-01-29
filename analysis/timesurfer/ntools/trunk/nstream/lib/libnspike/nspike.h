#ifndef __NSPIKE_H__
#define __NSPIKE_H__

#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdbool.h>
#include <stdarg.h>
#include "libnspike.h"
#include "cfulist.h"
#include "hw_defines.h"

// digital io status callback
// func(timestamp, ushort[nports], user_data_ptr);
typedef void (*dio_status_funcptr)(u32, unsigned short[NSPIKE_MAX_DIO_PORTS], void*);

struct NSpike {
    struct sockaddr_in master_address;
	int dioports;
    MasterDSP *master_dsp;
    MasterDSP *secondary_master_dsp;
    cfulist_t *auxdsps;
    DAQEngine *daq_engine;
    NetKeepAlive *net_keep_alive;
    int valid;
    Exception *ex;    
};

/*
 *
 * Setup and Configuration
 *
 *
 */

/* crease new nspike object */
NSpike *NSpike_create();

/* get/set ip address for nspike master dsp */
bool NSpike_set_master_address(NSpike *nsp, char *master_ip);
bool NSpike_get_master_address(NSpike *nsp, char *ip, int len);

/* get/set number of dioports installed in master dsp */
bool NSpike_set_dioport_count(NSpike *nsp, int dioports);
int NSpike_get_num_dioports(NSpike *nsp);

/* define an available auxdsp */
bool NSpike_add_auxdsp(NSpike *nsp, char *ip_address, int sample_rate);

/* define masterdsp(s) */
bool NSpike_add_masterdsp(NSpike *nsp, char *ip_address);

/* get number of aux dsps */
int NSpike_get_num_auxdsps(NSpike *nsp);

/* set channel configuration */
bool NSpike_set_channel(NSpike *nsp, Channel *ch);

/* get a channel by user channel number */
Channel *NSpike_get_channel_by_number(NSpike *nsp, int number);

/* apply changes to a given channel */
bool NSpike_apply_channel_changes(NSpike *nsp, Channel *ch);

/* initialize hardware (must set dioports and master ip first) */
bool NSpike_configure(NSpike *nsp);

/* get data socket */
int NSpike_get_data_socket(NSpike *nsp);

/* get messages socket */
int NSpike_get_messages_socket(NSpike *nsp);

/* register dio state change callback */
void NSpike_register_dio_state_change_callback(NSpike *nsp, dio_status_funcptr *callback, void *user_data_ptr);

/* get latest timestamp from nspike */
long int NSpike_get_latest_timestamp(NSpike *nsp);

/*
 *
 * Control Functions
 *
 */

/* raise digital line high (returns timestamp) */
int NSpike_raise_dio(NSpike *nsp, int digio_num);

/* set the dac to listen to a specific channel */
bool NSpike_set_dac_channel(NSpike *nsp, int channel_number, int dac_channel_number);

/* set the dac gain */
bool NSpike_set_dac_gain(NSpike *nsp, int channel_number, int gain);

/* start data acquisition (returns first timestamp or -1 on failure) */
int NSpike_start_acquire(NSpike *nsp);

/* stop data acquisition */
bool NSpike_stop_acquire(NSpike *nsp);

/* query acquire state */
bool NSpike_is_acquiring(NSpike *nsp);

/* get last error */
void NSpike_get_error(NSpike *nsp, char *msg, int len);

/* destroy nspike object and cleanup */
void NSpike_destroy(NSpike *nsp);

// private methods
//
bool NSpike__is_valid_object(NSpike *nsp);
void NSpike__set_error(NSpike *nsp, char *fmt, ...);

#endif
