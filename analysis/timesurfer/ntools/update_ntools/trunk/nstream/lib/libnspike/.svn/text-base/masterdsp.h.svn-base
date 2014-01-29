#ifndef __MASTERDSP_H__
#define __MASTERDSP_H__

#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdbool.h>
#include "hw_defines.h"
#include "libnspike.h"

struct MasterDSP {
    pthread_t reader_thread;
    struct sockaddr_in address;
    DSP *master_dsp;
    pthread_t dio_notify_thread;
    cfulist_t *inbound_message_queue;
    cfulist_t *dio_message_queue;
    pthread_mutex_t inbound_queue_mutex; 
    pthread_mutex_t dio_queue_mutex; 
    pthread_mutex_t dio_state_mutex; 
    pthread_cond_t dio_messageready;
    pthread_cond_t masterdsp_messageready;
    pthread_cond_t dio_synchronous_state_ready;
    dio_status_funcptr dio_status_callback;
    void *dio_status_callback_user_data;
    bool shutdown;   // flag for threads to terminate
    bool threads_started;
    int last_ping_time;

    Exception *ex;
    int dio_port_count;
    int dio_port_type[NSPIKE_MAX_DIO_PORTS];
    int dio_raised[NSPIKE_MAX_DIO_PORTS];
    unsigned short dio_statemachineptr[NSPIKE_DIO_N_STATE_MACHINES];  // the current pointer for each state machines */
    unsigned short dio_statemachinebuffer[NSPIKE_DIO_N_STATE_MACHINES];  // the start of the memory buffer to be written to for each state machine
    short dac_gain[NSPIKE_DAC_CHANNEL_COUNT];
    short dac_delay[NSPIKE_DAC_CHANNEL_COUNT];
    short dac_cutoff[NSPIKE_DAC_CHANNEL_COUNT];
    bool dac_mute[NSPIKE_DAC_CHANNEL_COUNT];
    Channel *dac_channels[NSPIKE_DAC_CHANNEL_COUNT];

    int last_dio_timestamp;
    unsigned short last_dio_state[NSPIKE_MAX_DIO_PORTS];
};


MasterDSP *MasterDSP_create();
bool MasterDSP_connect(MasterDSP *mdsp, char *address);
bool MasterDSP_get_address(MasterDSP *mdsp, char *ip, int len);
bool MasterDSP_initialize(MasterDSP *mdsp);
void MasterDSP_destroy(MasterDSP *mdsp);
bool MasterDSP_set_dac_gain(MasterDSP *mdsp, int channel, short gain);
bool MasterDSP_set_dac_cutoff(MasterDSP *mdsp, int channel, short cutoff);
bool MasterDSP_set_dac_delay(MasterDSP *mdsp, int channel, short delay);
bool MasterDSP_set_dac_mute(MasterDSP *mdsp, int channel, bool muted);
bool MasterDSP_set_dac_channel(MasterDSP *mdsp, Channel *ch, int dac_ch_num);
bool MasterDSP_configure_dac_channel_if_selected(MasterDSP *mdsp, Channel *ch);
bool MasterDSP_set_dio_mode(MasterDSP *mdsp, int diogroup, bool output);
bool MasterDSP_set_dio_port_count(MasterDSP *mdsp, int count);
bool MasterDSP_reset_clock(MasterDSP *mdsp);
void MasterDSP_register_dio_state_change_callback(MasterDSP *mdsp, dio_status_funcptr *callback, void *user_data_ptr);
int MasterDSP__write_dio_command(MasterDSP *mdsp, unsigned short *command, int len);
int MasterDSP__dio_next_state_machine(MasterDSP *mdsp);
bool MasterDSP_change_dio_out_state(MasterDSP *mdsp, int output, int raise);
void MasterDSP__reader_thread_main( void *void_mdsp_ptr );
int MasterDSP__read_message_from_queue(UDPSocket *udps, void *buffer, int size);
int MasterDSP__blockedread_message_from_queue(UDPSocket *udps, void *buffer, int size, struct timeval *tv);
void MasterDSP__dio_thread_main( void *void_mdsp_ptr );
int MasterDSP_wait_for_dio_change(MasterDSP *mdsp, int channel_num, int state);
void MasterDSP__ping(MasterDSP *mdsp);





#endif

