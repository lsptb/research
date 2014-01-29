#ifndef __DAQENGINE_H__
#define __DAQENGINE_H__

#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdbool.h>
#include "hw_defines.h"
#include "libnspike.h"

struct DAQEngine {
    pthread_t reader_thread;
    bool shutdown;
    bool threads_started;
    UDPSocket *data_socket;
    AuxDSP *auxdsps[NSPIKE_MAX_AUXDSPS];
    Exception *ex;
    int state;
    int auxdsp_count;
    int first_timestamp_for_first_dsp;
    int first_dsp_index;
    long int latest_timestamp;
    int daqengine_unix_socket;
    int clientapp_unix_socket;
};

#define DAQENGINE_ACQ_STOPPED 0
#define DAQENGINE_ACQ_STARTUP 1
#define DAQENGINE_ACQ_RUN     2

DAQEngine *DAQEngine_create();
bool DAQEngine_initialize(DAQEngine *daq);
int DAQEngine_start_acquire(DAQEngine *daq);
bool DAQEngine_stop_acquire(DAQEngine *daq);
bool DAQEngine_is_acquiring(DAQEngine *daq);
int DAQEngine_get_client_socket(DAQEngine *daq);
long int DAQEngine_get_latest_timestamp(DAQEngine *daq);
void DAQEngine_register_auxdsp(DAQEngine *daq, AuxDSP *adsp);
void DAQEngine_destroy(DAQEngine *daq);
void DAQEngine__reader_thread_main(void *void_daq_ptr);
void DAQEngine__process_packet_in_startup(DAQEngine *daq, int pipeid, int timestamp, unsigned short *read_data_ptr, int data_length);
void DAQEngine__process_packet(DAQEngine *daq, int pipeid, int timestamp, unsigned short *read_data_ptr, int data_length);

#endif

