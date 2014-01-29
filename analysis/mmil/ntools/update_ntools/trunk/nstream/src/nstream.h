#ifndef __NSTREAM_H__
#define __NSTREAM_H__
#include <stdbool.h> 

#define NSTREAM_NUM_COMEDI_CHANNELS 8 
#define NSTREAM_NUM_NSPIKE_CHANNELS 256 
#define NSTREAM_CLOCK_PFI NI_EXT_PFI(7)
#define NSTREAM_START_PFI NI_EXT_PFI(0) 
#define NSTREAM_COMEDI_DEVICE "/dev/comedi0"

#define NSTREAM_STATE_VARS_SHM_KEY 1234
#define NSTREAM_BEHAVIOR_SHM_KEY     5678
#define NSTREAM_NEURAL_SHM_KEY       8765

#define NSTREAM_SR 30000
#define NSTREAM_SS 2
#define NSTREAM_DURATION 10
#define NSTREAM_CIRCULAR_BUFFER_SIZE (NSTREAM_SR * NSTREAM_SS * NSTREAM_DURATION)
#define NSTREAM_COMEDI_BUFFER_SIZE   (NSTREAM_SR * NSTREAM_SS * NSTREAM_DURATION * NSTREAM_NUM_COMEDI_CHANNELS) 
#define NSTREAM_COMEDI_NUM_SAMPLES   (NSTREAM_SR * NSTREAM_DURATION * NSTREAM_NUM_COMEDI_CHANNELS) 
#define NSTREAM_NSPIKE_NUM_SAMPLES   (NSTREAM_SR * NSTREAM_DURATION * NSTREAM_NUM_NSPIKE_CHANNELS) 
#define NSTREAM_NSPIKE_BUFFER_SIZE   (NSTREAM_SR * NSTREAM_SS * NSTREAM_DURATION * NSTREAM_NUM_NSPIKE_CHANNELS) 
#define NSTREAM_ANALOG_STATE_NUM_SAMPLES 100
#define NSTREAM_ANALOG_STATE_DIMENSION 3
#define NSTREAM_DIO_STATE_NUM_ELEMENTS   100
#define NSTREAM_SPIKE_DETECT_DEFAULT_POS_THRESH 500
#define NSTREAM_SPIKE_DETECT_DEFAULT_NEG_THRESH -500
#define NSTREAM_SPIKE_DETECT_DEFAULT_REFRACTORY_PERIOD 30

#define NSTREAM_SPIKE_PRESAMPLE  5
#define NSTREAM_SPIKE_POSTSAMPLE  25
#define NSTREAM_SPIKE_TOTAL_SAMPLES (NSTREAM_SPIKE_PRESAMPLE + NSTREAM_SPIKE_POSTSAMPLE)
#define NSTREAM_SPIKE_NUM_ELEMENTS 100

#define NSTREAM_DATAGLOVE_NUM_SAMPLES  5000 

typedef struct spike_event
{
    int timestamp;
    short waveform[NSTREAM_SPIKE_TOTAL_SAMPLES];
    int cell_id;    
} spike_event;

typedef struct dio_event
{
    long int nspike_global_timestamp;
    unsigned short state[4];
} dio_event;

typedef struct dataglove_event
{
    long int timestamp;
    float thumb;
    float index;
    float middle;
    float ring;
    float little;
    int gesture;
} dataglove_event;

typedef struct shared_state_vars
{
    long int time;
    long int comedi_time;
    long int nspike_time;
    long int latest_timestamp;
    long int recording_start_time;  /* offset from nspike time as reported by nstream */
    long int analog_state_buffer[NSTREAM_ANALOG_STATE_DIMENSION][NSTREAM_ANALOG_STATE_NUM_SAMPLES];
    int analog_state_index;
    dio_event dio_message_buffer[NSTREAM_DIO_STATE_NUM_ELEMENTS];
    int dio_message_index;
    dataglove_event dataglove_buffer[NSTREAM_DATAGLOVE_NUM_SAMPLES];
    int dataglove_buffer_index;
    unsigned short comedi_buffer[NSTREAM_COMEDI_NUM_SAMPLES];
    short nspike_buffer[NSTREAM_NSPIKE_NUM_SAMPLES];
    unsigned short *comedi_buffer_write_ptr;
    short *nspike_buffer_write_ptr;
    unsigned short *comedi_buffer_end_ptr;
    short *nspike_buffer_end_ptr;
    spike_event spike_event_buffer[NSTREAM_NUM_NSPIKE_CHANNELS][NSTREAM_SPIKE_NUM_ELEMENTS];
    int spike_index[NSTREAM_NUM_NSPIKE_CHANNELS]; 
    short spike_detect_positive_threshold[NSTREAM_NUM_NSPIKE_CHANNELS];
    short spike_detect_negative_threshold[NSTREAM_NUM_NSPIKE_CHANNELS];
    int spike_detect_refractory_period;
} shared_state_vars;

#define NSTREAM_CONTROL_SOCKET_NUMBER 4010

void lprintf(const char *fmt, ...);

#endif
