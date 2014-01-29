#ifndef __CHANNEL_H__
#define __CHANNEL_H__

#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdbool.h>
#include "hw_defines.h"
#include "libnspike.h"

struct Channel {
    int number;
    short hw_number;
    int sample_rate;
    int dsp_slot_num;
    unsigned short highpass_filter;
    unsigned short lowpass_filter;
    short calibration_factor;
    Channel *reference_channel;
    AuxDSP *auxdsp;
    Exception *ex;
};


Channel *Channel_create();

void Channel_set_auxdsp(Channel *ch, AuxDSP *adsp);
AuxDSP *Channel_get_auxdsp(Channel *ch);

void Channel_set_number(Channel *ch, int number);
int Channel_get_number(Channel *ch);

void Channel_set_dsp_slot_num(Channel *ch, int number);
int Channel_get_dsp_slot_num(Channel *ch);

void Channel_set_hw_number(Channel *ch, int number);
int Channel_get_hw_number(Channel *ch);

bool Channel_set_sample_rate(Channel *ch, int sample_rate);
int Channel_get_sample_rate(Channel *ch);

void Channel_set_calibration_factor(Channel *ch, short factor);
short Channel_get_calibration_factor(Channel *ch);

/* defaults to ground, to set ground, use NULL */
void Channel_set_reference_channel(Channel *ch, Channel *ref_ch);
Channel *Channel_get_reference_channel(Channel *ch);

void Channel_set_highpass_filter(Channel *ch, short value);
short Channel_get_highpass_filter(Channel *ch);

void Channel_set_lowpass_filter(Channel *ch, short value);
short Channel_get_lowpass_filter(Channel *ch);

/* channel to configure, dsp to configure it on, channel slot on that dsp */
bool Channel_configure_on_dsp(Channel *ch, DSP *dsp, int dsp_slot_num);

/* configure the channel on it's home dsp */
bool Channel_configure_on_home_auxdsp(Channel *ch);

void Channel_destroy(Channel *ch);

void Channel__get_low_filter_coeff(Channel *ch, unsigned short *data);
void Channel__get_high_filter_coeff(Channel *ch, unsigned short *data);

#endif

