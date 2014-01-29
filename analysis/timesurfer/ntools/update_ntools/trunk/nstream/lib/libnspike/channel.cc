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
#include <math.h>
#include <limits.h>
#include "channel.h"
#include "dsp.h"

Channel *Channel_create()
{
    Channel *ch;
    int i;

    if ((ch = (Channel*)malloc(sizeof(Channel))) == NULL)
        return false;

    memset(ch, 0, sizeof(Channel));

    ch->highpass_filter        = NSPIKE_DEFAULT_HIGH_FILTER;
    ch->lowpass_filter         = NSPIKE_DEFAULT_LOW_FILTER;
    ch->calibration_factor = 16384;
    ch->sample_rate        = 30000;
    ch->ex                 = Exception_create("");

    return ch;
}

void Channel_set_auxdsp(Channel *ch, AuxDSP *adsp)
{
    ch->auxdsp = adsp;
}

AuxDSP *Channel_get_auxdsp(Channel *ch)
{
    return ch->auxdsp;
}

void Channel_set_number(Channel *ch, int number)
{
    ch->number = number;
}

int Channel_get_number(Channel *ch)
{
    return ch->number;
}

void Channel_set_hw_number(Channel *ch, int number)
{
    ch->hw_number = number;
}

int Channel_get_hw_number(Channel *ch)
{
    return ch->hw_number;
}

void Channel_set_dsp_slot_num(Channel *ch, int number)
{
    ch->dsp_slot_num = number;
}

int Channel_get_dsp_slot_num(Channel *ch)
{
    return ch->dsp_slot_num;
}

void Channel_set_calibration_factor(Channel *ch, short factor)
{
    ch->calibration_factor = factor;
}

short Channel_get_calibration_factor(Channel *ch)
{
    return ch->calibration_factor;
}

void Channel_set_reference_channel(Channel *ch, Channel *ref_ch)
{
    ch->reference_channel = ref_ch;
}

Channel *Channel_get_reference_channel(Channel *ch)
{
    return ch->reference_channel;
}

void Channel_set_highpass_filter(Channel *ch, short value)
{
    ch->highpass_filter = value;
}

short Channel_get_highpass_filter(Channel *ch)
{
    return ch->highpass_filter;
}

void Channel_set_lowpass_filter(Channel *ch, short value)
{
    ch->lowpass_filter = value;
}

short Channel_get_lowpass_filter(Channel *ch)
{
    return ch->lowpass_filter;
}

bool Channel_set_sample_rate(Channel *ch, int sample_rate)
{
    ch->sample_rate = sample_rate;
    return true;
}

int Channel_get_sample_rate(Channel *ch)
{
    return ch->sample_rate;
}

bool Channel_configure_on_home_auxdsp(Channel *ch)
{
    return Channel_configure_on_dsp(ch, ch->auxdsp->dsp, Channel_get_dsp_slot_num(ch));
}

bool Channel_configure_on_dsp(Channel *ch, DSP *dsp, int dsp_slot_num)
{
    unsigned short         data[24];
    int 		   i;

    // Send message to DSP
    /* create the array for programming the channel */
    data[0] = ch->hw_number;
    data[1] = ch->calibration_factor;

    if (ch->reference_channel != NULL) {
        data[2] = Channel_get_hw_number(ch->reference_channel);
        data[3] = -1 * Channel_get_calibration_factor(ch->reference_channel);
    }
    else {
        /* use channel 0 with a gain of 0 to reference to ground */
        data[2] = 0;
        data[3] = 0;
    }

    // printf("channel=%d dsp_ind=%d hw_number=%d ref=%d cal=%d\n", ch->number, ch->dsp_slot_num, ch->hw_number, data[2], data[3]);
    Channel__get_high_filter_coeff(ch, data+4);
    Channel__get_low_filter_coeff(ch, data+14);

    if (!DSP_write_data(dsp, NSPIKE_CHANNEL_SETTINGS_ADDR + NSPIKE_DSP_CHAN_ADDR_INC * 
            dsp_slot_num, 24, data))
    {
        Exception_propagate(dsp->ex, ch->ex);
        return false;
    }
    
    return true;
}


void Channel_destroy(Channel *ch)
{
    Exception_destroy(ch->ex);
    free(ch);
}

void Channel__get_high_filter_coeff(Channel *ch, unsigned short *data)
{
    /* data are as follows:
    * data[0]                 LPFa0lowval        
    * data[1]                 LPFa0highval        
    * data[2]                 LPFa1lowval        
    * data[3]                 LPFa1highval        
    * data[4]                 LPFa0lowval        
    * data[5]                 LPFa0highval        
    * data[6]                 LPFb2lowval        
    * data[7]                 LPFb2highval        
    * data[8]                 LPFb1lowval        
    * data[9]                 LPFb1highval        
    */

    double         d;                //damping factor
    double        alpha;
    double	  alphaout, alphaoutx2;
    double        beta;
    double	  nbetaout;
    double        gamma;
    double	  gammaout;
    double        fsample;        // sampling rate
    double        thetac;                // sampling rate
    double	  tmpmax;
    unsigned short        LPFa0loval;
    unsigned short         LPFa0hival;        
    unsigned short         LPFa1loval;        
    unsigned short         LPFa1hival;        
    unsigned short         LPFb2loval;        
    unsigned short         LPFb2hival;        
    unsigned short         LPFb1loval;        
    unsigned short         LPFb1hival;        
    short lowpass;

    lowpass = ch->lowpass_filter;

    if (lowpass > 10000) {
        memset(data, 0, 10 * sizeof(unsigned short));
        data[5] = 16384;
        return;
    }

    // set the sampling rate 
    fsample = NSPIKE_DSP_BASE_SAMP_RATE;

    // get the cutoff frequency in radians
    thetac = NSPIKE_TWOPI * (double) lowpass / fsample;

    // dampling factor.  2.0 is for butterworth filters
    d = 2.0;


    // alpha, beta, and gamma are the three filter coefficients 
    // derived from theta_c
    beta = 0.5 * (1 - 0.5 * d * sin(thetac))/(1 + 0.5 * d * sin(thetac));

    gamma = (0.5 + beta) * cos(thetac);

    alpha = (0.5 + beta - gamma)/4;
    
    /* convert the coefficients to signed longs */
    nbetaout = trunc(-beta * 32768 * 65536);
    tmpmax = pow(2,32);
    if (nbetaout < 0) {
	nbetaout = tmpmax + nbetaout;
    }
    LPFb2hival = (unsigned short) floor(nbetaout / 65536);
    LPFb2loval = (unsigned short) (((u32) nbetaout) % 65536);

    gammaout = (int) trunc(gamma * 32768 * 65536);
    if (gammaout < 0) {
	gammaout = tmpmax + gammaout;
    }
    LPFb1hival = (unsigned short) floor(gammaout / 65536);
    LPFb1loval = (unsigned short) (((u32) gammaout) % 65536);

    alphaout = (int) trunc(alpha * 32768 * 65536);
    if (alphaout < 0) {
	alphaout = tmpmax - (u32) alphaout;
    }
    LPFa0hival = (unsigned short) floor(alphaout / 65536);
    LPFa0loval = (unsigned short) (((u32) alphaout) % 65536);

    alphaoutx2 = (int) trunc(alpha * UINT_MAX);
    if (alphaoutx2 < 0) {
	alphaoutx2 = UINT_MAX + alphaoutx2;
    }
    LPFa1hival = (unsigned short) floor(alphaoutx2 / 65536);
    LPFa1loval = (unsigned short) (((u32) alphaoutx2) % 65536);



    data[0] = LPFa0loval;
    data[1] = LPFa0hival;
    data[2] = LPFa1loval;
    data[3] = LPFa1hival;
    data[4] = LPFa0loval;
    data[5] = LPFa0hival;
    data[6] = LPFb2loval;
    data[7] = LPFb2hival;
    data[8] = LPFb1loval;
    data[9] = LPFb1hival;
    return;
}


void Channel__get_low_filter_coeff(Channel *ch, unsigned short *data)
{
    /* data are as follows:
    * data[0]                 HPFa0lowval        
    * data[1]                 HPFa0highval        
    * data[2]                 HPFa1lowval        
    * data[3]                 HPFa1highval        
    * data[4]                 HPFa0lowval        
    * data[5]                 HPFa0highval        
    * data[6]                 HPFb2lowval        
    * data[7]                 HPFb2highval        
    * data[8]                 HPFb1lowval        
    * data[9]                 HPFb1highval        
    */

    double         d;                //damping factor
    double        alpha;
    double        alphaout, alphaoutx2;
    double        beta;
    double        nbetaout;
    double        gamma;
    double        gammaout;
    double        fsample;        // sampling rate
    double        thetac;                // sampling rate
    double 	  tmpmax;
    unsigned short        HPFa0loval;
    unsigned short         HPFa0hival;        
    unsigned short         HPFa1loval;        
    unsigned short         HPFa1hival;        
    unsigned short         HPFb2loval;        
    unsigned short         HPFb2hival;        
    unsigned short         HPFb1loval;        
    unsigned short         HPFb1hival;        
    short highpass;

    highpass = ch->highpass_filter;

    if (highpass < 2) {
        memset(data, 0, 10 * sizeof(unsigned short));
        data[5] = 16384;
        return;
    }
    // set the sampling rate 
    fsample = NSPIKE_DSP_BASE_SAMP_RATE;

    // get the cutoff frequency in radians
    thetac = NSPIKE_TWOPI * (double) highpass / fsample;

    // dampling factor.  2.0 is for butterworth filters
    d = 2.0;


    // alpha, beta, and gamma are the three filter coefficients 
    // derived from theta_c
    beta = 0.5 * (1 - 0.5 * d * sin(thetac))/(1 + 0.5 * d * sin(thetac));

    gamma = (0.5 + beta) * cos(thetac);

    alpha = (0.5 + beta + gamma)/4;
    
    tmpmax = pow(2,32);
    /* convert the coefficients to signed longs */
    nbetaout = trunc(-beta * 32768 * 65536);
    if (nbetaout < 0) {
	nbetaout = tmpmax + nbetaout;
    }
    HPFb2hival = (unsigned short) floor(nbetaout / 65536);
    HPFb2loval = (unsigned short) (((u32) nbetaout) % 65536);

    gammaout = (int) trunc(gamma * 32768 * 65536);
    if (gammaout < 0) {
	gammaout = tmpmax + gammaout;
    }
    HPFb1hival = (unsigned short) floor(gammaout / 65536);
    HPFb1loval = (unsigned short) (((u32) gammaout) % 65536);

    alphaout = (int) trunc(alpha * 32768 * 65536);
    if (alphaout < 0) {
	alphaout = tmpmax + alphaout;
    }
    HPFa0hival = (unsigned short) floor(alphaout / 65536);
    HPFa0loval = (unsigned short) (((u32) alphaout) % 65536);

    alphaoutx2 = (int) trunc(-alpha * UINT_MAX);
    if (alphaoutx2 < 0) {
	alphaoutx2 = UINT_MAX + alphaoutx2;
    }
    HPFa1hival = (unsigned short) floor(alphaoutx2 / 65536);
    HPFa1loval = (unsigned short) (((u32) alphaoutx2) % 65536);


    data[0] = HPFa0loval;
    data[1] = HPFa0hival;
    data[2] = HPFa1loval;
    data[3] = HPFa1hival;
    data[4] = HPFa0loval;
    data[5] = HPFa0hival;
    data[6] = HPFb2loval;
    data[7] = HPFb2hival;
    data[8] = HPFb1loval;
    data[9] = HPFb1hival;

    return;
}

