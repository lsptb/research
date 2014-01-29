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
#include <sys/time.h>
#include <unistd.h>
#include "masterdsp.h"
#include "udpsocket.h"


MasterDSP *MasterDSP_create()
{
    MasterDSP *mdsp;
    int i;
    struct timeval now;

    if ((mdsp = (MasterDSP*)malloc(sizeof(MasterDSP))) == NULL)
        return false;

    memset(mdsp, 0, sizeof(MasterDSP));

    mdsp->ex = Exception_create("Master DSP: ");

    for (i = 0; i < NSPIKE_DAC_CHANNEL_COUNT; i++) {
	    mdsp->dac_gain[i]   = NSPIKE_DEFAULT_DAC_GAIN;
        mdsp->dac_cutoff[i] = NSPIKE_DEFAULT_DAC_CUTOFF;
        mdsp->dac_delay[i]  = NSPIKE_DEFAULT_DAC_DELAY;
    }

    mdsp->shutdown = false;
    mdsp->threads_started = false;
    
    pthread_mutex_init(&mdsp->inbound_queue_mutex, NULL);
    pthread_mutex_init(&mdsp->dio_queue_mutex, NULL);
    pthread_mutex_init(&mdsp->dio_state_mutex, NULL);
    pthread_cond_init(&mdsp->dio_messageready, NULL);
    pthread_cond_init(&mdsp->dio_synchronous_state_ready, NULL);
    pthread_cond_init(&mdsp->masterdsp_messageready, NULL);

    gettimeofday(&now, NULL);
    mdsp->last_ping_time = now.tv_sec;

    return mdsp;
}

bool MasterDSP_connect(MasterDSP *mdsp, char *address)
{
    UDPSocket *socket;
    socket = UDPSocket_create();
    if (!UDPSocket_connect(socket, address, NSPIKE_MESSAGE_PORT))
    {
        Exception_propagate(socket->ex, mdsp->ex);
        UDPSocket_destroy(socket);
        return false;
    }

    UDPSocket_set_read_func(socket, (udpsocket_read_func_ptr_type*)&MasterDSP__read_message_from_queue, mdsp);
    UDPSocket_set_blockedread_func(socket, (udpsocket_blockedread_func_ptr_type*)&MasterDSP__blockedread_message_from_queue, mdsp);

    mdsp->master_dsp = DSP_create(0, socket);

    return true;
}

bool MasterDSP_initialize(MasterDSP *mdsp)
{
    int i;
    unsigned short data[24];
   
    mdsp->inbound_message_queue = cfulist_new();
    mdsp->dio_message_queue = cfulist_new();
  
    /* start reader thread */
    pthread_create(&mdsp->reader_thread, NULL, (void*(*)(void*))&MasterDSP__reader_thread_main, (void*)mdsp);

    /* start dio notifier thread */
    pthread_create(&mdsp->dio_notify_thread, NULL, (void*(*)(void*))&MasterDSP__dio_thread_main, (void*)mdsp);
 
    mdsp->threads_started = true;
    
    //
    // configure the DACs
    for (i = 0; i < NSPIKE_DAC_CHANNEL_COUNT; i++) {
	    if(!MasterDSP_set_dac_gain(mdsp, i, mdsp->dac_gain[i])) 
            return false;
        if(!MasterDSP_set_dac_cutoff(mdsp, i, mdsp->dac_cutoff[i]))
            return false;
        if(!MasterDSP_set_dac_delay(mdsp, i, mdsp->dac_delay[i]))
            return false;
    }

    //
    // configure digital io

    /* setup dio statemachine pointers */
	mdsp->dio_statemachinebuffer[0] = NSPIKE_DIO_STATE0_BUFFER_START + 1;
	mdsp->dio_statemachinebuffer[1] = NSPIKE_DIO_STATE1_BUFFER_START + 1;
	mdsp->dio_statemachinebuffer[2] = NSPIKE_DIO_STATE2_BUFFER_START + 1;
	mdsp->dio_statemachinebuffer[3] = NSPIKE_DIO_STATE3_BUFFER_START + 1;
	mdsp->dio_statemachineptr[0] = NSPIKE_DIO_STATE0_PTR;
	mdsp->dio_statemachineptr[1] = NSPIKE_DIO_STATE1_PTR;
	mdsp->dio_statemachineptr[2] = NSPIKE_DIO_STATE2_PTR;
	mdsp->dio_statemachineptr[3] = NSPIKE_DIO_STATE3_PTR;
    
	/* first set up the ports for input or output */
	data[0] = 0;
	for (i = 0; i < mdsp->dio_port_count; i++) {
	    if (mdsp->dio_port_type[i] == NSPIKE_DIO_OUTPUT_PORT) {
		data[0] |= 1 << i;
	    }
	}
	if (!DSP_write_data(mdsp->master_dsp, NSPIKE_DIO_OUT_ENABLE, 1, data)) 
        return false;

	/* Set all of the output port bits to be 0 */ 
	data[0] = 0;
	for (i = 0; i < mdsp->dio_port_count; i++) {
	    if (mdsp->dio_port_type[i] == NSPIKE_DIO_OUTPUT_PORT) {
		    if (!DSP_write_data(mdsp->master_dsp, NSPIKE_DIO_OUT_1 + i, 1, data))
                return false;
	    }
	}
	/* Set the mask for all of the ports to be 1's so that all changes are
	 * reported on the used ports, and set the mask to be all zeron if there are 
	 * unused ports*/
	for (i = 0; i < NSPIKE_MAX_DIO_PORTS; i++) {
	    if (i < mdsp->dio_port_count) {
		data[0] = 0xFFFF;
	    }
	    else {
		data[0] = 0;
	    }
	    if (!DSP_write_data(mdsp->master_dsp, NSPIKE_DIO_IN_1_MASK + i, 1, data))
            return false;
	}

	/* Set the debounce time for the inputs to be 0 ms  */
	data[0] = 20;
	if (!DSP_write_data(mdsp->master_dsp, NSPIKE_DIO_IN_DEBOUNCE_TIME, 1, data))
        return false;

	/* Enable the state machine */
	data[0] = 1;
	if (!DSP_write_data(mdsp->master_dsp, NSPIKE_DIO_STATE_ENABLE, 1, data))
        return false;

    return true;
}

bool MasterDSP_get_address(MasterDSP *mdsp, char *ip, int len)
{
    if (!DSP_get_address(mdsp->master_dsp, ip, len))
    {
        Exception_propagate(mdsp->master_dsp->ex, mdsp->ex);
        return false;
    }

    return true;
}

bool MasterDSP_set_dio_port_count(MasterDSP *mdsp, int count)
{
    if (!(-1 < count < NSPIKE_MAX_DIO_PORTS))
    {
        Exception_set(mdsp->ex, "MasterDSP_set_dio_port_count: invalid number of ports %d", count);
        return false;
    }

    mdsp->dio_port_count = count;
    return true;
}

bool MasterDSP_set_dio_mode(MasterDSP *mdsp, int diogroup, bool output)
{
    if (!(-1 < diogroup < mdsp->dio_port_count))
    {
        Exception_set(mdsp->ex, "MasterDSP_set_dio_mode: invalid dioport %d", diogroup);
        return false;
    }

    if (output == true)
        mdsp->dio_port_type[diogroup] = NSPIKE_DIO_OUTPUT_PORT;
    else
        mdsp->dio_port_type[diogroup] = 0;

    return true;
}

bool MasterDSP_reset_clock(MasterDSP *mdsp)
{
    unsigned short data[1];

    data[0] = 0;
    if (!DSP_write_data(mdsp->master_dsp, NSPIKE_SYNC_CONTROL_ADDR, 1, data))
    {
        Exception_propagate_with_text(mdsp->master_dsp->ex, mdsp->ex, "reset clock failed: ");
        return false;
    }

    data[0] = NSPIKE_BIT0 | NSPIKE_BIT1 | NSPIKE_BIT2;
    if (!DSP_write_data(mdsp->master_dsp, NSPIKE_SYNC_CONTROL_ADDR, 1, data))
    {
        Exception_propagate_with_text(mdsp->master_dsp->ex, mdsp->ex, "reset clock failed: ");
        return false;
    }

    sleep(3);

    return true;
}

void MasterDSP_destroy(MasterDSP *mdsp)
{
    if (mdsp->threads_started)
    {
        mdsp->shutdown = true;
        
        pthread_join(mdsp->reader_thread, NULL);
        pthread_join(mdsp->dio_notify_thread, NULL);
    }
    
    cfulist_destroy(mdsp->inbound_message_queue);
    cfulist_destroy(mdsp->dio_message_queue);

    if (mdsp->master_dsp)
        DSP_destroy(mdsp->master_dsp);

    pthread_mutex_destroy(&mdsp->inbound_queue_mutex);
    pthread_mutex_destroy(&mdsp->dio_queue_mutex);
    pthread_mutex_destroy(&mdsp->dio_state_mutex);
    pthread_cond_destroy(&mdsp->dio_messageready);
    pthread_cond_destroy(&mdsp->dio_synchronous_state_ready);
    pthread_cond_destroy(&mdsp->masterdsp_messageready);

    Exception_destroy(mdsp->ex);

    free(mdsp);
}

bool MasterDSP_set_dac_gain(MasterDSP *mdsp, int channel, short gain)
{
    unsigned short address;
    short converted_gain = 0;

    if (channel == 0) {
	    address = NSPIKE_DAC_0_GAIN_ADDR;
    }
    else if (channel == 1) {
	    address = NSPIKE_DAC_1_GAIN_ADDR;
    }
    else
    {
        Exception_set(mdsp->ex, "dac_gain: channel must be 0 or 1");
        return false;
    }

    /* the gains from 0 to 128 should be mapped to 0 to 32767 */
    converted_gain = (gain < 128) ? (gain * 128) : 32767;
    if(!DSP_write_data(mdsp->master_dsp, address, 1, (unsigned short *) &converted_gain))
        return false;
   
    mdsp->dac_gain[channel] = gain; 

    return true;
}


bool MasterDSP_set_dac_cutoff(MasterDSP *mdsp, int channel, short cutoff)
{
    unsigned short address1;
    unsigned short address2;

    if (channel == 0) {
	address1 = NSPIKE_DAC_0_POS_THRESH_ADDR;
	address2 = NSPIKE_DAC_0_NEG_THRESH_ADDR;
    }
    else {
	address1 = NSPIKE_DAC_1_POS_THRESH_ADDR;
	address2 = NSPIKE_DAC_1_NEG_THRESH_ADDR;
    }
    mdsp->dac_cutoff[channel] = cutoff;
    if (!DSP_write_data(mdsp->master_dsp, address1, 1, (unsigned short *) &cutoff))
        return false;

    cutoff = cutoff * -1;
    if (!DSP_write_data(mdsp->master_dsp, address2, 1, (unsigned short *) &cutoff))
        return false;

    return true;
}

bool MasterDSP_set_dac_delay(MasterDSP *mdsp, int channel, short delay)
{
    unsigned short address;

    if (channel == 0) {
	address = NSPIKE_DAC_0_DELAY_ADDR;
    }
    else {
	address = NSPIKE_DAC_1_DELAY_ADDR;
    }
    mdsp->dac_delay[channel] = delay;
    /* map the delay in ms to the number of samples */
    delay = delay * 30;
    if (!DSP_write_data(mdsp->master_dsp, address, 1, (unsigned short *) &delay))
        return false;

    return true;
}

bool MasterDSP_set_dac_mute(MasterDSP *mdsp, int channel, bool muted)
{
    unsigned short gain = 0;
    unsigned short address;
    if (channel == 0)
        address = NSPIKE_DAC_0_GAIN_ADDR;
    else if (channel == 1)
        address = NSPIKE_DAC_1_GAIN_ADDR;
    else
    {
        Exception_set(mdsp->ex, "dac_mute: channel must be 0 or 1");
        return false;
    }
        
    if (muted) {
	    /* set the gain to 0 */
	    if (!DSP_write_data(mdsp->master_dsp, address, 1, &gain))
            return false;
    }
    else {
        if(!MasterDSP_set_dac_gain(mdsp, channel, mdsp->dac_gain[channel]))
            return false;
    }

    mdsp->dac_mute[channel] = muted;
    return true;
}

bool MasterDSP_set_dac_channel(MasterDSP *mdsp, Channel *ch, int dac_ch_num)
{
    if (!(-1 < dac_ch_num < NSPIKE_DAC_CHANNEL_COUNT))
    {
        Exception_set(mdsp->ex, "set_dac_channel: channel must be 0 or 1");
        return false;
    }

    mdsp->dac_channels[dac_ch_num] = ch;
    return MasterDSP_configure_dac_channel_if_selected(mdsp, ch);
}



/* this function will send configuration to the masterdsp for a given channel
 * if it is selected as a dac channel.  if it is selected and configuration
 * succeeds, it will return true.  if the channel is not selected, it will
 * return true.  it only returns false on error when the channel is selected.
 */
bool MasterDSP_configure_dac_channel_if_selected(MasterDSP *mdsp, Channel *ch)
{
    int i;

    for (i = 0; i < NSPIKE_DAC_CHANNEL_COUNT; i++)
    {
        if (mdsp->dac_channels[i] == ch)
        {
            if (!Channel_configure_on_dsp(ch, mdsp->master_dsp, i))
            {
                return false;
            }
        }
    }

    return true;
}


int MasterDSP_trigger_dio_out(MasterDSP *mdsp, int output, int time_in_us)
{
    unsigned short command[NSPIKE_DIO_MAX_COMMAND_LEN];
    int tmptime, len;

    tmptime = time_in_us;


    /* to trigger an output we send in a seqence of commands to the current 
     * state machine */
    command[0] = NSPIKE_DIO_S_SET_OUTPUT_HIGH | output;
    tmptime = time_in_us;
    len = 1;
    /* string together a set of wait commands to wait for the specified
     * interval */
    while (tmptime > 10000) {
	command[len++] = NSPIKE_DIO_S_WAIT | (10000 * NSPIKE_SAMP_TO_TIMESTAMP);
	tmptime -= 10000;
    }
    command[len++] = NSPIKE_DIO_S_WAIT | (tmptime * NSPIKE_SAMP_TO_TIMESTAMP);
    /* now turn off the output */
    command[len++] = NSPIKE_DIO_S_SET_OUTPUT_LOW | output;
    /* send the command off to the DSP */
    if (!MasterDSP__write_dio_command(mdsp, command, len)) { 
        Exception_set(mdsp->ex, "Error writing DIO command to master DSP");
	return 0;
    }

    // mdsp->dio_raised[output] = 0;

    //sprintf(tmpstring, "Output triggered on port %d for %f ms", output, 
	//    (float) time_in_us / 10.0);
    // DisplayStatusMessage(tmpstring);
    return 1;
}

bool MasterDSP_change_dio_out_state(MasterDSP *mdsp, int output, int raise)
{
    unsigned short command[NSPIKE_DIO_MAX_COMMAND_LEN];

    /* to trigger an output we send in a seqence of commands to the current 
     * state machine */
    if (raise) {
	command[0] = NSPIKE_DIO_S_SET_OUTPUT_HIGH | output;
    }
    else {
	command[0] = NSPIKE_DIO_S_SET_OUTPUT_LOW | output;
    }
    /* send the command off to the DSP */
    if (!MasterDSP__write_dio_command(mdsp, command, 1)) { 
        Exception_set(mdsp->ex, "Error writing digital IO command to master DSP");
	    return false;
    }
/*    if (raise) {
	sprintf(tmpstring, "Output on port %d raised", output);
    }
    else {
	sprintf(tmpstring, "Output on port %d lowered", output);
    }
    */
    // mdsp->dio_raised[output] = raise;
    // DisplayStatusMessage(tmpstring);
    return true;
}


void MasterDSP_register_dio_state_change_callback(MasterDSP *mdsp, dio_status_funcptr *callback, void *user_data_ptr)
{
    mdsp->dio_status_callback = (dio_status_funcptr)callback;
    mdsp->dio_status_callback_user_data = user_data_ptr;
}


int MasterDSP__write_dio_command(MasterDSP *mdsp, unsigned short *command, int len)
{
    int s;
    
    /* check to see that the length is in bounds */
    if (len > NSPIKE_DIO_MAX_COMMAND_LEN) {
	Exception_set(mdsp->ex, "Error: command length %d too long to write to digital IO state machine. Must be <= %d", len, NSPIKE_DIO_STATE_SIZE);
	return 0;
    }

    /* we first need to find a free state machine */
    if ((s = MasterDSP__dio_next_state_machine(mdsp)) == -1) {
	Exception_set(mdsp->ex, "Error: all state machines busy, cannot program digital IO");
	return 0;
    }
    /* we have to append a jump statement to the command to tell the state
     * machine to jump back to the first instruction, which is "wait forever,"
     * when it completes the command */
    command[len++] = NSPIKE_DIO_S_JUMP_ABS | 0;

    /* Now write the command to the DSP */
    if (!DSP_write_data(mdsp->master_dsp, mdsp->dio_statemachinebuffer[s], len, command)) {
        Exception_set(mdsp->ex, "Error writing digital IO command to master DSP");
	    return 0;
    }
    /* Finally, we reset the state machine pointer to the beginning of the 
     * new command, which is at offset 1 from the beginning of the buffer */
    command[0] = 1;
    if (!DSP_write_data(mdsp->master_dsp, mdsp->dio_statemachineptr[s], 1, command)) {
        Exception_set(mdsp->ex, "Error writing digital IO state machine pointer to master DSP");
	    return 0;
    }
    return 1;
}


int MasterDSP__dio_next_state_machine(MasterDSP *mdsp)
    /* returns the index of the next available state machine or -1 if 
     * all are busy */
{
    unsigned short status[1];
    if (DSP_read_data(mdsp->master_dsp, NSPIKE_DIO_STATE_AVAILABLE, 1, status) != 1) {
        Exception_set(mdsp->ex, "Error reading Master DSP digital IO state pointer");
	    return NSPIKE_DIO_STATE_MACHINES_BUSY;
    }
    if ((short) status[0] != NSPIKE_DIO_STATE_MACHINES_BUSY) {
	return (int) status[0];
    }
    else {
	return -1;
    }
}

int MasterDSP__read_message_from_queue(UDPSocket *udps, void *buffer, int size)
{
    size_t readsize;       
    unsigned short *data;
    MasterDSP *mdsp;
    mdsp = (MasterDSP*)UDPSocket_get_read_userptr(udps);

    pthread_mutex_lock(&mdsp->inbound_queue_mutex);

    if (cfulist_num_entries(mdsp->inbound_message_queue) > 0)
    {
        cfulist_pop_data(mdsp->inbound_message_queue, (void**)&data, &readsize);
        if (readsize > size)
           readsize = size; 
        memcpy(buffer, data, readsize);
    }
    else
        readsize = -1;

    pthread_mutex_unlock(&mdsp->inbound_queue_mutex);

    ///printf("popping %d byte message off queue\n", readsize);
    return readsize;
}

int MasterDSP__blockedread_message_from_queue(UDPSocket *udps, void *buffer, int size, struct timeval *tv)
{
    size_t readsize;       
    unsigned short *data;
    MasterDSP *mdsp;
    mdsp = (MasterDSP*)UDPSocket_get_read_userptr(udps);
    struct timespec ts;
    struct timeval tv2;
    
    gettimeofday(&tv2, NULL);
    ts.tv_sec = tv2.tv_sec;
    ts.tv_nsec = tv2.tv_usec*1000;
    ts.tv_sec += tv->tv_sec;
    ts.tv_nsec += tv->tv_usec*1000;

    if (ts.tv_nsec >= 1000000000) {
        ts.tv_sec += ts.tv_nsec / 1000000000;
        ts.tv_nsec %= 1000000000;
    }

    pthread_mutex_lock(&mdsp->inbound_queue_mutex);

    if (cfulist_num_entries(mdsp->inbound_message_queue) < 1)
        pthread_cond_timedwait(&mdsp->masterdsp_messageready, &mdsp->inbound_queue_mutex, &ts);

    if (cfulist_num_entries(mdsp->inbound_message_queue) > 0)
    {
        cfulist_pop_data(mdsp->inbound_message_queue, (void**)&data, &readsize);
        if (readsize > size)
           readsize = size; 
        memcpy(buffer, data, readsize);
    }
    else
    {
        readsize = -1;
    }
    
    //printf("popping %d byte message off queue\n", readsize);
    pthread_mutex_unlock(&mdsp->inbound_queue_mutex);
    return readsize;
}

/* all functions below here live exclusively in the reader thread */

void MasterDSP__reader_thread_main( void *void_mdsp_ptr )
{
     MasterDSP *mdsp;
     mdsp = (MasterDSP*) void_mdsp_ptr;

     fd_set rfds;
     struct timeval timeout;
     unsigned short readbuf[NSPIKE_MAX_PACKET_SIZE];
     unsigned short *message;
     int readsize;
     char tmp[256];
     struct timeval tv;

     for(;;)
     {
//         printf(">>>>>>>READER LOOP\n\n");
         if (mdsp->shutdown)
         {
             pthread_exit(NULL);
         }

         timeout.tv_usec = 0; // 50ms
         timeout.tv_sec  = 1;
         if (-1 == (readsize = UDPSocket__native_blockedread(mdsp->master_dsp->control_socket, readbuf, NSPIKE_MAX_PACKET_SIZE * sizeof(short), &timeout)))
         {
             Exception_get(mdsp->master_dsp->control_socket->ex, tmp, sizeof(tmp));
             Exception_printf("udp read error: %s\n", tmp);
             continue;
         }
         
         if (readsize == 0)
             continue;
         
         message = (short unsigned int*)malloc(readsize);
         memcpy(message, readbuf, readsize);
         /* the master dsp is supposed to set a specific pipeid for dio messages
          * but when i ran the sniffer i didn't see that.  loren checks the packet size
          * so that's what we're going to do here.... */
         if (readsize == NSPIKE_DIO_MESSAGE_SIZE)
         {
            pthread_mutex_lock(&mdsp->dio_queue_mutex);
            cfulist_push_data(mdsp->dio_message_queue, (void *)message, readsize);
            // printf("got dio message packet: queue len %d\n", cfulist_num_entries(mdsp->dio_message_queue));
            pthread_cond_signal(&mdsp->dio_messageready);
            pthread_mutex_unlock(&mdsp->dio_queue_mutex);
         }
         else
         {
            pthread_mutex_lock(&mdsp->inbound_queue_mutex);
            cfulist_push_data(mdsp->inbound_message_queue, (void *)message, readsize);
            pthread_cond_signal(&mdsp->masterdsp_messageready);
            // printf("got message packet: queue len %d\n", cfulist_num_entries(mdsp->inbound_message_queue));
            pthread_mutex_unlock(&mdsp->inbound_queue_mutex);
         }
  
     }
}

int MasterDSP_wait_for_dio_change(MasterDSP *mdsp, int channel_num, int state)
{
    int timestamp = -1;
    int i; 
    for (i = 0; i < 20; i++)
    {
        pthread_mutex_lock(&mdsp->dio_state_mutex);
        if (mdsp->last_dio_timestamp != 0)
        {
            timestamp = mdsp->last_dio_timestamp;
            mdsp->last_dio_timestamp = 0;
            pthread_mutex_unlock(&mdsp->dio_state_mutex);
            break;
        }
        pthread_mutex_unlock(&mdsp->dio_state_mutex);
        usleep(250000);
    }

    return timestamp; 
}

void MasterDSP__dio_thread_main( void *void_mdsp_ptr )
{
    MasterDSP *mdsp;
    mdsp = (MasterDSP*) void_mdsp_ptr;
    size_t readsize;       
    unsigned short *data;
    int i;
    u32 timestamp;
    unsigned short dio_state[NSPIKE_MAX_DIO_PORTS];
    struct timeval tv;
    struct timespec ts;

    tv.tv_sec = 1;
    tv.tv_usec = 0;


    for (;;)
    {
        pthread_mutex_lock(&mdsp->dio_queue_mutex);
        MasterDSP__ping(mdsp);
        if (mdsp->shutdown)
        {
           pthread_exit(NULL);
        }
        gettimeofday(&tv, NULL);  
        tv.tv_sec += 1;

        ts.tv_sec = tv.tv_sec;
        ts.tv_nsec = tv.tv_usec*1000;

        if (pthread_cond_timedwait(&mdsp->dio_messageready, &mdsp->dio_queue_mutex, &ts) == 0)
        {
            // printf("dio thread loop!\n");

            if (cfulist_num_entries(mdsp->dio_message_queue) > 0)
            {
                cfulist_pop_data(mdsp->dio_message_queue, (void**)&data, &readsize);
                pthread_mutex_unlock(&mdsp->dio_queue_mutex);

                if (mdsp->dio_status_callback)
                {
                    DSP__byte_swap(data, NSPIKE_DIO_MESSAGE_SIZE / sizeof(unsigned short));
                    timestamp = (data[1] + (data[2] << 16));
                    for (i = 0; i < NSPIKE_MAX_DIO_PORTS; i++) {
                        dio_state[i] = data[i+3];
                    }
                    (mdsp->dio_status_callback)(timestamp, dio_state, mdsp->dio_status_callback_user_data);
                }
               
                /* update the last dio timestamp */ 
                pthread_mutex_lock(&mdsp->dio_state_mutex);
                mdsp->last_dio_timestamp = timestamp;
                memcpy(mdsp->last_dio_state, dio_state, sizeof(mdsp->last_dio_state));
                pthread_mutex_unlock(&mdsp->dio_state_mutex);
            }
            else
                pthread_mutex_unlock(&mdsp->dio_queue_mutex);
            
        }
        else
        {
            //printf("dio wakeup!\n");
            pthread_mutex_unlock(&mdsp->dio_queue_mutex);
        }
    }
}

/* XXX HACK HACK SLASH SLASH:  "Ping" the MasterDSP by setting it's gain so that the switch
* doesn't forget where we are */
void MasterDSP__ping(MasterDSP *mdsp)
{
    int gain_for_ping = 0;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    if (tv.tv_sec - mdsp->last_ping_time > NSPIKE_DSP_ECHO_INTERVAL)
    {
        if (mdsp->dac_mute[0])
            gain_for_ping = 0;
        else
            gain_for_ping = mdsp->dac_gain[0];
        
        if(!MasterDSP_set_dac_gain(mdsp, 0, gain_for_ping))
        {
            Exception_printf("WARNING: unable to set dac gain (ping) for network/switch keepalive\n");
        }

        // printf("*********** PING!! \n");
        mdsp->last_ping_time = tv.tv_sec;
    }  

}
