// #define _GNU_SOURCE
#include <stdio.h>
#include <stdarg.h>
#include <comedilib.h>
#include "libnspike.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <sys/select.h>
#include <arpa/inet.h>
#include <sys/un.h>
#include <sys/ipc.h> 
#include <sys/sem.h>
#include <sys/shm.h>
#include <fcntl.h>
#include <math.h>
#include <linux/fadvise.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include "nstream.h"
#include "fglove.h"
#include "dataglovedaq.h"
#include "nstreamclient.h"


#define WRITER_NOT_INITIALIZED 0
#define WRITER_STOPPED         1
#define WRITER_START           2
#define WRITER_SYNC            3
#define WRITER_STOP            4
#define WRITER_ACTIVE          5

int comedi_start_acquire();
void writer_thread_main(void *foo);

NSpike *nsp;
unsigned int channel_list[NSTREAM_NUM_COMEDI_CHANNELS];
comedi_t *device_handle;
comedi_cmd cmd;
DataGloveDAQ *dataglove_daq;

pthread_t writer_thread;

/* writer thread state variables */
int writer_state = WRITER_NOT_INITIALIZED;
char nspike_low_filename[8192];
char nspike_high_filename[8192];

char comedi_filename[8192];
char dio_filename[8192];
char log_filename[8192];
char dataglove_filename[8192];
unsigned short *writer_comedi_ptr;
short *writer_nspike_ptr;
unsigned long recording_start_time = 0;
int nspike_writer_decimation_factor = 1;
int comedi_writer_decimation_factor = 3;
int comedi_num_channels_to_write = NSTREAM_NUM_COMEDI_CHANNELS;
int nspike_num_channels_to_write_low = NSTREAM_NUM_NSPIKE_CHANNELS;
int nspike_num_channels_to_write_high = NSTREAM_NUM_NSPIKE_CHANNELS;
int nspike_high_channels_offset = 0;

/* acquisition active flag */
bool acquire_active = false;
bool nstream_shutdown = false;
bool verbose_state = true;
char start_error[512];

/* acquisition start timestamp */
long int start_timestamp;

/* input file descriptors */
int comedi_fd = 0;
int nspike_fd = 0;
int control_fd = 0;
int messages_fd = 0;
int maxreadfd_all = 0;

/* ring buffer for comedi data */
RingBuffer *comedi_rb;

/* output file pointers */
FILE *nspike_low_out;
FILE *nspike_high_out;
FILE *comedi_out;
FILE *dio_out;
FILE *log_out;
FILE *dataglove_out;

void process_control_message(int sock);
bool start_acquire();
bool stop_acquire();
bool start_save();
bool stop_save();
bool set_dac_gain(int dac_chan, int dac_gain);
bool set_dac_channel(int dac_chan, int analog_chan);
bool set_verbose(bool enable);
bool set_channel_filters(int data_channel, int high_pass, int low_pass);
bool comedi_init();
bool nspike_init();
int comedi_start_acquire();
int get_control_socket();
int interpolate_3x(void *buffer, RingBuffer *rb, int buf_len);
shared_state_vars *shared_state_vars_init(void);
long int write_shmem(char *circ_buf, char **circ_buf_write_ptr, char *circ_buf_end_ptr, char *in_buf, int in_len);
void handle_dio_data(int timestamp, unsigned short data[4], void *x);
bool check_shmem_ptrs_time();
int set_real_time_priority(void);


struct shared_state_vars *shared_state;

void lprintf(const char *fmt, ...)
{
    va_list args;
    char string[8192];
 
    va_start(args, fmt);
    vsnprintf(string, 8192, fmt, args);
    va_end(args);
    
    if (writer_state == WRITER_ACTIVE)
    {
        fprintf(log_out, "nstream: %s", string);
        fflush(log_out);
    }

    printf("nstream: %s", string);
}

void print_timestamp()
{
  char buffer[30];
  struct timeval tv;

  time_t curtime;

  gettimeofday(&tv, NULL); 
  curtime=tv.tv_sec;

  strftime(buffer,30,"%m-%d-%Y  %T.",localtime(&curtime));
  lprintf("%s%ld\n",buffer,tv.tv_usec);
}

void exit_clean(int code)
{
    if (writer_state != WRITER_NOT_INITIALIZED)
    {
        if (writer_state != WRITER_STOPPED)
            writer_state = WRITER_STOP;
        sleep(1);
        nstream_shutdown = true;
        pthread_join(writer_thread, NULL);
    }
    exit(code);
}


main(int argc, char *argv[])
{
    int i;
    char str[100];
    fd_set readfds;
    int maxreadfd     = 0;

    dataglove_daq = NULL;
   
    set_real_time_priority();

    /* initialization */
    lprintf("initializing shared memory...\n");
    if ((shared_state = shared_state_vars_init()) == NULL)
    {
        lprintf("failed, exiting.\n");
        exit_clean(-1);
    }

    lprintf("initializing libnspike...\n");
    if (!nspike_init())
    {
        lprintf("failed, exiting.\n");
        exit_clean(-1);
    }
    lprintf("libnspike initialized.\n");
    
    lprintf("initializing comedi...\n");
    if (!comedi_init())
    {
        lprintf("failed, exiting.\n");
        exit_clean(-1);
    }
    lprintf("comedi initialized.\n");

     
    /* get fd for unix domain control socket */
    lprintf("creating control socket...\n");
    control_fd = get_control_socket();

    maxreadfd_all = ((comedi_fd > nspike_fd) ? comedi_fd : nspike_fd);
    maxreadfd_all = comedi_fd;
    maxreadfd_all = ((maxreadfd_all > messages_fd) ? maxreadfd_all : messages_fd);
    maxreadfd_all = ((maxreadfd_all > control_fd) ? maxreadfd_all : control_fd) + 1;
#if 0 
    /* startup acquisition */ 
    if (!start_acquire())
    {
        lprintf("unable to start acquisition, exiting\n");
        exit_clean(-1);
    }
#endif
    /* start writer thread */
    writer_state = WRITER_STOPPED;
    pthread_create(&writer_thread, NULL, (void* (*)(void*))&writer_thread_main, NULL);

    /* ----------- main loop ---------------- */
    
    char buffer[65536];
    int bytes_read = 0;
    comedi_rb = RingBuffer_create(65536*2);
    struct timeval tv, tv2;

    // int l = 0;

    while(1)
    {
        FD_ZERO(&readfds);
        FD_SET(control_fd, &readfds);
        FD_SET(messages_fd, &readfds);
        if (acquire_active == true)
        {
            FD_SET(comedi_fd, &readfds);
            FD_SET(nspike_fd, &readfds);
            maxreadfd = maxreadfd_all;
        }
        else 
            maxreadfd = ((control_fd > messages_fd) ? control_fd : messages_fd) + 1;
        
        select(maxreadfd, &readfds, NULL, NULL, NULL);

        shared_state->latest_timestamp = NSpike_get_latest_timestamp(nsp) - start_timestamp;
 
        if (FD_ISSET(comedi_fd,&readfds))
        {
            bytes_read = read(comedi_fd, buffer, sizeof(buffer) - 1);
            /* we have to interpolate on even channel boundaries 
             * so put comedi data in a ring buffer and let the
             * interpolation function pull out the size it needs */
            if (bytes_read > 0)
            {
                if (!RingBuffer_push_data(comedi_rb, (char*)&buffer, bytes_read))
                {
                    lprintf("comedi ringbuffer overflow, bailing!\n");
                    stop_acquire();
                    exit_clean(-1);
                }
                memset(buffer, 0, sizeof(buffer));
                bytes_read = interpolate_3x(&buffer, comedi_rb, sizeof(buffer));
                
                write_shmem((char*)shared_state->comedi_buffer, 
                            (char**)&shared_state->comedi_buffer_write_ptr,
                            (char*)shared_state->comedi_buffer_end_ptr, 
                            buffer, bytes_read);
                
                shared_state->comedi_time += bytes_read / NSTREAM_NUM_COMEDI_CHANNELS / NSTREAM_SS;
            }
            if (bytes_read == -1)
            { 
                lprintf("comedi error, bailing!\n");
                stop_acquire();
                exit_clean(-1);
            }
                
        }
        if (FD_ISSET(nspike_fd, &readfds))
        {
            bytes_read = read(nspike_fd, buffer, sizeof(buffer));
            write_shmem((char*)shared_state->nspike_buffer, 
                        (char**)&shared_state->nspike_buffer_write_ptr, 
                        (char*)shared_state->nspike_buffer_end_ptr, 
                        buffer, bytes_read);
            shared_state->nspike_time += (bytes_read / NSTREAM_NUM_NSPIKE_CHANNELS / NSTREAM_SS);
        }
        if (FD_ISSET(control_fd, &readfds))
        {
            process_control_message(control_fd);
        }
        if (FD_ISSET(messages_fd, &readfds))
        {
            bytes_read = read(messages_fd, buffer, sizeof(buffer));
            buffer[bytes_read] = NULL;
            lprintf("%s",buffer);
        }
        
        /*
         * handle data glove data
         */
        if (dataglove_daq != NULL)
        {
            if (dataglove_daq->ReadData(shared_state->latest_timestamp, 
                                        &shared_state->dataglove_buffer[shared_state->dataglove_buffer_index]))
            {
                if (writer_state == WRITER_ACTIVE)
                {
                    fprintf(dataglove_out, "%u,%0.03f,%0.03f,%0.03f,%0.03f,%0.03f,%d\n",
                        shared_state->dataglove_buffer[shared_state->dataglove_buffer_index].timestamp - recording_start_time,
                        shared_state->dataglove_buffer[shared_state->dataglove_buffer_index].thumb,
                        shared_state->dataglove_buffer[shared_state->dataglove_buffer_index].index,
                        shared_state->dataglove_buffer[shared_state->dataglove_buffer_index].middle,
                        shared_state->dataglove_buffer[shared_state->dataglove_buffer_index].ring,
                        shared_state->dataglove_buffer[shared_state->dataglove_buffer_index].little,
                        shared_state->dataglove_buffer[shared_state->dataglove_buffer_index].gesture);
                }

                if (shared_state->dataglove_buffer_index + 1 == NSTREAM_DATAGLOVE_NUM_SAMPLES)
                    shared_state->dataglove_buffer_index = 0;
                else
                    shared_state->dataglove_buffer_index++;
            }
        }
        //lprintf("nspike - comedi write time: %d\n", shared_state->nspike_time - shared_state->comedi_time);
        //lprintf("nspike time: %d\tcomedi_time: %d\tdifference: %d\n", shared_state->nspike_time, shared_state->comedi_time, shared_state->nspike_time - shared_state->comedi_time);
//        check_shmem_ptrs_time();
//        writer_thread_main();
    }
}

long int write_shmem(char *circ_buf, char **circ_buf_write_ptr, char *circ_buf_end_ptr, char *in_buf, int in_len)
{
    int len_to_end = 0;
    int bytes_written = 0;

    if (*circ_buf_write_ptr + in_len > circ_buf_end_ptr)
    {
        len_to_end = circ_buf_end_ptr - *circ_buf_write_ptr;
        memcpy(*circ_buf_write_ptr, in_buf, len_to_end);
        bytes_written = len_to_end;
        *circ_buf_write_ptr = circ_buf;
    }

    memcpy(*circ_buf_write_ptr, in_buf + len_to_end, in_len - len_to_end);
    *circ_buf_write_ptr += (in_len - len_to_end);
    
    bytes_written += in_len - len_to_end;

    return bytes_written;
}
    
int interpolate_3x(void *buffer, RingBuffer *rb, int buf_len)
{
    unsigned short tmpbuf[65536];
    unsigned short *write_ptr = &tmpbuf[0];
    unsigned short *read_ptr = (short unsigned int *)buffer;
    int i;
    
    int shorts_in_buffer = RingBuffer_data_available(rb) / sizeof(unsigned short);
    if (shorts_in_buffer * 2 * 3 > 65536)
        shorts_in_buffer = 10922;

    int shorts_on_chan_boundary = shorts_in_buffer - (shorts_in_buffer % NSTREAM_NUM_COMEDI_CHANNELS);
    if ((i = RingBuffer_get_data(rb, (char*)buffer, shorts_on_chan_boundary * sizeof(unsigned short))) !=
            shorts_on_chan_boundary * sizeof(unsigned short))
        lprintf("ERROR: comedi data misaligned in ringbuffer, behavior data screwed.\n");

    for (i = 0; i < shorts_on_chan_boundary / NSTREAM_NUM_COMEDI_CHANNELS; i++)
    {
        memcpy(write_ptr, read_ptr, NSTREAM_NUM_COMEDI_CHANNELS * sizeof(unsigned short));
        memcpy(write_ptr + NSTREAM_NUM_COMEDI_CHANNELS, read_ptr, NSTREAM_NUM_COMEDI_CHANNELS * sizeof(unsigned short));
        memcpy(write_ptr + (2 * NSTREAM_NUM_COMEDI_CHANNELS), read_ptr, NSTREAM_NUM_COMEDI_CHANNELS * sizeof(unsigned short));
        read_ptr += NSTREAM_NUM_COMEDI_CHANNELS;
        write_ptr += NSTREAM_NUM_COMEDI_CHANNELS * 3;
    }

    memcpy(buffer, &tmpbuf, shorts_on_chan_boundary * 3 * sizeof(unsigned short));

    return shorts_on_chan_boundary * 3 * sizeof(unsigned short);
} 

void send_control_reply(bool success, char *msg, struct sockaddr_in *client)
{
    control_response cr;

    cr.success = success;
    if (!success)
        strncpy(cr.error_message, msg, sizeof(cr.error_message));
  
    if (sendto(control_fd, &cr, sizeof(control_response), 0, (sockaddr*)client, sizeof(struct sockaddr_in)) != sizeof(control_response))
        lprintf("error sending control reply");
}
 
/* ---------- parse and process control messages --------------- */
/* this could really stand to be updated to pass back meaningful errors
 * from libnspike or actual nstream state, but for now this will work
 * well enough */
void process_control_message(int sock)
{
    control_message message;
    struct sockaddr_in client_addr;
    int bytes_read = 0;
    char error_msg[256];
    Channel *ch;
  
    memset(&message, 0, sizeof(message));
    memset(&client_addr, 0, sizeof(client_addr));
    socklen_t len = sizeof(client_addr);

    bytes_read = recvfrom(sock, &message, sizeof(message), 0, (sockaddr*)&client_addr, &len);
    if (bytes_read < 1)
    {
        lprintf("error reading control message\n");
        return;
    }

    if (message.message_type == 'S')
    {
       if (writer_state == WRITER_ACTIVE)
       {
           send_control_reply(false, "save already active", &client_addr);
           lprintf("save already active\n");
           return;
       }

       if (!acquire_active)
       {
           send_control_reply(false, "can't start save when acquisition not active", &client_addr);
           lprintf("can't start save when acquisition not active\n");
           return;
       }
       
       lprintf("starting save...\n");
       
       sprintf(nspike_low_filename, "%s/%s.low.nspike.dat", message.recording_path, 
               message.recording_filename_root);
       sprintf(nspike_high_filename, "%s/%s.high.nspike.dat", message.recording_path, 
               message.recording_filename_root);
       sprintf(comedi_filename, "%s/%s.comedi.dat", message.recording_path, 
               message.recording_filename_root);
       sprintf(dio_filename, "%s/%s.dio.txt", message.recording_path, 
               message.recording_filename_root);
       sprintf(log_filename, "%s/%s.log.txt", message.recording_path, 
               message.recording_filename_root);
       sprintf(dataglove_filename, "%s/%s.glove.txt", message.recording_path, 
               message.recording_filename_root);

       lprintf("nspike low decimation_factor:     %d\n", message.decimation_factor);
       lprintf("effective low sample rate: %dHz\n", NSTREAM_SR/message.decimation_factor);
       lprintf("number of comedi channels to write:     %d\n", message.comedi_num_channels_to_write);
       lprintf("number of nspike channels to write low:     %d\n", message.nspike_num_channels_to_write_low);
       lprintf("number of nspike channels to write high:     %d\n", message.nspike_num_channels_to_write_high);
       lprintf("nspike high channels offset:     %d\n", message.nspike_high_channels_offset);
       lprintf("comedi_file:           %s\n", comedi_filename);
       lprintf("nspike_low_file:       %s\n", nspike_low_filename);
       lprintf("nspike_high_file:      %s\n", nspike_high_filename);
       lprintf("dio_file:              %s\n", dio_filename);
       lprintf("log_file:              %s\n", log_filename);
       if (dataglove_daq != NULL)
           lprintf("glove_file:            %s\n", dataglove_filename);

       nspike_writer_decimation_factor = message.decimation_factor;
       comedi_num_channels_to_write = message.comedi_num_channels_to_write;
       nspike_num_channels_to_write_low = message.nspike_num_channels_to_write_low;
       nspike_num_channels_to_write_high = message.nspike_num_channels_to_write_high;
       nspike_high_channels_offset = message.nspike_high_channels_offset;
       start_save();
       send_control_reply(true, NULL, &client_addr);
    }
    else if (message.message_type == 'T')
    {
       if (writer_state != WRITER_ACTIVE)
       {
           send_control_reply(false, "can't stop inactive save", &client_addr);
           lprintf("can't stop inactive save\n");
           return;
       }
       lprintf("stopping save\n");
       stop_save();
       send_control_reply(true, NULL, &client_addr);
    }
    else if (message.message_type == 'G')
    {
        if (!set_dac_gain(message.dac_channel, message.dac_gain))
           send_control_reply(false, "failed.  check log.", &client_addr);
        else
           send_control_reply(true, NULL, &client_addr);
    }
    else if (message.message_type == 'L')
    {
        if(!set_dac_channel(message.dac_channel, message.data_channel))
           send_control_reply(false, "failed.  check log.", &client_addr);
        else
           send_control_reply(true, NULL, &client_addr);
    }
    else if (message.message_type == 'F')
    {
        if(!set_channel_filters(message.data_channel, message.high_pass, message.low_pass))
           send_control_reply(false, "failed.  check log.", &client_addr);
        else
           send_control_reply(true, NULL, &client_addr);
    }
    else if (message.message_type == 'V')
    {
        if(!set_verbose(message.enable))
           send_control_reply(false, "failed.  check log.", &client_addr);
        else
           send_control_reply(true, NULL, &client_addr);
    }    
    else if (message.message_type == 'X')
    {
       lprintf("got exit message\n");
       send_control_reply(true, NULL, &client_addr);
       stop_acquire();
       lprintf("clean exit\n");
       exit_clean(1);
    }
    else if (message.message_type == 'C')
    {
        if (acquire_active)
        {
           send_control_reply(false, "can't configure channels with acquisition active.", &client_addr);
           return;
        }

        if (message.hw_channel < 0 || message.hw_channel > 255)
        {
           send_control_reply(false, "hardware channel must be 0 < number < 255", &client_addr);
           return;
        }


        lprintf("mapping software channel %d to hardware channel %d\n",message.sw_channel,message.hw_channel);
        /* we have two systems, each which see channels 1-127, so decrement here if we're over.
         * also: technically nspike is a 127 channel system, so replicate channel 0 on 128 */
        if (message.hw_channel > 127)
            message.hw_channel -= 128;
        if (message.hw_channel == 127)
            message.hw_channel = 0;

        ch = Channel_create(); 
        Channel_set_number(ch, message.sw_channel);
        Channel_set_hw_number(ch, message.hw_channel);
        Channel_set_sample_rate(ch, NSTREAM_SR);
        if (!NSpike_set_channel(nsp, ch))
        {
           Channel_destroy(ch);
           NSpike_get_error(nsp, (char*)&error_msg, 256);
           send_control_reply(false, (char*)&error_msg, &client_addr);
           return;
        }
        send_control_reply(true, NULL, &client_addr);
    }
    else if (message.message_type == 'Q')
    {
        if (acquire_active)
        {
           send_control_reply(false, "acquisition already active.", &client_addr);
           return;
        }

        if (start_acquire())
            send_control_reply(true, NULL, &client_addr);
        else
           send_control_reply(false, (char*)&start_error, &client_addr);

    }
    else if (message.message_type == 'A')
    {
        if (acquire_active)
        {
           send_control_reply(false, "can't change config after acquisition started.", &client_addr);
           return;
        }

        if (!NSpike_add_auxdsp(nsp, (char*)&message.ip_address,30000))
        {
           NSpike_get_error(nsp, (char*)&error_msg, 256);
           send_control_reply(false, (char*)&error_msg, &client_addr);
           return;
        }
        
        send_control_reply(true, NULL, &client_addr);
    }
    else if (message.message_type == 'M')
    {
        if (acquire_active)
        {
           send_control_reply(false, "can't change config after acquisition started.", &client_addr);
           return;
        }

        if (!NSpike_add_masterdsp(nsp, (char*)&message.ip_address))
        {
           NSpike_get_error(nsp, (char*)&error_msg, 256);
           send_control_reply(false, (char*)&error_msg, &client_addr);
           return;
        }
        
        send_control_reply(true, NULL, &client_addr);
    }
    else if (message.message_type == 'D')
    {
        if (acquire_active)
        {
           send_control_reply(false, "can't change config after acquisition started.", &client_addr);
           return;
        }
        
        if (dataglove_daq != NULL)
            delete dataglove_daq;

        dataglove_daq = new DataGloveDAQ();
        if (!dataglove_daq->Initialize((char*)&message.serial_port))
        {
            delete dataglove_daq;
            send_control_reply(false, "failed to initialize data glove.", &client_addr);
            return;
        }
        dataglove_daq->SetDelay(message.delay);

        send_control_reply(true, NULL, &client_addr);
    } 
    else
    {
        lprintf("warning: unknown control message type\n");
        send_control_reply(false, "unparsable control message", &client_addr);
    }
}

bool start_acquire()
{
    unsigned long trigger_timestamp   = 0;
    unsigned long comedi_pad_samples  = 0;
    int fd = 0;
    short *zeros = NULL;
    
    start_timestamp     = 0;

    NSpike_register_dio_state_change_callback(nsp, (dio_status_funcptr*)&handle_dio_data, NULL);

    lprintf("configuring nspike hardware...\n");
    if (!NSpike_configure(nsp))
    {
        NSpike_get_error(nsp, (char*)&start_error, sizeof(start_error));
        lprintf("nspike_configure() failed.  check hardware.\n");
        return false;
    }

    
    lprintf("setting up comedi acquisition\n");
    if ((comedi_fd = comedi_start_acquire()) < 0)
    {
        sprintf((char*)&start_error, "unable to start comedi acquire\n");
        lprintf("unable to start comedi acquire\n");
        return false;
    }

    nspike_fd = NSpike_get_data_socket(nsp);
    fcntl(nspike_fd,F_SETFL,O_NONBLOCK);
    maxreadfd_all = ((nspike_fd > maxreadfd_all) ? (nspike_fd + 1) : maxreadfd_all);
    
    lprintf("starting nspike acquisition\n"); 
    if ((start_timestamp = NSpike_start_acquire(nsp)) == -1)
    {
       NSpike_get_error(nsp, (char*)&start_error, sizeof(start_error));
       lprintf("nspike_start_acquire() failed\n"); 
       return false;
    }
   
    usleep(50000);
    
    lprintf("raising nspike trigger digital line\n"); 
    if((trigger_timestamp = NSpike_raise_dio(nsp, 0)) == -1)
    { 
        NSpike_get_error(nsp, (char*)&start_error, sizeof(start_error));
        lprintf("failed to get dio raise timestamp\n");
        return false;
    }
    
    comedi_pad_samples = (trigger_timestamp - start_timestamp) * NSTREAM_NUM_COMEDI_CHANNELS; 
    
    lprintf("trigger raised, trigger_timestamp = %d, start_timestamp = %d diff = %d\n", 
            trigger_timestamp, start_timestamp, trigger_timestamp - start_timestamp);

    /* add offset to shared memory */
    memset(shared_state->comedi_buffer, 0, comedi_pad_samples);
    // shared_state->comedi_buffer_write_ptr = &shared_state->comedi_buffer[0];
    
    shared_state->comedi_time = trigger_timestamp - start_timestamp;
    shared_state->comedi_buffer_write_ptr = &shared_state->comedi_buffer[0] + comedi_pad_samples;
    
    shared_state->nspike_buffer_write_ptr = &shared_state->nspike_buffer[0];
    
    shared_state->comedi_buffer_end_ptr = &shared_state->comedi_buffer[0] + NSTREAM_COMEDI_NUM_SAMPLES;
    shared_state->nspike_buffer_end_ptr = &shared_state->nspike_buffer[0] + NSTREAM_NSPIKE_NUM_SAMPLES;
    
    lprintf("acquisition started.\n");
    acquire_active = true;
    return true;
}

bool stop_acquire()
{
    char buffer[1024];
    lprintf("stopping comedi\n");
    if(comedi_cancel(device_handle, comedi_find_subdevice_by_type(device_handle,COMEDI_SUBD_AI,0)) != 0) 
    {
        lprintf("failed to stop comedi acquisition\n"); 
        return false;
    }

    lprintf("stopping nspike\n");
    if (!NSpike_stop_acquire(nsp))
    {
        lprintf("failed to stop nspike acquisition\n"); 
        return false; 
    }

   
    if (acquire_active == true) 
    {
        usleep(100000);
        /* drain fds */
        while(read(comedi_fd, buffer, sizeof(buffer)) > 0) { }
        while(read(nspike_fd, buffer, sizeof(buffer)) > 0) { }
    }

    RingBuffer_reset(comedi_rb);

    acquire_active = false; 

    return true; 
}


bool comedi_init()
{
    device_handle = comedi_open(NSTREAM_COMEDI_DEVICE);
    if(!device_handle)
    {
        fprintf(stderr, "comedi: Could not open device: %s.\n",NSTREAM_COMEDI_DEVICE);
        return false;
    }
    

    return true;
}

bool nspike_init()
{
    int i = 0;
    
    nsp = NSpike_create();
    messages_fd = NSpike_get_messages_socket(nsp);
    // NSpike_set_master_address(nsp, "10.1.2.10");
    NSpike_set_dioport_count(nsp, 32);
/*  NSpike_add_auxdsp(nsp, "10.1.2.11", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.12", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.13", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.14", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.15", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.16", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.17", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.18", 30000); 
    NSpike_add_auxdsp(nsp, "10.1.2.19", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.20", 30000); 
    NSpike_add_auxdsp(nsp, "10.1.2.21", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.22", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.23", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.24", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.25", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.26", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.27", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.28", 30000);
*/
    Channel *ch;
    int hw_num;
    i = 1;
    
    int group;
    int channel;
    int division;
    int twice;
#if 0
    for (twice = 0; twice < 2; twice++)
    {
        for (group = 0; group < 4; group++)
        { 
            for (division = 0; division < 8; division++)
            {
                for (channel = 0; channel < 4; channel++)
                {
                    hw_num = (division * 16) + (group * 4) + channel;
                   
                    if (hw_num > 127)
                       hw_num -= 128;
                    if (hw_num == 127)
                        hw_num = 0;
                    ch = Channel_create(); 
                    Channel_set_number(ch, i);
                    Channel_set_hw_number(ch, hw_num);
                    Channel_set_sample_rate(ch, 30000);
                    NSpike_set_channel(nsp, ch);
 
                    // lprintf("channel %d hw_num %d\n",i,hw_num);
                    i++;
                }
            }
        }
    }  
#endif

#if 0 
    for (i = 0; i < 256; i++)
    {
        hw_num = i;

        ch = Channel_create();
        Channel_set_number(ch, i);

        if (i > 127)
           hw_num -= 128;

        if (hw_num == 127)
            hw_num = 60;
        
        Channel_set_hw_number(ch, hw_num);

        Channel_set_sample_rate(ch, 30000);
        NSpike_set_channel(nsp, ch);
        lprintf("Added channel %d hw_num %d\n", i, hw_num);
    }

    if (!NSpike_configure(nsp))
    {
        lprintf("nspike_configure() failed.  check hardware.\n");
        return false;
    }
#endif

//    NSpike_register_dio_state_change_callback(nsp, (dio_status_funcptr*)&handle_dio_data, NULL);
    return true;
}

int comedi_start_acquire()
{
    char tmpstring[200];
    int i = 0;

    cmd.subdev = comedi_find_subdevice_by_type(device_handle,COMEDI_SUBD_AI,0);
    if (cmd.subdev < 0)
    {
        fprintf(stderr, "comedi: No analog input subdevice found.\n");
        return -1;
    }

    cmd.flags = 0;

    cmd.start_src      = TRIG_EXT;
    cmd.start_arg      = NSTREAM_START_PFI | CR_EDGE;
    
    cmd.scan_begin_src = TRIG_EXT;
    cmd.scan_begin_arg = NSTREAM_CLOCK_PFI | CR_EDGE;

    cmd.convert_src    = TRIG_TIMER;
    cmd.convert_arg    = 0;

    cmd.scan_end_src   = TRIG_COUNT;
    cmd.scan_end_arg   = NSTREAM_NUM_COMEDI_CHANNELS;

    cmd.stop_src       = TRIG_NONE;
    cmd.stop_arg       = 0;

    /* Set input range for all channels to [-5,5] */
    for (i = 0; i < NSTREAM_NUM_COMEDI_CHANNELS; i++)
        channel_list[i]=CR_PACK(i,1,AREF_DIFF);

    /* Set input range for channels 0,1 to [-0.5,0.5] for sound input levels. */
    channel_list[0]=CR_PACK(0,2,AREF_DIFF);  
    channel_list[1]=CR_PACK(1,2,AREF_DIFF);  

    cmd.chanlist = channel_list;
    cmd.chanlist_len = NSTREAM_NUM_COMEDI_CHANNELS;

    fcntl(comedi_fileno(device_handle),F_SETFL,O_NONBLOCK);
    
    if(comedi_command_test(device_handle,&cmd) < 0)
    {
        fprintf(stderr, "comedi: Error testing comedi command, errno=%d\n",errno);
        return -1;
    }

    if(comedi_command(device_handle,&cmd) < 0)
    {
        fprintf(stderr,"comedi: Error executing comedi command, errno=%d\n",errno);
        return -1;
    }

    return comedi_fileno(device_handle);
}

int get_control_socket()
{
    int sock;
    int len;
    struct sockaddr_in address;
    memset(&address, 0, sizeof(address)); 
    
    address.sin_port        = htons(NSTREAM_CONTROL_SOCKET_NUMBER);
    address.sin_family      = AF_INET;
    /* security for now */
    address.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
   
    if ((sock = socket(PF_INET, SOCK_DGRAM, 0)) < 0)
    {
        perror("can't create control socket");
        exit_clean(1);
    }

    if (bind(sock, (struct sockaddr *)&address, sizeof(address)) < 0)
    {
        perror("can't bind control socket");
        exit_clean(1);
    }
    
    fcntl(sock,F_SETFL,O_NONBLOCK);
    return sock;
}
   
bool set_dac_channel(int dac_chan, int analog_chan)
{
    lprintf("setting dac channel %d to listen to data channel %d\n", dac_chan, analog_chan);
    if (!NSpike_set_dac_channel(nsp, analog_chan, dac_chan))
    {
        lprintf("failed\n");
        return false;
    }
    return true;
}

bool set_verbose(bool enable)
{
    lprintf("setting verbose state to %d\n",enable);
    verbose_state = enable;
    return true;
}

bool set_channel_filters(int data_channel, int high_pass, int low_pass)
{
    Channel *ch;
    lprintf("setting filters for channel %d highpass=%d lowpass=%d\n", data_channel, high_pass, low_pass);

    ch = NSpike_get_channel_by_number(nsp, data_channel);
    if (ch == NULL)
    {
        lprintf("unknown channel %d, failed\n", data_channel);
        return false;
    }

    Channel_set_highpass_filter(ch, high_pass);
    Channel_set_lowpass_filter(ch, low_pass);

    if (!NSpike_apply_channel_changes(nsp, ch))
    {
        lprintf("unable to apply filters for channel %d, failed\n", data_channel);
        return false;
    }

    return true;
}

bool set_dac_gain(int dac_chan, int dac_gain)
{
    lprintf("setting dac channel %d gain to %d\n", dac_chan, dac_gain);
    if (!NSpike_set_dac_gain(nsp, dac_chan, dac_gain))
    {
        lprintf("failed\n");
        return false; 
    }
    
    return true;
}

shared_state_vars *shared_state_vars_init(void)
{
  int shm_id = 0;
  int i;
  void *shared_memory = NULL;
  struct shared_state_vars *share;

  shm_id = shmget((key_t)NSTREAM_STATE_VARS_SHM_KEY, sizeof(struct shared_state_vars), 0666 | IPC_CREAT);

  /* Check share memory */
  if(shm_id == -1) { perror("shmget function failed"); exit_clean(-1); }

  /* Make shared memory accessible by program */
  shared_memory = (struct shared_state_vars*)shmat(shm_id,(void *)0, 0);
  if(shared_memory == (void *)-1) { fprintf(stderr, "shmat function failed\n"); }

  /* Make pointer to the share memory */
  share = (struct shared_state_vars *)shared_memory;

  memset(share, 0, sizeof(struct shared_state_vars));

   /* set reasonable defaults for spike extraction */
   for (i=0; i<NSTREAM_NUM_NSPIKE_CHANNELS;i++)
    { 
        share->spike_detect_positive_threshold[i] = NSTREAM_SPIKE_DETECT_DEFAULT_POS_THRESH;
        share->spike_detect_negative_threshold[i] = NSTREAM_SPIKE_DETECT_DEFAULT_NEG_THRESH;
    }

    share->spike_detect_refractory_period = NSTREAM_SPIKE_DETECT_DEFAULT_REFRACTORY_PERIOD;

  return share;
}

int set_real_time_priority(void)
{
    struct sched_param schp;
    /*
     * set the process to real-time privs
     */
    memset(&schp, 0, sizeof(schp));
    schp.sched_priority = sched_get_priority_max(SCHED_FIFO) - 1;

    if (sched_setscheduler(0, SCHED_FIFO, &schp) != 0) {
            perror("sched_setscheduler");
            return -1;
    }

     return 0;

}


void handle_dio_data(int timestamp, unsigned short data[4], void *x)
{
    int pnum, i;
    char dio_text[256] = "";
    char dio_text_for_file[256] = "";

    if (acquire_active == true)
    {
        /* 
         * write to shm 
         */
        shared_state->dio_message_buffer[shared_state->dio_message_index].nspike_global_timestamp = timestamp - start_timestamp;
        for (pnum = 0; pnum < 2; pnum++)
        {
            shared_state->dio_message_buffer[shared_state->dio_message_index].state[pnum] = data[pnum];
        }

        if (shared_state->dio_message_index + 1 == NSTREAM_DIO_STATE_NUM_ELEMENTS)
            shared_state->dio_message_index = 0;
        else
            shared_state->dio_message_index++;

        /*
         *  format for screen and disk 
         */
        for (pnum = 0; pnum < 2; pnum++) {
            for (i = 0; i < 16; i++) {
                if (data[pnum] & (1 << 15-i)) {
                    strcat(dio_text, "1");
                    strcat(dio_text_for_file, "1");
                }
                else {
                    strcat(dio_text_for_file, "0");
                    strcat(dio_text, "0");
                }
                if ((i+1)%4 == 0)
                    strcat(dio_text, " ");

                strcat(dio_text_for_file, " ");
            }
        }
        // if we're saving adjust the timestamp to the start of the recording
        if (writer_state == WRITER_ACTIVE)
        {
            fprintf(dio_out, "%d %s\n",timestamp - (recording_start_time + start_timestamp), dio_text_for_file);
            
            /* quick fix for dio message loss bug when recordings not
             * not properly stopped.  XXX may cause performance issues
             * with chatty dio */
            fflush(dio_out);

            if (verbose_state)
                lprintf("dio state change: save_adj_ts=%d %s\n",timestamp - (recording_start_time + start_timestamp), dio_text);
        }
        else
            if (verbose_state)
                lprintf("dio state change: ts=%d %s\n",timestamp - start_timestamp, dio_text);

    }
}

bool start_save()
{
    writer_state = WRITER_START;
    return true;
}


bool stop_save()
{
    writer_state = WRITER_STOP;
    return true;
}

int comedi_samples_written = 0;
int nspike_samples_written = 0;

bool check_shmem_ptrs_time()
{
    unsigned short *comedi_ptr;
    short *nspike_ptr;

    comedi_ptr = (unsigned short*)&shared_state->comedi_buffer;
    comedi_ptr += ((shared_state->comedi_time * NSTREAM_NUM_COMEDI_CHANNELS) % NSTREAM_COMEDI_NUM_SAMPLES);
    
    nspike_ptr = (short*)&shared_state->nspike_buffer;
    nspike_ptr += ((shared_state->nspike_time * NSTREAM_NUM_NSPIKE_CHANNELS) % NSTREAM_NSPIKE_NUM_SAMPLES);


    if (nspike_ptr != shared_state->nspike_buffer_write_ptr)
    {
        lprintf(">>>>>>>>>>>>>>> nspike ptr/time mismatch %u %u \n", (unsigned long)nspike_ptr, (unsigned long)shared_state->nspike_buffer_write_ptr);
    }

    if (comedi_ptr != shared_state->comedi_buffer_write_ptr)
    {
        lprintf(">>>>>>>>>>>>>>> comedi ptr/time mismatch\n");
    }

    return true;
}

unsigned long last_write_time = 0;
int write_chunk_size = 14000;
unsigned short *comedi_write_to_ptr = NULL;
short *nspike_write_to_ptr = NULL;

void writer_thread_main(void *foo)
{
    unsigned long int min_time = 0;
    unsigned long int max_time = 0;
    int i = 0;

    for (;;)
    {
        usleep(250000);
        if (writer_state == WRITER_STOPPED)
        {
            if (nstream_shutdown == true)
                exit(NULL);
        }
        if (writer_state == WRITER_START)
        {
            /* open files */
            nspike_low_out = fopen(nspike_low_filename,"w");
            nspike_high_out = fopen(nspike_high_filename,"w");
            posix_fadvise(fileno(nspike_low_out), 0, 0, POSIX_FADV_DONTNEED); 
            posix_fadvise(fileno(nspike_high_out), 0, 0, POSIX_FADV_DONTNEED); 
            setvbuf(nspike_low_out, NULL, _IOFBF, 1048600000); 
            setvbuf(nspike_high_out, NULL, _IOFBF, 1048600000); 
            comedi_out = fopen(comedi_filename,"w");
            setvbuf(comedi_out, NULL, _IOFBF, 1024000); 
            dio_out = fopen(dio_filename, "w");
            log_out = fopen(log_filename, "w");
            if (dataglove_daq != NULL) 
                dataglove_out = fopen(dataglove_filename, "w");
            comedi_samples_written = 0;
            nspike_samples_written = 0;
            writer_state = WRITER_SYNC;
        }
        if (writer_state == WRITER_SYNC)
        {
            /* set initial write pointers to point to the sample pointed to by
             * the buffer that is currently lagging behind the other */

            /* initialize both to the beginning of their respective buffers */
            writer_comedi_ptr = (unsigned short*)&shared_state->comedi_buffer;
            writer_nspike_ptr = (short *)&shared_state->nspike_buffer;

            /* set start time to the buffer that is lagging behind the most */        
            last_write_time = (shared_state->comedi_time > shared_state->nspike_time) ? shared_state->nspike_time : shared_state->comedi_time;
            recording_start_time = last_write_time;
            shared_state->recording_start_time = recording_start_time;
            lprintf("recording_start_time: %ums\n",recording_start_time);
            
            /* set pointers accordingly */
            writer_comedi_ptr += ((last_write_time * NSTREAM_NUM_COMEDI_CHANNELS) % NSTREAM_COMEDI_NUM_SAMPLES);
            writer_nspike_ptr += ((last_write_time * NSTREAM_NUM_NSPIKE_CHANNELS) % NSTREAM_NSPIKE_NUM_SAMPLES);

            writer_state = WRITER_ACTIVE; 
        }
        if (writer_state == WRITER_ACTIVE)
        {
           
            //printf("min_time: %d\tmax_time: %d\tlast_write_time: %d\n",min_time,max_time,last_write_time); 
            /* get time for buffer that is lagging the most */
            min_time = (shared_state->comedi_time > shared_state->nspike_time) ? shared_state->nspike_time : shared_state->comedi_time;

            /* get time for buffer that is leading the most */
            max_time = (shared_state->comedi_time < shared_state->nspike_time) ? shared_state->nspike_time : shared_state->comedi_time;

            /* if buffer has popped, bail */
            if (max_time - last_write_time > (NSTREAM_NSPIKE_NUM_SAMPLES / NSTREAM_NUM_NSPIKE_CHANNELS))
            {
               /* god this is ugly */
               for (;;)
               {    
                  lprintf("*** ERROR: writer fell behind, ring buffer overflow, restart nstream!\n");
                  sleep(1);
               }
            }

            if (min_time - last_write_time > write_chunk_size)
            {
                comedi_write_to_ptr = (unsigned short*)&shared_state->comedi_buffer;
                comedi_write_to_ptr += ((min_time * NSTREAM_NUM_COMEDI_CHANNELS) % (NSTREAM_COMEDI_NUM_SAMPLES));
                nspike_write_to_ptr = (short*)&shared_state->nspike_buffer;
                nspike_write_to_ptr += ((min_time * NSTREAM_NUM_NSPIKE_CHANNELS) % (NSTREAM_NSPIKE_NUM_SAMPLES));

                /* write comedi data */
                /* if the writer has wrapped, write to the end of the buffer and wrap */
                if (comedi_write_to_ptr < writer_comedi_ptr)
                {
                    /* write out the end of the buffer and wrap */
                    posix_fadvise(fileno(comedi_out), 0, 0, POSIX_FADV_DONTNEED);
                    while (writer_comedi_ptr < shared_state->comedi_buffer_end_ptr)
                    {
		    /*  Write comedi data from the first channel to comedi_num_channels_to_write */
                        fwrite(writer_comedi_ptr,2,comedi_num_channels_to_write,comedi_out);
                        /*  Skip total comedi channel count to decimate and downsample to 10kHz */
                        for (i = 0; i < comedi_writer_decimation_factor; i++) 
                            writer_comedi_ptr += NSTREAM_NUM_COMEDI_CHANNELS;
                    }
                    /* reset to start of buffer */
                    writer_comedi_ptr = shared_state->comedi_buffer;
                } 

                /* if there's data to write, catch up to the writer */
                if (comedi_write_to_ptr > writer_comedi_ptr)
                {
                    posix_fadvise(fileno(comedi_out), 0, 0, POSIX_FADV_DONTNEED);
                    while (writer_comedi_ptr < comedi_write_to_ptr)
                    {
		    /*  Write comedi data from the first channel to comedi_num_channels_to_write */
                        fwrite(writer_comedi_ptr,2,comedi_num_channels_to_write,comedi_out);
                        /*  Skip total comedi channel count to decimate and downsample to 10kHz */
                        for (i = 0; i < comedi_writer_decimation_factor; i++) 
                            writer_comedi_ptr += NSTREAM_NUM_COMEDI_CHANNELS;
                    }
                }
                
                /* write nspike data */
                /* if the writer has wrapped, write to the end of the buffer and wrap */
                if (nspike_write_to_ptr < writer_nspike_ptr)
                {
                    posix_fadvise(fileno(nspike_low_out), 0, 0, POSIX_FADV_DONTNEED);
                    posix_fadvise(fileno(nspike_high_out), 0, 0, POSIX_FADV_DONTNEED);
                    /* write out the end of the buffer and wrap */
                    while (writer_nspike_ptr < shared_state->nspike_buffer_end_ptr)
                    { 
		    /*  Write nspike data from the first channel to nspike_num_channels_to_write_low */
                        fwrite(writer_nspike_ptr,2,nspike_num_channels_to_write_low,nspike_low_out);
                        for (i = 0; i < nspike_writer_decimation_factor; i++) 
                        {
		    /*  Write nspike data for the nspike_high_channels_offset to nspike_num_channels_to_write_high */
                            fwrite(writer_nspike_ptr + nspike_high_channels_offset,2,nspike_num_channels_to_write_high,nspike_high_out);
                            writer_nspike_ptr += NSTREAM_NUM_NSPIKE_CHANNELS;
                        }
                    }
                    /* reset to start of buffer */
                    writer_nspike_ptr = shared_state->nspike_buffer;
                } 
                /* if there's data to write, catch up to the writer */
                if (nspike_write_to_ptr > writer_nspike_ptr)
                {
                    posix_fadvise(fileno(nspike_low_out), 0, 0, POSIX_FADV_DONTNEED);
                    posix_fadvise(fileno(nspike_high_out), 0, 0, POSIX_FADV_DONTNEED);
                    while (writer_nspike_ptr < nspike_write_to_ptr)
                    { 
		    /*  Write nspike data from the first channel to nspike_num_channels_to_write_low */
                        fwrite(writer_nspike_ptr,2,nspike_num_channels_to_write_low,nspike_low_out);
                        for (i = 0; i < nspike_writer_decimation_factor; i++) 
                        {
		    /*  Write nspike data for the nspike_high_channels_offset to nspike_num_channels_to_write_high */
                            fwrite(writer_nspike_ptr + nspike_high_channels_offset,2,nspike_num_channels_to_write_high,nspike_high_out);
                            writer_nspike_ptr += NSTREAM_NUM_NSPIKE_CHANNELS;
                        }
                    }
                }
                last_write_time = min_time;
            }
        }
        if (writer_state == WRITER_STOP)
        {
            if (writer_state != WRITER_STOPPED)
            {
                // lprintf("save stopped nspike_time = %d, comedi_time = %d\n",shared_state->nspike_time,shared_state->comedi_time);
                fclose(comedi_out);
                fclose(nspike_low_out);
                fclose(nspike_high_out);
                fclose(dio_out);
                fclose(log_out);
                if (dataglove_daq != NULL)
                    fclose(dataglove_out);
                max_time = 0;
                min_time = 0;
                recording_start_time = 0;
                shared_state->recording_start_time = 0;
            }
            writer_state = WRITER_STOPPED;
        }
    }
}


