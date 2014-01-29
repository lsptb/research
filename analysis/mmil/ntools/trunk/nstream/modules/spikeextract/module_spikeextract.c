#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>

#include "../../src/nstream.h"

#define NSPIKE_MU_CHANNELS 128


#define DELAY 5 
#define BUFFER_SIZE (30*DELAY*NSTREAM_NUM_NSPIKE_CHANNELS)
#define POSTBUFFER_SIZE (NSTREAM_SPIKE_TOTAL_SAMPLES*NSTREAM_NUM_NSPIKE_CHANNELS)
#define POSSPIKETHRESH 700 
#define NEGSPIKETHRESH -10000
#define REFRACTORY_PERIOD 10000 

#define CIRC_SIZE (30*DELAY+NSTREAM_SPIKE_TOTAL_SAMPLES)

int number_of_samples = 30;
struct shared_state_vars *share;

main() 
{
  int shm_id = 0;
  int i = 0, j=0;
  unsigned long time = 0, timesince = 0;
  void *shared_memory = (void *) 0;
  short circ[CIRC_SIZE];
  int last_spike_time[NSPIKE_MU_CHANNELS];
  short buf[BUFFER_SIZE], postbuf[POSTBUFFER_SIZE], prebuf[POSTBUFFER_SIZE];
  unsigned long new_time = 0;
  int write_pos = 0, postwrite_pos = 0;
  unsigned long SpikeTime;

  struct timespec ts;

  for (i=0;i<POSTBUFFER_SIZE;i++) { postbuf[i] = 0; }
  for (i=0;i<CIRC_SIZE;i++) { circ[i] = 0; }
  for (i=0;i<NSPIKE_MU_CHANNELS;i++) { last_spike_time[i] = 0; }
  shm_id = shmget((key_t)NSTREAM_STATE_VARS_SHM_KEY, sizeof(struct shared_state_vars), 0666);
	
  /* Check share memory */
  if(shm_id == -1) { perror("shmget failed: "); exit(EXIT_FAILURE); }
	
  /* Make shared memory accessible by program */
  shared_memory = (struct shared_state_vars*)shmat(shm_id,(void *)0, 0);
  if(shared_memory == (void *)-1)
  { 
    fprintf(stderr, "shmat failed in %s\n", __FILE__);
    return; 
  } 
  /*  Make pointer to the share memory */
  share = (struct shared_state_vars *)shared_memory;
      
  time = share->nspike_time;
	printf("Time is %d\n",time);
  ts.tv_sec = 1;
  ts.tv_nsec = 1000000;
    //nanosleep(&ts,NULL);
    // usleep(20);

    for (;;)
  {
    ts.tv_sec = 0;
    ts.tv_nsec = 1000000;
    nanosleep(&ts,NULL);
      //printf("%d\t%d\t%d\n",time, new_time, new_time-time);
      //usleep(1);
      new_time = share->nspike_time;
      //printf("%d\t%d\t%d\n",time, new_time, new_time-time);
      if (new_time - time > 30*DELAY)
      {
        //printf("lag: %d\tnspike_time: %d\ttime: %d\n",new_time-time,new_time,time);
      /* ==================  READ DATA FROM SHM ============================== */
        write_pos = (time*NSTREAM_NUM_NSPIKE_CHANNELS)%(NSTREAM_NSPIKE_NUM_SAMPLES);
        write_pos -= 1;
        if (write_pos < 0)
            write_pos = NSTREAM_NSPIKE_NUM_SAMPLES + write_pos;
        postwrite_pos = write_pos;
           
        for(i=0; i < BUFFER_SIZE; i++)
        { 
         buf[BUFFER_SIZE - 1 - i] = share->nspike_buffer[write_pos];
         write_pos -=1; 	 
         if(write_pos < 0) { 	 
           write_pos = NSTREAM_NSPIKE_NUM_SAMPLES + write_pos; 
         } 
        } 

        memmove(&prebuf,&postbuf,POSTBUFFER_SIZE*sizeof(short));

        for(i=0; i < POSTBUFFER_SIZE; i++)
        {
         postbuf[POSTBUFFER_SIZE - 1 - i] = share->nspike_buffer[postwrite_pos];
         postwrite_pos -=1;
         if(postwrite_pos < 0) {
           postwrite_pos = NSTREAM_NSPIKE_NUM_SAMPLES + postwrite_pos;
         }
        }

        for (i=NSPIKE_MU_CHANNELS;i<NSTREAM_NUM_NSPIKE_CHANNELS;i++)
        {
            for(j=0;j<NSTREAM_SPIKE_TOTAL_SAMPLES;j++)
            {
                circ[j] = prebuf[i + NSTREAM_NUM_NSPIKE_CHANNELS*j];
            }
            for(j=0;j<DELAY*30;j++)
            {
                circ[NSTREAM_SPIKE_TOTAL_SAMPLES + j] = buf[i + NSTREAM_NUM_NSPIKE_CHANNELS*j];
            }

            /* ================== PROCESS IT ============================== */
            for(j=NSTREAM_SPIKE_PRESAMPLE;j<CIRC_SIZE-NSTREAM_SPIKE_POSTSAMPLE;j++)
            {
                if ((circ[j] < share->spike_detect_negative_threshold[i] || 
                    circ[j] > share->spike_detect_positive_threshold[i]))
                {
                   /* printf("time from last spike: %d\n", time - NSTREAM_SPIKE_TOTAL_SAMPLES + j - last_spike_time[i]); */
                    timesince = time - NSTREAM_SPIKE_TOTAL_SAMPLES + j - last_spike_time[i - NSPIKE_MU_CHANNELS];
                    if (timesince > share->spike_detect_refractory_period)
                    {
                        printf("Time since last spike: %d\n",timesince);
                        SpikeTime =  time - NSTREAM_SPIKE_TOTAL_SAMPLES + j;
                        printf("Spike detected at %d on channel %d\n",SpikeTime, i);
                        update_spike_buffer(SpikeTime, i, &circ[j-NSTREAM_SPIKE_PRESAMPLE]);
                        last_spike_time[i- NSPIKE_MU_CHANNELS] = SpikeTime;
                    
                    }
                }
            }
        } 
        time = new_time;
      } 
      else 
      { 
      }
  }

}

void update_spike_buffer(int timestamp, int channel, short *wave)
{
    int index = share->spike_index[channel];
    int i; 
    share->spike_event_buffer[channel][index].timestamp = timestamp;
    share->spike_event_buffer[channel][index].cell_id   = 1;
    /*memcpy(&share->spike_event_buffer[channel][index].waveform, 
           wave, NSTREAM_SPIKE_TOTAL_SAMPLES *sizeof(short));
    */
    
    for (i = 0; i < NSTREAM_SPIKE_TOTAL_SAMPLES; i++)
        share->spike_event_buffer[channel][index].waveform[i] = wave[i];
 
    share->spike_index[channel]++;
    if (index == NSTREAM_SPIKE_NUM_ELEMENTS)
        share->spike_index[channel] = 0;

    printf("Channel %d: Index: %d\n",channel,share->spike_index[channel]);
}

/*------------------------------------------------------------*/

