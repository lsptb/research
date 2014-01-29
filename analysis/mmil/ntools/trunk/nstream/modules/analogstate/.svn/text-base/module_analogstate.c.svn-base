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

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
void update_analog_state_buffer(unsigned long trial, unsigned long digital_state, unsigned long time);
int state_convert(int analog_state);
void bubbleSort(int a[], int n);
int median(int *circ);
int min(int *circ);
int max(int *circ);

#define DELAY 20
#define BUFFER_SIZE 30*DELAY
int number_of_samples = 30;
struct shared_state_vars *share;

main() 
{
  int shm_id = 0;
  int i = 0, j=0;
  unsigned long time = 0;
  void *shared_memory = (void *) 0;
  int circ[number_of_samples];
  int buf[BUFFER_SIZE];
  int analog_state, digital_state;
  int max_state, min_state;
  int old_state = 0, trial = 0;
  unsigned long new_time = 0;
  int write_pos = 0;

  struct timespec ts;

  for (i=0;i<number_of_samples;i++) { circ[i] = 0; }
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
      
  time = share->comedi_time;
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
      new_time = share->comedi_time;
      //printf("%d\t%d\t%d\n",time, new_time, new_time-time);
      if (new_time - time > BUFFER_SIZE)
      {
        //printf("lag: %d\tcomedi_time: %d\ttime: %d\n",new_time-time,new_time,time);
      /* ==================  READ DATA FROM SHM ============================== */
        write_pos = (time*NSTREAM_NUM_COMEDI_CHANNELS)%(NSTREAM_COMEDI_NUM_SAMPLES);
        write_pos -= 8;
        if (write_pos < 0)
            write_pos = NSTREAM_COMEDI_NUM_SAMPLES + write_pos;
           
        for(i=0; i < BUFFER_SIZE; i++)
        { 
         buf[BUFFER_SIZE - 1 - i] = share->comedi_buffer[write_pos];
         write_pos -=8; 	 
         if(write_pos < 0) { 	 
           write_pos = NSTREAM_COMEDI_NUM_SAMPLES + write_pos; 
         } 
        } 

        for (i=0;i<BUFFER_SIZE;i+=number_of_samples)
        {
            for(j=0;j<number_of_samples;j++)
            {
                circ[j] = buf[j + i];
            }

            /* ================== PROCESS IT ============================== */
            //for (i=0; i<30; i++) { printf("%d\t",circ[i]); }
            //printf("\n================================\n");
            bubbleSort(&circ,number_of_samples);
            //for (i=0; i<30; i++) { printf("%d\t",circ[i]); }
            //printf("\n********************************"); 
            analog_state = circ[number_of_samples/2];
            min_state = circ[0];
            max_state = circ[number_of_samples-1];

            digital_state = state_convert(analog_state);
         
            //printf("Analog state: %d\n", analog_state); 
            //printf("Digital state: %d\n", digital_state); 
            //printf("Min and max state: %d\t%d\n",min_state,max_state);
            if(max_state-min_state > 30){digital_state = old_state;}

            switch(digital_state)
            { //------------------ transition ----------------------//
              case 0:  
                if (digital_state != old_state)  
                { 
                  trial ++;
                  printf("%d\t%d\t%d\t%d\t%d\n", trial, digital_state, time);
                  update_analog_state_buffer(trial, digital_state, time);
                  
                } 
                old_state = digital_state;  
                break;
              case 1 ... 63 :                                            
                if (digital_state != old_state)
                {
                  printf("%d\t%d\t%d\t%d\t%d\n", trial, digital_state, time);
                  update_analog_state_buffer(trial, digital_state, time);
                }
                old_state = digital_state;
                break;
            } // end of switch stmt */
            time = time + number_of_samples;
        } 
      } 
      else 
      { 
    //    printf("*** lag: %d\tcomedi_time: %d\ttime: %d\n",new_time-time,new_time,time);
    //    printf("******************%d\t%d\t%d\n",time, new_time, new_time-time);
    //    printf("caughtup!\n");
      }
  }

}

void update_analog_state_buffer(unsigned long trial, unsigned long digital_state, unsigned long time)
{
    share->analog_state_buffer[0][share->analog_state_index] = trial;
    share->analog_state_buffer[1][share->analog_state_index] = digital_state;
    share->analog_state_buffer[2][share->analog_state_index] = time;

    if (share->analog_state_index + 1 == NSTREAM_ANALOG_STATE_NUM_SAMPLES)
        share->analog_state_index = 0;
    else
        share->analog_state_index++;
}

/*------------------------------------------------------------*/
/* returns median of state buffer */
int median(int *circ)
{
  int out;

  bubbleSort(circ, number_of_samples);
  out = circ[number_of_samples/2];
  return out;
}
/*------------------------------------------------------------*/
/* returns max of state buffer */
int max(int *circ)
{
  int out;

  bubbleSort(circ, number_of_samples);
  out = circ[number_of_samples-1];
  return out;
}
/*------------------------------------------------------------*/
/* returns min of state buffer */
int min(int *circ)
{
  int out;

  bubbleSort(circ, number_of_samples);
  out = circ[0];
  return out;
}

/*------------------------------------------------------------*/
/* It sorts in non-decreasing order the first N positions of A. */
/* It uses the bubble sort method. */
void bubbleSort(int a[], int n)
{
  int lcv;
  int limit = n-1;
  int temp;
  int lastChange;

  while (limit)
  {
    lastChange = 0;
    for (lcv=0;lcv<limit;lcv++)
      /* Notice that the values in positions LIMIT+1 .. N are in
       * their final position, i.e. they are sorted right */
        if (a[lcv]>a[lcv+1]) {
          temp = a[lcv];
          a[lcv] = a[lcv+1];
          a[lcv+1] = temp;
          lastChange = lcv;
        }
    limit = lastChange;
  }
}
/*---------------------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* converts analog state to digital state */
int state_convert(int analog_state)
{
  float del, offset, m;
  int out = 0;
  float Range = 65536;
  float Cardinality = 64;

  del = Range/Cardinality;
  offset = del/2;

  out = round((analog_state-offset)/(del));
  //printf("*********************%f\t%f\t%d\n",del,offset,out);
  return out;
}

