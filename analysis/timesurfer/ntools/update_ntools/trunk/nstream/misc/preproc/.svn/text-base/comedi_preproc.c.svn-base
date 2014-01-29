/*---------------------------------------------------------------------------
 * broker_child.c - Child process that is forked from broker. This process
 * will read stdin in this case pipe from broker and calculate status and
 * time and write it into share memory circula buffer.
 * IPC communications used 1.) Semaphores 2.) Share Memory 3.) Circular Buffer
 *---------------------------------------------------------------------------*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/*---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
/*---------------------------------------------------------------------------*/
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include "myshare.h"

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

/*---------------------------------------------------------------------------*/
int median(int *circ);
int min(int *circ);
int max(int *circ);
void mt(int trial, int stat, int time);
int state_convert(int analog_state);
void bubbleSort(int a[], int n);
/*---------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{

  int i, c;
  int eye_fd, touch_fd, state_fd, joystick_fd, displayraw_fd;
  FILE *state_fp;

  /*------- Circular Buffer ------------*/
  ushort *buf;
  int k, m, n, val;
  float stim_time;
  
  /*--------- Preprocess Vars. ---------*/
  int circ[HIST];
  /* char achlist[CH-RAW_CHANNEL_START]; */
  int analog_state, digital_state;
  int max_state, min_state;
  int old_state=0, trial=0, time=0, j=0, prev_stim=0;
  long int   data_len = HIST*CH*SS;
 
  /*------ set output filenames --------------------
  char name[10];
  sscanf(&argv[1], "%s", name); 
  printf("Name is %s \n", name);*/

  /*-----  Initialization -------------------------*/
  /* memset(achlist,'\0',sizeof(achlist)); */

  /*------ Parse input options ---------------------*/
  while ((i=getopt(argc,argv,"c:f")) !=EOF)
    {
      switch(i)
	{
	  /*	case 'c':
	  strcpy(achlist,optarg);	  break;
	default:
	exit(1); */
	} /* end switch statement */
    }  /* end while statement */
  
  /*       for (c = RAW_CHANNEL_START; c < CH; c++)
      {
        if (strncmp(&achlist[c-RAW_CHANNEL_START],"1",1)==0){
	printf("Printing channel %d\n", c-RAW_CHANNEL_START);
	}
      }
  */  
  /*------ intialize eye position output file ------*/
  eye_fd = creat(EYE_FILE, 0666);
  if(eye_fd == -1) { fprintf(stderr, "Open %s failure\n", EYE_FILE);  exit(EXIT_FAILURE);}

  /*------ intialize touch position output file ------*/
  touch_fd = creat(TOUCH_FILE, 0666);
  if(touch_fd == -1) { fprintf(stderr, "Open %s failure\n", TOUCH_FILE);  exit(EXIT_FAILURE);}

  /*------ intialize joystick data output file ------*/
  joystick_fd = creat(JOYSTICK_FILE, 0666);
  if(joystick_fd == -1) { fprintf(stderr, "Open %s failure\n", JOYSTICK_FILE);  exit(EXIT_FAILURE);}

  /*------ intialize raw state output file ------*/
  state_fd = creat(RAWSTATE_FILE, 0666);
  if(state_fd == -1) { fprintf(stderr, "Open %s failure\n", RAWSTATE_FILE);  exit(EXIT_FAILURE);}

  /*------ intialize display raw data output file ------*/
  displayraw_fd = creat(DISPLAYRAW_FILE, 0666);
  if(displayraw_fd == -1) { fprintf(stderr, "Open %s failure\n", DISPLAYRAW_FILE);  exit(EXIT_FAILURE);}

  /*------ intialize state output file ------*/
  state_fp = fopen(EV_FILE, "w");
  if(state_fp == NULL) { fprintf(stderr, "Open %s failure\n", EV_FILE);  exit(EXIT_FAILURE);}
  
  buf=malloc(data_len);
   
  while ((k=read(0, buf, data_len)) != 0) 
  {
    i = 0;
    time++;
    /* fprintf(stderr,"%d %d\n",HIST,sizeof(circ)); */
    for (n = 0; n < k/SS; n+=CH)
    {
      circ[i] = buf[n+STATE_CHANNEL]; // put in buffer
      write(state_fd,&buf[n+STATE_CHANNEL],SS); // write to state file
      write(displayraw_fd,&buf[n+DISPLAYRAW_CHANNEL],SS); // write to display raw file

      if (i==0) 
      {
	    write(eye_fd,&buf[n+EYE_CHANNEL_HOR],SS);
	    write(eye_fd,&buf[n+EYE_CHANNEL_VER],SS);
	    write(touch_fd,&buf[n+TOUCH_CHANNEL_HOR],SS);
	    write(touch_fd,&buf[n+TOUCH_CHANNEL_VER],SS);
        write(joystick_fd,&buf[n+JOYSTICK_CHANNEL_HOR],SS);
        write(joystick_fd,&buf[n+JOYSTICK_CHANNEL_VER],SS);
      }
      i++;   
    }
   
    analog_state = median(circ);
    min_state = min(circ);
    max_state = max(circ);
    digital_state = state_convert(analog_state);
 
    if(max_state-min_state > 30){digital_state = old_state;}
    //    printf("Digital state: %d\n", digital_state); 

    switch(digital_state)
    { //------------------ transition ----------------------//
      case 0:  if (old_state != 0)  { 
        trial ++;
	fprintf(state_fp, "%d\t%d\t%d\n", trial, digital_state, time);
      } old_state = 0;  break;
      case 1 ... 63 :                                            
        if (digital_state != old_state)
        {
	  fprintf(state_fp, "%d\t%d\t%d\n", trial, digital_state, time);
        }
      old_state = digital_state;
      break;
    } // end of switch stmt
  }  /* ======================== END OF WHILE LOOP ==============================*/

  return EXIT_SUCCESS;
}
/*------------------------------------------------------------*/
/* write data for each transition */
void mt(int trial, int stat, int time)
{
  /*  printf("Trial: %d\tState: %d\tTime: %d\t\n", trial, stat, time); */
  fprintf(stdout, "%d\t%d\t%d\n", trial, stat, time);
}
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
  return out;
}
/*------------------------------------------------------------*/
/* returns median of state buffer */
int median(int *circ)
{
  int out;

  bubbleSort(circ, HIST);
  out = circ[HIST/2];
  return out;
}
/*------------------------------------------------------------*/
/* returns max of state buffer */
int max(int *circ)
{
  int out;

  bubbleSort(circ, HIST);
  out = circ[HIST-1];
  return out;
}
/*------------------------------------------------------------*/
/* returns min of state buffer */
int min(int *circ)
{
  int out;

  bubbleSort(circ, HIST);
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






