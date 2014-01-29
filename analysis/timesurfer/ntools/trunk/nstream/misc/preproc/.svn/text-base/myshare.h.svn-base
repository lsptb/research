/*----------------------- Share Memory Declaration ------------------------------*/
#include <sys/shm.h>
#define SS 2                  /* sample size */
#define PERIOD 75000
#define VERBOSE 0		/* 0 = not verbose, 1 = verbose */
/*----------------------- Broker Child path -------------------------------------*/
/*----------------------- Broker Child path -------------------------------------*/
#define COMEDI_DEVICE "/dev/comedi0"
/*----------------------- Channel Information -----------------------------------*/
#define CH     8		/* total number of channels recording from         */
#define STATE_CHANNEL 0	/* which channel of 0..CH-1 to look for state data */
#define EYE_CHANNEL_HOR 3	/* which channels to look for eye position data    */
#define EYE_CHANNEL_VER 4	/* which channels to look for eye position data    */
#define TOUCH_CHANNEL_HOR 1	/* which channels to look for hand position data    */
#define TOUCH_CHANNEL_VER 2	/* which channels to look for hand position data    */
#define JOYSTICK_CHANNEL_HOR 5	/* which channel to look for joystick position data    */
#define JOYSTICK_CHANNEL_VER 6	/* which channel to look for joystick position data */
#define DISPLAYRAW_CHANNEL 7	/* which channel to look for joystick position data */
#define JOYSTICK_OFFSET 6 /*  needed for comedi shared memory mex functions for reading backwards */

/* ---------------------- Sampling Rates ----------------------------------------*/
#define SR 30000	   /* Sampling Rate                                   */
#define EYE_SR 1000	       /* Sampling rate of eye position                   */
#define TOUCH_SR 1000	       /* Sampling rate of hand position                   */
#define JOY_SR 1000    /* Sampling rate of joystick position */
#define DISPLAY_SR SR          /* Sampling rate of display data        */         
                                                                                
#define DUR 100			/* Size of circular buffer in seconds              */
#define HIST SR/EYE_SR 		/* history buffer size in elements                 */
#define DT SR/EYE_SR 		/* size of stim buffer size in elements */
#define STIM_HIST 2             /* history buffer size in elements                 */

#define SHARE_SEM 5215        /* The Key to create sempahore                     */
#define SHARE_KEY 6023        /* The Key to create share memory                  */
/*-------------------------------------------------------------------------------*/
#define RAW1_FILE "file.raw1.dat"
#define RAW2_FILE "file.raw2.dat"
#define RAW3_FILE "file.raw3.dat"
#define RAW4_FILE "file.raw4.dat"
#define EYE_FILE "file.eye.dat"
#define TOUCH_FILE "file.hnd.dat"
#define JOYSTICK_FILE "file.joy.dat"
#define RAWSTATE_FILE "file.state.dat"
#define DISPLAYRAW_FILE "file.display.dat"
#define EV_FILE "file.ev.txt"
/*-------------------------------------------------------------------------------*/
#define STATE_OUTER 3  		/* Width of state circular buffer (trial,state,time) */
#define STATE_INNER 1000  	/* Length of state circular buffer                 */

#define EYE_OUTER 2    		/* Width of eye position circular buffer (x,y)     */
#define EYE_INNER DUR*EYE_SR  /* Length of eye position circular buffer          */
#define TOUCH_OUTER 2    		/* Width of hand position circular buffer (x,y)     */
#define TOUCH_INNER DUR*TOUCH_SR  /* Length of hand position circular buffer          */
#define JOY_OUTER 2           /* Width of hand position circular buffer (x,y)     */
#define JOY_INNER DUR*JOY_SR  /* Length of hand position circular buffer          */
#define DISPLAY_OUTER 2           /* Width of display circular buffer (trial,time) */
#define DISPLAY_INNER 1000        /* Length of display circular buffer                 */
#define DISPLAYRAW_OUTER 1
#define DISPLAYRAW_INNER DUR*DISPLAY_SR  /* Length of display circular buffer                   */


/*-------------------------------------------------------------------------------*/
int shm_id;
typedef struct shared_use_st 
{
	short int file_state;
	char file_name[128];
	long  int  time;
	unsigned short int  state_pos;
	unsigned long  int  eye_pos;
	unsigned long  int  touch_pos;
	unsigned long  int  joy_pos;
    unsigned long  int  displayraw_pos;    
    unsigned long  int  state[STATE_OUTER][STATE_INNER];
	unsigned short int  eye[EYE_OUTER][EYE_INNER];
	unsigned short int  touch[TOUCH_OUTER][TOUCH_INNER];
	unsigned short int  joy[JOY_OUTER][JOY_INNER];
    unsigned long  int  displayraw[DISPLAYRAW_OUTER][DISPLAYRAW_INNER];
    unsigned short int  check;
}shared_use_st;
/*----------------------------------------------------------------------*/
