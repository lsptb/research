#ifndef __HWDEFINES_H__
#define __HWDEFINES_H__

typedef unsigned int u32;

/* misc defines, not neccessarily hw related */
#define NSPIKE_TWOPI   6.28318530717959
#define NSPIKE_UNIX_DATA_SOCKET_NAME "libnspike-data"
#define NSPIKE_MAX_AUXDSPS 32
#define NSPIKE_AUXDSP_DATA_BUFFER_SIZE 9280000 // 
#define NSPIKE_BIT0    0x0001
#define NSPIKE_BIT1    0x0002
#define NSPIKE_BIT2    0x0004
#define NSPIKE_SOCKET_RCVBUF_SIZE 33554428
#define NSPIKE_DSP_ECHO_PORT 7
#define NSPIKE_DSP_ECHO_INTERVAL 60

/* end misc */

#define NSPIKE_MAX_CHAN_PER_DSP	16
#define NSPIKE_DSP_BASE_SAMP_RATE	((int) 30000)
#define NSPIKE_DEFAULT_LOW_FILTER		12000
#define NSPIKE_DEFAULT_HIGH_FILTER	    1

#define NSPIKE_DAC_CHANNEL_COUNT    2   // number of dac channels available from the hw
#define NSPIKE_MAX_DAC_GAIN	128
#define NSPIKE_MAX_DAC_CUTOFF	32767
#define NSPIKE_MAX_DAC_DELAY	500	// in ms
#define NSPIKE_DEFAULT_DAC_GAIN		100  
#define NSPIKE_DEFAULT_DAC_CUTOFF	0
#define NSPIKE_DEFAULT_DAC_DELAY		0
#define NSPIKE_MAX_DSP_DIO_PORTS	4  // Note that this could cause problems if there are only 2 ports (32 channels) of digitial IO

#define NSPIKE_MESSAGE_PORT 4001
#define NSPIKE_DATA_PORT    4003
#define NSPIKE_BOGUS_ECHO_PORT 9999
#define NSPIKE_SAMP_TO_TIMESTAMP	3   //divide by 3 to turn 30KHz timestamps into 100 usec timestamps

// Hardware Addresses
//

/* commands to read and write data */
#define NSPIKE_SHORT_READ		0x80     // read a few bytes
#define NSPIKE_SHORT_WRITE		0x82     // write a few bytes
#define NSPIKE_MAX_SHORT_WRITE_SIZE	1000	// we can write at most 1000 bytes at a time
#define NSPIKE_SERIAL_WRITE_READ	0x85     // write bytes and then read them back
#define NSPIKE_RESET_DEVICE		0x84	 // reset the device as if it had been power cycled

#define NSPIKE_DEVICE_CONTROL_ADDR	0xFFFB   // control device
#define NSPIKE_DSP_SRAM		0x4000  // base location for dsp memory. This is passed in as the high word of addresses when writing to the DSP 
#define NSPIKE_MAX_PACKET_SIZE                500  // the maximum size of a data buffer from the DSPs in terms of the number of short integers in the packet
#define NSPIKE_DAC_0_GAIN_ADDR		0xD1  // gain for NSPIKE_DAC 0
#define NSPIKE_DAC_1_GAIN_ADDR		0xD2  // gain for NSPIKE_DAC 1
#define NSPIKE_BLOCKS_TO_SEND_ADDR	0xD2  // the number of packets to send

/* delays for analog output. This can be used to improve the synchronization of
 * the visual and auditory traces */
#define NSPIKE_DAC_0_DELAY		0  // the delay, in samples, of analog output 0
#define NSPIKE_DAC_1_DELAY		0  // the delay, in samples, of analog output 1


/* DSP addresses */
#define	NSPIKE_DSP_CODE_REV_ADDR	0x13  // The current revision number
#define NSPIKE_RESYNC_COUNT_ADDR	0xC3  // test for synchronization error (not functional)
#define NSPIKE_DAC_0_DELAY_ADDR	0xCC  // delay, in samples, for audio output 0
#define NSPIKE_DAC_1_DELAY_ADDR	0xCD  // delay, in samples, for audio output 1
#define NSPIKE_SYNC_CONTROL_ADDR	0xCF  // control synchronization betwen dsps

/* Main DSP addresses */
#define NSPIKE_DAC_0_GAIN_ADDR		0xD1  // gain for NSPIKE_DAC 0
#define NSPIKE_DAC_1_GAIN_ADDR		0xD2  // gain for NSPIKE_DAC 1
#define NSPIKE_DAC_0_POS_THRESH_ADDR	0xD3  // the positive threshold for the diode emulation for NSPIKE_DAC 0
#define NSPIKE_DAC_0_NEG_THRESH_ADDR	0xD4  // the negative threshold for the diode emulation for NSPIKE_DAC 0
#define NSPIKE_DAC_1_POS_THRESH_ADDR	0xD5  // the positive threshold for the diode emulation for NSPIKE_DAC 1
#define NSPIKE_DAC_1_NEG_THRESH_ADDR	0xD6  // the negative threshold for the diode emulation for NSPIKE_DAC 1

/* Main DSP Digital IO address */
#define NSPIKE_MAX_DIO_PORTS       4
#define NSPIKE_DIO_OUTPUT_PORT     1
#define NSPIKE_DIO_OUT_ENABLE		0xD9   	// bit 0 = 1 -> bits 1-16 output
					// bit 1 = 1 -> bits 17-32 output
					// bit 2 = 1 -> bits 33-48 output
					// bit 3 = 1 -> bits 49-64 output
					// value of 0 = input
#define NSPIKE_DIO_OUT_1		0xDA	// control output bits 1-16
#define NSPIKE_DIO_OUT_2		0xDB	// control output bits 17-32
#define NSPIKE_DIO_OUT_3		0xDC	// control output bits 33-48
#define NSPIKE_DIO_OUT_4		0xDD	// control output bits 49-64
#define NSPIKE_DIO_IN_1		0xDE	// state of input bits 1-16
#define NSPIKE_DIO_IN_2		0xDF	// state of input bits 17-32
#define NSPIKE_DIO_IN_3		0xE0	// state of input bits 33-48
#define NSPIKE_DIO_IN_4		0xE1	// state of input bits 33-48
#define NSPIKE_DIO_N_STATE_MACHINES	4	// the total number of state machines
#define NSPIKE_DIO_STATE_SIZE		62	// instructions per state machine. The first instruction is always left at "wait forever", and the final instruction is always a jump to instruction 0, so there are 64 - 2 = 62 instructions available for programming 
#define NSPIKE_DIO_STATE_AVAILABLE	0xE7	// the index of a currently available state machine
#define NSPIKE_DIO_STATE_MACHINES_BUSY -1		// the result of reading from the NSPIKE_DIO_S_AVAILABLE pointer if no state machines are free
#define NSPIKE_DIO_STATE_ENABLE	0xED	// if bit 0 = 1, run the state machine
#define NSPIKE_DIO_STATE0_PTR		0xF0	// Pointer to current instruction for state machine 0
#define NSPIKE_DIO_STATE1_PTR		0xF3	// Pointer to current instruction for state machine 1
#define NSPIKE_DIO_STATE2_PTR		0xF6	// Pointer to current instruction for state machine 2
#define NSPIKE_DIO_STATE3_PTR		0xF9	// Pointer to current instruction for state machine 3
#define NSPIKE_DIO_STATE0_BUFFER_START	0x0200 // the start offset for the state pointer for state machine 0
#define NSPIKE_DIO_STATE1_BUFFER_START	0x0240 // the start offset for the state pointer for state machine 1
#define NSPIKE_DIO_STATE2_BUFFER_START	0x0280 // the start offset for the state pointer for state machine 2
#define NSPIKE_DIO_STATE3_BUFFER_START	0x02C0 // the start offset for the state pointer for state machine 3
#define NSPIKE_DIO_MESSAGE_SIZE	14	// 7 * sizeof(unsigned short) bytes per digital IO message packet

#define NSPIKE_DIO_IN_1_MASK		0xE2	// Mask for input 1
#define NSPIKE_DIO_IN_2_MASK		0xE3	// Mask for input 2
#define NSPIKE_DIO_IN_3_MASK		0xE4	// Mask for input 3
#define NSPIKE_DIO_IN_4_MASK		0xE5	// Mask for input 4
#define NSPIKE_DIO_IN_DEBOUNCE_TIME	0xE8	// debounce time for inputs
#define NSPIKE_DIO_PIPE_ID		0xEE	// ID for packets containing changed state information

#define NSPIKE_DIO_MAX_COMMAND_LEN	NSPIKE_DIO_STATE_SIZE // the maximum length of a dio command
#define NSPIKE_DIO_S_HALT		0xffff  // halt the state machine
#define NSPIKE_DIO_S_WAIT		0x8000  // wait for time in bits 1-15
#define NSPIKE_DIO_S_WAIT_TIME		0x7000  // wait for absolute time given in argument
#define NSPIKE_DIO_S_WAIT_MASKED_INPUT	0x6000  // wait for masked input to change
#define NSPIKE_DIO_S_WAIT_INPUT_HIGH	0x4100  // wait for input in bits 1-6 to be high
#define NSPIKE_DIO_S_WAIT_INPUT_LOW	0x4000  // wait for input in bits 1-6 to be low
#define NSPIKE_DIO_S_SET_OUTPUT_HIGH	0x3100  // set output in bits 1-6 to be high
#define NSPIKE_DIO_S_SET_OUTPUT_LOW	0x3000  // set output in bits 1-6 to be low
#define NSPIKE_DIO_S_SET_PORT		0x2000  // set port in bits 1-2 to argument
#define NSPIKE_DIO_S_JUMP_ABS		0x0100  // jump to absolute address in bits 1-8
#define NSPIKE_DIO_S_JUMP_REL		0x0000  // jump to bits 1-8 instructions relative to current address

/* aux DSP addresses */
#define NSPIKE_PIPE_ID_ADDR		0xD0  // the ID # for each DSP
#define NSPIKE_BLOCKS_TO_SEND_ADDR	0xD2  // the number of packets to send
#define NSPIKE_NUM_CHAN_ADDR		0xD3  // the number of channels on this dsp
#define NSPIKE_NUM_SAMPLES_ADDR	0xD4  // the number of samples per packet
#define NSPIKE_DECIMATION_ADDR		0xD5  // the decimation factor
#define NSPIKE_SAMPLE_COUNT_LOW_ADDR   0xD8  // control synchronization betwen dsps
#define NSPIKE_DIO_IN_DEBOUNCE_TIME	0xE8	// debounce time for inputs
#define NSPIKE_SAMPLE_COUNT_HIGH_ADDR  0xD9  // control synchronization betwen dsps

/* Main and aux DSP addresses for channel information */
#define NSPIKE_CHANNEL_SETTINGS_ADDR	0x180 // the start of the channel settings
#define NSPIKE_DSP_CHAN_ADDR_INC	24    // the number of words per channel

#endif

