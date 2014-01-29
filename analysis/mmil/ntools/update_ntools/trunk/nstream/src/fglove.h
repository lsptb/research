/*--------------------------------------------------------------------------*/
// fglove.h
//
// 5DT Data Glove driver SDK
// Version 2.1
//
// Copyright (C) 2000-2005, 5DT <Fifth Dimension Technologies>
// http://www.5dt.com/
/*--------------------------------------------------------------------------*/
#ifndef _FGLOVE_H_
#define _FGLOVE_H_
/*--------------------------------------------------------------------------*/

#include <termios.h> // port stuff
#include <pthread.h> // posix thread stuff

// For old C++ compilers ..
#ifndef bool
#define bool int
#define false (0)
#define true (!false)
#endif

#define FGLOVEAPI



#ifndef LPVOID
#define LPVOID void*
#endif

#ifndef HANDLE
#define HANDLE void*
#endif


#define USER_DEBUG 0	// level 0 = user debug
#define DEV_DEBUG 1 	// level 1 = development debug

// Uncomment the debug level you want -
#define DEBUG_LEVEL USER_DEBUG
//#define DEBUG_LEVEL DEV_DEBUG

#define VENDOR_ID_5DT 0x5d70

// Glove type / capability defines for DG14 Ultra
#define DG14U_R			0x00
#define DG14U_L			0x01
#define DG14U_RW		0x02
#define DG14U_LW		0x03
#define DG14U_RW_SLAVE	0x06
#define DG14U_LW_SLAVE	0x07
#define DG5U_R			0x10
#define DG5U_L			0x11
#define DG5U_RW			0x12
#define DG5U_LW			0x13
#define DG5U_RW_SLAVE	0x16
#define DG5U_LW_SLAVE	0x17
#define DG14U_INFO		0x80		// Mask for info packet
#define MASK_WIRELESS	0x02


// The packet size for DG14 Ultra
#define PACKET_SIZE_DG14U	29
// Serial Number length
#define LENGTH_SERIAL_NUMBER 11
// Length of date strings
#define LENGTH_DATE 8


/*--------------------------------------------------------------------------*/
/*---------------------------------- enums ---------------------------------*/
enum EfdGloveHand
{
	FD_HAND_LEFT,   // left-handed glove
	FD_HAND_RIGHT   // right-handed glove
};

enum EfdGloveTypes
{
	FD_GLOVENONE,   // no glove
	FD_GLOVE5U,		// DG5 Ultra serial
	FD_GLOVE5UW,	// DG5 Ultra serial, wireless
	FD_GLOVE5U_USB,	// DG5 Ultra USB
	FD_GLOVE7,      // 7-sensor
	FD_GLOVE7W,     // 7-sensor, wireless
	FD_GLOVE16,     // 16-sensor
	FD_GLOVE16W,    // 16-sensor, wireless
	FD_GLOVE14U,	// DG14 Ultra serial
	FD_GLOVE14UW,	// DG14 Ultra serial, wireless
	FD_GLOVE14U_USB	// DG14 Ultra USB
};

enum EfdSensors
{
	FD_THUMBNEAR=0,
	FD_THUMBFAR,
	FD_THUMBINDEX,
	FD_INDEXNEAR,
	FD_INDEXFAR,
	FD_INDEXMIDDLE,
	FD_MIDDLENEAR,
	FD_MIDDLEFAR,
	FD_MIDDLERING,
	FD_RINGNEAR,
	FD_RINGFAR,
	FD_RINGLITTLE,
	FD_LITTLENEAR,
	FD_LITTLEFAR,
	FD_THUMBPALM,
	FD_WRISTBEND,
	FD_PITCH,
	FD_ROLL
};

/*--------------------------------------------------------------------------*/

typedef struct
{

	int				m_nGloveHand;				// Left or right handed glove
	int				m_nGloveType;				// Type of currently connected glove
	int				m_nSensors;					// Number of sensors (18 for now, regardless of glove type)
	unsigned short	*m_pSensorRaw;				// Raw (A/D) sensor values
	unsigned short	*m_pUpperCal, *m_pLowerCal;	// Upper & lower calibration values
	float			*m_pSensorMax;				// Max scaled sensor value
	float			*m_pSensorScaled;			// Scaled sensor values
	float			*m_pUpperThr, *m_pLowerThr;	// Upper & lower thresholds for gesture recognition
	int				m_nGestures;				// Number of gesture that are currently defined
	char			**m_ppGestTable;			// Gesture definition table
	int				m_nGestIndex;				// Internal gesture counter
	int				m_aGestBuffer[4];			// Internal gesture buffer
	int				m_nCurrGesture;				// Current gesture being performed
	pthread_t		m_hThread;					// Posix thread handle
	int				m_idComDev;					// File handle for device
	struct termios	m_tty;						// port info
	struct termios	m_ttyBefore;				// port info
	char			m_buf[16384];				// read-ahead buffer
	int				m_iBufBegin;				// read-ahead buffer info
	int				m_iBufEnd;					// read-ahead buffer info
	unsigned char	m_aInfo[32];				// Glove info data block
	bool			m_bNewData;					// New data available indicator
	void			(*m_pCallbackFunc);			// New data arrived callback function
	LPVOID			m_pCallbackParam;			// Callback function parameters
	float			m_fPacketTime;				// Last packet period or 1/Instant packet rate
	int				m_iPacketRate;				// Packet rate over a second period
	bool			m_bPauseUpdateThread;				// Are we updating firmware
	void			*m_pCompanion;				// Wireless companion glove
	unsigned short		m_usVersion;				// USB Version number
	bool			m_bAutoCalibrate;			// Auto calibration on/off
	// Wireless stuff
	bool			m_bIsCompanion;				// IS this glove a companion ?
	bool			m_bHasCompanion;			// Does this glove own a companion wireless glove?
	bool			m_bOpenCompanion;			// Is the companion glove open?
	int			m_nPort;				// The COM port number that this glove is on - for wireless use only

} fdGlove;

/*--------------------------------------------------------------------------*/
FGLOVEAPI fdGlove *fdOpen(char *pPort);
FGLOVEAPI int   fdClose(fdGlove *pFG);
FGLOVEAPI int   fdGetGloveHand(fdGlove *pFG);
FGLOVEAPI int   fdGetGloveType(fdGlove *pFG);
FGLOVEAPI int   fdGetNumSensors(fdGlove *pFG);
FGLOVEAPI void  fdGetSensorRawAll(fdGlove *pFG, unsigned short *pData);
FGLOVEAPI unsigned short fdGetSensorRaw(fdGlove *pFG, int nSensor);
FGLOVEAPI void  fdSetSensorRawAll(fdGlove *pFG, unsigned short *pData);
FGLOVEAPI void  fdSetSensorRaw(fdGlove *pFG, int nSensor, unsigned short nRaw);
FGLOVEAPI void  fdGetSensorScaledAll(fdGlove *pFG, float *pData);
FGLOVEAPI float fdGetSensorScaled(fdGlove *pFG, int nSensor);
FGLOVEAPI int   fdGetNumGestures(fdGlove *pFG);
FGLOVEAPI int   fdGetGesture(fdGlove *pFG);
FGLOVEAPI void  fdGetCalibrationAll(fdGlove *pFG, unsigned short *pUpper, unsigned short *pLower);
FGLOVEAPI void  fdGetCalibration(fdGlove *pFG, int nSensor, unsigned short *pUpper, unsigned short *pLower);
FGLOVEAPI void  fdSetCalibrationAll(fdGlove *pFG, unsigned short *pUpper, unsigned short *pLower);
FGLOVEAPI void  fdSetCalibration(fdGlove *pFG, int nSensor, unsigned short nUpper, unsigned short nLower);
FGLOVEAPI void  fdResetCalibration(fdGlove *pFG, int nSensor);
FGLOVEAPI void  fdResetCalibration(fdGlove *pFG);
FGLOVEAPI void  fdGetSensorMaxAll(fdGlove *pFG, float *pMax);
FGLOVEAPI float fdGetSensorMax(fdGlove *pFG, int nSensor);
FGLOVEAPI void  fdSetSensorMaxAll(fdGlove *pFG, float *pMax);
FGLOVEAPI void  fdSetSensorMax(fdGlove *pFG, int nSensor, float fMax);
FGLOVEAPI void  fdGetThresholdAll(fdGlove *pFG, float *pUpper, float *pLower);
FGLOVEAPI void  fdGetThreshold(fdGlove *pFG, int nSensor, float *pUpper, float *pLower);
FGLOVEAPI void  fdSetThresholdAll(fdGlove *pFG, float *pUpper, float *pLower);
FGLOVEAPI void  fdSetThreshold(fdGlove *pFG, int nSensor, float fUpper, float fLower);
FGLOVEAPI void  fdGetGloveInfo(fdGlove *pFG, unsigned char *pData);
FGLOVEAPI void  fdGetDriverInfo(fdGlove *pFG, unsigned char *pData);
FGLOVEAPI void	fdSetCallback(fdGlove *pFG,void *pFunc,LPVOID param);
FGLOVEAPI int	fdGetPacketRate(fdGlove *pFG);
FGLOVEAPI bool	fdNewData(fdGlove *pFG);
FGLOVEAPI int	fdGetFWVersionMajor(fdGlove *pFG);
FGLOVEAPI int	fdGetFWVersionMinor(fdGlove *pFG);
FGLOVEAPI bool	fdGetAutoCalibrate(fdGlove *pFG);
FGLOVEAPI bool	fdSetAutoCalibrate(fdGlove *pFG, bool bAutoCalibrate);
FGLOVEAPI bool	fdSaveCalibration(fdGlove *pFG,const char *pFileName);
FGLOVEAPI bool	fdLoadCalibration(fdGlove *pFG,const char *pFileName);

// Global defines for wireless gloves
static bool		g_bSecondWireless[10];			// Static value to indicate presence of second glove
static void		*g_ppWireless[10];			// Other wireless glove pointer
// Additional functions
FGLOVEAPI fdGlove *fdOpen(char *pPort,int nGloveHand);
FGLOVEAPI float fdGetPacketTime(fdGlove *pFG);
FGLOVEAPI bool	fdCheckWireless(char *pPort);
FGLOVEAPI bool	fdGetSerialNumber(fdGlove *pFG, char *pData);
// Info functions
FGLOVEAPI bool	fdGetDOB(fdGlove *pFG, char *pData);
FGLOVEAPI bool	fdGetDOS(fdGlove *pFG, char *pData);

/*--------------------------------------------------------------------------*/
#endif // #ifndef _FGLOVE_H_
/*--------------------------------------------------------------------------*/
