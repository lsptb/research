/*--------------------------------------------------------------------------*/
// A simple cross-platform console application to test the glove
//
// WIN32 must be defined when compiling for Windows.
// For Visual C++ this is normally already defined.
//
// Copyright (C) 2000, 5DT <Fifth Dimension Technologies>
/*--------------------------------------------------------------------------*/
#include "fglove.h"
#include "dataglovedaq.h"
#include "nstream.h"
#include <stdio.h>
#include <string.h>
#include <unistd.h>  // for usleep

DataGloveDAQ::DataGloveDAQ()
{
    glove_ptr_      = NULL;
    last_timestamp_ = 0;
    delay_          = 600;      // default to 50hz sample rate (30000/50)
}

/*--------------------------------------------------------------------------*/
bool DataGloveDAQ::Initialize(char *szPort)
{
	int      glovetype = FD_GLOVENONE;
	int      i;

	// szPort = "/dev/ttyS5";

	// Initialize glove
	lprintf( "initializing data glove on %s\n", szPort );
	if (NULL == (glove_ptr_ = fdOpen( szPort )))
	{
		lprintf( "failed.\n" );
		return false;
	}

	//---------------------------------------------------------------------------------------
	
	char *szType = "?";
	glovetype = fdGetGloveType(glove_ptr_);
	switch (glovetype) {
	case FD_GLOVENONE:    szType = "None"; break;
	case FD_GLOVE5U:      szType = "Data Glove 5 Ultra"; break;
	case FD_GLOVE5UW:     szType = "Data Glove 5 Ultra W"; break;
	case FD_GLOVE5U_USB:  szType = "Data Glove 5 Ultra USB"; break;
	case FD_GLOVE7:       szType = "Data Glove 5"; break;
	case FD_GLOVE7W:      szType = "Data Glove 5W"; break;
	case FD_GLOVE16:      szType = "Data Glove 16"; break;
	case FD_GLOVE16W:     szType = "Data Glove 16W"; break;
	case FD_GLOVE14U:     szType = "DG14 Ultra serial"; break;
	case FD_GLOVE14UW:    szType = "DG14 Ultra W"; break;
	case FD_GLOVE14U_USB: szType = "DG14 Ultra USB"; break;
	}
	
	lprintf( "data glove type: %s\n", szType );
	lprintf( "data glove handedness: %s\n", fdGetGloveHand(glove_ptr_)==FD_HAND_RIGHT?"right":"left" );
    return true;
}
	
	//---------------------------------------------------------------------------------------
	// Now continuously display the sensor data

bool DataGloveDAQ::ReadData(long int timestamp, dataglove_event *dest_buffer)
{	
	float gloveA_scaled[18];
    
    if (timestamp - last_timestamp_ < delay_)
        return false;

    last_timestamp_ = timestamp;
    fdGetSensorScaledAll(glove_ptr_, gloveA_scaled);

/* 
    printf("A:%0.1f %0.1f||%0.1f||%0.1f %0.1f||%0.1f||%0.1f %0.1f||%0.1f||%0.1f %0.1f||%0.1f||%0.1f %0.1f",
    gloveA_scaled[FD_THUMBNEAR],
    gloveA_scaled[FD_THUMBFAR],
    gloveA_scaled[FD_THUMBINDEX],
    gloveA_scaled[FD_INDEXNEAR],
    gloveA_scaled[FD_INDEXFAR],
    gloveA_scaled[FD_INDEXMIDDLE],
    gloveA_scaled[FD_MIDDLENEAR],
    gloveA_scaled[FD_MIDDLEFAR],
    gloveA_scaled[FD_MIDDLERING],
    gloveA_scaled[FD_RINGNEAR],
    gloveA_scaled[FD_RINGFAR],
    gloveA_scaled[FD_RINGLITTLE],
    gloveA_scaled[FD_LITTLENEAR],
    gloveA_scaled[FD_LITTLEFAR]);
*/
    
    // printf(" >> %d\n", fdGetGesture(glove_ptr_));
    dest_buffer->timestamp = timestamp;
    dest_buffer->thumb     = gloveA_scaled[FD_THUMBNEAR];
    dest_buffer->index     = gloveA_scaled[FD_INDEXNEAR];
    dest_buffer->middle    = gloveA_scaled[FD_MIDDLENEAR];
    dest_buffer->ring      = gloveA_scaled[FD_RINGNEAR];
    dest_buffer->little    = gloveA_scaled[FD_LITTLENEAR];
    dest_buffer->gesture   = fdGetGesture(glove_ptr_); 

    return true; 
}

void DataGloveDAQ::SetDelay(int delay)
{
    delay_ = delay;
}

DataGloveDAQ::~DataGloveDAQ()
{
	// Close glove
    if (glove_ptr_ != NULL)
    {
        lprintf( "closing data glove...\n" );
        fdClose( glove_ptr_ );
        lprintf( "data glove closed.\n" );
    }
}
