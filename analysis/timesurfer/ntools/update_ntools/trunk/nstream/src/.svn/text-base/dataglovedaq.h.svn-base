#ifndef __DATAGLOVEDAQ_H__
#define __DATAGLOVEDAQ_H__
#include <stdbool.h>
#include "nstream.h"
#include "fglove.h"

class DataGloveDAQ
{
    public:
    DataGloveDAQ();
    ~DataGloveDAQ();
    bool Initialize(char *szPort);
    bool ReadData(long int timestamp, dataglove_event *dest_buffer);
    void SetDelay(int delay);
    private:
    fdGlove *glove_ptr_;
    long int last_timestamp_;
    int delay_;
};

#endif
