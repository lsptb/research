#ifndef __EDFCONVERT_H__
#define __EDFCONVERT_H__

/* defaults for main header */
#define EH_DEFAULT_VERSION "0" 
#define EH_DEFAULT_START_DATE "01.01.70"
#define EH_DEFAULT_START_TIME "00.00.00" 

/* defaults for signal header */
#define EH_S_DEFAULT_TRANSDUCER	"Active electrode"
#define EH_S_DEFAULT_PHYS_DIM "uV"
#define EH_S_DEFAULT_PHYS_MIN "-3276"
#define EH_S_DEFAULT_PHYS_MAX "3276" 
#define EH_S_DEFAULT_DIG_MIN "-32768" 
#define EH_S_DEFAULT_DIG_MAX "32767" 
#define EH_S_DEFAULT_PREFILTER "No Filter" 
#define EH_S_DEFAULT_RESERVED "" 

#define DATA_RECORD_LEN		.1  /* duration of data records in seconds */

/* end user tweakable parameters */

#define EH_VERSION_LEN              8
#define EH_LOCAL_PATIENT_ID_LEN     80 
#define EH_LOCAL_RECORDING_ID_LEN   80 
#define EH_START_DATE_LEN           8 
#define EH_START_TIME_LEN           8 
#define EH_HEADER_BYTES_LEN         8 
#define EH_RESERVED_LEN             44 
#define EH_NUM_DATA_RECORDS_LEN     8 
#define EH_DATA_RECORD_DURATION_LEN 8 
#define EH_NUM_SIGNALS_LEN          4 

#define EH_MAX_BUFFER_LEN           80

struct edf_header {
	char version[EH_VERSION_LEN];
	char local_patent_id[EH_LOCAL_PATIENT_ID_LEN];
	char local_recording_id[EH_LOCAL_RECORDING_ID_LEN];
	char start_date[EH_START_DATE_LEN];
	char start_time[EH_START_TIME_LEN];
	char header_bytes[EH_HEADER_BYTES_LEN];
	char reserved[EH_RESERVED_LEN];
	char num_data_records[EH_NUM_DATA_RECORDS_LEN];
	char data_record_duration[EH_DATA_RECORD_DURATION_LEN];
	char num_signals[EH_NUM_SIGNALS_LEN];
};

typedef struct edf_header edf_header_t;


#define EH_S_LABEL_LEN       16
#define EH_S_TRANSDUCER_LEN  80
#define EH_S_PHYS_DIM_LEN    8
#define EH_S_PHYS_MIN_LEN    8
#define EH_S_PHYS_MAX_LEN    8
#define EH_S_DIG_MIN_LEN     8 
#define EH_S_DIG_MAX_LEN     8 
#define EH_S_PREFILTER_LEN   80
#define EH_S_NUM_SAMPLES_LEN 8 
#define EH_S_RESERVED_LEN    32 

#define EH_S_TOTAL_LEN (EH_S_LABEL_LEN + EH_S_TRANSDUCER_LEN + \
EH_S_PHYS_DIM_LEN + EH_S_PHYS_MIN_LEN + EH_S_PHYS_MAX_LEN + EH_S_DIG_MIN_LEN \
+ EH_S_DIG_MAX_LEN + EH_S_PREFILTER_LEN + EH_S_NUM_SAMPLES_LEN \
+ EH_S_RESERVED_LEN)

struct edf_signal_header {
	char *buffer;
	int length;
};

typedef struct edf_signal_header edf_signal_header_t;

char *pad_copy(char *dst, char *src, int len);
edf_header_t *build_edf_header(xml_header_t *x_hdr, int signal_header_len);
edf_signal_header_t *build_signal_header(xml_header_t *x_hdr);
bool copy_frame(FILE *in, FILE *out, int samplecount, int num_chan);
 
#endif
