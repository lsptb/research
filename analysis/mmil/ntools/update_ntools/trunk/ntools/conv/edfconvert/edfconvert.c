#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hdrxml.h"
#include "edfconvert.h"

int main(int argc, char *argv[])
{
	xml_header_t *x_hdr;
	char *new_filename;

	if (argc != 2) {
		fprintf(stderr, "Usage: %s hdrxmlfile\n", argv[0]);
		return EXIT_FAILURE;
	}

	fprintf(stdout, "edfconvert $Rev: 119 $ \n");
	fprintf(stdout, "Parsing input file: \n");
	if ((x_hdr = parse_hdrxml(argv[1])) == NULL) {
		fprintf(stderr, "Unable to parse hdrxml file %s\n", argv[1]);
		return EXIT_FAILURE;
	}
	
	FILE *in;
	FILE *out;
	if ((in = fopen(x_hdr->fname, "r")) == NULL) {
		perror("Unable to open source data file");
		return EXIT_FAILURE;
	}

	edf_signal_header_t *es_hdr = build_signal_header(x_hdr);
	edf_header_t *e_hdr = build_edf_header(x_hdr, es_hdr->length);

	new_filename = (char*)malloc(strlen(x_hdr->fname));

	if (!new_filename) {
		perror("malloc");
		return EXIT_FAILURE;
	}
		
	strncpy(new_filename, x_hdr->fname, strlen(x_hdr->fname));
	strncpy(new_filename + (strlen(new_filename) - 4), ".edf", 4);
	
	if ((out = fopen(new_filename, "w")) == NULL) {
		perror("Unable to open destination data file");
		return EXIT_FAILURE;
	}

	if (fwrite(e_hdr, 1, sizeof(edf_header_t), out) != sizeof(edf_header_t)) {
		perror("Error writing EDF header");
		return EXIT_FAILURE;
	}
	
	if (fwrite(es_hdr->buffer, 1, es_hdr->length, out) != es_hdr->length) {
		perror("Error writing EDF signal headers");
		return EXIT_FAILURE;
	}
	

	int frame = 0;
	int frames_to_copy = (int)(floor(x_hdr->duration) / DATA_RECORD_LEN);

	fprintf(stdout, "output file: %s\n", new_filename);
	fprintf(stdout, "%d frames to convert\n", frames_to_copy);
	for (frame = 0; frame < frames_to_copy; frame++) {
		if (!copy_frame(in, out, (x_hdr->sfreq * DATA_RECORD_LEN), 
		x_hdr->num_channels)) {
			fprintf(stderr, "Failed\n");
			return EXIT_FAILURE;
		}
		if ((frame % 100) == 0)
			fprintf(stdout, "%d frames converted.\n", frame);
	}

	fprintf(stdout, "%d frames converted total.  Done.\n",frame);
	fclose(in);
	fclose(out);	
	return EXIT_SUCCESS;	

}

bool copy_frame(FILE *in, FILE *out, int samplecount, int num_chan)
{
	short *in_buf;
	short *in_buf_rp;
	short *out_buf;
	short *out_buf_wp;

	if ((in_buf = (short*)malloc(samplecount*sizeof(short)*num_chan)) == NULL) {
		fprintf(stderr, "malloc() failed for in_buf\n");
		return false;
	}

	if ((out_buf = (short*)malloc(samplecount*sizeof(short)*num_chan)) == NULL) {
		fprintf(stderr, "malloc() failed for out_buf\n");
		return false;
	}

	if (fread(in_buf, sizeof(short), samplecount*num_chan, in) != 
		samplecount*num_chan) {
		perror("fread in_buf");
		return false;
	}

	memset(out_buf, 0, sizeof(out_buf));
	
	in_buf_rp = in_buf;
	out_buf_wp = out_buf;

	int chan;
	for (chan = 0; chan < num_chan; chan++) {
		in_buf_rp = in_buf + chan;
		do {
			*out_buf_wp = *in_buf_rp;
			
			out_buf_wp++;	
			in_buf_rp += num_chan;
		} while (in_buf_rp < (in_buf + (samplecount*num_chan)));
	}


	if (fwrite(out_buf, sizeof(short), samplecount*num_chan, out) != 
		samplecount*num_chan) {
		perror("fwrite out_buf");
		return false;
	}

	free(in_buf);
	free(out_buf);

	return true;
}
	
edf_header_t *build_edf_header(xml_header_t *x_hdr, int signal_header_len) {
	edf_header_t *e_hdr;

	if ((e_hdr = (edf_header_t*)malloc(sizeof(edf_header_t))) == NULL)
		return NULL;

	char tmp[EH_MAX_BUFFER_LEN];
	memset(&tmp, 0, sizeof(tmp));

	memset(e_hdr, 0, sizeof(edf_header_t));
	pad_copy(e_hdr->version, EH_DEFAULT_VERSION, EH_VERSION_LEN); 
	pad_copy(e_hdr->local_patent_id, x_hdr->fname, EH_LOCAL_PATIENT_ID_LEN); 
	pad_copy(e_hdr->local_recording_id, x_hdr->fname, EH_LOCAL_RECORDING_ID_LEN); 
	pad_copy(e_hdr->start_date, EH_DEFAULT_START_DATE, EH_START_DATE_LEN); 
	pad_copy(e_hdr->start_time, EH_DEFAULT_START_TIME, EH_START_TIME_LEN);

	snprintf(tmp, EH_MAX_BUFFER_LEN, "%d", 
		(int)(signal_header_len + sizeof(edf_header_t)));
	pad_copy(e_hdr->header_bytes, tmp, EH_HEADER_BYTES_LEN);

	pad_copy(e_hdr->reserved, "", EH_RESERVED_LEN);


	snprintf(tmp, EH_MAX_BUFFER_LEN, "%d", 
		(int)(floor(x_hdr->duration) / DATA_RECORD_LEN));
	pad_copy(e_hdr->num_data_records, tmp, EH_NUM_DATA_RECORDS_LEN);

	snprintf(tmp, EH_MAX_BUFFER_LEN, "%.02f", DATA_RECORD_LEN);
	pad_copy(e_hdr->data_record_duration, tmp, EH_DATA_RECORD_DURATION_LEN);
	
	snprintf(tmp, EH_MAX_BUFFER_LEN, "%d", x_hdr->num_channels);
	pad_copy(e_hdr->num_signals, tmp, EH_NUM_SIGNALS_LEN);

	return e_hdr;
}

edf_signal_header_t *build_signal_header(xml_header_t *x_hdr)
{
	char *buf;
	char *buf_ptr;
	edf_signal_header_t *es_hdr;

	char tmp[EH_MAX_BUFFER_LEN];
	memset(&tmp, 0, sizeof(tmp));
	
	int i;

	if ((buf = (char*)malloc(EH_S_TOTAL_LEN * x_hdr->num_channels)) == NULL) 
		return NULL;

	if ((es_hdr = (edf_signal_header_t*)malloc(sizeof(edf_signal_header_t))) 
            == NULL) {
		free(buf);
		return NULL;
	}

	es_hdr->buffer = buf;
	es_hdr->length = (EH_S_TOTAL_LEN * x_hdr->num_channels);

	buf_ptr = buf;

	/* lets keep these here for now, even though it's ugly, so we have the
 	 * the option of tweaking these on a per channel basis */

	for (i = 0; i < x_hdr->num_channels; i++) {
		buf_ptr = pad_copy(buf_ptr, 
			get_channel_name(x_hdr->channel_list, i+1), EH_S_LABEL_LEN);
	}

	for (i = 0; i < x_hdr->num_channels; i++) {
		buf_ptr = pad_copy(buf_ptr, 
			EH_S_DEFAULT_TRANSDUCER, EH_S_TRANSDUCER_LEN);	
	}

	for (i = 0; i < x_hdr->num_channels; i++) {
		buf_ptr = pad_copy(buf_ptr, EH_S_DEFAULT_PHYS_DIM, EH_S_PHYS_DIM_LEN);	
	}

	for (i = 0; i < x_hdr->num_channels; i++) {
		buf_ptr = pad_copy(buf_ptr, EH_S_DEFAULT_PHYS_MIN, EH_S_PHYS_MIN_LEN);	
	}

	for (i = 0; i < x_hdr->num_channels; i++) {
		buf_ptr = pad_copy(buf_ptr, EH_S_DEFAULT_PHYS_MAX, EH_S_PHYS_MAX_LEN);	
	}

	for (i = 0; i < x_hdr->num_channels; i++) {
		buf_ptr = pad_copy(buf_ptr, EH_S_DEFAULT_DIG_MIN, EH_S_DIG_MIN_LEN);	
	}

	for (i = 0; i < x_hdr->num_channels; i++) {
		buf_ptr = pad_copy(buf_ptr, EH_S_DEFAULT_DIG_MAX, EH_S_DIG_MAX_LEN);	
	}

	for (i = 0; i < x_hdr->num_channels; i++) {
		buf_ptr = pad_copy(buf_ptr, EH_S_DEFAULT_PREFILTER, 
			EH_S_PREFILTER_LEN);	
	}

	snprintf(tmp, EH_MAX_BUFFER_LEN, "%.f", x_hdr->sfreq * DATA_RECORD_LEN);
	for (i = 0; i < x_hdr->num_channels; i++) {
		buf_ptr = pad_copy(buf_ptr, tmp, EH_S_NUM_SAMPLES_LEN);
	}

	for (i = 0; i < x_hdr->num_channels; i++) {
		buf_ptr = pad_copy(buf_ptr, EH_S_DEFAULT_RESERVED, EH_S_RESERVED_LEN);
	}

	return es_hdr;
}


char *pad_copy(char *dst, char *src, int len)
{
	int srclen = strlen(src);
	int i = 0;

	if (srclen > len)
		return NULL;

	strncpy(dst, src, srclen);
	dst += srclen;
	
	for (i = 0; i < (len - srclen); i++)
		*dst++ = ' ';

	return dst;
}	
