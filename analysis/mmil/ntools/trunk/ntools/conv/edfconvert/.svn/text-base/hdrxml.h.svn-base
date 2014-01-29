#ifndef __HDRXML_H__
#define __HDRXML_H__

#include <stdbool.h>

#include <libxml/tree.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#define MAX_NAME_LEN 512

struct channel {
	int num;
	char name[MAX_NAME_LEN];
	struct channel *next;
};
	
typedef struct channel channel_t;

enum duration_unit_type {
	DUR_UNIT_UNKNOWN,
	DUR_UNIT_SECONDS
};

enum data_type { 
	DATA_UNKNOWN,
	DATA_INT16
};

struct xml_header {
	char fname[MAX_NAME_LEN];
	int sfreq;
	int num_channels;
	int num_samples;
	enum data_type data_type;
	float duration;
	enum duration_unit_type duration_units;
	channel_t *channel_list;
};


typedef struct xml_header xml_header_t;

void print_hdr(xml_header_t *hdr);
xml_header_t *parse_hdrxml(char *filename);
bool parse_fname(xml_header_t *hdr, xmlNode *node);
bool parse_node(xml_header_t *hdr, xmlDoc *doc, xmlNode *node);
int parse_int(xmlDoc *doc, xmlNode *node);
channel_t *parse_channel_list(xmlDoc *doc, xmlNode *node);
bool parse_string(xmlDoc *doc, xmlNode *node, char *dst, int len);
channel_t *make_channel(int num, char *name);
channel_t *append(channel_t *list, channel_t *item);
char *get_channel_name(channel_t *list, int num);
bool parse_stringprop(char *dst, xmlNode *node, char *name, int len);
bool parse_duration(xml_header_t *hdr, xmlDoc *doc, xmlNode *node);
bool parse_data_type(xml_header_t *hdr, xmlDoc *doc, xmlNode *node);
int parse_float(xmlDoc *doc, xmlNode *node);
#endif
