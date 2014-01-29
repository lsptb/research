#include <stdbool.h>
#include <string.h>

#include <libxml/tree.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#include "hdrxml.h"

static void
print_element_names(xmlNode * a_node)
{
	xmlNode *cur_node = NULL;

	for (cur_node = a_node; cur_node; cur_node = cur_node->next) {
		if (cur_node->type == XML_ELEMENT_NODE) {
//			printf("node type: Element, name: %s\n", cur_node->name);
		}

	        print_element_names(cur_node->children);
    	}
}


xml_header_t *parse_hdrxml(char *filename)
{
	xmlDoc *doc = NULL;
	xmlNode *root_element = NULL;

	xml_header_t *hdr = (xml_header_t*)malloc(sizeof(xml_header_t));
	
	LIBXML_TEST_VERSION

	if ((doc = xmlReadFile(filename, NULL, 0)) == NULL) {
		fprintf(stderr, 
			"Unable to parse XML header file: %s\n", filename);
		return NULL;
	}

	root_element = xmlDocGetRootElement(doc);
	
	if (root_element == NULL) {
		fprintf(stderr, "XML document empty\n");
		xmlFreeDoc(doc);
		return NULL;
	}

	if (xmlStrcmp(root_element->name, (const xmlChar *)"hdr") != 0) {
		fprintf(stderr, "root node is not <hdr>\n");
		xmlFreeDoc(doc);
		return NULL;
	}

	
	if (!parse_stringprop(hdr->fname, root_element, "fname", MAX_NAME_LEN)) {
		fprintf(stderr, "unable to parse fname\n");
		return NULL;
	}

	xmlNode *cur_node = NULL;

	for (cur_node = root_element->xmlChildrenNode; cur_node; 
		cur_node = cur_node->next) {
		if (cur_node->type == XML_ELEMENT_NODE) {
			parse_node(hdr, doc, cur_node);
		}

    	}

	print_hdr(hdr);		

	xmlFreeDoc(doc);
	xmlCleanupParser();

	return hdr;
}

void print_hdr(xml_header_t *hdr)
{
	fprintf(stdout, "fname:\t\t%s\n",hdr->fname);
	fprintf(stdout, "sfreq:\t\t%d\n",hdr->sfreq);
	fprintf(stdout, "num_channels:\t%d\n",hdr->num_channels);
	fprintf(stdout, "num_samples:\t%d\n",hdr->num_samples);
	fprintf(stdout, "data_type:\t%d\n",hdr->data_type);
	fprintf(stdout, "duration:\t%f\n",hdr->duration);
	fprintf(stdout, "dur_units:\t%d\n",hdr->duration_units);

//	channel_t *cur;
//	for (cur = hdr->channel_list; cur; cur = cur->next) {
//		printf("channel %d: name=%s\n",cur->num, cur->name);
//	}
}

bool parse_stringprop(char *dst, xmlNode *node, char *name, int len)
{
	xmlChar *prop;

	if ((prop = xmlGetProp(node, (xmlChar *)name)) == NULL)
		return false;
		
	strncpy(dst, (char*)prop, len);
	xmlFree(prop);

	return true;
}

int parse_chnum(xmlNode *node)
{
	xmlChar *fname;
	int num;

	if ((fname = xmlGetProp(node, (const xmlChar *)"num")) == NULL)
		return false;
		
	num = atoi((char*)fname);
	xmlFree(fname);

	return num;
}

bool parse_node(xml_header_t *hdr, xmlDoc *doc, xmlNode *node)
{
	if (xmlStrcmp(node->name, (const xmlChar *) "sfreq") == 0)
		hdr->sfreq = parse_int(doc, node);
	else if (xmlStrcmp(node->name, (const xmlChar *)"num_channels") == 0)
		hdr->num_channels = parse_int(doc, node);
	else if (xmlStrcmp(node->name, (const xmlChar *)"num_samples") == 0)
		hdr->num_samples = parse_int(doc, node);
	else if (xmlStrcmp(node->name, (const xmlChar *)"channel_names") == 0)
		hdr->channel_list = parse_channel_list(doc, node->xmlChildrenNode);
	else if (xmlStrcmp(node->name, (const xmlChar *)"duration") == 0)
		parse_duration(hdr, doc, node);
	else if (xmlStrcmp(node->name, (const xmlChar *)"data_type") == 0)
		parse_data_type(hdr, doc, node);

	return true;
}

bool parse_duration(xml_header_t *hdr, xmlDoc *doc, xmlNode *node)
{
	char duration_unit[MAX_NAME_LEN];

	if (!parse_stringprop(duration_unit, node, "units", MAX_NAME_LEN)) {
		fprintf(stderr, "unable to parse duration units\n");
		return false;
	}

	if (strncmp(duration_unit, "seconds", MAX_NAME_LEN) != 0) {
		fprintf(stderr, "unknown duration units\n");
		return false;
	}

	hdr->duration_units = DUR_UNIT_SECONDS;

	hdr->duration = parse_float(doc, node);

	return true;
}
	

bool parse_data_type(xml_header_t *hdr, xmlDoc *doc, xmlNode *node)
{
	char data_type[MAX_NAME_LEN];

	parse_string(doc, node, (char *)&data_type, MAX_NAME_LEN);

	if (strncmp(data_type, "int16", MAX_NAME_LEN) != 0) {
		fprintf(stderr, "unknown data type\n");
		return false;
	}

	hdr->data_type = DATA_INT16;

	return true;
}


channel_t *parse_channel_list(xmlDoc *doc, xmlNode *node)
{
	xmlNode *cur_node;
	channel_t *channel_list_head = NULL;
	int num;
	char name[MAX_NAME_LEN];
	
	for (cur_node = node; cur_node; cur_node = cur_node->next) {
		if (xmlStrcmp(cur_node->name, (const xmlChar *) "channel_name") == 0) {
			parse_string(doc, cur_node, (char *)&name, MAX_NAME_LEN);
			num = parse_chnum(cur_node);	 
			channel_list_head = append(channel_list_head,
				make_channel(num, name));
		}
								
	}

	return channel_list_head;
}

char *get_channel_name(channel_t *list, int num)
{
	channel_t *c;
	
	for (c = list; c; c = c->next) { 
		if (c->num == num)
			return c->name;	
	}

	return NULL;
}

	
channel_t *append(channel_t *list, channel_t *item)
{
	if (list == NULL)
		return item;

	channel_t *c;
	
	for (c = list; c->next; c = c->next) { }
	c->next = item;

	return list;
}

channel_t *make_channel(int num, char *name)
{
	channel_t *c;
	
	if ((c = (channel_t*)malloc(sizeof(channel_t))) == NULL)
		return NULL;

	memset(c, 0, sizeof(channel_t));
	
	c->num = num;
	strncpy(c->name, name, MAX_NAME_LEN);
	return c;
}


int parse_int(xmlDoc *doc, xmlNode *node)
{
	xmlChar *text;
	int val;
	text = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
	val = atoi((char*)text);
	xmlFree(text);

	return val;
}

int parse_float(xmlDoc *doc, xmlNode *node)
{
	xmlChar *text;
	float val;
	text = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
	val = atof((char*)text);
	xmlFree(text);

	return val;
}


bool parse_string(xmlDoc *doc, xmlNode *node, char *dst, int len)
{
	xmlChar *text;
	
	text = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
	strncpy(dst, (char *)text, len);
	xmlFree(text);

	return true;
}

