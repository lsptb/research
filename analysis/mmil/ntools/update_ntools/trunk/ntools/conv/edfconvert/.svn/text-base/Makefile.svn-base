#
# Makefile for edfconvert 
#

CC=gcc
CFLAGS=-g -O2 -Wall -Werror `xml2-config --cflags`
LDFLAGS=-g `xml2-config --libs`
OBJS=edfconvert.o hdrxml.o

all: edfconvert 

.c.s:
	$(CC) $(CFLAGS) $(CPPFLAGS) -S $<

edfconvert: edfconvert.o hdrxml.o
	$(CC) $(LDFLAGS) $(OBJS) -o $@

clean:
	$(RM) edfconvert *.o *~
