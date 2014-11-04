# Makefile for faerie
# Written by yiduan (2013-07-10)

AllObjects = EntityExtraction.o HeapBinary.o Util.o array.o

CC = g++
CPPFLAGS = -Wall -O3
LDFLAGS = -lrt

faerie :$(AllObjects)
	$(CC) -o faerie $(LDFLAGS) $(AllObjects)
	
HeapBinary.o : HeapBinary.cpp HeapBinary.h Util.o
	$(CC) -c -o $@ $(CPPFLAGS) $<

Util.o : Util.cpp Util.h array.o
	$(CC) -c -o $@ $(CPPFLAGS) $<

array.o : array.cc array.h
	$(CC) -c -o $@ $(CPPFLAGS) $<

.PHONY: clean
clean :
	rm *.o
