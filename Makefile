# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=source/main.cpp source/BDFReader.cpp source/Display.cpp source/DSPStatistic.cpp source/ECG.cpp source/NoiseGen.cpp source/edflib/edflib.c
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=psyphibio-proc
IDIR =-I./include -I./include/edflib -I./include/gtest
OBJDIR = obj
LDIR = -L./lib
LIBS= -lgtest -lgtest_main
OUTDIR = bin/lin/

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(IDIR) -o $(OUTDIR)$@

.cpp.o:
	$(CC) $(CFLAGS) $(IDIR) $< -o $@

psyphibio-proc_test: $(OBJECTS) source/main_test.o
	$(CC) $(LDFLAGS) $(OBJECTS) source/main_test.o $(IDIR) $(LDIR) $(LIBS) -o $(OUTDIR)$@
	
main_test.o: main.cpp
	$(CC) $(CFLAGS) $(IDIR) source/main_test.cpp

clean:
	rm -rf source/*o hello

