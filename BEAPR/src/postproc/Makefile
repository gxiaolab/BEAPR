CC=g++
CFLAGS=-c 
LDFLAGS= 
SOURCES= TOK.cpp main.cpp  TextIO.cpp GenomeSegment.cpp Prep.cpp SNP.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=postproc
all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@  

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm -rf *.bak chk *~  *.o $(EXECUTABLE)
install:
	mv $(EXECUTABLE) ../../bin/.

