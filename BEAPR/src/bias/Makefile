CC=g++
CFLAGS=-c 
LDFLAGS= 
SOURCES= TOK.cpp main.cpp  TextIO.cpp GenomeSegment.cpp Prep.cpp SNP.cpp Normalizer.cpp CrossLinker.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=bias
all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@  

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm -rf *.bak chk *~  *.o $(EXECUTABLE)
install:
	mv $(EXECUTABLE) ../../bin/.

