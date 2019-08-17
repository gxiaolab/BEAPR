CC=g++
CFLAGS=-c 
LDFLAGS= 
SOURCES= TOK.cpp main.cpp  TextIO.cpp Motif.cpp 
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

