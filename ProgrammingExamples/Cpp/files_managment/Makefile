CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS =		readFile.o

LIBS = 

TARGET =	readFile

%.o: %.c
	$(COMPILE) -c $< -o $@

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
