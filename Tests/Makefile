CXXFLAGS =	-g -p -Wall -fmessage-length=0

OBJS =		Tests.o forces.o pdb_utils.o

LIBS =		-lgsl -lgslcblas -lm

TARGET =	Tests

%.o: %.c
	$(COMPILE) -c $< -o $@

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
