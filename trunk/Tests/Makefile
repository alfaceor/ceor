CXXFLAGS = -g -pg -Wall -fmessage-length=0

OBJS =		Tests.o forces.o pdb_utils.o

LIBS =		-lgsl -lgslcblas -lm

TARGET =	Tests

%.o: %.c
	$(COMPILE) -pg -c $< -o $@

$(TARGET):	$(OBJS)
	$(CXX) -pg -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
