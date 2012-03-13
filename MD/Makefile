CXXFLAGS =	-g -Wall -fmessage-length=0

OBJS =		MD.o Monomer.o Conformation.o util_functions.o

LIBS = -lgsl -lgslcblas -lm

TARGET =	MD

%.o: %.c
	$(COMPILE) -c $< -o $@

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
