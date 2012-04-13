TARGET=main

OBJS=main.o

LIBS=-lgsl -lgslcblas -lm

%.o: %.c
	$(COMPILE) -c $< -o $@

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all: $(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)

