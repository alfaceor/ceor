TARGET=main

OBJS=main.o

%.o: %.c
	$(COMPILE) -c $< -o $@

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)
