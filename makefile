FC=ifort
CFLAGS=-O2

SRCS=$(wildcard *.F90)
OBJS=$(SRCS:.F90=.o)
TARGET=DE

$(TARGET): $(OBJS)
	$(FC) $(CFLAGS) $^ -o $@
%.o: %.F90
	$(FC) $(CFLAGS) -c $^
clean:
	rm $(OBJS) $(TARGET)
