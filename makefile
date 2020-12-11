FC=ifort
CFLAGS=-O2

SRCS=$(wildcard *.F90)
OBJS=$(SRCS:.F90=.o)
TARGET=DE
MODa=f_cmb.mod f_gw.mod lamb.mod
MODb=fiducial.mod typedef.mod
MODc=constants.mod imsl.mod

$(TARGET): $(OBJS)
	$(FC) $(CFLAGS) $^ -o $@

DE.o: DE.F90 $(MODa)
	$(FC) $(CFLAGS) -c DE.F90
f_cmb.o f_cmb.mod: f_cmb.F90 $(MODb)
	$(FC) $(CFLAGS) -c f_cmb.F90
f_gw.o f_gw.mod: f_gw.F90 $(MODb)
	$(FC) $(CFLAGS) -c f_gw.F90
lamb.o lamb.mod: lamb.F90 $(MODb)
	$(FC) $(CFLAGS) -c lamb.F90
fiducial.o fiducial.mod: fiducial.F90 $(MODc)
	$(FC) $(CFLAGS) -c fiducial.F90
constants.o constants.mod: constants.F90
	$(FC) $(CFLAGS) -c constants.F90
typedef.o typedef.mod: typedef.F90
	$(FC) $(CFLAGS) -c typedef.F90
imsl.o imsl.mod: imsl.F90
	$(FC) $(CFLAGS) -c imsl.F90
clean:
	rm $(OBJS) $(TARGET) $(MODa) $(MODb) $(MODc)
