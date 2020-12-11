ifort -c imsl.F90
ifort -c fiducial.F90
ifort -c f_gw.F90
ifort -c lamb.F90
ifort -c f_cmb.F90
ifort typedef.o imsl.o fiducial.o f_gw.o lamb.o f_cmb.o DE.F90 -o DE

