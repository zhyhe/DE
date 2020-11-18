module constants     !常数作为module
        implicit none
        real(8)              :: N_GW=1000.
        real(8),   parameter :: PI=3.141592653589793_8
        complex(8),parameter :: cj=(0.,1.)
        real(8),   parameter :: c0=299792458._8!m/s
        real(8),   parameter :: Msun=1.989D30!kg
        real(8),   parameter :: GM3=1.547451289752685D-5 !PI*G*Msun/c0^3
        real(8),   parameter :: G=6.67259D-11!Nm^2/kg^2
        real(8),   parameter :: pc=3.0856777581467192D16!m
        real(8),   parameter :: F0=4396.99667411512_8 !c0^3/[6^(3/2)*PI*Msun*G]
        !*************************************************!
        real(8),   parameter :: w0=-1,wa=0.,Omega_k=0.
        real(8),   parameter :: Omega_m=0.2736_8,Omega_de=0.7264_8
        real(8),   parameter :: h0=0.705_8
        !*************************************************!
        real(8),   parameter :: deg2nat=1.745329300562541D-2 !pi/180
        real(8),   parameter :: ET(4)=(/43.54_8,10.42_8,19.48_8,60._8/)!(varphi,lambda,gama,zeta)
end module constants

