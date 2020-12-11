module F_GW
        use fiducial
        use typedef
        implicit none
        !private Hub,F_dC,simp,d_C
contains
        !******************    光度距离    ******************!
        real(8) function Hub(z,P) !Hubble parameter
                implicit none
                real(8),intent(in) :: z
                type(para)         :: P
                real(8)            :: E,Om_de
                E=(1+z)**(3*(1+P.w0+P.wa))*Exp(-3*P.wa*z/(1+z))
                Om_de=1-P.Omega_m-P.Omega_k
                Hub=P.h0*1D2*(P.Omega_m*(1+z)**3+P.Omega_k*(1+z)**2+Om_de*E)**0.5!km/s/Mpc
        end function Hub

        real(8) function F_dC(z,P) !共动距离的积分量
                implicit none
                real(8),intent(in)    :: z
                type(para),intent(in) :: P
                F_dC=c0*1D-3/Hub(z,P)
        end function F_dC

       !   特定的simpson积分，专门为了计算上面这个函数的积分   !
        real(8) function simp(P,a,b,width)
                implicit none
                real(8),intent(in)    :: a,b,width
                type(para),intent(in) :: P
                integer               :: i,n
                real(8)               :: width_
                if(b-a<width)then
                        simp=0.
                else
                        n=(b-a)/width
                        if(mod(n,2)==0) n=n+1
                        width_=(b-a)/(n-1)
                        simp=F_dC(a,P)+F_dC(b,P)
                        do i=2,n-1
                        if(mod(i,2)==0)then
                                simp=simp+4*F_dC(a+(i-1)*width_,P)
                        else
                                simp=simp+2*F_dC(a+(i-1)*width_,P)
                        end if
                        enddo
                        simp=simp*width_/3
                endif
        end function simp

        real(8) function d_C(z,P)!共动距离(Mpc)
                implicit none
                real(8),intent(in) :: z
                type(para)         :: P
                real(8)            :: dz=1D-2
                d_C=simp(P,0._8,z,dz)
        end function d_C

        real(8) function d_L(z,P) !光度距离(Mpc)
                implicit none
                real(8),intent(in)    :: z
                type(para),intent(in) :: P
                real(8)               :: e=1D-8
                real(8)               :: msk
                if(P.Omega_k<-e)then
                        msk=P.h0*1D5*sqrt(-P.Omega_k)/c0
                        d_L=(1+z)/msk*sin(msk*d_C(z,P))
                elseif(P.Omega_k>e)then
                        msk=P.h0*1D5*sqrt(P.Omega_k)/c0
                        d_L=(1+z)/msk*sinh(msk*d_C(z,P))
                else
                        d_L=(1+z)*d_C(z,P)
                endif
        end function d_L


        !***********    光度距离d_L的偏微分dd_L    ***********!
        function dd_L(z,P)
                implicit none
                real(8),intent(in)    :: z
                type(para),intent(in) :: P
                type(para) P1,P2
                real(8) :: h
                real(8) :: dd_L(6)
                
                h=1D-7
                P1=P
                P2=P
                P1.w0=P1.w0-h
                P2.w0=P2.w0+h
                dd_L(1)=(d_L(z,P2)-d_L(z,P1))/(2*h*d_L(z,P))
                
                h=1D-7
                P1=P
                P2=P
                P1.wa=P1.wa-h
                P2.wa=P2.wa+h
                dd_L(2)=(d_L(z,P2)-d_L(z,P1))/(2*h*d_L(z,P))
                
                h=1D-3
                P1=P
                P2=P
                P1.Omega_m=P1.Omega_m-h
                P2.Omega_m=P2.Omega_m+h
                dd_L(3)=(d_L(z,P2)-d_L(z,P1))/(2*h*d_L(z,P))
                
                h=1D-3
                P1=P
                P2=P
                P1.Omega_k=P1.Omega_k-h
                P2.Omega_k=P2.Omega_k+h
                dd_L(4)=(d_L(z,P2)-d_L(z,P1))/(2*h*d_L(z,P))
                
                h=1D-4
                P1=P
                P2=P
                P1.h0=P1.h0-h
                P2.h0=P2.h0+h
                dd_L(5)=(d_L(z,P2)-d_L(z,P1))/(2*h*d_L(z,P))
                
                h=1D-7
                dd_L(6)=(d_L(z+h,P)-d_L(z-h,P))/(2*h*d_L(z,P))
        end function dd_L

        !**********        Fisher Matrix        **********!
        subroutine Fisher(P)
                implicit none
                integer,parameter :: m=5
                type(para),intent(in) :: P
                real(8) :: Fish(m,m),iFish(m,m)
                real(8) :: fis(2,2),ifis(2,2)
                real(8) :: dLs(6)
                real(8) :: B,C
                real(8) :: a,z,dL,dz,ddL
                integer :: i,j,k,l,n=10000
                Fish=0
                open(unit=30,file='/home/zhyhe/workspace/DE.data/New_SNRall_ET2CE10000.dat')
                read(30,*)
                do l=1,19

                read(30,*) B,C
                Fish=0

                do k=1,n
                read(30,*) a,z,dL,ddL,dz
        !if (a>10) then
                dLs=dd_L(z,P)
                do j=1,m
                        do i=1,m
                                Fish(i,j)=Fish(i,j)+dLs(i)*dLs(j)/(ddL**2+(dz*dLs(6))**2)
                        enddo
                enddo
        !endif
                enddo
                Fish=210000*Fish/n
                iFish=inverse(Fish,m)
                fis=dob(dob(dob(Fish,5,5),4,4),3,3)
                ifis=inverse(fis,2)
                print*,''
                print "(2F7.2,2F14.9,A8,F14.9)",B,C,sqrt(iFis(1,1)),sqrt(iFis(2,2)),'FoM=',1/sqrt(Det(ifis,2))
                print "(5F14.9)",sqrt(iFish(1,1)),sqrt(iFish(2,2)),sqrt(iFish(3,3)),sqrt(iFish(4,4)),sqrt(iFish(5,5))

                enddo
                close(30)
        end subroutine Fisher

end module F_GW

        




