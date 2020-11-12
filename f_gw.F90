module F_GW
        use fiducial
        use typedef
        implicit none
        private A_12,r2,Hub,F_dC,simp,d_C,Fish_simp
contains
        !********    fit of ax+bx^2+cx^3    ********!
        subroutine fit_poly3(X,Y,m,Rot)
                implicit none
                integer,intent(in) :: m
                real(8),intent(in) :: X(m),Y(m)
                integer            :: i,j,k
                real(8)            :: coe(3,3),re(3),Rot(3) !存放系数矩阵与值
                coe=0.
                re=0.
                do i=1,3
                  do j=1,i
                    do k=1,m
                     coe(i,j)=coe(i,j)+X(k)**(i+j)
                    enddo
                  enddo
                enddo
                do i=1,3
                  do j=i+1,3
                    coe(i,j)=coe(j,i) 
                  enddo
                enddo
                do i=1,3
                  do k=1,m
                    re(i)=re(i)+Y(k)*X(k)**i
                  enddo
                enddo
                Rot=Root(coe,re,3)
        end subroutine fit_poly3

        real(8) function A_12(z) !函数A^(-1/2)
                implicit none
                real(8),intent(in) :: z
                integer            :: i
                real(8)            :: Rot(3)
                open(unit=10,file='coe.txt')
                do i=1,3
                        read(10,*) Rot(i)
                enddo
                close(10)
                !A_12=0.1449_8*z-0.0118_8*z**2+0.0012_8*z**3
                A_12=Rot(1)*z+Rot(2)*z**2+Rot(3)*z**3
        end function A_12

        !********    evolution of the source rate    ********!
        real(8) function r2(z)
                implicit none
                real(8),intent(in) :: z
                if(z<=1.)then
                        r2=1+2*z
                elseif(z<5)then
                        r2=(15-3*z)/4
                else
                        r2=0.
                endif
        endfunction r2

        !*************    number distribution    ************!
        real(8) function f(k,z)
                implicit none
                real(8),intent(in) :: z
                integer,intent(in) :: k
                select case(k)
                case(1)  !uniform distribution
                        f=1.089826289902180D-3*4*PI*(dC(z))**2/((1+z)*H(z))
                case(2)  !nonuniform distribution
                        f=4.298333357962829D-4*r2(z)*4*PI*(dC(z))**2/((1+z)*H(z))
                case default
                        print *,'ERROR'
                        stop
                endselect
        end function f

        !******************    光度距离    ******************!
        real(8) function Hub(z,P) !Hubble parameter
                implicit none
                real(8),intent(in) :: z
                type(para)         :: P
                real(8)            :: E,Om_de
                E=(1+z)**(3*(1+P.w0+P.wa))*Exp(-3*P.wa*z/(1+z))
                Om_de=1-P.Omega_m-P.Omega_k
                Hub=P.h0*100*(P.Omega_m*(1+z)**3+P.Omega_k*(1+z)**2+Om_de*E)**0.5!km/s/Mpc
        end function Hub

        real(8) function F_dC(z,P) !共动距离的积分量
                implicit none
                real(8),intent(in)    :: z
                type(para),intent(in) :: P
                F_dC=c0*1D-3/Hub(z,P)
        end function F_dC

       !   特定的simpson积分，专门为了计算上面这个函数的积分   !
        real(8) function simp(P,a,b,width)
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
                real(8)            :: dz=1D-3
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
        real(8) function dd_L(z,P,j)
                implicit none
                integer,intent(in) :: j
                real(8),intent(in) :: z
                type(para),intent(in) :: P
                type(para)         :: P1,P2
                real(8)            :: h
                P1=P
                P2=P

                select case(j)
                case(1)
                        h=1D-3
                        P1.w0=P1.w0-h
                        P2.w0=P2.w0+h
                        dd_L=(d_L(z,P2)-d_L(z,P1))/(2*h*d_L(z,P))
                case(2)
                        h=1D-2
                        P1.wa=P1.wa-h
                        P2.wa=P2.wa+h
                        dd_L=(d_L(z,P2)-d_L(z,P1))/(2*h*d_L(z,P))
                case(3)
                        h=1D-3
                        P1.Omega_m=P1.Omega_m-h
                        P2.Omega_m=P2.Omega_m+h
                        dd_L=(d_L(z,P2)-d_L(z,P1))/(2*h*d_L(z,P))
                case(4)
                        h=1D-3
                        P1.Omega_k=P1.Omega_k-h
                        P2.Omega_k=P2.Omega_k+h
                        dd_L=(d_L(z,P2)-d_L(z,P1))/(2*h*d_L(z,P))
                case(5)
                        h=1D-4
                        P1.h0=P1.h0-h
                        P2.h0=P2.h0+h
                        dd_L=(d_L(z,P2)-d_L(z,P1))/(2*h*d_L(z,P))
                case default
                        print *,'dd_L函数整数参数有误'
                        stop
                endselect
        end function dd_L

        !**********        Fisher Matrix        **********!
        real(8) function Fish_simp(P,a,b,width,i,j,k) !为下面的Fisher函数作准备
                implicit none
                real(8),intent(in)    :: a,b,width
                type(para),intent(in) :: P
                integer,intent(in)    :: i,j,k !分别是i,j导数及fk
                integer               :: l,n
                real(8)               :: width_,x
                n=(b-a)/width
                if(mod(n,2)==0) n=n+1
                width_=(b-a)/(n-1)
                Fish_simp=dd_L(a,P,i)*dd_L(a,P,j)*f(k,a)/(A_12(a))**2+dd_L(b,P,i)*dd_L(b,P,j)*f(k,b)/(A_12(b))**2
                do l=2,n-1
                x=a+(l-1)*width_
                if(mod(l,2)==0)then
                        Fish_simp=Fish_simp+4*dd_L(x,P,i)*dd_L(x,P,j)*f(k,x)/(A_12(x))**2
                else
                        Fish_simp=Fish_simp+2*dd_L(x,P,i)*dd_L(x,P,j)*f(k,x)/(A_12(x))**2
                end if
                enddo
                Fish_simp=Fish_simp*width_/3
        end function Fish_simp

        function Fisher(P,k)  !返回Fisher Matrix F_GW
                type(para),intent(in) :: P
                real(8) :: Fisher(5,5)
                integer,intent(in) :: k
                integer :: i,j
                real(8) :: dz=1D-2
                do j=1,5
                        do i=1,5
                                Fisher(i,j)=Fish_simp(P,1D-2,2._8,dz,i,j,k)
                        enddo
                enddo
        end function Fisher

end module F_GW

        




