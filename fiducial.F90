module fiducial
        use constants
        use IMSL
        implicit none
        private F_dC   !这个函数是为了计算dC来的，不外传

                
        !################# 计算过程中一些最基本的函数 #################!
        contains
        
                !****************** Hubble parameter ******************!
                real(8) function H(z)
                        implicit none
                        real(8),intent(in) :: z
                        H=h0*1D2*(Omega_m*(1+z)**3+Omega_de)**0.5!km/s/Mpc
                end function H                
        
                !******************     共动距离     ******************!
                real(8) function F_dC(z)
                        implicit none
                        real(8),intent(in) :: z
                        F_dC=c0*1D-3/H(z)
                end function F_dC
                real(8) function dC(z)!Mpc
                        implicit none
                        real(8),intent(in) :: z
                        real(8)            :: dz=1D-3
                        dC=simpson(F_dC,0._8,z,dz)
                end function dC

                !******************     光度距离     ******************!
                real(8) function dL(z)!Mpc
                        implicit none
                        real(8),intent(in) :: z
                        dL=(1+z)*dC(z)
                end function dL

                !********    产生lbound到ubound之间的随机数    ********!
                real(8) function random(lbound,ubound)
                        implicit none
                        real(8),intent(in) :: lbound,ubound
                        real(8)            :: len,t
                        len=ubound-lbound
                        call random_number(t)
                        random=lbound+len*t
                end function random

                subroutine rand(theta,phi,angle_lim) !产生球面上的随机数,angle_lim表示theta的最大可取值
                        implicit none
                        real(8) :: x,y,z,r
                        real(8) :: theta,phi,angle_lim
                        10 call random_number(x)
                        call random_number(y)
                        call random_number(z)
                        x=x-5D-1
                        y=y-5D-1
                        z=z-5D-1
                        r=sqrt(x**2+y**2+z**2)
                        if(r .gt. 5D-1 .or. r .le. 1D-1 .or. z .le. cos(angle_lim*pi/180)/2) goto 10
                        theta=-asin(z/r)+pi/2
                        call random_number(phi)
                        phi=2*pi*phi
                end subroutine rand
                
end module fiducial




