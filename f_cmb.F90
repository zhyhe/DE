module F_CMB
        use fiducial
        use typedef
        implicit none
        real(8),parameter :: tau0=10

contains
        real(8) function Kappa_dot(tau)
                real(8),intent(in) :: tau
                real(8),parameter  :: ne=1,xe=1,sigmaT=1
                Kappa_dot=ne*xe*sigmaT*tau
        end function Kappa_dot

        real(8) function Kappa(tau)
                real(8),intent(in) :: tau
                real(8)            :: d=1D-2
                Kappa=simpson(Kappa_dot,tau,tau0,d)
        end function Kappa

        real(8) function g_(tau)
                real(8),intent(in) :: tau
                g_=Kappa_dot(tau)*Exp(-Kappa(tau))
        end function g_

        subroutine Fisher_CMB()
                implicit none
                integer,parameter :: n=5
                integer :: i
                real(8) :: Fis(n,n),iFis(n,n)
                real(8),parameter :: Fis_cmb(n,n)=\
                (/41430.3_8,  11508.5_8,  287229._8, -678690._8,  339923._8,\
                  11508.5_8,  3202.04_8,  79741.5_8, -190373._8,  94349.8_8,\
                  287229._8,  79741.5_8, 2198540._8,-4656630._8, 2518130._8,\
                 -678690._8, -190373._8,-4656630._8,13691200._8,-5480910._8,\
                  339923._8,  94349.8_8, 2518130._8,-5480910._8, 2915870._8/)


                open(unit=26,file='/home/zhyhe/workspace/DE.data/Fis.txt')
                do i=1,5
                read(26,'(5F25.15)'),Fis(:,i)
                enddo
                close(26)

                Fis=Fis+Fis_cmb
                iFis=inverse(Fis,n)
                print "(5F14.9)",sqrt(iFis(1,1)),sqrt(iFis(2,2)),sqrt(iFis(3,3)),sqrt(iFis(4,4)),sqrt(iFis(5,5))
        end subroutine Fisher_CMB




end module F_CMB





