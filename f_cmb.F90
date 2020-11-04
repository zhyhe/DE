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



end module F_CMB





