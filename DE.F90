Program DE
        use Lamb
        use f_gw
        implicit none
        integer :: i,j
        integer,parameter :: n=7
        real(8) :: z=2D-1
        type(para) :: P,P0
        type(source) :: S,S0
 
        real(8) :: dee

        real(8) :: x(20),A(20),Rot(3)
        real(8) :: Fis_U(5,5),Fis_N(5,5),iFis(5,5)
        real(8) :: fis4(4,4),fis3(3,3),ifis3(3,3),fis2(2,2),ifis2(2,2)
        real(8) :: z_p,dw_p,rho2
        call random_seed()
        P0=para(w0,wa,Omega_m,Omega_k,h0)
        A=A_()

end
