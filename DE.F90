Program DE
        use Lamb
        use f_gw
        use f_cmb
        implicit none
        integer :: i,j
        type(para) :: P,P0
        type(source) :: S,S0
 

        real(8) :: Fis_U(5,5),Fis_N(5,5),iFis(5,5)
        real(8) :: fis4(4,4),fis3(3,3),ifis3(3,3),fis2(2,2),ifis2(2,2)
        real(8) :: z_p,dw_p,rho2
        call random_seed()

        P0=para(w0,wa,Omega_m,Omega_k,h0)
        print*,d_L(0.01_8,P0)
        call Fisher_S()
        Fis_U=Fisher(P0)
        iFis=inverse(Fis_U,5)
        print *,' '
        print*,'Fisher Matrix of sum of sources:'
        print "(5E25.16)",Fis_U
        !print *,' '
        !print "(5F25.16)",iFis
        print *,' '
        print "(5F12.5)",sqrt(iFis(1,1)),sqrt(iFis(2,2)),sqrt(iFis(3,3)),sqrt(iFis(4,4)),sqrt(iFis(5,5))

        fis4=del(Fis_U,5,5,5)
        fis3=del(fis4,4,4,4)
        fis2=del(fis3,3,3,3)
        ifis2=inverse(fis2,2)
        print *,' '
        print*,'w0 & wa'
        !print "(2E25.16)",fis2
        !print *,' '
        !print "(2F25.16)",ifis2
        !print *,' '
        !print "(2F12.4)",sqrt(iFis2(1,1)),sqrt(iFis2(2,2))
        print "(2F12.4)",sqrt(iFis2(1,1)),sqrt(iFis2(2,2))
        print *,' '
        print *,'FoM'
        print*,1/sqrt(Det(ifis2,2))
        print *,' '
        z_p=-1/(1+iFis2(2,2)/iFis2(1,2))
        rho2=iFis(1,2)**2/(iFis(1,1)*iFis(2,2))
        dw_p=sqrt(iFis2(1,1)*(1-rho2))
        print *,' '
        print*,'z_p & dw_p'
        print "(2F12.4)",z_p,dw_p
        !print *,Kappa(z)
        fis4=del(Fis_U,5,1,1)
        fis3=del(fis4,4,1,1)
        ifis3=inverse(fis3,3)
        print *,' '
        print*,'Omega_m, Omega_k & h0'
        !print "(3E25.16)",fis3
        !print *,' '
        !print "(3F25.16)",ifis3
        !print *,' '
        print "(3F12.5)",sqrt(iFis3(1,1)),sqrt(iFis3(2,2)),sqrt(iFis3(3,3))


end
