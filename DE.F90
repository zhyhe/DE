Program DE
        use Lamb
        use f_gw
        use f_cmb
        implicit none
        integer :: i,j
        type(para) :: P,P0
        type(source) :: S,S0

        real(8) :: z_p,dw_p,rho2
        real(8) z
        call random_seed()

        z=0.1_8
        print*,Omega_m

        P0=para(w0,wa,Omega_m,Omega_k,h0)
        !print*,'ddddddddddL',dL(z),d_L(z,P0)
        !print*,h0,P0.h0

        !call Fisher_S()
        call Fisher(P0)
        !Fis_U=Fisher(P0)
        !iFis=inverse(Fis_U,5)
        !print *,' '
        !print*,'Fisher Matrix of sum of sources:'
        !print "(5E25.16)",Fis_U
        !print *,' '
        !print "(5F25.16)",iFis
        !print *,' '
        !print "(5F12.5)",sqrt(iFis(1,1)),sqrt(iFis(2,2)),sqrt(iFis(3,3)),sqrt(iFis(4,4)),sqrt(iFis(5,5))

        !print "(2E25.16)",fis2
        !print *,' '
        !print "(2F25.16)",ifis2
        !print *,' '
        !print "(2F12.4)",sqrt(iFis2(1,1)),sqrt(iFis2(2,2))


end
