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
        real(8) :: z_p,dw_p
        call random_seed()
        P0=para(w0,wa,Omega_m,Omega_k,h0)
        A=A_()
        open(unit=13,file='A.txt')
        do i=1,20
        write(13,*),A(i)
        enddo
        close(13)


        !*******************************
        !print *,d_L(z,P0)
        !print *,dL(z)

        !open(unit=11,file='appen/S_h.txt')
        !do i=1,1500
        !        write(11,*) sqrt(S_h(dfloat(i)))
        !enddoa

        !********************************
        do i=1,20
        x(i)=1D-1*i
        enddo
        call fit_poly3(x,A,20,Rot)

        open(unit=14,file='coe.txt')
        do i=1,3
        write(14,*),Rot(i)
        enddo
        close(14)



       









        Fis_U=Fisher(P0,1)
        Fis_N=Fisher(P0,2)
        iFis=inverse(Fis_U,5)
        print *,' '
        print*,'Uniform Distribution'
        print "(5E25.16)",Fis_U
        !print *,' '
        !print "(5F25.16)",iFis
        print *,' '
        print "(5F12.5)",sqrt(iFis(1,1)),sqrt(iFis(2,2)),sqrt(iFis(3,3)),sqrt(iFis(4,4)),sqrt(iFis(5,5))
       
        iFis=inverse(Fis_N,5)
        print *,' '
        print*,'Nonuniform Distribution'
        print "(5E25.16)",Fis_N
        !print *,' '
        !print "(5F25.16)",iFis
        print *,' '
        print "(5F12.5)",sqrt(iFis(1,1)),sqrt(iFis(2,2)),sqrt(iFis(3,3)),sqrt(iFis(4,4)),sqrt(iFis(5,5))
        
        fis4=del(Fis_U,5,5,5)
        fis3=del(fis4,4,4,4)
        fis2=del(fis3,3,3,3)
        ifis2=inverse(fis2,2)
        print *,' '
        print*,'Uniform Distribution of w0 & wa'
        !print "(2E25.16)",fis2
        !print *,' '
        !print "(2F25.16)",ifis2
        !print *,' '
        print "(2F12.4)",sqrt(iFis2(1,1)),sqrt(iFis2(2,2))
        !z_p=-1/(1+sqrt(iFis2(2,2))/(iFis2(1,2)*sqrt(iFis2(1,1))))
        !dw_p=sqrt(iFis2(1,1)*(1-iFis2(1,2)**2))
        !print*,'z_p & dw_p'
        !print "(2F12.4)",z_p,dw_p
        !print *,Kappa(z)
        fis4=del(Fis_U,5,1,1)
        fis3=del(fis4,4,1,1)
        ifis3=inverse(fis3,3)
        print *,' '
        print*,'Uniform Distribution of Omega_m, Omega_k & h0'
        !print "(3E25.16)",fis3
        !print *,' '
        !print "(3F25.16)",ifis3
        !print *,' '
        print "(3F12.5)",sqrt(iFis3(1,1)),sqrt(iFis3(2,2)),sqrt(iFis3(3,3))        
        
        fis4=del(Fis_N,5,5,5)
        fis3=del(fis4,4,4,4)
        fis2=del(fis3,3,3,3)
        ifis2=inverse(fis2,2)
        print *,' '
        print*,'Nonuniform Distribution of w0 & wa'
        !print "(2E25.16)",fis2
        !print *,' '
        !print "(2F25.16)",ifis2
        !print *,' '
        print "(2F12.4)",sqrt(iFis2(1,1)),sqrt(iFis2(2,2))        
        !z_p=-1/(1+sqrt(iFis2(2,2))/(iFis2(1,2)*sqrt(iFis2(1,1))))
        !dw_p=sqrt(iFis2(1,1)*(1-iFis2(1,2)**2))
        !print*,'z_p & dw_p'
        !!print "(2F12.4)",z_p,dw_p

end
