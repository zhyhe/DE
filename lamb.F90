module Lamb
        use fiducial
        use typedef
        implicit none
        integer,parameter :: dimen=6

        !$$$$$$$$$$$$$$$     计算Fisher Matrix Lambda     $$$$$$$$$$$$$$$!
contains
        !*******************     计算H(f)     *******************!
        subroutine Pat_Func(S,Det,F1,Fx)
                implicit none
                type(source) :: S
                integer,intent(in) :: Det
                real(8) :: F1,Fx,at,bt
                real(8) :: Salpha,Svarphi,Sdelta
                real(8) :: Slambdar,Sgama,Szeta,Spsi
                if( Det/=1 .and. Det/=2 .and. Det/=3 ) then
                        print *,'探测器参数有问题，请检查程序'
                        stop
                end if
                !       计算变换矩阵
                Salpha=S.alpha
                Svarphi=S.varphi
                Sdelta=S.delta
                Spsi=ET(1)*deg2nat
                Slambdar=ET(2)*deg2nat
                Sgama=ET(3)*deg2nat-2*pi/3*(Det-1)
                Szeta=ET(4)*deg2nat

                !       计算文章中的a，b，以及两个响应函数
                at = 1./16*sin(2*Sgama)*(3-cos(2*Spsi))*(3-cos(2*Sdelta))*cos(2*(Salpha-Slambdar))
                at = at-1./4*cos(2*Sgama)*sin(Spsi)*(3-cos(2*Sdelta))*sin(2*(Salpha-Slambdar))
                at = at+1./4*sin(2*Sgama)*sin(2*Spsi)*sin(2*Sdelta)*cos(Salpha-Slambdar)
                at = at-1./2*cos(2*Sgama)*cos(Spsi)*sin(2*Sdelta)*sin(Salpha-Slambdar)
                at = at+3./4*sin(2*Sgama)*cos(Spsi)*cos(Spsi)*cos(Sdelta)*cos(Sdelta)

                bt = cos(2*Sgama)*sin(Spsi)*sin(Sdelta)*cos(2*(Salpha-Slambdar))
                bt = bt+1./4*sin(2*Sgama)*(3-cos(2*Spsi))*sin(Sdelta)*sin(2*(Salpha-Slambdar))
                bt = bt+cos(2*Sgama)*cos(Spsi)*cos(Sdelta)*cos(Salpha-Slambdar)
                bt = bt+1./2*sin(2*Sgama)*sin(2*Spsi)*cos(Sdelta)*sin(Salpha-Slambdar)

                F1 = sin(Szeta)*(at*cos(2*Svarphi)+bt*sin(2*Svarphi))
                Fx = sin(Szeta)*(bt*cos(2*Svarphi)-at*sin(2*Svarphi))
        end subroutine Pat_Func



        REAL(8) FUNCTION phi_(cos_iota,F_plus,F_cross) !函数,计算delta20,结果是对的
                real(8),intent(in) :: cos_iota,F_plus,F_cross
                phi_=atan(-(2*cos_iota*F_cross)/((1+cos_iota**2)*F_plus))
        END FUNCTION phi_
        REAL(8) FUNCTION psi_(f,eta,M)  !code书写正确,没有错误
                integer            :: i
                real(8),intent(in) :: f,eta,M
                real(8) :: psi(0:7)     !数组元素输入正确
                real(8) :: piMf
                real(8) :: lambda,alpha
                piMf=GM3*M*f
                psi(0)=1.
                psi(1)=0.
                psi(2)=3715._8/756+55._8/9*eta
                psi(3)=-16*PI
                psi(4)=15293365._8/508032+27145._8/504*eta+3085._8/72*eta*eta
                psi(5)=PI*(38645._8/756-65._8/9*eta)*(1+3*log(2*6._8**(3./2)*piMf))   !log前面加了一个3,log后面加了一个2
                psi(6)=11583231236531._8/4694215680._8-640._8/3*PI**2-6848._8/21*0.577_8+(-15737765635._8/3048192._8+\
                        2255._8/12*PI**2)*eta+76055._8/1728*eta**2-127825._8/1296*eta**3-6848._8/63*log(64*2*piMf)!piMf前面加了个２
                psi(7)=PI*(77096675._8/254016+378515._8/1512*eta-74045._8/756*eta**2)
                psi_=0.
                do i=0,7
                        psi_=psi_+(3./(256*(2*piMf)**(5._8/3)*eta))*psi(i)*(2*piMf)**(dfloat(i)/3) !(分母上加了(2piMf)^(5/3))(no)
                        !psi_=psi_+(3./(256*eta))*psi(i)*(2*piMf)**(dfloat(i)/3) !分母上加了(2piMf)^(5/3)
                end do
        END FUNCTION psi_

        COMPLEX(8) FUNCTION H_(f,S,Det) !函数,计算h的傅里叶变换H(f),这里的Det表示探测器,只能取1,2,3
                implicit none
                integer,intent(in)      :: Det
                real(8),intent(in)      :: f
                type(source),intent(in) :: S
                real(8)         :: NF_plus,NF_cross
                real(8)         :: A,M
                M=S.M_c/(S.eta**(3._8/5))
                call Pat_Func(S,Det,NF_plus,NF_cross)
                !NF_plus=F_plus(S,Det)
                !NF_cross=F_cross(S,Det)
                A=(G**(5._8/6)*c0**(-3._8/2)/(S.d_L*pc*1D6))*sqrt(NF_plus**2*(1+S.cos_iota**2)**2+\
                  NF_cross**2*4*S.cos_iota**2)*sqrt(5*PI/96)*PI**(-7._8/6)*(S.M_c*Msun)**(5._8/6)
                H_=A*f**(-7._8/6)*exp(cj*(2*PI*f*S.t_c-PI/4-2*S.psi_c+2*psi_(f/2,S.eta,M)-phi_(S.cos_iota,NF_plus,NF_cross)))
                !H_=A*f**(-7._8/6)*exp(cj*(2*PI*f*S.t_c-PI/4-2*S.psi_c-phi_(S.cos_iota,NF_plus,NF_cross)))
        END FUNCTION H_

        !*******************     计算H(f)的偏导     *******************!
        FUNCTION D_H(f,S,Det)
                integer,intent(in)  ::  Det
                real(8),intent(in)  ::  f
                type(source),intent(in):: S
                type(source) S1,S2
                complex(8):: D_H(6)
                integer      i
                real(8)      h
                S1=S
                S2=S
                h=1D-8
                S1.M_c=S.M_c-h !M_c=1.4
                S2.M_c=S.M_c+h
                D_H(1)=(H_(f,S2,Det)-H_(f,S1,Det))*S.M_c/(2*h)

                S1=S
                S2=S
                h=1D-5
                S1.eta=S.eta-h !eta=0.25
                S2.eta=S.eta+h
                D_H(2)=(H_(f,S2,Det)-H_(f,S1,Det))*S.eta/(2*h)

                !S1=S
                !S2=S
                !h=1D-8
                !S1.t_c=S.t_c-h !t_c=0
                !S2.t_c=S.t_c+h
                !D_H=(H_(f,S2,Det)-H_(f,S1,Det))/(2*h)
                D_H(3)=2*PI*f*cj*H_(f,S,Det)

                !S1=S
                !S2=S
                !h=1D-15
                !S1.psi_c=S.psi_c-h !psi_c=0
                !S2.psi_c=S.psi_c+h
                !D_H=(H_(f,S2,Det)-H_(f,S1,Det))/(2*h)
                D_H(4)=-2*cj*H_(f,S,Det)
                
                S1=S
                S2=S
                h=1D-3
                S1.varphi=S.varphi-h
                S2.varphi=S.varphi+h
                D_H(5)=(H_(f,S2,Det)-H_(f,S1,Det))/(2*h)

                !S1=S
                !S2=S
                !h=1D-4
                !S1.d_L=S.d_L-h !d_L*=1000
                !S2.d_L=S.d_L+h
                !D_H=(H_(f,S2,Det)-H_(f,S1,Det))*S.d_L/(2*h)
                D_H(6)=-H_(f,S,Det)
        END FUNCTION D_H

        !***************        Fisher Matrix of a Source        ***************!
        real(8) function S_h(f) !函数，计算ET的PSD,为了计算内积
                implicit none
                real(8),parameter :: f0=200.
                real(8),parameter :: S0=1.449D-52
                real(8),parameter :: p1=-4.05_8,p2=-0.69_8
                real(8),parameter :: a1=185.62_8,a2=232.56_8
                real(8),parameter :: b1=31.18_8,b2=-64.72_8,b3=52.24_8,b4=-42.16_8,b5=10.17_8,b6=11.53_8
                real(8),parameter :: c1=13.58_8,c2=-36.46_8,c3=18.56_8,c4=27.43_8
                real(8),intent(in):: f
                real(8)           :: x
                x=f/f0
                S_h=S0*(x**p1+a1*x**p2+a2*((1+b1*x+b2*x**2+b3*x**3+b4*x**4+b5*x**5+b6*x**6)/(1+c1*x+c2*x**2+c3*x**3+c4*x**4)))
        end Function S_h
        
        real(8) function S_h1(f) !函数，计算ET的PSD,为了计算内积
                implicit none
                real(8),intent(in) :: f
                real(8) :: x,y,t
                open(unit=15,file='ETB.txt')
                read(15,*) x,y
                do while(f>x)
                read(15,*) x,y
                enddo
                S_h1=y**2
                close(15)
        end Function S_h1

        function lambda(S)  !返回Fisher Matrix Lambda
                implicit none
                type(source),intent(in):: S      !自变量是个源
                real(8)      :: lambda(dimen,dimen)      !这个函数的返回值是个数组
                integer      :: i,j,Det,l
                integer,parameter:: n=41        !积分次数,只能取奇数
                real(8)          :: def,ef_lower=0.
                real(8)          :: f_upper,ef_upper,M
                real(8)          :: lamb_simp(6,6,3)!3表示3个探测器
                complex(8)       :: dH(6,n) !6表示fisher matrix的参数数目，n表示积分步数
                real(8)          :: Sh(n),x(n)     !Sh(n)表示每一步的噪声
                M=S.M_c/S.eta**(3._8/5)
                f_upper=F0/M
                ef_upper=log(f_upper)

                def=(ef_upper-ef_lower)/(n-1) !这里是计算步长
                lambda=0.

                if(ef_upper-ef_lower<def)then !这里写得看起来很复杂，是为了算得更快
                        lambda=0.             !先循环f，再循环Det,i,j
                else                          !所以用dH来保存D_H(f,S,Det)
                        do Det=1,3            !lamb_simp(i,j,Det)保存三个探测器分别的fisher matrix
                        do i=1,n
                        x(i)=exp(ef_lower+(i-1)*def)
                        Sh(i)=S_h(x(i))
                        dH(:,i)=D_H(x(i),S,Det)
                        enddo
                        do j=1,dimen
                        do i=j,dimen
                        lamb_simp(i,j,Det)=4*x(1)*(real(dH(i,1))*real(dH(j,1))+aimag(dH(i,1))*aimag(dH(j,1)))/Sh(1)+\
                        4*x(n)*(real(dH(i,n))*real(dH(j,n))+aimag(dH(i,n))*aimag(dH(j,n)))/Sh(n)
                        do l=2,n-1
                        if(mod(l,2)==0)then
                                lamb_simp(i,j,Det)=lamb_simp(i,j,Det)+16*x(l)*(real(dH(i,l))*real(dH(j,l))+aimag(dH(i,l))*aimag(dH(j,l)))/Sh(l)
                        else
                                lamb_simp(i,j,Det)=lamb_simp(i,j,Det)+8*x(l)*(real(dH(i,l))*real(dH(j,l))+aimag(dH(i,l))*aimag(dH(j,l)))/Sh(l)
                        end if
                        enddo
                        lamb_simp(i,j,Det)=lamb_simp(i,j,Det)*def/3
                        enddo
                        enddo
                        lambda=lambda+lamb_simp(:,:,Det)
                        enddo
                        
                        do j=1,dimen
                        do i=1,j-1
                        lambda(i,j)=lambda(j,i)
                        enddo
                        enddo
                endif
        end function lambda

        !**********************************************************************************************!
        subroutine Fisher_S()
                implicit none
                type(source) :: S
                integer      :: i,j,k,l
                integer,parameter :: n=10000
                real(8) :: A(n)
                real(8),parameter :: m1=1.35_8,m2=1.35_8
                real(8)          :: lamda(dimen,dimen),M,M_c,ilamb(dimen,dimen)
                real(8) :: SNR,iota
                M=m1+m2
                S.eta=0.25_8
                M_c=M*S.eta**6D-1
                S.t_c=0.
                S.psi_c=pi/4
                A=0.
                open(unit=10,file='~/workspace/DE.data/New_SNRall_ET2CE7200000_1to1.dat')
                do i=1,2
                read (10,*)
                enddo
                open(unit=11,file='~/workspace/DE.data/dd_L.txt')
                do i=1,10000
                        read(10,*) SNR,S.z,S.d_L,S.alpha,S.delta,S.varphi,iota
                        !print*,S
                        !print*,iota
                        S.alpha=deg2nat*S.alpha
                        S.delta=deg2nat*S.delta
                        S.varphi=deg2nat*S.varphi
                        S.cos_iota=cos(deg2nat*iota)
                        S.d_L=1000*S.d_L
                        S.M_c=(1+S.z)*M_c

                        lamda=lambda(S)
                        ilamb=inverse(lamda,dimen)
                        !print '(6E25.15)',lamda
                        !print *,' '
                        !print '(6E25.15)',matmul(lamda,ilamb)
                        A(i)=ilamb(6,6)+0.0025*S.z**2
                        !print'(F6.4,F21.15)',S.z,sqrt(A(i))
                        write(11,'(F6.4,E25.15,F13.4)')S.z,A(i),rho(S)
                        !write(*,'(F6.4,2F25.15)')S.z,rho(S),Sqrt(A(i))
                        if(mod(i,1000)==0) print*,i
                enddo
                close(11)
                close(10)
        end subroutine Fisher_S

        !*******************        Combined SNR        *******************!
        real(8) function rho(S)  !返回Fisher Matrix Lambda
                implicit none
                type(source),intent(in):: S      !自变量是个源
                integer      :: i,Det,l
                integer,parameter:: n=41        !积分次数,只能取奇数
                real(8)          :: def,ef_lower=0.
                real(8)          :: f_upper,ef_upper,M
                real(8)          :: a,b,x,rho_simp
                M=S.M_c/S.eta**(3._8/5)
                f_upper=F0/M
                ef_upper=log(f_upper)
                def=(ef_upper-ef_lower)/(n-1) !这里是计算步长
                a=exp(ef_lower)
                b=f_upper
                rho=0.

                if(ef_upper-ef_lower<def)then
                        rho=0.
                else
                        do Det=1,3
                        rho_simp=4*a*ABS(H_(a,S,Det))**2/S_h(a)+4*b*ABS(H_(b,S,Det))**2/S_h(b)
                        do l=2,n-1
                        x=exp(ef_lower+(l-1)*def)
                        if(mod(l,2)==0)then
                                rho_simp=rho_simp+16*x*ABS(H_(x,S,Det))**2/S_h(x)
                        else
                                rho_simp=rho_simp+8*x*ABS(H_(x,S,Det))**2/S_h(x)
                        end if
                        enddo
                        rho_simp=rho_simp*def/3
                        rho=rho+rho_simp
                        enddo
                endif
                rho=sqrt(rho)
        end function rho
        !*******************        Combined SNR        *******************!
end module Lamb
