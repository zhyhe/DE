module Lamb
        use fiducial
        use typedef
        implicit none
        integer,parameter :: dimen=7

        !$$$$$$$$$$$$$$$     计算Fisher Matrix Lambda     $$$$$$$$$$$$$$$!
contains
        !*******************     计算H(f)     *******************!
        real(8) Function F_plus(S,Det)    !函数,计算F+,结果是对的
                implicit none
                integer,intent(in)      :: Det   !表示探测器
                type(source),intent(in) :: S
                if( Det/=1 .and. Det/=2 .and. Det/=3 ) then
                        print *,'探测器参数有问题，请检查程序'
                        stop
                end if
                F_plus=sqrt(3._8)/2*(0.5*(1+cos(S.theta)**2)*cos(2*(S.phi+2*PI/3*(Det-1)))*cos(2*S.psi)-cos(S.theta)*sin(2*(S.phi+2*PI/3*(Det-1)))*sin(2*S.psi))
        end Function F_plus
        real(8) Function F_cross(S,Det)   !函数,计算Fx,结果是对的
                implicit none
                integer,intent(in)      :: Det
                type(source),intent(in) :: S
                if( Det/=1 .and. Det/=2 .and. Det/=3 ) then
                        print *,'探测器参数有问题，请检查程序'
                        stop
                end if
                F_cross=sqrt(3._8)/2*(0.5*(1+cos(S.theta)**2)*cos(2*(S.phi+2*PI/3*(Det-1)))*sin(2*S.psi)+cos(S.theta)*sin(2*(S.phi+2*PI/3*(Det-1)))*cos(2*S.psi))
        end Function F_cross

        REAL(8) FUNCTION phi_(cos_iota,F_plus,F_cross) !函数,计算phi20,结果是对的
                real(8),intent(in) :: cos_iota,F_plus,F_cross
                phi_=atan(-(2*cos_iota*F_cross)/((1+cos_iota**2)*F_plus))
        END FUNCTION phi_
        REAL(8) FUNCTION psi_(f,eta,M)  !code书写正确,没有错误
                integer            :: i
                real(8),intent(in) :: f,eta,M
                real(8) :: psi(0:7)     !数组元素输入正确
                real(8) :: piMf
                real(8) :: lambda,theta
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
                NF_plus=F_plus(S,Det)
                !NF_plus=-0.126538817540157_8
                NF_cross=F_cross(S,Det)
                !NF_cross=-0.740002081832980_8
                A=(G**(5._8/6)*c0**(-3._8/2)/(S.d_L*pc*1D6))*sqrt(NF_plus**2*(1+S.cos_iota**2)**2+\
                  NF_cross**2*4*S.cos_iota**2)*sqrt(5*PI/96)*PI**(-7._8/6)*(S.M_c*Msun)**(5._8/6)
                H_=A*f**(-7._8/6)*exp(cj*(2*PI*f*S.t0-PI/4-2*S.PHI0+2*psi_(f/2,S.eta,M)-phi_(S.cos_iota,NF_plus,NF_cross)))
                !H_=A*f**(-7._8/6)*exp(cj*(2*PI*f*S.t0-PI/4-2*S.PHI0-phi_(S.cos_iota,NF_plus,NF_cross)))
        END FUNCTION H_

        !*******************     计算H(f)的偏导     *******************!
         COMPLEX(8) FUNCTION D_H(f,S,Det,j)
                integer,intent(in):: Det,j !j代表对第几个参数进行求导
                real(8),intent(in):: f
                real(8)           :: h
                type(source),intent(in) :: S
                type(source)            :: S1,S2
                integer :: i
                S1=S
                S2=S
                select case(j)
                case(1)
                        h=1D-8
                        S1.M_c=S.M_c-h !M_c=1.4
                        S2.M_c=S.M_c+h
                        D_H=(H_(f,S2,Det)-H_(f,S1,Det))*S.M_c/(2*h)
                case(2)
                        h=1D-5
                        S1.eta=S.eta-h !eta=0.25
                        S2.eta=S.eta+h
                        D_H=(H_(f,S2,Det)-H_(f,S1,Det))*S.eta/(2*h)
                case(3)
                        !h=1D-8
                        !S1.t0=S.t0-h !t0=0
                        !S2.t0=S.t0+h
                        !D_H=(H_(f,S2,Det)-H_(f,S1,Det))/(2*h)
                        D_H=2*PI*f*cj*H_(f,S,Det)
                case(4)
                        !h=1D-15
                        !S1.PHI0=S.PHI0-h !PHI0=0
                        !S2.PHI0=S.PHI0+h
                        !D_H=(H_(f,S2,Det)-H_(f,S1,Det))/(2*h)
                        D_H=-2*cj*H_(f,S,Det)
                case(5)
                        h=1D-1
                        S1.cos_iota=S.cos_iota-h !cos_iota=1
                        S2.cos_iota=S.cos_iota+h
                        D_H=(H_(f,S2,Det)-H_(f,S1,Det))/(2*h)
                case(6)
                        h=1D-3
                        S1.psi=S.psi-h 
                        S2.psi=S.psi+h
                        D_H=(H_(f,S2,Det)-H_(f,S1,Det))/(2*h)
                case(7)
                        !h=1D-4
                        !S1.d_L=S.d_L-h !d_L*=1000
                        !S2.d_L=S.d_L+h
                        !D_H=(H_(f,S2,Det)-H_(f,S1,Det))*S.d_L/(2*h)
                        D_H=-H_(f,S,Det)
                case default
                        print *,'求导的参数有误，请检查程序'
                        stop
                end select
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



        real(8) function lamb_simp(S,low,up,width,i,j,Det) !为下面的lambda函数作准备
                implicit none
                real(8),intent(in)    :: low,up,width
                type(source),intent(in) :: S
                integer,intent(in)    :: i,j,Det !分别是i,j导数及Det
                integer               :: l,n
                real(8)               :: width_,x
                real(8)               :: a,b
                if(up-low<width)then
                        lamb_simp=0.
                else
                        n=(up-low)/width
                        if(mod(n,2)==0) n=n+1
                        width_=(up-low)/(n-1)
                        a=exp(low)
                        b=exp(up)
                        lamb_simp=4*a*(real(D_H(a,S,Det,i))*real(D_H(a,S,Det,j))+aimag(D_H(a,S,Det,i))*aimag(D_H(a,S,Det,j)))/S_h(a)+\
                        4*b*(real(D_H(b,S,Det,i))*real(D_H(b,S,Det,j))+aimag(D_H(b,S,Det,i))*aimag(D_H(b,S,Det,j)))/S_h(b)
                        do l=2,n-1
                        x=exp(a+(l-1)*width_)
                        if(mod(l,2)==0)then
                                lamb_simp=lamb_simp+16*x*(real(D_H(x,S,Det,i))*real(D_H(x,S,Det,j))+aimag(D_H(x,S,Det,i))*aimag(D_H(x,S,Det,j)))/S_h(x)
                        else
                                lamb_simp=lamb_simp+8*x*(real(D_H(x,S,Det,i))*real(D_H(x,S,Det,j))+aimag(D_H(x,S,Det,i))*aimag(D_H(x,S,Det,j)))/S_h(x)
                        end if
                        enddo
                        lamb_simp=lamb_simp*width_/3
                endif
        end function lamb_simp

        real(8) function lamb_norm(S,low,up,width,i,j,Det) !为下面的lambda函数作准备
                implicit none
                real(8),intent(in)    :: low,up,width
                type(source),intent(in) :: S
                integer,intent(in)    :: i,j,Det !分别是i,j导数及Det
                integer               :: l,n
                real(8)               :: x
                real(8)               :: a
                if(up-low<width)then
                        lamb_norm=0.
                else
                        a=exp(low)
                        n=(up-low)/width
                        lamb_norm=0.
                        do l=1,n
                        x=exp(a+l*width)
                                lamb_norm=lamb_norm+4*x*width*(real(D_H(x,S,Det,i))*real(D_H(x,S,Det,j))+aimag(D_H(x,S,Det,i))*aimag(D_H(x,S,Det,j)))/S_h(x)
                        enddo
                endif
        end function lamb_norm

        function lambda(S)  !返回Fisher Matrix Lambda
                implicit none
                type(source),intent(in):: S      !自变量是个源
                real(8)      :: lambda(dimen,dimen)      !这个函数的返回值是个数组
                integer      :: i,j,k
                real(8),parameter:: df=5D-1,ef_lower=0.
                real(8)          :: f_upper,ef_upper,M
                M=S.M_c/S.eta**(3._8/5)
                f_upper=F0/M
                ef_upper=log(f_upper)

                do j=1,dimen
                        do i=j,dimen
                                lambda(i,j)=0.
                                do k=1,3
                                        !lambda(i,j)=lambda(i,j)+lamb_simp(S,ef_lower,ef_upper,df,i,j,k)
                                        lambda(i,j)=lambda(i,j)+lamb_norm(S,ef_lower,ef_upper,df,i,j,k)
                                enddo
                        end do
                end do
                do j=1,dimen
                        do i=1,j-1
                                lambda(i,j)=lambda(j,i)
                        enddo
                enddo
        end function lambda

        !**********************************************************************************************!
        function A_()
                implicit none
                real(8) :: A_(20)
                type(source) :: S
                integer      :: i,j,k,l
                integer,parameter :: n=1
                real(8),parameter :: m1=1.4_8,m2=1.4_8
                real(8)          :: lamda(dimen,dimen),M,M_c,z,ilamb(dimen,dimen)
                real(8) :: l6(6,6),il6(6,6)
                real(8) :: x,iota


                M=m1+m2
                S.eta=0.25_8
                M_c=M*S.eta**6D-1
                S.t0=0.
                S.PHI0=pi/4
                A_=0.
                do i=1,1
                        z=i*1D-1
                        print'(A2,F3.1)','z=',z
                        S.d_L=dL(z)
                        S.M_c=(1+z)*M_c
                        do j=1,n
                                call rand(S.theta,S.phi,180._8)
                                call rand(iota,S.psi,20._8)
                                S.cos_iota=cos(iota)
                                print*,iota
                                print*,S

                                lamda=lambda(S)
                                l6=del(lamda,7,5,5)
                                il6=inverse(l6,6)
                                print '(6D25.15)',l6
                                !ilamb=inverse(lamda,dimen)
                                !print*,sqrt(ilamb(7,7))
                                print *,' '
                                print '(6E25.15)',il6
                                !print *,' '
                                A_(i)=A_(i)+1/(il6(6,6)+0.0025*z**2)
                        enddo
                        A_(i)=1/sqrt(A_(i)/n)
                        print'(A12,F15.7)','delta(lndL)=',A_(i)
                enddo
        end function A_

        !*******************        Combined SNR        *******************!
        real(8) function rho_simp(S,a,b,width,Det) !为下面的rho函数作准备
                implicit none
                real(8),intent(in)    :: a,b,width
                type(source),intent(in) :: S
                integer,intent(in)    :: Det
                integer               :: l,n
                real(8)               :: width_,x
                if(b-a<width)then
                        rho_simp=0.
                else
                        n=(b-a)/width
                        if(mod(n,2)==0) n=n+1
                        width_=(b-a)/(n-1)
                        rho_simp=4*ABS(H_(a,S,Det))**2/S_h(a)+4*ABS(H_(b,S,Det))**2/S_h(b)
                        do l=2,n-1
                        x=a+(l-1)*width_
                        if(mod(l,2)==0)then
                                rho_simp=rho_simp+16*ABS(H_(x,S,Det))**2/S_h(x)
                        else
                                rho_simp=rho_simp+8*ABS(H_(x,S,Det))**2/S_h(x)
                        end if
                        enddo
                        rho_simp=rho_simp*width_/3
                endif
        end function rho_simp
        real(8) function rho(S)
                implicit none
                integer :: i
                type(source),intent(in) :: S
                real(8),parameter       :: df=1D-1,f_lower=1.
                real(8)                 :: f_upper,M
                M=S.M_c/S.eta**(3._8/5)
                f_upper=F0/M
                rho=0.
                do i=1,3
                rho=rho+rho_simp(S,f_lower,f_upper,df,i)
                enddo
                rho=sqrt(rho)
        end function rho
        !*******************        Combined SNR        *******************!
        




end module Lamb
