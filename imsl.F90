module IMSL  !@@@@@@@@@@    自定义函数库    @@@@@@@@@@!
        implicit none
contains

        !**********    n阶矩阵A的行列式    **********!
        recursive function det(A,n) result(D)
                implicit none
                integer :: n,k
                real(8) :: A(n,n),B(n-1,n-1)
                real(8) :: D
                if(n>1) then
                        D=0.
                        do k=1,n
                        if(k==1)then
                                B=A(2:n,2:n)
                        elseif(k<n)then
                                B(1:k-1,:)=A(1:k-1,2:n)
                                B(k:n-1,:)=A(k+1:n,2:n)
                        else
                                B=A(1:n-1,2:n)
                        endif
                        D=D+(-1)**(1+k)*A(k,1)*det(B,n-1)
                        enddo
                else
                        D=A(1,1)
                endif
        end function det

        function del(A,n,row,col) result(B)
                implicit none
                integer :: n,row,col
                real(8) :: A(n,n),B(n-1,n-1)
                real(8) :: C(n-1,n)
                if(n<2)then
                        print *,'ERROR'
                endif

                if(row==1)then
                        C=A(2:n,:)
                elseif(row<n)then
                        C(1:row-1,:)=A(1:row-1,:)
                        C(row:n-1,:)=A(row+1:n,:)
                else
                        C=A(1:n-1,:)
                endif
                if(col==1)then
                        B=C(:,2:n)
                elseif(col<n)then
                        B(:,1:col-1)=C(:,1:col-1)
                        B(:,col:n-1)=C(:,col+1:n)
                else
                        B=C(:,1:n-1)
                endif
        end function del

        !**********    n阶矩阵A的逆矩阵    **********!
        function inv(A,n)
                implicit none
                integer :: n,i,j
                real(8) :: A(n,n)
                real(8) :: inv(n,n),b(n-1,n-1)
                real(8) :: De
                De=det(A,n)
                print *,'de=====',De
                do i=1,n
                        do j=1,n
                                b=del(A,n,i,j)
                                inv(j,i)=(-1)**(i+j)*det(b,n-1)/De
                        end do
                end do
        end function inv

        !**********    n阶矩阵A的LU分解    **********!
        subroutine lu(A,n)
                implicit none
                integer,intent(in) :: n
                real(8) :: A(n,n),m
                integer :: i,j,k
                do i=2,n
                        do j=i,n
                                A(j,i-1)=A(j,i-1)/A(i-1,i-1)
                                m=A(j,i-1)
                                do k=i,n
                                A(j,k)=A(j,k)-m*A(i-1,k)
                                enddo
                        enddo
                enddo
        end subroutine lu

        !**********    矩阵A的LU分解法求逆    **********!
        function inverse(A,n)
                implicit none
                integer,intent(in) :: n
                real(8),intent(in) :: A(n,n)
                integer :: i,j,k
                real(8) :: B(n,n),s,de
                real(8) :: inverse(n,n),L(n,n),U(n,n)
                B=A
                call lu(B,n)
                de=1.
                do i=1,n
                        de=de*B(i,i)
                enddo
                do i=1,n
                        L(i,i)=1.
                        U(i,i)=1/B(i,i)
                        do k=i-1,1,-1
                                s=0.
                                do j=k+1,i
                                        s=s+B(k,j)*U(j,i)
                                enddo
                                L(k,i)=0.
                                U(k,i)=-s/B(k,k)
                        enddo
                        do k=i+1,n
                                L(k,i)=0.
                                do j=i,k-1
                                        L(k,i)=L(k,i)-B(k,j)*L(j,i)
                                enddo
                                U(k,i)=0.
                        enddo
                enddo
                inverse=matmul(U,L)
        end function inverse

        !**********    分块迭代法计算一个实对称矩阵A的逆矩阵    **********!
        function inver(A,n)
                implicit none
                integer,intent(in) :: n
                real(8),intent(in) :: A(n,n)
                real(8) :: inver(n,n)
                real(8) :: b(n-1),p
                integer :: t,i,j
                inver=0.
                inver(1,1)=1/A(1,1)
                do t=2,n
                        b=0.
                        p=A(t,t)
                        do i=1,t-1
                                do j=1,t-1
                                        b(i)=b(i)-inver(i,j)*A(j,t)
                                enddo
                                p=p+A(t,i)*b(i)
                        enddo

                        do i=1,t-1
                                do j=1,t-1
                                        inver(i,j)=inver(i,j)+b(i)*b(j)/p
                                enddo
                        enddo
                        do i=1,t
                                inver(i,t)=b(i)/p
                                inver(t,i)=inver(i,t)
                        enddo
                        inver(t,t)=1/p
                enddo
        end function inver

        !**********    计算一个函数的积分    **********!
        real(8) function simpson(func,a,b,width) !simpson积分，精度很高
                real(8),intent(in) :: a,b,width
                real(8),external   :: func
                integer            :: i,n
                real(8)            :: width_
                if(b-a<width)then
                        simpson=0.
                else
                        n=(b-a)/width
                        if(mod(n,2)==0) n=n+1
                        width_=(b-a)/(n-1)
                        simpson=func(a)+func(b)
                        do i=2,n-1
                        if(mod(i,2)==0)then
                                simpson=simpson+4*func(a+(i-1)*width_)
                        else
                                simpson=simpson+2*func(a+(i-1)*width_)
                        end if
                        enddo
                        simpson=simpson*width_/3
                endif
        end function simpson


        !**********    计算一个函数的微分    **********!
        real(8) function Deriv(f,x,h)
                real(8),intent(in) :: x,h
                real(8),external   :: f
                print *,f(x-h),f(x+h),f(x)
                Deriv=(f(x+h)-f(x-h))/(2*h)
        end function Deriv

        !***************  球贝塞尔函数  ****************!
        recursive real(8) function bj(l,x) result(ans)
                integer,intent(in) :: l
                real(8),intent(in) :: x
                if(l>0)then
                        ans=(2*l-1)/x*bj(l-1,x)-bj(l-2,x)
                elseif(l<-1)then
                        ans=(2*l+3)/x*bj(l+1,x)-bj(l+2,x)
                elseif(l==0)then
                        ans=sin(x)/x
                else
                        ans=cos(x)/x
                endif
        end function bj

        !***************  勒让德多项式  ****************!
        recursive real(8) function legP(n,x) result(ans)
                integer,intent(in) :: n
                real(8),intent(in) :: x
                if(n>1)then
                        ans=(2*n-1)*x/n*legP(n-1,x)-(n-1)*legP(n-2,x)/n
                elseif(n==1)then
                        ans=x
                elseif(n==0)then
                        ans=1
                else
                        print *,'ERROR'
                        return
                endif
        end function legP

        !***************  方程组求根  ****************!
        function Root(coe,y,n)
                implicit none
                integer,intent(in) :: n             !n是矩阵阶数
                real(8),intent(in) :: coe(n,n),y(n) !coe是系数矩阵，y是函数值
                integer            :: i
                real(8)            :: ch(n,n),den   !denominator分母
                real(8)            :: Root(n)
                den=det(coe,n)
                do i=1,n
                ch=coe
                ch(:,i)=y
                Root(i)=det(ch,n)/den
                enddo
        end function Root

        
        !***************  最小二乘法拟合  ****************!



               


end module IMSL

