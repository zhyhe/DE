module typedef      
        implicit none
        !$$$$$$$$$$     数据类型source,存放源的参数     $$$$$$$$$$!
        type source
                real(8) :: M_c,eta,t0,PHI0,cos_iota,psi,d_L,theta,phi
        end type source

        !$$$$$$$$$$     数据类型para,存放宇宙学参数     $$$$$$$$$$!
        type para
                real(8) :: w0,wa,Omega_m,Omega_k,h0
        end type para
end module typedef

