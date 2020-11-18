module typedef      
        implicit none
        !$$$$$$$$$$     数据类型source,存放源的参数     $$$$$$$$$$!
        type source
                real(8) :: alpha,delta,varphi,cos_iota,M_c,eta,t_c,psi_c,d_L,z
        end type source

        !$$$$$$$$$$     数据类型para,存放宇宙学参数     $$$$$$$$$$!
        type para
                real(8) :: w0,wa,Omega_m,Omega_k,h0
        end type para
end module typedef

