pro S_h
	set_plot,'PS'
	device,filename='S_h.ps'
	openr,11,'./S_h.txt'
	srr=dindgen(1500)
        brr=dblarr(1500)
        readf,11,brr
        close,11
        plot,srr,brr,thick=3,color='000000'XL,background='FFFFFF'XL,xrange=[1,2000],/xlog,/ylog
end
