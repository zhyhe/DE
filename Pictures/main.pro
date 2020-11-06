pro main
	set_plot,'PS'
	device,filename='main.ps'
	device,decomposed=1
	openr,10,'./DE.txt'
	srr=(dindgen(20)+1)/10
	arr=dblarr(20)
	brr=dblarr(20)
	readf,10,arr
	readf,10,brr
	close,10
	plot,srr,arr,thick=3,color='000000'XL,background='FFFFFF'XL
	oplot,srr,brr,thick=3,color='0000FF'XL
	XYOUTS,1,250,'-:Uniform distribution',color='000000'XL,charsize=1.5,font=1,charthick=1.5
	XYOUTS,1,180,'-:Nonuniform distribution',color='0000FF'XL,charsize=1.5,font=1,charthick=1.5


end
