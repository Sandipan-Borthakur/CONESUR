mode='blue'
;mode='red'
if mode eq 'red' then begin
  dir='Reduced_red'
  cd,dir
  f=file_search('*[0-9].ech') 
  for ifile=0,n_elements(f)-1 do begin
    rdech,e,f[ifile],/NOCONT
    e.sig=e.sig*sqrt(e.cont)
    restore,'harps_red.ord_default.sav'
    blzcoef=blzcoef>1.
    ee=make_cont(e,blzcoef,colrange=col_range,param=[10,8.d5,1.d7,1.],/SCALING,/PLOT, $
                 yr=[0,10],wrange=[[5577,5577.2],[5885,5902],[6299.9,6300.2]])
    nmout=strmid(f[ifile],0,strpos(f[ifile],'.ech'))+'c.ech'
    rdech,e,f[ifile],/RAW
    wdech,nmout,ee.head,ee.spec,sig=ee.sig,cont=ee.cont,wave=e.wave $
         ,orders=e.orders,/OVERWRITE
  endfor
  cd,'..'
  
  
endif else if mode eq 'blue' then begin
  dir='Reduced_blue'
  cd,dir
  f=file_search('*[0-9].ech')
  for ifile=0,n_elements(f)-1 do begin
    rdech,e,f[ifile],/NOCONT
    e.sig=e.sig*sqrt(e.cont)
    restore,'harps_blue.ord_default.sav'
    
    ee=make_cont(e,blzcoef,colrange=col_range,param=[10,1.d5,1.d8,1.],/SCALING,/PLOT, $
                 yr=[0,10],order_scales=order_scales,weights=weights)
    nmout=strmid(f[ifile],0,strpos(f[ifile],'.ech'))+'c.ech'
    rdech,e,f[ifile],/RAW
    wdech,nmout,ee.head,ee.spec,sig=ee.sig,cont=ee.cont,wave=e.wave $
         ,orders=e.orders,/OVERWRITE
  endfor
  cd,'..'
endif
end
