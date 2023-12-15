function make_cont, ech, blaze, TEMPLATE=template, DEBUG=debug $
                  , COLRANGE=col_range, YRANGE=yrange, PARAM=param $
                  , WRANGE=wrange, POLY=ipoly, PLOT=plot, SCALING=scaling $
                  , SPLICE_ONLY=splice_only,WEIGHTS=weights $
                  , ORDER_SCALES=order_scales, FILENAME=filename
;
  if(n_params() lt 2) then begin
    print,'Usage: ech_new=make_cont(echfile,blaze $'
    print,'               [,TEMPLATE=template[,DEBUG=debug $'
    print,'               [,/PLOT[,YRANGE=yrange[,PARAM=param $'
    print,'               [,COLRANGE=colrange[,WRANGE=wrange]]]]]]])'
    print,' ech       - ech-structure containing the raw spectrum and the wavelength solution.'
    print,' blaze     - a 2D array conting the blaze function.'
    print,' template  - a [2,*] array with a template continuum normalized spectrum'
    print,'             with the wavelength in the first column.'
    print,' inst_mode - a string that is used to select the default set of parameters.'
    print,' yrange    - sets vertical plot range for the plot.'
    print,' param     - an array of 4 parameters (see the code).'
    print,' colrange  - an array in form [2,nord] giving startin,ending columns'
    print,' wrange    - an array in form [[Wl1,Wr1],[Wl2,Wr2],...[WlK,WrK]] giving'
    print,'             left and right boundaries to excluded from the fit. Usefull'
    print,'             when dealing with emission lines or strong cosmics.'
    print,' filename  - optional label for plots'
    print,' /POLY     - use polynomial instead of optimal filtering.'
    print,' /SCALING  - scale spectral orders when splicing (recommended).'
    print,' /PLOT     - displays interative adjustments the continuum fit.'
    print,' /DEBUG    - allows examining every order and prevents making'
    print,'             modifications in the ECH file. This flag also sets the PLOT flag.'
    print,' /SPLICE_ONLY - return spliced spectrum, wavelength scale, blaze functions and'
    print,'             uncertainties in a single structure.'
    print,'The output ech structure contains spliced spectrum subdivided to the spectral'
    print,'intervals of the additional orders and the corresponding continuum normalization.'
    return,0
  endif

  if(size(ech, /TYPE) eq 8) then begin
    e = ech
    i=where(tag_names(e) eq 'SIG', has_sig)
    has_sig=(has_sig gt 0)?1B:0B
  endif else if(size(ech, /TYPE) eq 7) then begin
    rdech, e, ech, /NOCONT
    i=where(tag_names(e) eq 'SIG', has_sig)
    has_sig=(has_sig gt 0)?1B:0B
  endif else begin
    print,'MAKE_CONT: first parameter must be a valid ech structure or a file'
    print,'          with such structure.'
    return,-1
  endelse

  e_orig=e

  if(keyword_set(filename)) then fname=', File:'+filename else fname=''

  if(keyword_set(debug)) then plot=1

  i=where(tag_names(e) eq 'WAVE', ni)
  if(ni le 0) then begin
    print,'MAKE_CONT: the ech structure does not contain the wavelength scale!'
    return,-2
  endif 

  nord=n_elements(e.spec[0,*]) ; Number of orders
  npix=n_elements(e.spec[*,0])

  if(keyword_set(col_range)) then begin
    colrange=col_range
  endif else begin
    colrange=intarr(2,nord)
    colrange[0,*]=0
    colrange[1,*]=npix-1
  endelse

  if(keyword_set(template)) then begin
    template_wave=reform(template[0,*]) ; template wavelength scale
    template_sp  =reform(template[1,*]) ; template spectrum
;    i=uniq(template_wave)
;    template_wave=template_wave[i]
;    template_sp  =template_sp  [i]
  endif
  
  target=sxpar(e.head,'OBJECT')

  w=e.wave

  if(npix ne n_elements(blaze[*,0])) then begin
    print,'Inconsistent row length in blaze and spectrum arrays:'
    help,blaze,e.spec
    return,-3
  endif

  if(nord ne n_elements(blaze[0,*])) then begin
    print,'MAKE_CONT: inconsistent column length in blaze and spectrum arrays:'
    help,blaze,e.spec
    return,-4
  endif

  wmin=min(w[colrange[0,     0]:colrange[1,     0],     0]) $
      <min(w[colrange[0,nord-1]:colrange[1,nord-1],nord-1])
  wmax=max(w[colrange[0,     0]:colrange[1,     0],     0]) $
      >max(w[colrange[0,nord-1]:colrange[1,nord-1],nord-1])

  if(n_elements(param) ge 4) then begin
    par_0=param[0]  ; Number of iterations
    par_1=param[1]  ; Smoothness for initial guess
    par_2=1.d-4
    par_3=param[2]  ; Smoothness  for the polishing touch
    par_4=0.01*(1.d-2<1.d0/sqrt(median(e.spec)))
    par_5=param[3]  ; Vertical scaling
  endif else begin
    par_0=10
    par_1=5.d5
    par_2=1.d-4
    par_3=5.d6
    par_4=0.01*(1.d-2<1.d0/sqrt(median(e.spec)))
    par_5=1.d0
  endelse

  b=blaze>1.                             ; Avoid 0 blaze function
  for iord=0,nord-1 do begin
    i0=colrange[0,iord]
    i1=colrange[1,iord]
    b[i0:i1,iord]=middle(b[i0:i1,iord],1.)
  endfor

  if(keyword_set(yrange)) then yr=yrange; else yr=[0,2*median(sB)]

; Create an equispaced wavelength grid covering the whole thing
  dwave=abs(w[npix/2,nord/2]-w[npix/2-1,nord/2])*0.5d0
  nwave=ceil((wmax-wmin)/dwave)+1L
  wave=(wmax-wmin)*dindgen(nwave+1L)/nwave+wmin

  for iord=1,nord-1 do begin                 ; Mark gaps between non-overlapping orders
    miss=where(wave ge w[colrange[1,iord-1],iord-1] $
           and wave le w[colrange[0,iord  ],iord  ], nmiss)
    if(nmiss gt 0) then begin
      if(keyword_set(wrange)) then begin
        wrange=[[wrange],[minmax(wave[miss])]]
      endif else begin
        wrange=[[minmax(wave[miss])]]
      endelse
    endif
  endfor

  e.cont=b

;====================================================================================
  splice_ech,e,wsort,ss,bb,index,COLRANGE=colrange,SIG=sig,SCALING=scaling $ ; Splice
            ,DEBUG=debug,YRANGE=yr,ORDER_SCALES=order_scales,WEIGHTS=weights $
            ,WRANGE=wrange

  if keyword_set(splice_only) then begin
    out={wave:wsort, spec:ss, cont:bb, unc:sig}
    return, out
  endif

  j=uniq(wsort)
  wsort_u=wsort[j]

  if(keyword_set(template)) then begin
    tmp=bezier_interp(template_wave-wmin, template_sp $ ; Interpolate throughput calibrated on the Sun
       ,bezier_init  (template_wave-wmin, template_sp,/DOUBLE),wsort_u-wmin,/DOUBLE)
  endif else begin
    tmp=1
  endelse

  sB=ss[j]/bb[j]/tmp                          ; Divide by the blaze
;  sB=clean_cosmic(wsort_u-wmin,sB,PLOT=plot,NHIST=n_histogram,SMOOTH1=smooth1,SMOOTH2=smooth2)

  weight=middle(sB,0.5,x=wsort_u-wmin)
  nwgt=n_elements(weight)
  weight=weight/middle(weight,3.*par_1) $
        +[0.,2.*weight[1:nwgt-2]-weight[0:nwgt-3]-weight[2:nwgt-1],0]
  weight=weight>0.

  weight=bezier_interp(wsort_u-wmin,weight $  ; Interpolate blaze on an equispaced grid
        ,bezier_init  (wsort_u-wmin,weight,/DOUBLE),wave-wmin,/DOUBLE)
  weight=weight/max(weight)
; Incorporate user-rejected regions
  if(keyword_set(wrange)) then begin     ; If exclusion regions are given, set weights to 0
    for i=0,n_elements(wrange)/2-1 do begin
      weight[where(wave ge wrange[0,i] and wave le wrange[1,i])]=0.
    endfor
  endif

;  if(has_sig) then sig=sig/sqrt(bb)
  spec_B2=bezier_init  (wsort_u-wmin,sB,/DOUBLE)
  ss_B   =bezier_interp(wsort_u-wmin,sB,spec_B2,wave-wmin,/DOUBLE)

  bbb=middle(bb,1.)

  niter=0
  cont_B=1.d0
iter:
;======================================================================================
; All the same except for the variable names
  if(keyword_set(wrange)) then begin     ; If exclusion regions are given, set weights to 0
    for i=0,n_elements(wrange)/2-1 do begin
      ss_B[where(wave ge wrange[0,i] and wave le wrange[1,i])]=0.
    endfor
  endif

  c=ss_B/cont_B

  if(keyword_set(ipoly)) then begin
    if(keyword_set(template)) then begin
      for ii=0,par_0 do c=middle(c,(ipoly gt 2)?ipoly:15,eps=par_2,WEIGHT=weight,/POLY)>c
      c=middle(c,(ipoly gt 2)?ipoly:15,eps=par_4,WEIGHT=weight,/POLY)*cont_B
    endif else begin
      for ii=0,par_0 do c=top(c,(ipoly gt 2)?ipoly:15,eps=par_2,WEIGHT=weight,/POLY)>c
      c=top(c,(ipoly gt 2)?ipoly:15,eps=par_4,WEIGHT=weight,/POLY)*cont_B
    endelse
  endif else begin
    if(keyword_set(template)) then begin
      for ii=0,par_0 do c=middle(c,par_1,eps=par_2,WEIGHT=weight,LAM2=par_3)>c
      c=middle(c,par_1,eps=par_4,WEIGHT=weight,LAM2=par_3)*cont_B
    endif else begin
      for ii=0,par_0 do c=top(c,par_1,eps=par_2,WEIGHT=weight,LAM2=par_3)>c
      c=top(c,par_1,eps=par_4,WEIGHT=weight,LAM2=par_3)*cont_B
    endelse
  endelse

  cont_B=c*par_5
  cont_B=middle(cont_B,1.)

  if(keyword_set(plot)) then begin
    !p.region=0
    !p.position=0
    !p.multi=0
    plot,wsort_u,sB,xs=3,xtit='Wavelength',tit=target+' Iteration:'+strtrim(niter,2) $
        +fname,yr=yr,psym=3,/NODATA;,/YLOG
    oplot,wsort_u,sB,col=c24(3),psym=3
    oplot,wave,cont_B,thick=1,psym=3
  endif

  if(keyword_set(debug)) then begin
    answ=get_kbrd(1)
  endif else if(keyword_set(plot)) then begin
    wait,1
  endif

; Interpolate the continuum back on every sp. order
  cont_B2=bezier_init(wave-wmin,cont_B,/DOUBLE)
  cont=bezier_interp(wave-wmin,cont_B,cont_B2,wsort-wmin,/DOUBLE)

  for iord=0,nord-1 do begin
    i1=colrange[0,iord]
    i2=colrange[1,iord]
    i=where(index eq iord)
    e.wave[i1:i2,iord]=wsort[i]
    e.spec[i1:i2,iord]=ss[i]
    if(has_sig) then e.sig[i1:i2,iord]=sig[i]
    e.cont[*,iord]=1.
    e.cont[i1:i2,iord]=cont[i]*bbb[i]

; Inspect every sp. order
    if(keyword_set(plot)) then begin
      !p.multi=[0,1,2]
      d=max(e.spec[i1:i2,iord],j)
      d=d/b[j+i1,iord]
      cc=replicate(1.,npix)
;      cc[i1:i2]=e.cont[i1:i2,iord]*b[i1:i2,iord]/bbb[i]
;      cc[i1:i2]=middle(cc[i1:i2],1.d0)
      cc[i1:i2]=e.cont[i1:i2,iord]
      plot,e.wave[*,iord],cc,xs=1,xr=minmax(w[*,iord]) $
          ,ys=3,yr=[0,max(b[i1:i2,iord]*d)>max(e.cont[*,iord])>max(e.spec[i1:i2,iord])] $
          ,tit='Order #'+strtrim(iord+1,2) $
          +' Iteration: '+strtrim(niter+1,2)+' of '+strtrim(round(par_0),2) $
          +fname 
;      oplot,w[i1:i2,iord],b[i1:i2,iord]*d,col=c24(4)
;      oplot,e.wave[*,iord],cc,col=c24(3)
;      oplot,ech.wave[i1:i2,iord],ech.spec[i1:i2,iord]
      oplot,w[i1:i2,iord],b[i1:i2,iord]*d,col=c24(4)
      oplot,e.wave[*,iord],e.cont[*,iord],col=c24(3)
      oplot,e.wave[i1:i2,iord],e.spec[i1:i2,iord]
      legend,['Extracted spectrum','Blaze function','Derived continuum'] $
                                                   ,col=c24([1,4,3]),line=[0,0,0]
      plot,e.wave[i1:i2,iord],e.spec[i1:i2,iord]/e.cont[i1:i2,iord] $
                             ,xs=1,xr=minmax(w[*,iord]),yr=[0,1.1],xtit='Wavelength [Å]'
      oplot,!x.crange,[1,1],line=3,col=c24(3)
      if(keyword_set(debug)) then answ=get_kbrd(1) else wait,0.5
    endif
  endfor
  niter=niter+1
  if(niter lt par_0) then begin
    weight=(ss_b/cont_B)<(cont_B/(ss_b>1))
    goto,iter
  endif

;  for iord=0,nord-1 do begin
;    i1=colrange[0,iord]
;    i2=colrange[1,iord]
;    i=where(index eq iord)
;;    w=e_orig.spec[i1:i2,iord]/e.spec[i1:i2,iord]
;    w=weights[i1:i2,iord]/order_scales[iord]
;    e.cont[*,iord]=1.
;    e.cont[i1:i2,iord]=(cont[i]*bbb[i]*w)>1.
;    plot,e_orig.wave[*,iord],e_orig.spec[*,iord],xs=1
;    oplot,e_orig.wave[*,iord],e.cont[*,iord],col=c24(3),thick=3
;    plot,e_orig.wave[i1:i2,iord],e_orig.spec[i1:i2,iord]/e.cont[i1:i2,iord],xs=1,yr=[0,1.2]
;    oplot,!x.crange,[1,1],line=3
;    stop
;  endfor

  !p.multi=[0,1,2]
  plot,wsort,sB,xtit='Wavelength',tit=target+fname $
      ,xs=1,xr=minmax(e.wave),yr=yr,psym=3,/NODATA;,/YLOG
  oplot,wsort,sB,col=c24(3),psym=3
  oplot,wave,cont_B,thick=1,psym=3

  plot,e.wave,e.spec/e.cont,xs=1,yr=[0,1.1],psym=3,xtit='Wavelength [Å]'
  oplot,!x.crange,[1,1],line=3,col=c24(3)
  !p.multi=0

  return,e

end
