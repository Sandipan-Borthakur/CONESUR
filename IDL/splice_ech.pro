pro splice_ech,ech,wave,spec,blaz,index,sig=sig,ORDERS=orders,COLRANGE=colrange $
              ,WEIGHTS=weights,SCALING=scaling,ORDER_SCALES=order_scales $
              ,DEBUG=debug,YRANGE=yr,WRANGE=wrange
; splice_ech takes an ech-structure which must include the blaze functions
; as ech.cont, and the wavelength solution as ech.wave and an optional
; list of orders to splice (e.g. with polarization data odd and even "orders"
; must be spliced separately). splice_ech returns three 1D arrays: wavelength,
; spectrum and the blaze.

; Sanity check
  i=where(tag_names(ech) eq 'SPEC', has_spec)
  if(has_spec eq 0) then begin
    print,'Ilegitimate ech structure: must include SPEC tag'
    stop
  endif
  i=where(tag_names(ech) eq 'WAVE', has_wave)
  if(has_wave eq 0) then begin
    print,'splice_ech expects the wavelength solution to be included to the ech structure'
    stop
  endif
  i=where(tag_names(ech) eq 'SIG', has_sig)
  has_sig=(has_sig gt 0)?1B:0B
  i=where(tag_names(ech) eq 'CONT', has_cont)
  if(has_cont eq 0) then begin
    print,'splice_ech expects the blaze functions to be included to the ech structure as CONT'
    stop
  endif

  npix=n_elements(ech.spec[*,0]) ; Order length in pixels
  nord=n_elements(ech.spec[0,*]) ; Number of sp. orders

  weights=dblarr(npix,nord)+1
  if(keyword_set(colrange)) then begin
    sz=size(colrange)
    if(sz[0] ne 2) then begin
      print,'COLRANGE should match the number of spectral orders'
      help,colrange
      stop
    endif
    if(sz[2] ne nord) then begin
      print,'COLRANGE should match the number of spectral orders'
      help,colrange
      stop
    endif
    colr=colrange
  endif else begin
    colr=intarr(2,nord)
    colr[1,*]=npix-1
  endelse

  if(keyword_set(orders)) then begin   ; If a subset of orders was specified 
    sp=ech.spec[*,orders]              ; select only these orders
    ww=ech.wave[*,orders]
    bb=ech.cont[*,orders]
    if(has_sig) then unc=ech.sig[*,orders]
    npix=n_elements(sp[*,0]) ; Order length in pixels
    nord=n_elements(sp[0,*]) ; Number of sp. orders
    weights=dblarr(npix,nord)+1
  endif else begin
    sp=ech.spec
    ww=ech.wave
    if(keyword_set(scaling)) then begin
      bb=ech.cont>1.
      for iord=0,nord-1 do begin
        i0=colr[0,iord]
        i1=colr[1,iord]
;        scale=mean(top(ech.spec[i0:i1,iord],1.d5))/mean(top(ech.cont[i0:i1,iord],1.d5))
        scale=median(ech.spec[i0:i1,iord])/median(median(ech.cont[i0:i1,iord],5))
        bb[i0:i1,iord]=median(ech.cont[i0:i1,iord],5)*scale
        sp[i0:i1,iord]=ech.spec[i0:i1,iord]
      endfor
    endif else begin
      bb=ech.cont>1.
      for iord=0,nord-1 do begin
        i0=colr[0,iord]
        i1=colr[1,iord]
        bb[i0:i1,iord]=median(ech.cont[i0:i1,iord],5)
      endfor
    endelse
    if(has_sig) then unc=ech.sig
  endelse

;==================================================================================
  order_scales=replicate(1.,nord)
  order_overlap=replicate(-1,6,nord)
  dummy=max(median(sp/(bb>0.1),dimension=1),iord0); Order with largest signal
  beg1=colr[0,iord0]
  end1=colr[1,iord0]
  w1=ww[beg1:end1,iord0]
  s1=sp[beg1:end1,iord0]
  b1=bb[beg1:end1,iord0]
  if(has_sig) then sig1=unc[beg1:end1,iord0]
  if(iord0 gt 0) then begin                          ; Loop backward from iord0
    for iord=iord0-1,0,-1 do begin
      beg0=beg1                                      ; Shift current order to previous
      end0=end1
      w0=w1
      s0=s1
      b0=b1
      if(has_sig) then sig0=sig1
      beg1=colr[0,iord]                              ; New current order
      end1=colr[1,iord]
      w1=ww[beg1:end1,iord]
      s1=sp[beg1:end1,iord]
      b1=bb[beg1:end1,iord]
      if(has_sig) then sig1=unc[beg1:end1,iord]

      i0=where(w0 ge min(w1) and w0 le max(w1), ni0) ; Overlap within the previous order
      ii0=indgen(ni0)
      if(keyword_set(wrange)) then begin     ; If exclusion regions are given check the overlaps
        for iwrange=0,n_elements(wrange)/2-1 do begin
         iii0=where(w0[i0[ii0]] lt wrange[0,iwrange] or w0[i0[ii0]] gt wrange[1,iwrange], nii0)
         if nii0 gt 0 then ii0=ii0[iii0] else ni0=0
        endfor
      endif

      i1=where(w1 ge min(w0) and w1 le max(w0), ni1) ; Overlap within the current order
      ii1=indgen(ni1)
      if(keyword_set(wrange)) then begin     ; If exclusion regions are given check the overlaps
        for iwrange=0,n_elements(wrange)/2-1 do begin
         iii1=where(w1[i1[ii1]] lt wrange[0,iwrange] or w1[i1[ii1]] gt wrange[1,iwrange], nii1)
         if nii1 gt 0 then ii1=ii1[iii1] else ni1=0
        endfor
      endif

      if(ni0 gt 0 and ni1 gt 0) then begin ; We have overlapping orders
        tmpS0=bezier_interp(w1,s1,bezier_init(w1,s1,/DOUBLE),w0[i0],/DOUBLE)
        tmpB0=bezier_interp(w1,b1,bezier_init(w1,b1,/DOUBLE),w0[i0],/DOUBLE)
        if(has_sig) then tmpU0=bezier_interp(w1,sig1,bezier_init(w1,sig1,/DOUBLE),w0[i0],/DOUBLE)
        tmpS1=bezier_interp(w0,s0,bezier_init(w0,s0,/DOUBLE),w1[i1],/DOUBLE)
        tmpB1=bezier_interp(w0,b0,bezier_init(w0,b0,/DOUBLE),w1[i1],/DOUBLE)
        if(has_sig) then tmpU1=bezier_interp(w0,sig0,bezier_init(w0,sig0,/DOUBLE),w1[i1],/DOUBLE)
        if(keyword_set(scaling)) then begin
          scl0=total(s0[i0[ii0]]/b0[i0[ii0]])/total(tmpS0[ii0]/tmpB0[ii0])
          scl1=total(s1[i1[ii1]]/b1[i1[ii1]])/total(tmpS1[ii1]/tmpB1[ii1])
          scl =sqrt(scl0/scl1)
          s1=s1*scl
          if has_sig then sig1=sig1*scl
          tmpS0=tmpS0*scl
          tmpU0=tmpU0*scl
          order_scales[iord]=scl
        endif
        if(min(ww[0,iord]) lt min(ww[0,iord+1])) then begin ; Current order (1) is bluer than the previous (0)
          wgt0=dindgen(ni0)/(ni0-1) & wgt1=1.d0-wgt0
          weights[i0,iord+1]=weights[i0,iord+1]*s0[i0]/(s0[i0]*wgt0+tmpS0*wgt1)
          if(beg0 gt 0     ) then weights[     0:beg0-1,iord+1]=0.d0
          s0[i0]=s0[i0]*wgt0+tmpS0*wgt1
          b0[i0]=b0[i0]*wgt0+tmpB0*wgt1
          if(has_sig) then sig0[i0]=sqrt(sig0[i0]*sig0[i0]*wgt0+tmpU0*tmpU0*wgt1)
          wgt1=dindgen(ni1)/(ni1-1) & wgt0=1.d0-wgt1
          weights[i1,iord  ]=weights[i1,iord  ]*s1[i1]/(s1[i1]*wgt0+tmpS1*wgt1)
          if(end1 lt npix-1) then weights[end1+1:npix-1,iord  ]=0.d0
          s1[i1]=s1[i1]*wgt0+tmpS1*wgt1
          b1[i1]=b1[i1]*wgt0+tmpB1*wgt1
          if(has_sig) then sig1[i1]=sqrt(sig1[i1]*sig1[i1]*wgt0+tmpU1*tmpU1*wgt1)
;plot,w1,s1,xs=3,xr=minmax([w0,w1]),psym=3&oplot,w0,s0,col=c24(3)
        endif else begin
          wgt1=dindgen(ni0)/(ni0-1) & wgt0=1.d0-wgt1
          weights[i0,iord+1]=weights[i0,iord+1]*s0[i0]/(s0[i0]*wgt0+tmpS0*wgt1)
          if(end0 lt npix-1) then weights[end0+1:npix-1,iord+1]=0.d0
          s0[i0]=s0[i0]*wgt0+tmpS0*wgt1
          b0[i0]=b0[i0]*wgt0+tmpB0*wgt1
          if(has_sig) then sig0[i0]=sqrt(sig0[i0]*sig0[i0]*wgt0+tmpU0*tmpU0*wgt1)
          wgt0=dindgen(ni1)/(ni1-1) & wgt1=1.d0-wgt0
          weights[i1,iord  ]=weights[i1,iord  ]*s1[i1]/(s1[i1]*wgt0+tmpS1*wgt1)
          if(beg1 gt 0     ) then weights[     0:beg1-1,iord  ]=0.d0
          s1[i1]=s1[i1]*wgt0+tmpS1*wgt1
          b1[i1]=b1[i1]*wgt0+tmpB1*wgt1
          if(has_sig) then sig1[i1]=sqrt(sig1[i1]*sig1[i1]*wgt0+tmpU1*tmpU1*wgt1)
        endelse
        if(keyword_set(debug)) then begin
          plot,ww,sp/bb,psym=3,xs=1,yr=yr,ys=1
          oplot,w0,s0/b0,col=c24(2)
          oplot,w1,s1/b1,col=c24(3)
          answ=get_kbrd(1)
        endif
      endif else begin
        xmid=(beg0-end1)/2.d0
        scl0=top(s0/b0,1,/POLY)
        b0=b0*scl0
        scl0=top(s0/b0,1,/POLY)
        scl0=poly_fit(w0,scl0,1,/DOUBLE)
        scl1=top(s1/b1,1,/POLY)
        scl1=poly_fit(w1,scl1,1,/DOUBLE)
        xx=dindgen(101)/100*(min(w0)-max(w1))+max(w1)
        if(keyword_set(debug)) then begin
          plot,ww,sp/bb,psym=3,xs=1,yr=yr,ys=1
          oplot,w1,s1/b1,col=c24(2)
          oplot,xx,poly(xx,scl1),col=c24(2),line=3
          oplot,w0,s0/b0,col=c24(3)
          oplot,xx,poly(xx,scl0),col=c24(3),line=3
          answ=get_kbrd(1)
        endif
        scl =total(scl0[1]*scl1[1]*xx*xx $
                  +scl0[1]*scl1[0]*xx    $
                  +scl1[1]*scl0[0]*xx    $
                  +scl1[0]*scl0[0]) $
            /total((scl1[1]*xx+scl1[0])^2)
        s1=s1*scl
        order_scales[iord]=scl
;        tmpS0=tmpS0*scl
      endelse
      sp[beg0:end0,iord+1]=s0>0.
      bb[beg0:end0,iord+1]=b0>1.
      if(has_sig) then unc[beg0:end0,iord+1]=sig0
      sp[beg1:end1,iord]=s1>0.
      bb[beg1:end1,iord]=b1>1.
      if(has_sig) then unc[beg1:end1,iord]=sig1
    endfor
  endif
;==================================================================================
  beg1=colr[0,iord0]
  end1=colr[1,iord0]
  w1=ww[beg1:end1,iord0]
  s1=sp[beg1:end1,iord0]
  b1=bb[beg1:end1,iord0]
  if(has_sig) then sig1=unc[beg1:end1,iord0]
  if(iord0 lt nord-1) then begin                  ; Loop forward from iord0
    for iord=iord0+1,nord-1 do begin
      beg0=beg1                                       ; Shift current order to previous
      end0=end1
      w0=w1
      s0=s1
      b0=b1
      if(has_sig) then sig0=sig1
      beg1=colr[0,iord]                             ; New current order
      end1=colr[1,iord]
      w1=ww[beg1:end1,iord]
      s1=sp[beg1:end1,iord]
      b1=bb[beg1:end1,iord]
      if(has_sig) then sig1=unc[beg1:end1,iord]

      i0=where(w0 ge min(w1) and w0 le max(w1), ni0) ; Overlap within the previous order
      ii0=indgen(ni0)
      if(keyword_set(wrange)) then begin     ; If exclusion regions are given check the overlaps
        for iwrange=0,n_elements(wrange)/2-1 do begin
         iii0=where(w0[i0[ii0]] lt wrange[0,iwrange] or w0[i0[ii0]] gt wrange[1,iwrange], nii0)
         if nii0 gt 0 then ii0=ii0[iii0] else ni0=0
        endfor
      endif

      i1=where(w1 ge min(w0) and w1 le max(w0), ni1) ; Overlap within the current order
      ii1=indgen(ni1)
      if(keyword_set(wrange)) then begin     ; If exclusion regions are given check the overlaps
        for iwrange=0,n_elements(wrange)/2-1 do begin
         iii1=where(w1[i1[ii1]] lt wrange[0,iwrange] or w1[i1[ii1]] gt wrange[1,iwrange], nii1)
         if nii1 gt 0 then ii1=ii1[iii1] else ni1=0
        endfor
      endif

      if(ni0 gt 0 and ni1 gt 0) then begin
        tmpS0=bezier_interp(w1,s1,bezier_init(w1,s1,/DOUBLE),w0[i0],/DOUBLE)
        tmpB0=bezier_interp(w1,b1,bezier_init(w1,b1,/DOUBLE),w0[i0],/DOUBLE)
        if(has_sig) then tmpU0=bezier_interp(w1,sig1,bezier_init(w1,sig1,/DOUBLE),w0[i0],/DOUBLE)
        tmpS1=bezier_interp(w0,s0,bezier_init(w0,s0,/DOUBLE),w1[i1],/DOUBLE)
        tmpB1=bezier_interp(w0,b0,bezier_init(w0,b0,/DOUBLE),w1[i1],/DOUBLE)
        if(has_sig) then tmpU1=bezier_interp(w0,sig0,bezier_init(w0,sig0,/DOUBLE),w1[i1],/DOUBLE)
        if(keyword_set(scaling)) then begin
          scl0=total(abs(s0[i0[ii0]]/b0[i0[ii0]]))/total(abs(tmpS0[ii0]/tmpB0[ii0]))
          scl1=total(abs(s1[i1[ii1]]/b1[i1[ii1]]))/total(abs(tmpS1[ii1]/tmpB1[ii1]))
          scl =sqrt(scl0/scl1)
          s1=s1*scl
          if has_sig then sig1=sig1*scl
          tmpS0=tmpS0*scl
          tmpU0=tmpU0*scl
          order_scales[iord]=scl
        endif
        if(min(ww[0,iord]) lt min(ww[0,iord-1])) then begin ; Current order (1) is redder than the previous (0)
          wgt0=dindgen(ni0)/(ni0-1) & wgt1=1.d0-wgt0
          weights[i0,iord-1]=weights[i0,iord-1]*s0[i0]/(s0[i0]*wgt0+tmpS0*wgt1)
          if(end0 lt npix-1) then weights[end0+1:npix-1,iord-1]=0.d0
          s0[i0]=s0[i0]*wgt0+tmpS0*wgt1
          b0[i0]=b0[i0]*wgt0+tmpB0*wgt1
          if(has_sig) then sig0[i0]=sqrt(sig0[i0]*sig0[i0]*wgt0+tmpU0*tmpU0*wgt1)
          wgt1=dindgen(ni1)/(ni1-1) & wgt0=1.d0-wgt1
          weights[i1,iord  ]=weights[i1,iord  ]*s1[i1]/(s1[i1]*wgt0+tmpS1*wgt1)
          if(beg1 gt 0     ) then weights[     0:beg1-1,iord  ]=0.d0
          s1[i1]=s1[i1]*wgt0+tmpS1*wgt1
          b1[i1]=b1[i1]*wgt0+tmpB1*wgt1
          if(has_sig) then sig1[i1]=sqrt(sig1[i1]*sig1[i1]*wgt0+tmpU1*tmpU1*wgt1)
        endif else begin
          wgt1=dindgen(ni0)/(ni0-1) & wgt0=1.d0-wgt1
          weights[i0,iord-1]=weights[i0,iord-1]*s0[i0]/(s0[i0]*wgt0+tmpS0*wgt1)
          if(beg0 gt 0     ) then weights[     0:beg0-1,iord-1]=0.d0
          s0[i0]=s0[i0]*wgt0+tmpS0*wgt1
          b0[i0]=b0[i0]*wgt0+tmpB0*wgt1
          if(has_sig) then sig0[i0]=sqrt(sig0[i0]*sig0[i0]*wgt0+tmpU0*tmpU0*wgt1)
          wgt0=dindgen(ni1)/(ni1-1) & wgt1=1.d0-wgt0
          weights[i1,iord  ]=weights[i1,iord  ]*s1[i1]/(s1[i1]*wgt0+tmpS1*wgt1)
          if(end1 lt npix-1) then weights[end1+1:npix-1,iord  ]=0.d0
          s1[i1]=s1[i1]*wgt0+tmpS1*wgt1
          b1[i1]=b1[i1]*wgt0+tmpB1*wgt1
          if(has_sig) then sig1[i1]=sqrt(sig1[i1]*sig1[i1]*wgt0+tmpU1*tmpU1*wgt1)
        endelse
        if(keyword_set(debug)) then begin
          plot,ww,sp/bb,psym=3,xs=1,yr=yr,ys=1
          oplot,w0,s0/b0,col=c24(2)
          oplot,w1,s1/b1,col=c24(3)
          answ=get_kbrd(1)
        endif
      endif else begin
        xmid=(beg0-end1)/2.d0
        scl0=top(s0/b0,1,/POLY)
        scl0=poly_fit(w0,scl0,1,/DOUBLE)
        scl1=top(s1/b1,1,/POLY)
        b1=b1*scl1
        scl1=top(s1/b1,1,/POLY)
        scl1=poly_fit(w1,scl1,1,/DOUBLE)
        xx=dindgen(101)/100*(min(w0)-max(w1))+max(w1)
        scl =total(scl0[1]*scl1[1]*xx*xx $
                  +scl0[1]*scl1[0]*xx    $
                  +scl1[1]*scl0[0]*xx    $
                  +scl1[0]*scl0[0]) $
            /total((scl1[1]*xx+scl1[0])^2)
        order_scales[iord]=scl
        if(keyword_set(debug)) then begin
          plot,ww,sp/bb,psym=3,xs=1,yr=yr,ys=1
          oplot,w0,s0/b0,col=c24(2)
          oplot,xx,poly(xx,scl0),col=c24(2),line=3
          oplot,w1,s1/b1,col=c24(3)
          oplot,xx,poly(xx,scl1),col=c24(3),line=3
          answ=get_kbrd(1)
        endif
        s1=s1*scl
;        tmpS0=tmpS0*scl
      endelse
      sp[beg0:end0,iord-1]=s0>0.
      bb[beg0:end0,iord-1]=b0>1.
      if(has_sig) then unc[beg0:end0,iord-1]=sig0
      sp[beg1:end1,iord]=s1>0.
      bb[beg1:end1,iord]=b1>1.
      if(has_sig) then unc[beg1:end1,iord]=sig1
    endfor
  endif
  if(keyword_set(debug)) then begin
    plot,ww,sp/bb,psym=3,xs=1,yr=yr,ys=1
  endif
;============================================================================
; Form the output arrays
; INDEX array will track association of pixels with the original CCD pixels

  for iord=0,nord-1 do begin
    i1=colr[0,iord]
    i2=colr[1,iord]
    if(iord eq 0) then begin
      wave=ww[i1:i2,iord]
      spec=sp[i1:i2,iord]
      blaz=bb[i1:i2,iord]
      if(has_sig) then sig=unc[i1:i2,iord]
      index=replicate(iord,i2-i1+1)
    endif else begin
      wave=[wave,ww[i1:i2,iord]]
      spec=[spec,sp[i1:i2,iord]]
      blaz=[blaz,bb[i1:i2,iord]]
      if(has_sig) then sig=[sig,unc[i1:i2,iord]]
      index=[index,replicate(iord,i2-i1+1)]
    endelse
  endfor

  isort=sort(wave)
  wave=wave[isort]
  spec=spec[isort]
  blaz=blaz[isort]
  index=index[isort]
  if(has_sig) then sig=sig[isort]

  return
end
