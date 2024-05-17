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
