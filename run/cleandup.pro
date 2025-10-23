; clean the duplicats in the matched peaks

PRO cleandup,sid,outputstring
  ;sid = id2str(id)
  readcol,'matchpeak_'+sid+'.dat',xpix,ypix,cx,cy,sn,rcore,m,z,ihal,match,$
    format='I,I,D,D,D,D,D,D,I,I',/silent
  ; first remove 
  id = where(match eq 1)
  outputstring=string(n_elements(id),n_elements(match),$
               format='("matching peaks:",I,"/",I)')
  nt = n_elements(id)

  openw,1,'matchcluster_'+sid+'.dat'
  form = '(I5,I5,F8.2,F8.2,F8.4,F5.2,E11.4,F6.3,I8,I2)'
  printf,1,"# xpix - ypix - cxpix - cypix - sn - rcore - m - z - ihal - matchID"
  for i=0,nt-1 do begin
    mid = where(xpix eq xpix[id[i]] and ypix eq ypix[id[i]])
    if n_elements(mid) eq 1 then begin
      if mid[0] eq -1 then stop
      if mid[0] ne id[i] then stop
      printf,1,string(xpix[mid],ypix[mid],cx[mid],cy[mid],sn[mid],$
                      rcore[mid],m[mid],z[mid],ihal[mid],match[mid],format=form)
    endif else begin
      max = max(m[mid],maxid)
      mid = mid[maxid]
      printf,1,string(xpix[mid],ypix[mid],cx[mid],cy[mid],sn[mid],$
                      rcore[mid],m[mid],z[mid],ihal[mid],match[mid],format=form)
    endelse
  endfor
  close,1
END
