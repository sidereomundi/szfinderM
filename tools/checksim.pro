; first round to check the sim
; focus on cluster and galaxy distribution

; input galaxy file
PRO inputgal,fname,xpix,ypix,z
  readcol,fname,isub,xpix,ypix,z,format='I, D, D, X, D'
END

PRO inputclu,fname,ihal,xpix,ypix,r,m,n
  readcol,fname,ihal,xpix,ypix,r,m,format='I, D, D, D, D',/SILENT,count=n
  if n eq 0 then return
  id = where(m ge 3e14,n)
  if n eq 0 then return
  xpix = 4096*xpix[id]
  ypix = 4096*ypix[id]
  r = 4096*r[id]
  m = m[id]
  ihal = ihal[id]
END

FUNCTION dist2d,xmin,xmax,ymin,ymax,cx,cy
  nx = xmax-xmin  
  ny = ymax-ymin
  dis = dblarr(nx+1,ny+1,/nozero)
  for i=0,nx do begin
    for j=0,ny do begin
       dis[i,j] = sqrt((xmin+i+0.5-cx)^2+(ymin+j+0.5-cy)^2)
    endfor
  endfor
  return,dis  
END

PRO extract_radial,img,xpix,ypix,rpix,rx,yy
; set up
  dim = size(img,/dimensions,/L64)
  xid = round(xpix,/L64)
  yid = round(ypix,/L64)
  rid = ceil(rpix,/L64)*2
; set up range
  xmin = xid - rid
  xmax = xid + rid
  ymin = yid - rid
  ymax = yid + rid
  if xmin lt 0 then xmin = 0
  if ymin lt 0 then ymin = 0
  if xmax ge dim[0] then xmax = dim[0]-1
  if ymax ge dim[1] then ymax = dim[1]-1
; set image
  oimg = img[xmin:xmax,ymin:ymax]
; extract radial profile
; find peak with in 5 pixel of xpix,ypix
  xpixmin = floor(xpix-5)
  xpixmax = floor(xpix+5)
  ypixmin = floor(ypix-5)
  ypixmax = floor(ypix+5)
  if xpixmin lt 0 then xpixmin = 0
  if ypixmin lt 0 then ypixmin = 0
  if xpixmax ge dim[0] then xpixmax = dim[0]-1
  if ypixmax ge dim[1] then ypixmax = dim[1]-1

  maxsn = max(img[xpixmin:xpixmax,ypixmin:ypixmax],id)
  ypeak = id/(xpixmax-xpixmin+1)
  xpeak = id - ypeak*(xpixmax-xpixmin+1) +xpixmin+0.5
  ypeak = ypeak+ypixmin+0.5

;  rad = dist2d(xmin[0],xmax[0],ymin[0],ymax[0],xpix,ypix)
  rad = dist2d(xmin[0],xmax[0],ymin[0],ymax[0],xpeak,ypeak)
  rad = rad/rpix

  histxy,rad,oimg,rx,yy,binsize=0.05
  yy = yy/max(yy)
; plot 2D image
  oimg = alog(oimg)
  oimg = (oimg-min(oimg))/(max(oimg)-min(oimg))*254

  dis = sqrt((xpeak-xpix)^2+(ypeak-ypix)^2)/rpix
  tvim,oimg,xtitle=string(xpix)+'('+string(xmin,format='(I4)')+')',ytitle=string(ypix)+'('+string(ymin,format='(I4)')+')',title='peak offset (in rpix unit):'+string(dis)
; locate peak
  nx=xmax-xmin+1
  ny=ymax-ymin+1
  centery = findgen(ny)
  centerx = fltarr(ny)+xpeak-xmin
  oplot,centerx,centery,color=255
  centerx = findgen(nx)
  centery = fltarr(nx)+ypeak-ymin
  oplot,centerx,centery,color=255
; locate cluster
  centery = findgen(ny)
  centerx = fltarr(ny)+xpix-xmin
  oplot,centerx,centery
  centerx = findgen(nx)
  centery = fltarr(nx)+ypix-ymin
  oplot,centerx,centery
END

PRO readgeometry,path,z,id
   filename=path+'wmap.geometry.dat'
   readcol,filename,id,z,f='I,F',skipline=3
END

; check the sz profile
; also general sz map quick look
PRO checksimclu4
  window,0,xsize=900,ysize=500
  !P.multi=[0,2,1]
  path = './LightCones/TSZ_combined_flat_8.8_8.8/1234/'
  rx=[]
;openw,10,'peakcenter/single_profile_peakcenter.html'
openw,10,'clustercenter/single_profile_clustercenter.html'
printf,10,'<html><body>'
  readgeometry,path,z,zid
  for i=0,32 do begin
    id = string(i*4+16,format='(I03)')
    iz = where(id eq zid)

    print,'input ',id,z[iz]
    filename = path+'wmap.'+id+'.cluster.dat'
    inputclu,filename,ihal,xpix,ypix,r,m,n
    if n eq 0 then continue
    filename = path+'wmap.'+id+'.a.fits'
    img = mrdfits(filename,/silent,/dscale)
    for j=0,n-1 do begin
       extract_radial,img,xpix[j],ypix[j],r[j],tx,ty
       plot,tx,ty,xrange=[0,3],yrange=[1e-4,1],/ylog,$
           xtitle='r/r!Dvir',ytitle='Normalized y!Dsz',title='cluster '+string(j,format='(I)')+' @ z ='+string(z[iz],format='(F4.2)')
       rx = min(tx)
; overlay beta model
       rbeta = findgen(30)*0.1+mean(rx)
       for k=0,5 do begin
          rcore = 0.1
          beta = 0.8+k*0.04
          ybeta = (1+(rbeta/rcore)^2)^((1-3*beta)/2)
          ybeta = ybeta/max(ybeta)
          oplot,rbeta,ybeta,linestyle=2,color=k*50+3
       endfor
;       filename='clustercenter/single_profile_clustercenter'+id+'_'+string(j,format='(I01)')+'.png'
       filename='peakcenter/single_profile_peakcenter'+id+'_'+string(j,format='(I01)')+'.png'
       write_png,filename,tvrd(true=1)
       printf,10,'<a href="'+filename+'">'+filename+'</a><br>'
    endfor
 endfor
printf,10,'</body></html>'
close,10
END

; check the sz profile
; also general sz map quick look
PRO checksimclu5,img
  window,0,xsize=900,ysize=500
  !P.multi=[0,2,1]
  path = './LightCones/TSZ_combined_flat_8.8_8.8/1234/'
  rx=[]
openw,10,'lightcone_peakcenter/lightcone_profile_peakcenter.html'
;openw,10,'lightcone_clustercenter/lightcone_profile_clustercenter.html'
printf,10,'<html><body>'
  readgeometry,path,z,zid
  for i=0,32 do begin
    id = string(i*4+16,format='(I03)')
    iz = where(id eq zid)

    print,'input ',id,z[iz]
    filename = path+'wmap.'+id+'.cluster.dat'
    inputclu,filename,ihal,xpix,ypix,r,m,n
    if n eq 0 then continue
    for j=0,n-1 do begin
       extract_radial,img,xpix[j],ypix[j],r[j],tx,ty
       plot,tx,ty,xrange=[0,3],yrange=[1e-4,1],/ylog,$
           xtitle='r/r!Dvir',ytitle='Normalized y!Dsz',title='cluster '+string(j,format='(I)')+' @ z ='+string(z[iz],format='(F4.2)')
       rx = min(tx)
; overlay beta model
       rbeta = findgen(30)*0.1+mean(rx)
       for k=0,5 do begin
          rcore = 0.1
          beta = 0.8+k*0.04
          ybeta = (1+(rbeta/rcore)^2)^((1-3*beta)/2)
          ybeta = ybeta/max(ybeta)
          oplot,rbeta,ybeta,linestyle=2,color=k*50+3
       endfor
;       filename='lightcone_clustercenter/lightcone_profile_clustercenter'+id+'_'+string(j,format='(I01)')+'.png'
       filename='lightcone_peakcenter/lightcone_profile_peakcenter'+id+'_'+string(j,format='(I01)')+'.png'
       write_png,filename,tvrd(true=1)
       printf,10,'<a href="'+filename+'">'+filename+'</a><br>'
    endfor
 endfor
printf,10,'</body></html>'
close,10
END

PRO checksim
;  checksimclu4
img = mrdfits('LightCones/TSZ_combined_flat_8.8_8.8/1234/wmap.a.fits',/silent,/dscale)
checksimclu5,img
END
