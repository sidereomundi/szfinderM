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

; check galaxy redshift distribution
PRO checksimgal
  path = './LightCones/TSZ_combined_flat_8.8_8.8/2345/'
  zbin = []
  for i=0,32 do begin
    id = string(i*4+16,format='(I03)')
    print,'input ',id
    filename = path+'wmap.'+id+'.galaxies.dat'
    inputgal,filename,tx,ty,tz
    zbin = [[zbin],[minmax(tz)]]
 endfor
 plot,zbin[0,*],psym=-1
 oplot,zbin[1,*],psym=-1,color=220
END

; check the sz profile
; I general sz map quick look
PRO checksimclu1
  path = './LightCones/TSZ_combined_flat_8.8_8.8/1234/'
  for i=0,32 do begin
    id = string(i*4+16,format='(I03)')
    print,'input ',id
    filename = path+'wmap.'+id+'.cluster.dat'
    inputclu,filename,ihal,xpix,ypix,r,m,n
    if n eq 0 then continue
    filename = 'clu'+id+'.reg'
    openw,1,filename
      for j=0,n-1 do begin
        printf,1,'circle (',xpix[j],ypix[j],r[j],') # text = {',ihal[j],'}'
      endfor
    close,1
    command = 'ds9 -geometry 600x800 -scale log -cmap Heat '+path+'wmap.'+id+$
              '.a.fits -zoom to fit -regions '+filename+' -saveimage jpeg clu'+id+'.jpg'
    spawn,command
  endfor
END

; check the repeatness of clusters in the simulation
PRO checksimclur,simid
  if not keyword_set(simid) then simid=1234
  path = './LightCones/TSZ_combined_flat_13_13/'+string(simid,format='(I4)')+'/'
  tihal = []
  for i=5,32 do begin
    id = string(i*4+16,format='(I03)')
    filename = path+'wmap.'+id+'.cluster.dat'
    inputclu,filename,ihal,xpix,ypix,r,m,n
    if n eq 0 then continue
    tihal = [tihal,ihal]
  endfor
  histx,tihal,x,y,binsize=1
  plot,x,y,psym=10,xtitle='ihal',ytitle='count',title='repeatness of cluster in lightcone'
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

  rad = dist2d(xmin[0],xmax[0],ymin[0],ymax[0],xpix,ypix)
;  rad = dist2d(xmin[0],xmax[0],ymin[0],ymax[0],xpeak,ypeak)
  rad = rad/rpix
  histxy,rad,oimg,rx,yy,binsize=0.05
  yy = yy/max(yy)
END

; check the sz profile
; I general sz map quick look
PRO checksimclu2
  path = './LightCones/TSZ_combined_flat_8.8_8.8/1234/'
    plot,[0],/nodata,xrange=[0,3],yrange=[1e-4,1],/ylog,$
           xtitle='r/r!Dvir',ytitle='Normalized y!Dsz',title='sz profile'
    rx=[]
    yy=[]
  for i=0,32 do begin
    id = string(i*4+16,format='(I03)')
    print,'input ',id
    filename = path+'wmap.'+id+'.cluster.dat'
    inputclu,filename,ihal,xpix,ypix,r,m,n
    if n eq 0 then continue
    filename = path+'wmap.'+id+'.a.fits'
    img = mrdfits(filename,/silent,/dscale)

    for j=0,n-1 do begin
       extract_radial,img,xpix[j],ypix[j],r[j],tx,ty
       oplot,tx,ty,color=i*9
       rx = [[rx],min(tx)]
    endfor
 endfor
; overlay beta model
rbeta = findgen(30)*0.1+mean(rx)
for i=0,5 do begin
;rcore = i*0.02+0.05
rcore = 0.1
beta = 0.8+i*0.04
ybeta = (1+(rbeta/rcore)^2)^((1-3*beta)/2)
ybeta = ybeta/max(ybeta)
oplot,rbeta,ybeta,thick=2,linestyle=2
endfor
END

; check the sz profile
; I general sz map quick look
PRO checksimclu3,img
  path = './LightCones/TSZ_combined_flat_8.8_8.8/1234/'
    plot,[0],/nodata,xrange=[0,3],yrange=[1e-3,1],/ylog,$
           xtitle='r/r!Dvir',ytitle='Normalized y!Dsz',title='sz profile'
    rx=[]
    yy=[]
  for i=0,32 do begin
    id = string(i*4+16,format='(I03)')
    print,'input ',id
    filename = path+'wmap.'+id+'.cluster.dat'
    inputclu,filename,ihal,xpix,ypix,r,m,n
    if n eq 0 then continue
;    filename = path+'wmap.'+id+'.a.fits'
;    img = mrdfits(filename,/silent,/dscale)

    for j=0,n-1 do begin
       extract_radial,img,xpix[j],ypix[j],r[j],tx,ty
       oplot,tx,ty,color=i*9
       rx = [[rx],min(tx)]
    endfor
 endfor
; overlay beta model
rbeta = findgen(30)*0.1+mean(rx)
for i=0,5 do begin
;rcore = i*0.02+0.05
rcore = 0.1
beta = 0.8+i*0.04
ybeta = (1+(rbeta/rcore)^2)^((1-3*beta)/2)
ybeta = ybeta/max(ybeta)
oplot,rbeta,ybeta,thick=2
endfor

END

PRO checkclupos
END

PRO checksim
window,0,xsize=800,ysize=800
!P.multi=0
;   img = mrdfits('LightCones/TSZ_combined_flat_8.8_8.8/1234/wmap.a.fits',/silent,/dscale)
;   checksimclu3,img
  checksimclu2
END
