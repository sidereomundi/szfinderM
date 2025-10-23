;+
; simulate one mock SPT map set
; check the filtered maps
; output result
;-

PRO recordcommand,command,output,fileid
  printf,fileid,"########### "+command
  for i=0,n_elements(output)-1 do printf,fileid,output[i]
END

PRO simuone_alex,simid,rseed=rseed,backid=backid,log=log
; rseed: random seed
; backid: copy mock files to *back*
; log: save result to "log" file
  if not keyword_set(rseed) then rseed = 917
  if keyword_set(log) then openw,101,log
  ran1 = long(randomu(rseed)*1000000)
;  ran2 = long(randomu(rseed)*10000)
print,backid,ran1
; generate mock maps
;  command = 'mkoneszmap '+string(ran1)+' 152.4 18.0 beam150 '+path+$
;    'wmap.a.fits spt150.fits '+string(ran2)
;  spawn,command,output
;  if keyword_set(log) then recordcommand,command,output,101

;  ran2 = long(randomu(rseed)*10000)
;  command = 'mkoneszmap '+string(ran1)+' 95. 44. beam95 '+path+$
;    'wmap.a.fits spt95.fits '+string(ran2)
;  spawn,command,output
;  if keyword_set(log) then recordcommand,command,output,101
  ran2 = long(randomu(rseed+1)*1000000)
  command = './mkoneszmap '+string(ran1)+' 98.0780 44 beam95T '+backid+'_90.fits spt95.fits '+string(ran2)
  spawn,command,output
  if keyword_set(log) then recordcommand,command,output,101
  ran2 = long(randomu(rseed+2)*1000000)
  command = './mkoneszmap '+string(ran1)+' 154.477 18 beam150T '+backid+'_150.fits spt150.fits '+string(ran2)
  spawn,command,output
  if keyword_set(log) then recordcommand,command,output,101

; filtering
  command='./filtermapM spt'
  spawn,command,output
  if keyword_set(log) then recordcommand,command,output,101

; peak finder
  command='./szpeakfinder filteredM '+simid
  spawn,command,output
  if keyword_set(log) then recordcommand,command,output,101 

; peak compare
  command='./szpeakcmp '+simid
  spawn,command,output
  if keyword_set(log) then recordcommand,command,output,101 

; clean low mass match
  cleandup,simid,output
  if keyword_set(log) then recordcommand,'./cleandup',output,101

  if keyword_set(backid) then begin
      spawn,'mv spt150.fits spt150_'+backid+'.fits'
      spawn,'mv spt95.fits spt95_'+backid+'.fits'
      spawn,'mv szpeaks_'+simid+'.dat szpeaks_'+backid+'.dat'
      spawn,'mv matchpeak_'+simid+'.dat matchpeak_'+backid+'.dat'
      spawn,'mv matchcluster_'+simid+'.dat matchcluster_'+backid+'.dat'
      if keyword_set(log) then recordcommand,"moving to "+backid,"",101
  endif

  if keyword_set(log) then close,101
END
