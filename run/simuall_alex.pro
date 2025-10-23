PRO simuall_alex
  id = [1234,1327,4183,8372,9138]
  rseed=5614
  for i=0,4 do begin
    sid = string(id[i],format='(I04)')
    print,sid
    for j=0,0 do begin
        backid = sid;+'_'+string(j,format='(I02)')
        simuone_alex,sid,rseed=rseed,log=backid+'.log',backid=backid
    endfor
  endfor
END
