;+
; convert ID to 4-digit string
;-
FUNCTION id2str,id
  if (id lt 0) or (id gt 9999) then begin
      print, "!! wrong ID !!"
      stop
  endif
  return, string(id,format='(I04)')
END
