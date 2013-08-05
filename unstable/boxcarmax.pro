pro boxcarmax, wl, fl, pix, wl_out, ms_out, keep_nan = keep_nan
	nel = n_elements(fl)
	ms = dblarr(nel)
	for i=0, nel-1 do begin
		if ~finite(fl[i]) then ms[i] = !values.f_nan else begin
			temp = max(fl[0 > (i - pix) : (nel-1) < (i + pix)],/nan)
			ms[i] = temp
		endelse
	endfor
	if keyword_set(keep_nan) then begin
		wl_out = wl
		fl_out = ms
		return
	endif
	g = where(finite(ms))
	wl_out = wl[g]
	ms_out = ms[g]
end
