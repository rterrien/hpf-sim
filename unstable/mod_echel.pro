pro mod_echel, xs_in, ys_in, ws_in, xs_out, ys_out, ws_out
	ys = ys_in
	sy = size(ys,/dimen)
	for i=0, sy[1]-1 do begin
		ys[*,i] = mean(ys[*,i])
	endfor
	
	
	xs = xs_in
	sx = size(xs,/dimen)
	dx1 = xs[1,0] - xs[0,0]
;;	for i=1, sx[0]-1 do begin
;;		xs[i,0] = xs[0,0] + double(i) * dx1
;;	endfor
;;	
;;	for i=0, sx[0]-1 do xs[i,*] = xs[i,0]
	
	
	ws = ws_in
	sw = size(ws,/dimen)
;;	for i=0, sw[1] - 1 do begin
;;		dw = ws[1,i] - ws[0,i]
;;;;		dw = 1.1d * dw * 0.847
;;		for j=1, sw[0]-1 do ws[j,i] = ws[0,i] + double(j) * dw
;;	endfor
	
	xs_out = xs
	ys_out = ys
	ws_out = ws
end
	