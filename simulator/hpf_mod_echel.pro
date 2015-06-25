;+
; NAME:
;  hpf_mod_echel
;
; PURPOSE:
;
;  Alter the shape of the echellogram orders for diagnostic purposes
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  hpf_mod_echel, xs_in, ys_in, ws_in, xs_out, ys_out, ws_out
;
; INPUTS:
;
;	xs,ys,ws_in: The input x/y coordinates and wavelengths for points on the echellogram orders
;	
;	
; OUTPUTS:
;	
;	xs,ys,ws_out: The output x/y coordinates and wavelengths for points on the echellogram orders
;
; KEYWORD PARAMETERS:
;
;	ver: 0 = normal
;        1 = straight/horizontal
;        2 = straight/slanted
;        3 = curved but with approx half the swing in y
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-


pro hpf_mod_echel, xs_in, ys_in, ws_in, xs_out, ys_out, ws_out, ver=ver
	;ver = 0, return as normal
	;ver = 1, flat orders
	;ver = 2, linear y
	;ver = 3, curved with less curvature
	;ver = 4, flat orders + linear wavelengths
	ys = ys_in
	sy = size(ys,/dimen)
	if n_elements(ver) eq 0 then ver = 0
	ab = dindgen(sy[0])
	
	case ver of
	0: begin
		xs_out = xs_in
		ys_out = ys_in
		ws_out = ws_in
		return
	end
	1: begin
		for i=0, sy[1]-1 do begin
			ys[*,i] = mean(ys[*,i])
		endfor
		ws = ws_in
	end
	2: begin
		for i=0, sy[1]-1 do begin
			dy = ys[1:*,i] - ys[*,i]
			dy1 = dy[0] * .1d
			ys[*,i] = ys[0,i] + ab * dy1
		endfor
		ws = ws_in
	end
	3: begin
		for i=0, sy[1]-1 do begin
			ly = ys[0,i]
			for j=0, sy[0]-1 do begin
				my = (ly + ys[j,i])/2d
				ys[j,i] = my
			endfor
		endfor
		ws = ws_in
	end
	4: begin
		ws = ws_in
		sw = size(ws,/dimen)
		wx = dindgen(sw[0])
		for i=0, sy[1]-1 do begin
			ys[*,i] = mean(ys[*,i])
			dw = ws[1,i] - ws[0,i]
			wst = wx * dw + ws[0,i]
			ws[*,i] = wst
		endfor
	end
	endcase
	
	xs = xs_in
	sx = size(xs,/dimen)
	dx1 = xs[1,0] - xs[0,0]
;;	for i=1, sx[0]-1 do begin
;;		xs[i,0] = xs[0,0] + double(i) * dx1
;;	endfor
;;	
;;	for i=0, sx[0]-1 do xs[i,*] = xs[i,0]
	
	
;	ws = ws_in
;	sw = size(ws,/dimen)
;;	for i=0, sw[1] - 1 do begin
;;		dw = ws[1,i] - ws[0,i]
;;;;		dw = 1.1d * dw * 0.847
;;		for j=1, sw[0]-1 do ws[j,i] = ws[0,i] + double(j) * dw
;;	endfor
	
	xs_out = xs
	ys_out = ys
	ws_out = ws
end
	