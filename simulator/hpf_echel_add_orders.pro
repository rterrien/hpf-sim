;+
; NAME:
;  hpf_echel_add_orders
;
; PURPOSE:
;
;  Add orders to the original echellogram from Larry
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  out = hpf_echel_add_orders(model,blue_extra,red_extra)
;
; INPUTS:
;
;	model = optical model structure, with xs,ys,ws (n x norders) and orders (norders)
;
;	blue,red_extra = number of extra orders to add on each end
;	
;	
; OUTPUTS:
;	
;	out = modified model structure
;
; KEYWORD PARAMETERS:
;
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 04-28-14
;-


function hpf_echel_add_orders, model, blue_extra, red_extra
	
	max_order = max(model.orders,min=min_order)
	norders = n_elements(model.orders) + blue_extra + red_extra
	norders_orig = n_elements(model.orders)
	
	new_min_order = min_order - red_extra
	new_max_order = max_order + blue_extra
	orders_new_red = [indgen(red_extra) + new_min_order]
	orders_new_blue = [indgen(blue_extra) + max_order + 1]
	
	orders = model.orders
	ys = model.ys
	ws = model.ws
	
	;new_xs_blue = (model.xs[*,0],blue_extra)
	;new_xs_red = replicate(model.xs[*,0],red_extra)
	
	n_pts_per_order = (size(model.xs,/dimen))[0]
	
	new_xs = dblarr(n_pts_per_order,norders)
	new_ys = dblarr(n_pts_per_order,norders)
	new_ws = dblarr(n_pts_per_order,norders)
	new_orders = [orders_new_red,model.orders,orders_new_blue]
	
	new_ys_blue = dblarr(n_pts_per_order,blue_extra)
	new_ys_red = dblarr(n_pts_per_order,red_extra)
	new_ws_blue = dblarr(n_pts_per_order,blue_extra)
	new_ws_red = dblarr(n_pts_per_order,red_extra)
	
	for i=0, n_pts_per_order - 1 do begin
		start = [0d,0d]
		;ycs = mpfitfun('poly',1d/orders,ys[i,*],1d,start,yfit=ycs_fit)
		;wcs = mpfitfun('poly',1d/orders,ws[i,*],1d,start,yfit=wcs_fit)
		ycs = poly_fit(1d/orders,ys[i,*],1,/double)
		wcs = poly_fit(1d/orders,ws[i,*],1,/double)
		o1 = red_extra - 1
		o2 = red_extra + norders_orig - 1
		o3 = -1 * blue_extra
		new_xs[i,*] = model.xs[i,0]
		;new_ys[i,0:o1] = poly(1d/orders_new_red,ycs)
		;new_ys[i,o1 + 1:o2] = model.ys[i,*]
		;new_ys[i,o3: *] = poly(1d/orders_new_blue,ycs)
		new_ys[i,*] = poly(1d/new_orders,ycs)
		;new_ws[i,0:o1] = poly(1d/orders_new_red,wcs)
		;new_ws[i,o1+1:o2] = model.ws[i,*]
		;new_ws[i,o3:*] = poly(1d/orders_new_blue,wcs)
		new_ws[i,*] = poly(1d/new_orders,wcs)
	endfor
	
	model_out = {xs:new_xs, ys:new_ys, ws:new_ws, orders:new_orders}
	
	
	return,model_out
	
end
		

	