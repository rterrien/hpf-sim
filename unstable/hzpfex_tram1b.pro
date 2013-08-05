;; hzpfex_tram1.pro (and nonan.pro)
;; HPF DETECTOR SIMULATOR
;; This is the extraction code
;; It takes data from the image files, adds up with pre-defined orders
;; (profile fitting will be implemented in the future?)
;; 
;; CALLS
;; 
;; 
;; DIRECTLY MODIFIES
;;
;; CALLED BY
;; extract
;;
;; NOTES
;; 
;;
;; parameters:
;; infile - the input image file
;; inwlfile - the input wavelength image
;; outfile - output file
;; varfile - do not use
;; tellfile - telluric spectrum .sav file for a simple telluric correction

PRO nonan, array

;;Removes nan's in a 2-d array by averaging over their
;;neighbors
;;from chad's nonaninf.pro

  s=SIZE(array)
  
  junk = where(finite(array), complement=indx)

  ;;there is nothing to do
  IF indx[0] EQ -1 THEN RETURN

  FOR i=0,N_ELEMENTS(indx)-1 DO BEGIN
     xy=WHEREUNWRAP(array,indx[i])
     x1=(xy[0]-1) > 0
     x2=(xy[0]+1) < (s[1]-1)
     y1=(xy[1]-1) > 0
     y2=(xy[1]+1) < (s[2]-1)
     array[xy[0],xy[1]]=MEAN(array[x1:x2,y1:y2],/NaN)
  ENDFOR

END

function locmax, arr, width

;return an array of the same size as the input
;with original values where the mins are and 0's elsewhere

if n_elements(arr) lt width then message, 'arr is too small'
if (width mod 2) eq 0 then message, 'width must be odd'

ic = (width-1)/2
arrc = arr(ic:*)
b=bytarr(n_elements(arr)-width+1) + 1b
for i=1, ic do $
   b = b and (arrc ge arr(ic-i:*)) and (arrc ge arr(ic+i:*))
return, [arr(0:ic-1)*0, b*arrc, arr(0:ic-1)*0]

end 




PRO hzpfex_tram1b, infile, inwlfile, outfile, varfile = varfile, tellfile = tellfile, crimg = crimg, diag_output = diag_output, orders_lambdalow = orders_lambdalow, orders_lambdahigh = orders_lambdahigh, orders_gaps = orders_gaps, fiber_scale = fiber_scale, fiber_core_um = fiber_core_um, fiber_cladding_um = fiber_cladding_um, fiber_buffer_um = fiber_buffer_um,nfibers = nfibers, fiber_extra_sep_um = fiber_extra_sep_um, echellogram_file = echellogram_file, echellogram_wl_file = echellogram_wl_file

	if n_elements(echellogram_wl_file) eq 0 then echellogram_wl_file = 'hpf demag=2-0x f8-5 2012dec15 v10-1-wavelengths.dat'
	if n_elements(echellogram_file) eq 0 then echellogram_file = 'hpf demag=2-0x f8-5 2012dec15 v10-1-echelleogram.dat'
	
	pixel_size = 18d-3 ;mm
	;get the optical model parameters
	inp = dblarr(23,29)
	openr,1,echellogram_wl_file
	readf,1,inp
	close,1

	inp2 = dblarr(42,29)
	openr,1,echellogram_file
	readf,1,inp2
	close,1
	evens = indgen(21) * 2
	odds = indgen(21) * 2 + 1
	ys = inp2[odds,*]
	xs = inp2[evens,*]
	ws = inp[2:*,*]
	
	;ys = mod_echel(ys)
;;	mod_echel,xs,ys,ws,xs1,ys1,ws1
;;	xs = xs1
;;	ys = ys1
;;	ws = ws1
	
	buffer = 100 ; the extra pixels (/2) on each side for the big arrays which get filled in and convolved
	upfactor = 1d 
	
	ys = -1d * ys ;flip y around (stuart does this in his code)
	norders = (size(inp2))[2]
	
	warray=MAKE_ARRAY(2048, norders, /Double, Value=!Values.F_NAN)
	
	xs_pix = xs / pixel_size
	ys_pix = ys / pixel_size
	;xs,ys are the xy locations of the order sample points in mm
	;xs,ys_pix are the xy locations in pixels (with -1024 to 1023, not 0 to 2047)
	
	
	image=MRDFITS(infile, 0, hdr)


	;read in the stellar spectrum
	;input_wl = w * 1d3
	;input_fl = f
	
	extracted_fl = dblarr(2048,norders)
	extracted_wl = dblarr(2048,norders)

	
	
	xs_base = ((dindgen((2048+buffer)*upfactor) - (1024.+buffer/2)*upfactor)/upfactor)
	ws_base = dblarr((2048 + buffer)*upfactor,norders)
	;find left and right limits of each pixel for making the spectrum
	xs_base_left = xs_base - 0.5
	xs_base_right = xs_base + 0.5
	ws_base_left = dindgen((2048 + buffer) * upfactor,norders)
	ws_base_right = dindgen((2048 + buffer) * upfactor,norders)
	
	ys_base = dblarr((2048 + buffer) * upfactor,norders)
	xs_base_2d = dblarr((2048 + buffer) * upfactor,norders)

	for i=0, norders-1 do begin
		;for each order, there are 2048 x pixels
		;derive the wl for each pixel based on the x/wl map alone
		ws_base[*,i] = interpol(ws[*,i],xs_pix[*,i],xs_base,/spline)
		;same for left and right limits of each pixel
		ws_base_left[*,i] = interpol(ws[*,i],xs_pix[*,i],xs_base_left,/spline)
		ws_base_right[*,i] = interpol(ws[*,i],xs_pix[*,i],xs_base_right,/spline)
		;find the y for each x based on the x/y map alone
		ys_base[*,i] = interpol(ys_pix[*,i],xs_pix[*,i],xs_base,/spline)
		;replicate the x base array for convolving and plotting
		xs_base_2d[*,i] = xs_base
	endfor
	
	;create a kernel based on the fiber inputs
	;assusme a row of 2d gaussians for now
	;kernel size in pixels
	kernel_size_pix = fiber_core_um * fiber_scale
	fiber_extra_sep_pix = fiber_extra_sep_um * fiber_scale
	;separation between gaussians
	kernel_sep_pix = fiber_buffer_um * fiber_scale + fiber_extra_sep_pix

	;number of gaussians = nfibers
	y_size_kernel = (floor((nfibers + 1d) * kernel_sep_pix)); + kernel_size_pix) + kernel_sep_pix)
	x_size_kernel = (floor(3. * kernel_size_pix))*upfactor
	if x_size_kernel mod 2 eq 0 then x_size_kernel += 1.
	kernel = dblarr(x_size_kernel,y_size_kernel)
	x_center_kernel = floor(x_size_kernel / 2)

	extract_width = round(kernel_size_pix/2d + 1)

	;keep track of which orders are entirely contained on the chip and extract only these
	good_orders = intarr(norders)

	;stop
	
	;;;;;; temporary
	restore,inwlfile
	extracted_wl2 = warray
	;stop
	;;;;;;;;;;;;;;;;
	extracted_wl = extracted_wl2

	
	for i=0, norders-1 do begin
	
		;create the rectified arrays
		xs_rect = dblarr((2048+buffer)*upfactor,200)
		ys_rect = dblarr((2048+buffer)*upfactor,200)
		fs_rect = dblarr((2048+buffer)*upfactor,200)

		;fill in the xs with all the same value,
		;fill in the ys with sliding values based on the original
		for j=0, 200-1 do begin
			xs_rect[*,j] = xs_base
			ypos = j - 100d
			ys_rect[*,j] = ys_base[*,i] + ypos
		endfor

		;stop
		
		
		nx_rect = (size(ys_rect))[1]
		ny_rect = (size(ys_rect))[2]
		x_grid_rect_d = rebin(lindgen(nx_rect),nx_rect,ny_rect)
		y_grid_rect_d = rebin(1#lindgen(ny_rect),nx_rect,ny_rect)
	
		xs_d_offset = min(xs_rect)
		ys_d_offset = min(ys_rect)
		xs_rect_d = (xs_rect - xs_d_offset) * upfactor
		ys_rect_d = (ys_rect - ys_d_offset); * upfactor
		x_grid_rect = x_grid_rect_d + xs_d_offset
		y_grid_rect = y_grid_rect_d + ys_d_offset
	
		x_size_d = max(floor(xs_rect_d)) - min(floor(xs_rect_d))
		y_size_d = max(floor(ys_rect_d)) - min(floor(ys_rect_d))
		if x_size_d mod 2 eq 1 then x_size_d += 1
		if y_size_d mod 2 eq 1 then y_size_d += 1
	
		y_min = min(floor(ys_rect))
		y_max = max(floor(ys_rect))

		
		;stop
		
		
		;find the chunk of the detector that this corresponds to
		y0_det = min(y_grid_rect) + 1024; + buffer/2
		y1_det = y0_det + y_size_d - 1
		y0_subarr = 0
		y1_subarr = y_size_d - 1
		x0 = buffer/2
		x1 = x0 + 2047
		xd0 = 0
		xd1 = 2047
		
		mid = (100. + ys_d_offset + 1024.) < 100.

		good_orders[i] = 1
		
		if y0_det + (100. - y_size_kernel/3.) lt 0 then good_orders[i] = 0
		if y1_det - (100. - y_size_kernel/3.) gt 2047 then good_orders[i] = 0

		;stop
		if y0_det lt 0 then begin
			extra = abs(y0_det)
			y0_det = 0
			y0_subarr = extra
		endif
		if y1_det gt 2047 then begin
			extra = floor(max(y_grid_rect)+1024) - 2047
			y1_det = 2047
			y1_subarr = y_size_d - 1 - extra
		endif
		
		;stop
		
		
		print,'order ',i,' = ',good_orders[i]
		
		if y1_det lt 0 or y0_det gt 2047 then continue
		if good_orders[i] ne 1 then continue
		
		t1 = image[xd0:xd1 , y0_det:y1_det]
		

		t1a = dblarr(2047+buffer+1,y_size_d < round((abs(y0_det - y1_det) + 1)))
		t1a[buffer/2:2047+buffer/2,*] = t1
		;t2 =  warp_tri(xs_rect_d,ys_rect_d,x_grid_rect_d,y_grid_rect_d,t1,output_size = [x_size_d,y_size_d],/quintic)
		youts = ny_rect < (size(t1a,/dimen))[1]
		t2 = warp_tri(x_grid_rect_d,y_grid_rect_d,xs_rect_d,ys_rect_d,t1a,output_size = [nx_rect,ny_rect],/quintic)
		
		;t3 = warp_tri(xs_rect_d,ys_rect_d,x_grid_rect_d,y_grid_rect_d,t1a,output_size = [nx_rect,ny_rect])
	
		;extracted_wl[*,i] = ws_base[buffer/2:2047+buffer/2,i]
		;warray[*,i] = (rebin(ws_base,2048+buffer,norders))[buffer/2 : 2047 + buffer/2,i]
		
		


		
		for j=buffer/2, 2047+buffer/2 do begin ; loop over x pixels (recall that the overall array is "big")
			column = t2[j,mid-y_size_kernel/3:100+y_size_kernel/3]
			if n_elements(column) lt 3 then begin
				temptot = total(column) 
				;if temptot ne 0 then stop
			endif else begin
				colmaxs = locmax(column,3)
				maxlocs = where(colmaxs ne 0,nmaxs)
				;if nmaxs ne nfibers then stop
				temptot = 0d
				for k=0, nmaxs-1 do begin
					k1 = 0 > (maxlocs[k] - extract_width)
					k2 = (maxlocs[k] + extract_width) < (n_elements(column) - 1)
					temptot += total(column[k1:k2])
					;stop
				endfor
				;if j eq 200 then stop
			endelse
			extracted_fl[j-buffer/2,i] = temptot
		endfor
		
		wls = reform(extracted_wl[*,i])
		mp = (wls + wls[1:*])/2.
		xsi = mp[1:*] - mp
		xsi = [xsi[0],xsi,xsi[-1]]

		extracted_fl[*,i] /= xsi
		plot,extracted_wl[*,i],extracted_fl[*,i]
		
		;stop
		
	endfor
	
	orders = indgen(norders)
	good_orders = where(good_orders eq 1, n_good_orders)
	final_wl = [];dblarr(n_good_orders * 2048)
	final_fl = [];dblarr(n_good_orders * 2048)
	;now normalize and put into final wl and fl arrays
	for i=0, n_good_orders-1 do begin
		order = orders[good_orders[i]]
		subwave = reform(extracted_wl[*,order])
		print,minmax(subwave)
		subflux = reform(extracted_fl[*,order])
		subflux_orig = subflux ;for debugging
		nel = n_elements(subwave)
		boxcarmax,subwave,subflux,nel/50.,wlsm,ms
		wlfit = ((subwave - min(subwave))/ (max(subwave) - min(subwave)))*2.d - 1.d
		wlfitin = wlfit[where(subflux eq subflux)]
		fit = svdfit(wlfitin,ms,5,/legendre,/double)
		fl_nvals = double(rleg(wlfit,fit))
		subflux = subflux / fl_nvals
		if i eq 0 then begin
			final_wl = subwave
			final_fl = subflux
		endif else begin
			mean_pixel_size = mean(subwave[1:*]-subwave)
			gap_size = min(subwave) - max(extracted_wl[*,order-1])
			n_gap_pixels = floor(gap_size / mean_pixel_size)
			if n_gap_pixels lt 0 then begin
				gls = where(subwave gt max(final_wl))
				subwave = subwave[gls]
				subflux = subflux[gls]
				final_wl = [final_wl,subwave]
				final_fl = [final_fl,subflux]
			endif else begin
				if floor(n_gap_pixels) eq 0 then begin
					final_wl = [final_wl,subwave]
					final_fl = [final_fl,subflux]
				endif else begin
					gap_wl = dindgen(n_gap_pixels)/double(n_gap_pixels)*gap_size + max(extracted_wl[*,order-1])
					gap_fl = make_array(n_gap_pixels,/double,value=!values.f_nan)
					final_wl = [final_wl,gap_wl,subwave]
					final_fl = [final_fl,gap_fl,subflux]
				endelse
			endelse
		endelse
		

	endfor	
	order = sort(final_wl)
	final_wl = final_wl[order]
	final_fl = final_fl[order]
	final_wl /= 1d3
	mwrfits, [1d#final_wl, 1d#final_fl], outfile, hdr,/create
	diag_output = "test"
END          

	
	
	
	
	
;;	for i=0, norders-1 do begin
;;		y1s = floor(ys_base[*,i] - 0.3 * y_size_kernel + 1024)
;;		;y2s = floor(ys_base[*,i] + 0.45 * y_size_kernel + 1024)
;;		y2s = y1s + .6 * y_size_kernel
;;		if max(y2s) le 2047 and min(y1s) ge 0 then good_orders[i] = 1
;;		for j=buffer/2, 2047+buffer/2 do begin ; loop over x pixels (recall that the overall array is "big")
;;			;x i want is j
;;			;y limits are y1/2s
;;			extracted_wl[j-buffer/2,i] = ws_base[j,i]
;;			
;;			if y2s[j] lt 0 or y1s[j] gt 2047 then continue
;;			y1t = 0 > y1s[j]
;;			y2t = 2047 < y2s[j]
;;			
;;			column = image[j-buffer/2,y1t:y2t]
;;			if n_elements(column) lt 3 then begin
;;				temptot = total(column) 
;;				;if temptot ne 0 then stop
;;			endif else begin
;;				colmaxs = locmax(column,3)
;;				maxlocs = where(colmaxs ne 0,nmaxs)
;;				;if nmaxs ne nfibers then stop
;;				temptot = 0d
;;				for k=0, nmaxs-1 do begin
;;					k1 = 0 > (maxlocs[k] - extract_width)
;;					k2 = (maxlocs[k] + extract_width) < (n_elements(column) - 1)
;;					temptot += total(column[k1:k2])
;;					;stop
;;				endfor
;;				;if temptot ne 0 then stop
;;			endelse
;;			extracted_fl[j-buffer/2,i] = temptot
;;			
;;		endfor
;;		wls = reform(extracted_wl[*,i])
;;		mp = (wls + wls[1:*])/2.
;;		xsi = mp[1:*] - mp
;;		xsi = [xsi[0],xsi,xsi[-1]]
;;
;;		extracted_fl[*,i] /= xsi
;;		
;;;;		display,image
;;;;		oplot,y1s,co=cgcolor('green')
;;;;		oplot,y2s,co=cgcolor('green')
;;;;		oplot,ys_base[*,i] + 1024,co=cgcolor('red')
;;;;		stop
;;	endfor
;;	
;;	;;;;;; temporary
;;	restore,inwlfile
;;	extracted_wl = warray
;;	;stop
;;	;;;;;;;;;;;;;;;;
;;	
;;	orders = indgen(norders)
;;	good_orders = where(good_orders eq 1, n_good_orders)
;;	final_wl = [];dblarr(n_good_orders * 2048)
;;	final_fl = [];dblarr(n_good_orders * 2048)
;;	;now normalize and put into final wl and fl arrays
;;	for i=0, n_good_orders-1 do begin
;;		order = orders[good_orders[i]]
;;		subwave = reform(extracted_wl[*,order])
;;		print,minmax(subwave)
;;		subflux = reform(extracted_fl[*,order])
;;		subflux_orig = subflux ;for debugging
;;		nel = n_elements(subwave)
;;		boxcarmax,subwave,subflux,nel/50.,wlsm,ms
;;		wlfit = ((subwave - min(subwave))/ (max(subwave) - min(subwave)))*2.d - 1.d
;;		wlfitin = wlfit[where(subflux eq subflux)]
;;		fit = svdfit(wlfitin,ms,5,/legendre,/double)
;;		fl_nvals = double(rleg(wlfit,fit))
;;		subflux = subflux / fl_nvals
;;		if i eq 0 then begin
;;			final_wl = subwave
;;			final_fl = subflux
;;		endif else begin
;;			mean_pixel_size = mean(subwave[1:*]-subwave)
;;			gap_size = min(subwave) - max(extracted_wl[*,order-1])
;;			n_gap_pixels = floor(gap_size / mean_pixel_size)
;;			if n_gap_pixels lt 0 then begin
;;				gls = where(subwave gt max(final_wl))
;;				subwave = subwave[gls]
;;				subflux = subflux[gls]
;;				final_wl = [final_wl,subwave]
;;				final_fl = [final_fl,subflux]
;;			endif else begin
;;				if floor(n_gap_pixels) eq 0 then begin
;;					final_wl = [final_wl,subwave]
;;					final_fl = [final_fl,subflux]
;;				endif else begin
;;					gap_wl = dindgen(n_gap_pixels)/double(n_gap_pixels)*gap_size + max(extracted_wl[*,order-1])
;;					gap_fl = make_array(n_gap_pixels,/double,value=!values.f_nan)
;;					final_wl = [final_wl,gap_wl,subwave]
;;					final_fl = [final_fl,gap_fl,subflux]
;;				endelse
;;			endelse
;;		endelse
;;		
;;	;	stop
;;	endfor	
;;	order = sort(final_wl)
;;	final_wl = final_wl[order]
;;	final_fl = final_fl[order]
;;	final_wl /= 1d3
;;	mwrfits, [1d#final_wl, 1d#final_fl], outfile, hdr,/create
;;	diag_output = "test"
;;END          
