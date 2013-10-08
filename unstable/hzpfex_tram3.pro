;; hzpfex_tram3.pro (and nonan.pro)
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




PRO hzpfex_tram3, infile, inwlfile, outfile, varfile = varfile, tellfile = tellfile, crimg = crimg, diag_output = diag_output, orders_lambdalow = orders_lambdalow, orders_lambdahigh = orders_lambdahigh, orders_gaps = orders_gaps, fiber_scale = fiber_scale, fiber_core_um = fiber_core_um, fiber_cladding_um = fiber_cladding_um, fiber_buffer_um = fiber_buffer_um,nfibers = nfibers, fiber_extra_sep_um = fiber_extra_sep_um, optical_model = optical_model, orders_shape = orders_shape

	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;GET ECHELLOGRAM INFO (copied from tram_projection3)
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;	 ;set file defaults
;;	 if n_elements(echellogram_wl_file) eq 0 then echellogram_wl_file = 'hpf demag=2-0x f8-5 2012dec15 v10-1-wavelengths.dat'
;;	 if n_elements(echellogram_file) eq 0 then echellogram_file = 'hpf demag=2-0x f8-5 2012dec15 v10-1-echelleogram.dat'
;;	
;;	 pixel_size = 18d-3 ;mm
;;	 ;read in files
;;	 inp = dblarr(23,29)
;;	 openr,1,echellogram_wl_file
;;	 readf,1,inp
;;	 close,1
;;
;;	 inp2 = dblarr(42,29)
;;	 openr,1,echellogram_file
;;	 readf,1,inp2
;;	 close,1
;;	
;;	 ;process according to Barnes' format
;;	 evens = indgen(21) * 2
;;	 odds = indgen(21) * 2 + 1
;;	 ys = inp2[odds,*]
;;	 xs = inp2[evens,*]
;;	 ws = inp[2:*,*]
;;	
;;	 ;if requested, straighten out the orders (debugging tool to examine the effects of curvature)
;;	 if keyword_set(straight_orders) then begin
;;		 mod_echel, xs, ys, ws, xs1, ys1, ws1
;;		 xs_old = xs
;;		 ys_old = ys
;;		 ws_old = ws
;;		 xs = xs1
;;		 ys = ys1
;;		 ws = ws1
;;	 endif	
;;	
;;	 ys = -1d * ys ;flip y around (stuart does this in his code)
;;	
;;	 norders = (size(inp2))[2]

	pixel_size = 18d-3 ;mm

	case optical_model of
		'barnes_1212': optical_model_file = 'model_barnes_1212.sav'
		'ramsey_513': optical_model_file = 'model_ramsey_513.sav'
		else: optical_model_file = 'model_ramsey_513.sav'
	endcase
	
	restore,optical_model_file
	
	xs = model.xs
	ys = model.ys
	ws = model.ws
	
	if keyword_set(orders_shape) then begin
		mod_echel, xs, ys, ws, xs1, ys1, ws1, ver=orders_shape
		xs_old = xs
		ys_old = ys
		ws_old = ws
		xs = xs1
		ys = ys1
		ws = ws1
	endif	

	
	norders = n_elements(model.orders)

	;make the warray, although this is not *required* for extraction now	
	warray=MAKE_ARRAY(2048, norders, /Double, Value=!Values.F_NAN)
	
	
	;xs,ys are the xy locations of the order sample points in mm
	;xs,ys_pix are the xy locations in pixels (with -1024 to 1023, not 0 to 2047)
	xs_pix = xs / pixel_size
	ys_pix = ys / pixel_size



	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;READ IN THE IMAGE
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	image=MRDFITS(infile, 0, hdr)


	extracted_fl = dblarr(2048,norders)
	extracted_wl = dblarr(2048,norders)

	
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;PROCESS ECHELLOGRAM INFO (copied from tram_projection3)
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;_base arrays are interpolated versions of the Barnes input to match what is required
	;for the array
	
	buffer = 100
	upfactor = 1
	fs_base = dblarr((2048 + buffer)*upfactor,norders)
	fs_base2 = fs_base
	
	xs_base = ((dindgen((2048+buffer)*upfactor) - (1024.+buffer/2)*upfactor)/upfactor)
	ws_base = dblarr((2048 + buffer)*upfactor,norders)
	xs_base_nobuf = xs_base[buffer/2 : 2047+buffer/2]
	
	;find left and right limits of each pixel for binning the spectrum
	xs_base_left = xs_base - 0.5
	xs_base_right = xs_base + 0.5
	ws_base_left = dindgen((2048 + buffer) * upfactor,norders)
	ws_base_right = dindgen((2048 + buffer) * upfactor,norders)
	
	ys_base = dblarr((2048 + buffer) * upfactor,norders)
	xs_base_2d = dblarr((2048 + buffer) * upfactor,norders)

	;fill in the wavelengths that correspond to each bin and the edge of each bin
	;note that wavelength depends only on x-position
	for i=0, norders-1 do begin
		w1a = mpfitfun('poly',xs_pix[*,i],ws[*,i],1d,[0d,0d,0d],yfit=w1b,/quiet)
		w1c = poly(xs_base,w1a)
		ws_base[*,i] = w1c
		dw1 = (ws_base[1:*,i] - ws_base[*,i])/2d
		dex = interpol(dw1,xs_base[0:-2],xs_base[-1])
		dw2 = [dw1,dex]
		ws_base_right[*,i] = ws_base[*,i] + dw2
		ws_base_left[*,i] = ws_base[*,i] - dw2
		
		y1a = mpfitfun('poly',xs_pix[*,i],ys_pix[*,i],1d,[0d,0d,0d],yfit=y1b,/quiet)
		y1c = poly(xs_base,y1a)
		ys_base[*,i] = y1c
		;replicate the x base array for convolving and plotting
		xs_base_2d[*,i] = xs_base
	endfor
	
	ys_base_nobuf = ys_base[buffer/2 : 2047 + buffer/2,*]




	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;GET SOME KERNEL INFO
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;create a kernel based on the fiber inputs
	;assusme a row of 2d gaussians for now
	;kernel size in pixels
	kernel_size_pix = fiber_core_um * fiber_scale
	fiber_extra_sep_pix = fiber_extra_sep_um * fiber_scale
	
	;separation between gaussians
	kernel_sep_pix = fiber_buffer_um * fiber_scale + fiber_extra_sep_pix
	
	;number of gaussians = nfibers
	y_size_kernel = (floor((nfibers + 1d) * kernel_sep_pix))*upfactor; + kernel_size_pix) + kernel_sep_pix)
	x_size_kernel = (floor(3. * kernel_size_pix))*upfactor
	extract_width = round(kernel_size_pix/2d + 1)

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;EXTRACT THE SPECTRA
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;keep track of which orders are entirely contained on the chip and extract only these
	good_orders = intarr(norders)

	for i=0, norders-1 do begin

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;RECREATE THE ARRAYS USED FOR THE FLUX ARRAY WARPING (copied from tram_projection3)
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;create the rectified arrays
		;these map the rectified fs to the warped positions xs_rect, ys_rect
		xs_rect = dblarr((2048+buffer)*upfactor,200*upfactor)
		ys_rect = dblarr((2048+buffer)*upfactor,200*upfactor)
		fs_rect = dblarr((2048+buffer)*upfactor,200*upfactor)
;;		xs_rect = dblarr((2048)*upfactor,200*upfactor)
;;		ys_rect = dblarr((2048)*upfactor,200*upfactor)
;;		fs_rect = dblarr((2048)*upfactor,200*upfactor)
		
		for j=0, 200*upfactor-1 do begin
			xs_rect[*,j] = xs_base
			ypos = double(j)/double(upfactor) - 100d
			ys_rect[*,j] = ys_base[*,i] + ypos
		endfor


		;create indices for the rect arrays
	
		nx_rect = (size(fs_rect))[1] ;- buffer
		ny_rect = (size(fs_rect))[2]
		x_grid_rect_d = rebin(lindgen(nx_rect),nx_rect,ny_rect)
		y_grid_rect_d = rebin(1#lindgen(ny_rect),nx_rect,ny_rect)
	
		;find offset between indices (0,2048) and coords (-1024,1023)
		xs_d_offset = min(xs_rect) ;+buffer/2;- buffer/2
		ys_d_offset = min(ys_rect) ;+buffer/2
	
		;reposition the rect arrays to (0,2047)
		xs_rect_d = (xs_rect - xs_d_offset) * upfactor 
		ys_rect_d = (ys_rect - ys_d_offset) * upfactor 
	
		;reposition the grid to (-1024,1023)
		x_grid_rect = x_grid_rect_d + xs_d_offset ;- buffer/2
		y_grid_rect = y_grid_rect_d + ys_d_offset
	
		;how big must the output be
		x_size_d = max(floor(xs_rect_d)) - min(floor(xs_rect_d))
		y_size_d = max(floor(ys_rect_d)) - min(floor(ys_rect_d))
		if x_size_d mod 2 eq 1 then x_size_d += 1
		if y_size_d mod 2 eq 1 then y_size_d += 1
		
		;x_size_d = 2048d
		;y_size_d = 200d
	
		y_min = min(floor(ys_rect))
		y_max = max(floor(ys_rect))

		;find the chunk of the detector that this corresponds to
		;y locations on detector
		y0_det = min(y_grid_rect) + 1024; + buffer/2
		y1_det = y0_det + y_size_d/upfactor - 1
		;y locations in the sub-array
		y0_subarr = 0
		y1_subarr = y_size_d/upfactor - 1
		;x locations in sub-array and detector always the same
		x0 = buffer/2
		x1 = x0 + 2047
		xd0 = 0
		xd1 = 2047
	
		;if there is any that spills off the end of the detector, adjust the two coordiantes
		extra = 0d
		y1 = 100d - 0.3 * y_size_kernel 
		y2 = 100d + 0.3 * y_size_kernel 

		if y0_det lt 0 then begin
			extra = abs(y0_det)
			y0_det = 0
			y0_subarr = extra
			y1 = 100d - 0.3 * y_size_kernel - extra
			y2 = 100d + 0.3 * y_size_kernel - extra
		endif
		if y1_det gt 2047 then begin
			;extra = floor(max(y_grid_rect)+1024) - 2047
			extra = floor(y1_det - 2047)
			y1_det = 2047
			y1_subarr = y_size_d/upfactor - 1 - extra/upfactor
		endif


		y1_check = min(ys_base[*,i]) - y_size_kernel/2d - 1d + 1023d
		y2_check = max(ys_base[*,i]) + y_size_kernel/2d + 1d + 1023d
		if y1_check lt 0 or y2_check gt 2047 then begin
			good_orders[i] = 0
			continue
		endif
		;if nothing is on the detector skip this part
		if y1_det lt 0 or y0_det gt 2047 then begin
			good_orders[i] = 0
			continue
		endif
		good_orders[i] = 1
	
		;fill in
;;		t1 = fs_det[xd0:xd1 , y0_det:y1_det]
;;		t2 = t1 + tmp_conv_rebin[x0:x1, y0_subarr:y1_subarr]
;;		fs_det[xd0:xd1 , y0_det:y1_det] = t2
		
		tmp = dblarr(2148,(floor(y1_det) - floor(y0_det))+1)
		tmp[buffer/2:2047+buffer/2,*] = image[xd0:xd1, y0_det:y1_det]
		;tmp = image[xd0:xd1, y0_det:y1_det]
;;		x_grid_rect_d_nobuf = x_grid_rect_d[buffer/2:2047+buffer/2, y0_subarr:y1_subarr]
;;		y_grid_rect_d_nobuf = y_grid_rect_d[buffer/2:2047+buffer/2, y0_subarr:y1_subarr]
;;		xs_rect_d_nobuf = xs_rect_d[buffer/2:2047+buffer/2, y0_subarr:y1_subarr]
;;		ys_rect_d_nobuf = ys_rect_d[buffer/2:2047+buffer/2, y0_subarr:y1_subarr]
		;xs_rect_d -= 50.
		
		
		tmp_warp = warp_tri(x_grid_rect_d,y_grid_rect_d,xs_rect_d,ys_rect_d,tmp)
		;remove below and uncomment above
		;tmp_warp = tmp
		
		;tmp_warp = warp_tri(x_grid_rect_d_nobuf,y_grid_rect_d_nobuf,xs_rect_d_nobuf,ys_rect_d_nobuf,tmp)
		
		tmp_warp_sub = tmp_warp[buffer/2:2047+buffer/2,*]
	
		;stop

		;y1s and y2s are lower and upper limits of extraction regions
		;y1s = floor(ys_base[*,i] - 0.3 * y_size_kernel + 1024)
		;y2s = y1s + .6 * y_size_kernel
		
		;stop
		
		;if entire region is on the array then it is good
		;if max(y2s) le 2047 and min(y1s) ge 0 then good_orders[i] = 1
		
		;iterate over pixels in the extracted array
		;extracted_wl is created with buffer limits so only fill in the non-buffer region
		;;;;;; temporary
		restore,inwlfile
		extracted_wl = warray ;* 1d3
		;stop
		;;;;;;;;;;;;;;;;

		for j=buffer/2, 2047+buffer/2 do begin ; loop over x pixels (recall that the overall array is "big")
			;x i want is j
			;y limits are y1/2s
			;extracted_wl[j-buffer/2,i] = ws_base[j,i]  ;change is this back ryan 7-22-13
			
			;if there's nothing then continue, otherwise pull out what i can
			;if y2s[j] lt 0 or y1s[j] gt 2047 then continue
			;y1t = 0 > y1s[j]
			;y2t = 2047 < y2s[j]
			
			;pull the entire column
			column = tmp_warp_sub[j-buffer/2,y1:y2]
			
			;if there's  only a few then just take the total (this doesn't matter since these are not "good" orders)
			if n_elements(column) lt 3 then begin
				temptot = total(column) 
				;if temptot ne 0 then stop
			endif else begin
			
				;find the max locations
				colmaxs = locmax(column,3)
				maxlocs = where(colmaxs ne 0,nmaxs)
				temptot = 0d
				
				;pull out the flux within extract_width of the max locations
				for k=0, nmaxs-1 do begin
					k1 = 0 > (maxlocs[k] - extract_width)
					k2 = (maxlocs[k] + extract_width) < (n_elements(column) - 1)
					temptot += total(column[k1:k2])
					;stop
				endfor
			endelse
			;store
			extracted_fl[j-buffer/2,i] = temptot
			;stop
		endfor
		
		;convert to flux/wavelength
		wls = reform(extracted_wl[*,i])
		mp = (wls + wls[1:*])/2.
		xsi = mp[1:*] - mp
		xsi = [xsi[0],xsi,xsi[-1]]

		extracted_fl[*,i] /= xsi
		;stop
;;		display,image
;;		oplot,y1s,co=cgcolor('green')
;;		oplot,y2s,co=cgcolor('green')
;;		oplot,ys_base[*,i] + 1024,co=cgcolor('red')
;;		stop
	endfor
	
	
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;NORMALIZE THE ORDERS AND FILL IN TO FINAL ARRAY
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	orders = indgen(norders)
	good_orders = where(good_orders eq 1, n_good_orders)
	
	;now normalize and put into final wl and fl arrays
	for i=0, n_good_orders-1 do begin
		order = orders[good_orders[i]]
		
		;wavelengths for this order
		subwave = reform(extracted_wl[*,order])
		subflux = reform(extracted_fl[*,order])
		subflux_orig = subflux ;for debugging
		;diag
		print,minmax(subwave)
		
		nel = n_elements(subwave)
		
		;find the sliding maximums
		boxcarmax,subwave,subflux,nel/50.,wlsm,ms
		
		;normalize using a legendre polynomial
		wlfit = ((subwave - min(subwave))/ (max(subwave) - min(subwave)))*2.d - 1.d ; wl array remapped to -1,1
		wlfitin = wlfit[where(subflux eq subflux)] ; find finite elements
		fit = svdfit(wlfitin,ms,5,/legendre,/double)
		fl_nvals = double(rleg(wlfit,fit))
		subflux = subflux / fl_nvals
		
		;if we're on the first order create the final arrays, otherwise append to them and fill in gaps
		if i eq 0 then begin
			final_wl = subwave
			final_fl = subflux
		endif else begin
			;figure out pixel size for filling the gap
			mean_pixel_size = mean(subwave[1:*]-subwave)
			
			;how many pixels needed to fill gap?
			gap_size = min(subwave) - max(extracted_wl[*,order-1])
			n_gap_pixels = floor(gap_size / mean_pixel_size)
			
			;fill them in
			if n_gap_pixels lt 0 then begin
				;if there's overlap, then take the non-redundant part of this order and fill in
				gls = where(subwave gt max(final_wl))
				subwave = subwave[gls]
				subflux = subflux[gls]
				final_wl = [final_wl,subwave]
				final_fl = [final_fl,subflux]
			endif else begin
				;if there's no overlap
				if floor(n_gap_pixels) eq 0 then begin
					;if there's no need for filler pixels
					final_wl = [final_wl,subwave]
					final_fl = [final_fl,subflux]
				endif else begin
					;if there is need for filler pixels
					gap_wl = dindgen(n_gap_pixels)/double(n_gap_pixels)*gap_size + max(extracted_wl[*,order-1])
					gap_fl = make_array(n_gap_pixels,/double,value=!values.f_nan)
					final_wl = [final_wl,gap_wl,subwave]
					final_fl = [final_fl,gap_fl,subflux]
				endelse
			endelse
		endelse

	endfor	
	
	;stop
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;STORE RESULT
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	order = sort(final_wl)
	final_wl = final_wl[order]
	final_fl = final_fl[order]
	final_wl /= 1d3
	mwrfits, [1d#final_wl, 1d#final_fl], outfile, hdr,/create
	diag_output = "test"
END          
