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




PRO hzpfex_tram2, infile, inwlfile, outfile, varfile = varfile, tellfile = tellfile, crimg = crimg, diag_output = diag_output, orders_lambdalow = orders_lambdalow, orders_lambdahigh = orders_lambdahigh, orders_gaps = orders_gaps, fiber_scale = fiber_scale, fiber_core_um = fiber_core_um, fiber_cladding_um = fiber_cladding_um, fiber_buffer_um = fiber_buffer_um,nfibers = nfibers, fiber_extra_sep_um = fiber_extra_sep_um, echellogram_file = echellogram_file, echellogram_wl_file = echellogram_wl_file, straight_orders = straight_orders

	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;GET ECHELLOGRAM INFO (copied from tram_projection3)
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;set file defaults
	if n_elements(echellogram_wl_file) eq 0 then echellogram_wl_file = 'hpf demag=2-0x f8-5 2012dec15 v10-1-wavelengths.dat'
	if n_elements(echellogram_file) eq 0 then echellogram_file = 'hpf demag=2-0x f8-5 2012dec15 v10-1-echelleogram.dat'
	
	pixel_size = 18d-3 ;mm
	;read in files
	inp = dblarr(23,29)
	openr,1,echellogram_wl_file
	readf,1,inp
	close,1

	inp2 = dblarr(42,29)
	openr,1,echellogram_file
	readf,1,inp2
	close,1
	
	;process according to Barnes' format
	evens = indgen(21) * 2
	odds = indgen(21) * 2 + 1
	ys = inp2[odds,*]
	xs = inp2[evens,*]
	ws = inp[2:*,*]
	
	;if requested, straighten out the orders (debugging tool to examine the effects of curvature)
	if keyword_set(straight_orders) then begin
		mod_echel, xs, ys, ws, xs1, ys1, ws1
		xs_old = xs
		ys_old = ys
		ws_old = ws
		xs = xs1
		ys = ys1
		ws = ws1
	endif	
	
	ys = -1d * ys ;flip y around (stuart does this in his code)
	
	norders = (size(inp2))[2]

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
		;y1s and y2s are lower and upper limits of extraction regions
		y1s = floor(ys_base[*,i] - 0.3 * y_size_kernel + 1024)
		y2s = y1s + .6 * y_size_kernel
		
		;if entire region is on the array then it is good
		if max(y2s) le 2047 and min(y1s) ge 0 then good_orders[i] = 1
		
		;iterate over pixels in the extracted array
		;extracted_wl is created with buffer limits so only fill in the non-buffer region
		for j=buffer/2, 2047+buffer/2 do begin ; loop over x pixels (recall that the overall array is "big")
			;x i want is j
			;y limits are y1/2s
			extracted_wl[j-buffer/2,i] = ws_base[j,i]
			
			;if there's nothing then continue, otherwise pull out what i can
			if y2s[j] lt 0 or y1s[j] gt 2047 then continue
			y1t = 0 > y1s[j]
			y2t = 2047 < y2s[j]
			
			;pull the entire column
			column = image[j-buffer/2,y1t:y2t]
			
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
				endfor
			endelse
			;store
			extracted_fl[j-buffer/2,i] = temptot
		endfor
		
		;convert to flux/wavelength
		wls = reform(extracted_wl[*,i])
		mp = (wls + wls[1:*])/2.
		xsi = mp[1:*] - mp
		xsi = [xsi[0],xsi,xsi[-1]]

		extracted_fl[*,i] /= xsi
		
;;		display,image
;;		oplot,y1s,co=cgcolor('green')
;;		oplot,y2s,co=cgcolor('green')
;;		oplot,ys_base[*,i] + 1024,co=cgcolor('red')
;;		stop
	endfor
	
	;;;;;; temporary
	restore,inwlfile
	extracted_wl = warray
	;stop
	;;;;;;;;;;;;;;;;
	
	
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
