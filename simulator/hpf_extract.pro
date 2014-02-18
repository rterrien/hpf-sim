;+
; NAME:
;  hpf_extract
;
; PURPOSE:
;
;  Extract the spectrum from a simulated detector exposure
;
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  hpf_extract, spec_params, optical_params, proj_params, det_params, infile, outfile
;
; INPUTS:
;
;	spec_params: structure defined in initialize_spec_params
;
;	optical_params: structure defined in initialize_optical_params
;
;	proj_params: structure defined by process_optical_model
;
;	det_params: structure defined by initialize_det_params
;
;	infile: fits file containing the exposure
;
; OUTPUTS:
;	
;	outfile: output fits file to contain extracted spectra
;
; KEYWORD PARAMETERS:
;
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-



pro hpf_extract, spec_params, optical_params, proj_params, det_params, infile, outfile

	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;GET ECHELLOGRAM INFO
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	xs = proj_params.xs
	ys = proj_params.ys
	ws = proj_params.ws
	xs_pix = proj_params.xs_pix
	ys_pix = proj_params.ys_pix
	fs_base = proj_params.fs_base_noup
	xs_base = proj_params.xs_base_noup
	ws_base = proj_params.ws_base_noup
	xs_base_left = proj_params.xs_base_left_noup
	xs_base_right = proj_params.xs_base_right_noup
	ws_base_left = proj_params.ws_base_left_noup
	ws_base_right = proj_params.ws_base_right_noup
	ys_base = proj_params.ys_base_noup
	xs_base_2d = proj_params.xs_base_2d_noup
	
	norders = (size(proj_params.xs))[2]
	nspec = n_elements(spec_params.spec_file)

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; FIBER and KERNEL PARAMETERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	fiber_core_um = optical_params.fiber_core_um
	fiber_buffer_um = optical_params.fiber_buffer_um
	fiber_scale = optical_params.fiber_scale
	fiber_extra_sep_um = optical_params.fiber_extra_sep_um

	;kernel size in pixels
	kernel_size_pix = fiber_core_um * fiber_scale
	fiber_extra_sep_pix = fiber_extra_sep_um * fiber_scale
	
	;separation between gaussians
	kernel_sep_pix = fiber_buffer_um * fiber_scale + fiber_extra_sep_pix
	y_size_kernel = fiber_buffer_um * fiber_scale + 2d*fiber_extra_sep_pix
	x_size_kernel = (floor(3. * kernel_size_pix))		

	extract_width = round(kernel_size_pix/2d + 1) ;+ 2d ;CHANGE THIS +2d back to 0 when done rt 10-9-13

	buffer = 100.
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;READ IN THE IMAGE
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	image=MRDFITS(infile, 0, hdr)


	extracted_fl = dblarr(2048,norders,nspec)
	extracted_wl = dblarr(2048,norders)

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;EXTRACT THE SPECTRA
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	extracted_wl = mrdfits(optical_params.wlimg_file)
	;extracted_wl = warray ;* 1d3

	;keep track of which orders are entirely contained on the chip and extract only these
	good_orders = intarr(norders)
	
	if nspec ne 3 then message,'CHECK TO MAKE SURE EXTRACTION WINDOW IS POSITIONED PROPERLY ON RECTIFIED ORDER'

	for i=0, norders-1 do begin
	
		print,i
		boxsi = floor(max(ys_base[*,i] + nspec * y_size_kernel/2d) - min(ys_base[*,i] - nspec * y_size_kernel/2d)) + 12.
		
		;find the chunk of the detector that this corresponds to
		;y locations on detector
		mid_offset = 0.
		case optical_params.orders_shape of
		0: offset = 12.
		1: offset = 2.
		2: offset = 8.
		3: offset = 8.
		endcase
		; with new linear model, there is some changing offset that creeps in, 
		; the following clumsily corrects for it RCT 2-13-14
		if optical_params.model_file eq 'support/model_020714_shift.sav' then begin
			;boxsi = 400.
			if i le 10 then offset = 60. - i * 2.
			if i gt 10 and i lt 15 then offset = 60 - i*1.75
			if i ge 15 and i lt 20 then offset = 33 - (15-i) *1.5
			if i eq 20 then offset = 23
			if i ge 20 then offset = 28 ;- (i-20.) * 1.
			;mid_offset = 58.
		endif
		
		xs_rect = dblarr(2048 + buffer, boxsi)
		ys_rect = dblarr(2048 + buffer, boxsi)
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;RECREATE THE ARRAYS USED FOR THE FLUX ARRAY WARPING
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		for j=0, boxsi-1 do begin
			xs_rect[*,j] = xs_base
			ypos = double(j) - boxsi/2.
			ys_rect[*,j] = ys_base[*,i] + ypos
		endfor

		;find offset between indices (0,2048) and coords (-1024,1023)
		xs_d_offset = min(xs_rect) ;+buffer/2;- buffer/2
		ys_d_offset = min(ys_rect) ;+buffer/2
	
		;reposition the rect arrays to (0,2047)
		xs_rect_d = (xs_rect - xs_d_offset)
		ys_rect_d = (ys_rect - ys_d_offset)
	
		y_min = min(floor(ys_rect))
		y_max = max(floor(ys_rect))

		y0_det = min(ys_rect) + 1024 + offset; + buffer/2
		y1_det = y0_det + boxsi - 1
		x0 = buffer/2
		x1 = x0 + 2047
		xd0 = 0
		xd1 = 2047
		
		if y0_det lt 0 or y0_det gt 2047 then begin
			good_orders[i] = 0
			continue
		endif
		if y1_det lt 0 or y1_det gt 2047 then begin
			good_orders[i] = 0
			continue
		endif

		good_orders[i] = 1
	
		;fill in
		
		tmp = dblarr(2148,(floor(y1_det) - floor(y0_det))+1)
		tmp[buffer/2:2047+buffer/2,*] = image[xd0:xd1, y0_det:y1_det]
		
		nn = y1_det - y0_det + 1
		
		
		gx = 1d
		gy = 1d
		tmp_warp = tri_surf(tmp,xs_rect_d,reverse(ys_rect_d,2),gs=[gx,gy],missing=0d)
		tmp_warp = reverse(tmp_warp,2)
		mid = (size(tmp_warp,/dimen))[1]/2 - 1 + mid_offset
		y1 = mid - .4 * y_size_kernel
		y2 = mid + .4 * y_size_kernel
				
		tmp_warp_sub = tmp_warp[buffer/2:2047+buffer/2,*]
		tmp_warp_sub_maxs = tmp_warp_sub
		
		yy1 = mid - (nspec/2d) * y_size_kernel
		yy2 = mid + (nspec/2d) * y_size_kernel

		tmp_warp_sub_maxs[*,0:yy1] = 0d
		tmp_warp_sub_maxs[*,yy2:*] = 0d
		
		median_tmp = median(tmp_warp_sub_maxs[where(tmp_warp_sub_maxs gt 0.)],/double)
		horiz_tot = total(tmp_warp_sub_maxs,1,/double)
		horiz_norm = total(tmp_warp_sub_maxs gt median_tmp/2d,1,/double)
		zl = where(horiz_norm eq 0,nzl)
		if nzl ne 0 then horiz_norm[zl] = 1d
		;low_locs = where(horiz_norm lt max(horiz_norm)/1d4,nlow_locs)
		;if nlow_locs ne 0 then horiz_norm[low_locs] = 0d
		horiz_collapse = horiz_tot * horiz_norm
		;horiz_collapse[0:yy1] = 0d
		;horiz_collapse[yy2:*] = 0d
		horiz_maxs = hpf_locmax(horiz_collapse,3)
		;filter maxs to exclude noise peaks, use median_image
		low_locs = where(horiz_maxs lt max(horiz_collapse)*1d-4,nlow_locs)
		if nlow_locs ne 0 then horiz_maxs[low_locs] = 0d
		maxlocs = where(horiz_maxs ne 0,nmaxs)
		
		;This deals with a few special cases:
		;If there are too many spurious maxima, take the highest:
		if nmaxs gt optical_params.nfibers then begin
			rollmax = horiz_maxs
			maxlocs = []
			for j=0, optical_params.nfibers-1 do begin
				ll = (where(rollmax eq max(rollmax)))[0]
				rollmax[ll] = 0
				maxlocs = [maxlocs,ll]
			endfor
			maxlocs = maxlocs[sort(maxlocs)]
			ds = maxlocs[1:*] - maxlocs
			if total(ds lt extract_width) ne 0 then message,'finding maxima too close together'
			nmaxs = optical_params.nfibers
		endif
		;If the middle fiber is off in a 3 fiber config, take the average position of the end fibers
		if nmaxs eq 2 and optical_params.nfibers eq 3 then begin
			if (maxlocs[1] - maxlocs[0]) gt 2*extract_width then begin
				mm = floor(mean(maxlocs[0]))
				maxlocs = [maxlocs[0],mm,maxlocs[1]]
			endif else begin
				message,'still bad'
			endelse
		endif else begin
			if nmaxs ne optical_params.nfibers then message,'found different nmaxs and nfibers'
		endelse
		;if i mod 5 eq 0 and i ge 20 then stop
		;display,alog(tmp_warp+1.)
		;hline,mid,color=cgcolor('red')
		;hline,[yy1,yy2],color=cgcolor('blue')
		;wait,1
		for j=buffer/2, 2047+buffer/2 do begin ; loop over x pixels (recall that the overall array is "big")
			;x i want is j
			;y limits are y1/2s
			;extracted_wl[j-buffer/2,i] = ws_base[j,i]  ;change is this back ryan 7-22-13
			
			;if there's nothing then continue, otherwise pull out what i can
			;if y2s[j] lt 0 or y1s[j] gt 2047 then continue
			;y1t = 0 > y1s[j]
			;y2t = 2047 < y2s[j]
			
			;pull the entire column
			;column = tmp_warp_sub[j-buffer/2,y1:y2]
			column = reform(tmp_warp_sub[j-buffer/2,*])
			
			
			for k=0, nmaxs-1 do begin
				k1 = 0 > (maxlocs[k] - extract_width)
				k2 = (maxlocs[k] + extract_width) < (n_elements(column) - 1)
				extracted_fl[j-buffer/2,i,k] = total(column[k1:k2],/double)
				
			endfor
		endfor
		;plot,column
		;vline,maxlocs,color=cgcolor('red')
		;vline,[maxlocs-extract_width,maxlocs+extract_width]
		;wait,1.

		;stop
		
		;convert to flux/wavelength
		wls = reform(extracted_wl[*,i])
		mp = (wls + wls[1:*])/2.
		xsi = mp[1:*] - mp
		xsi = [xsi[0],xsi,xsi[-1]] 
		for j=0, nspec-1 do begin
			if spec_params.output_per_wavelength[j] then extracted_fl[*,i,j] /= xsi
		endfor
		
		;if output_per_wavelength then begin
		;	for j=0, nspec-1 do extracted_fl[*,i,j] /= xsi
		;endif
		;stop
;;		display,image
;;		oplot,y1s,co=cgcolor('green')
;;		oplot,y2s,co=cgcolor('green')
;;		oplot,ys_base[*,i] + 1024,co=cgcolor('red')
;;		stop

	endfor
	
	
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;NORMALIZE THE STELLAR ORDERS AND FILL IN TO FINAL ARRAY
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	orders = indgen(norders)
	good_orders = where(good_orders eq 1, n_good_orders)
	
	final_fl_all = ptrarr(nspec,/allocate_heap)
	
	;need to reverse these to have them strung together properly?
	if optical_params.model_file eq 'support/model_020714_shift.sav' then begin
		good_orders = reverse(good_orders)
		extracted_wl = reverse(extracted_wl,2)
		extracted_fl = reverse(extracted_fl,2)
	endif
	
	for j=0, nspec-1 do begin
		;now normalize and put into final wl and fl arrays
		for i=0, n_good_orders-1 do begin
			order = orders[good_orders[i]]
			
			;wavelengths for this order
			subwave = reform(extracted_wl[*,order])
			subflux = reform(extracted_fl[*,order,j])
			subflux_orig = subflux ;for debugging
			;diag
			print,minmax(subwave)
			
			nel = n_elements(subwave)
		
		
			if spec_params.normalize_output[j] then begin
				;find the sliding maximums
				hpf_boxcarmax,subwave,subflux,nel/50.,wlsm,ms
	   
				;normalize using a legendre polynomial
				wlfit = ((subwave - min(subwave))/ (max(subwave) - min(subwave)))*2.d - 1.d ; wl array remapped to -1,1
				wlfitin = wlfit[where(subflux eq subflux)] ; find finite elements
				fit = svdfit(wlfitin,ms,5,/legendre,/double)
				fl_nvals = double(hpf_rleg(wlfit,fit))
				subflux = subflux / fl_nvals
			endif

		
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
		*final_fl_all[j] = final_fl	
	endfor	
	;stop
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;STORE RESULT
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	order = sort(final_wl)
	final_wl = final_wl[order]
	final = dblarr(nspec+1,n_elements(final_wl))
	final[0,*] = reform(final_wl)
	for i=0, nspec-1 do final[i+1,*] = reform((*final_fl_all[i])[order])
	mwrfits,final,outfile,hdr,/create
	

END          
