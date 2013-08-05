;; tram_projection1
;; HPF DETECTOR SIMULATOR
;; This routine does the projection of fibers onto the detector using the tramline method
;; from Stuart Barnes
;; 
;; CALLS
;; 
;; 
;; DIRECTLY MODIFIES
;;
;; CALLED BY
;; HZPFSIM_IMG
;;
;; NOTES
;; 
;;
;; parameters:
;; w - spectrum wavelengths in microns
;; f - spectrum fluxes
;; res - resolution
;; pixel_sampling - number of pixels per resolution element
;; wlimg - image of wavelengths
;; specimg - image of fluxes
;; warray - wavelength array (for extraction)
;; fiber_size - fiber size in pixels
;; fiber_gap - gap between fibers in pixels
;; fiber_fractions - dblarr(7) with fiber power fractions
;; fiber_core/cladding/buffer - diameters of elements of fiber in um
;; nfibers - number of fibers to put into the kernel
;; fiber_extra_sep_um - extra separation between fibers
;; optical model - choice of 'barnes_1212' or 'ramsey_513' orders, wavelengths (nm), and x/y positions in mm
;; slitwidth_um - width of the slit
;; straight_orders - keyword, if set it will set all y positions for each order to the mean of that order
;; upfactor - factor to upsample the warping/convolving array


pro tram_projection3, w, f, res, pixel_sampling, specimg, calw=calw, calf=calf, diag_out = diag_out, fiber_fractions = fiber_fractions, warray = warray, fiber_scale = fiber_scale, fiber_core_um = fiber_core_um, fiber_cladding_um = fiber_cladding_um, fiber_buffer_um = fiber_buffer_um, nfibers = nfibers, fiber_extra_sep_um = fiber_extra_sep_um, optical_model = optical_model, slitwidth_um = slitwidth_um, straight_orders = straight_orders, upfactor = upfactor

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;PARAMETERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	buffer = 100 ; the extra pixels (/2) on each side for the big arrays which get filled in and convolved
	if n_elements(upfactor) eq 0 then upfactor = 1
	;upfactor = 2d ; the factor to upsample the array on which the warping/convolving is performed
	
	fs_det = dblarr(2048,2048) ; detector array

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;READ IN STELLAR SPECTRUM
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	input_wl = double(w * 1d3)
	input_fl = double(f)


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;GET ECHELLOGRAM INFO
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
;;	
;;	 ys = -1d * ys ;flip y around (stuart does this in his code)

	pixel_size = 18d-3 ;mm
	
	if n_elements(optical_model) eq 0 then optical_model = 'ramsey_513'
	case optical_model of
		'barnes_1212': optical_model_file = 'model_barnes_1212.sav'
		'ramsey_513': optical_model_file = 'model_ramsey_513.sav'
		else: optical_model_file = 'model_ramsey_513.sav'
	endcase
	
	restore,optical_model_file
	
	
	
	xs = model.xs
	ys = model.ys
	ws = model.ws
	
	if keyword_set(straight_orders) then begin
		mod_echel, xs, ys, ws, xs1, ys1, ws1
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
	;;PROCESS ECHELLOGRAM INFO
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;_base arrays are interpolated versions of the Barnes input to match what is required
	;for the array
	
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
	;;CREATE FLUX ARRAY
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	for i=0, norders-1 do begin
		;create highly sampled individual spectrum for integrating
		;have to upsample to get many pixels per output pixel
		;there is an option here to simply interpolate to the pixel values
		
		;define high sample wavelengths/pixel edges
		wl_order_highs = dindgen(1d6)/1d6 * (max(ws_base_right[*,i]) - min(ws_base_left[*,i]) + 10.) + min(ws_base_left[*,i]) - 5.
		gg = where(input_wl gt (min(ws_base_right[*,i]) - 1.) and input_wl lt (max(ws_base_left[*,i]) + 1.))
		;fl_order_highs = input_fl[gg]
		;wl_order_highs = input_wl[gg]
		;change above back 7-19-13

		
		midpts = (wl_order_highs[1:*] + wl_order_highs)/2.d
		xsis = midpts[1:*] - midpts
		xsis = [xsis[0],xsis,xsis[-1]]
		wl_order_highs_left = wl_order_highs - xsis/2d
		wl_order_highs_right = wl_order_highs + xsis/2d
		
		;are there any pixels outside of the input spectrum?
		outlocs = where((max(input_wl) lt wl_order_highs) or (min(input_wl) gt wl_order_highs),noutlocs)
		
		;interpolate to create highly sampled fluxes
		fl_order_highs = interpol(input_fl,input_wl,wl_order_highs);,/quadratic)

		;make sure the interpolation didn't extrapolate anywhere
		if noutlocs gt 0 then fl_order_highs[outlocs] = 0
		
		;find where each output pixel falls in the upsampled array and integrate that part
		left_inds = value_locate(wl_order_highs_left,ws_base_left[*,i]) + 1
		right_inds = value_locate(wl_order_highs_right,ws_base_right[*,i])

		
		
		for j=0, (size(fs_base))[1]-1 do begin
			;pixels entirely inside the output pixel
			totmid = total(fl_order_highs[left_inds[j]:right_inds[j]] * xsis[left_inds[j]:right_inds[j]],/double)
			;xx = wl_order_highs[left_inds[j]:right_inds[j]];xsis[left_inds[j]:right_inds[j]]
			;yy = fl_order_highs[left_inds[j]:right_inds[j]]
			;totmid = int_tabulated(xx,yy,/double,/sort)
			;left and right partial pixels
			li = left_inds[j]-1
			ri = right_inds[j]+1
			;fractions covered by the partial pixels
			lf = (wl_order_highs_left[left_inds[j]] - ws_base_left[j,i]) / xsis[li]
			rf = (ws_base_right[j,i] - wl_order_highs_right[right_inds[j]]) / xsis[ri]
			
			;sum up everything
			ltot = lf * fl_order_highs[li] * xsis[li]
			rtot = rf * fl_order_highs[ri] * xsis[ri]
			tot = totmid + rtot + ltot
			
			fs_base[j,i] = tot
			;if j eq 200 then stop
						
		endfor
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;FILL IN WARRAY
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		warray[*,i] = (rebin(ws_base,2048+buffer,norders))[buffer/2 : 2047 + buffer/2,i]
		
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;CREATE/FILL IN RECTIFIED ARRAYS
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
		;create the rectified arrays
		;these map the rectified fs to the warped positions xs_rect, ys_rect
		xs_rect = dblarr((2048+buffer)*upfactor,200*upfactor)
		ys_rect = dblarr((2048+buffer)*upfactor,200*upfactor)
		fs_rect = dblarr((2048+buffer)*upfactor,200*upfactor)

		;fill in the xs with all the same value,
		;fill in the ys with sliding values based on the original
		for j=0, 200*upfactor-1 do begin
			xs_rect[*,j] = xs_base
			ypos = double(j)/double(upfactor) - 100d
			ys_rect[*,j] = ys_base[*,i] + ypos
		endfor
		
		;put flux only in the very center row (so the tri_surf knows to bring flux back down right away)
		fs_rect[*,100*upfactor] = fs_base[*,i]
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;CREATE KERNEL
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
		
		;make sure width is odd
		if x_size_kernel mod 2 eq 0 then x_size_kernel += 1.
		kernel = dblarr(x_size_kernel,y_size_kernel)
		x_center_kernel = floor(x_size_kernel / 2)
		
		;fill in gaussians
		for j=0, nfibers-1 do begin
			ycenter = (j + 1.) * kernel_sep_pix * upfactor
			tk = psf_gaussian(npixel=[x_size_kernel,y_size_kernel],fwhm=[kernel_size_pix*upfactor,kernel_size_pix*upfactor],centroid=[x_center_kernel,ycenter],/double,/normalize)/double(nfibers)
			
			;if there is a slit, apply it
			if n_elements(slitwidth_um) ne 0 then begin
				slitwidth_pix = slitwidth_um * fiber_scale * upfactor
				slitkernel = psf_gaussian(npixel=x_size_kernel,fwhm=slitwidth_pix,centroid=x_center_kernel,/double,/normalize,ndimen=1)
				for k=0, y_size_kernel-1 do tk[*,k] = tk[*,k] * slitkernel
			endif
			
			kernel += tk
		endfor
		
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;SET UP FOR WARPING/CONVOLUTION
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		;create indices for the rect arrays
		
		nx_rect = (size(fs_rect))[1]
		ny_rect = (size(fs_rect))[2]
		x_grid_rect_d = rebin(lindgen(nx_rect),nx_rect,ny_rect)
		y_grid_rect_d = rebin(1#lindgen(ny_rect),nx_rect,ny_rect)
		
		;find offset between indices (0,2048) and coords (-1024,1023)
		xs_d_offset = min(xs_rect)
		ys_d_offset = min(ys_rect)
		
		;reposition the rect arrays to (0,2047)
		xs_rect_d = (xs_rect - xs_d_offset) * upfactor
		ys_rect_d = (ys_rect - ys_d_offset) * upfactor
		
		;reposition the grid to (-1024,1023)
		x_grid_rect = x_grid_rect_d + xs_d_offset
		y_grid_rect = y_grid_rect_d + ys_d_offset
		
		;how big must the output be
		x_size_d = max(floor(xs_rect_d)) - min(floor(xs_rect_d))
		y_size_d = max(floor(ys_rect_d)) - min(floor(ys_rect_d))
		if x_size_d mod 2 eq 1 then x_size_d += 1
		if y_size_d mod 2 eq 1 then y_size_d += 1
		
		y_min = min(floor(ys_rect))
		y_max = max(floor(ys_rect))


		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;PERFORM THE WARPING / CONVOLUTION
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		tmp = warp_tri(xs_rect_d,ys_rect_d,x_grid_rect_d,y_grid_rect_d,fs_rect,output_size = [x_size_d,y_size_d]);,/quintic)
		;remove below and uncomment above
		;tmp = fs_rect
		tmp_conv = convol(tmp,kernel,/center,/normalize,/edge_zero)
		;tmp_conv = tmp
		
		;rebin back down if necessary
		if upfactor ne 1 then tmp_conv_rebin = rebin(tmp_conv,x_size_d/upfactor,y_size_d/upfactor) else tmp_conv_rebin = tmp_conv
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;FILL IN TO DETECTOR ARRAY
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
		if y0_det lt 0 then begin
			extra = abs(y0_det)
			y0_det = 0
			y0_subarr = extra
		endif
		if y1_det gt 2047 then begin
			;extra = floor(max(y_grid_rect)+1024) - 2047
			extra = floor(y1_det - 2047)
			y1_det = 2047
			y1_subarr = y_size_d/upfactor - 1 - extra/upfactor
		endif
		
		;if nothing is on the detector skip this part
		if y1_det lt 0 or y0_det gt 2047 then continue
		
		;fill in
		t1 = fs_det[xd0:xd1 , y0_det:y1_det]
		t2 = t1 + tmp_conv_rebin[x0:x1, y0_subarr:y1_subarr]
		fs_det[xd0:xd1 , y0_det:y1_det] = t2
		
		print,'order ',i

	endfor
		
	specimg = fs_det
	
end