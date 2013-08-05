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

pro smooth_res_ryan, res, midwl, wl, fl, fsm
	;routine to do a simple convolution
	;with a gaussian of a fixed wavelength width
	;res = resolution
	;midwl = wavelength where resolution is converted to pixels
	;fsm = output flux
	
	;figure out which pixel has midwl
	a = abs(wl - midwl)
	al = where(a eq min(a))
	al = al[0]
	
	;what wavelength interval is this pixel
	pixwl = abs(wl[al] - wl[al-1])
	
	;what should be the wavelength interval of the gaussian (FWHM)
	dl = midwl/res
	;in pixels
	dlpix = dl/pixwl
	;convert to sigma from FWHM
	sigmap = dlpix / 2.35
	
	;diag
	print,'----------'
	print,'midwl',midwl
	print,'pixwl',pixwl
	print,'dl',dl
	print,'dlpix',dlpix
	
	;form kernel with 20X FWHM
	kernel_width = round(dlpix * 20)
	kernel_center = floor(kernel_width / 2)
	if kernel_width mod 2 eq 0 then kernel_width += 1
	print, 'kernel width ',kernel_width
	gx = dindgen(kernel_width)
	kernel = gaussian(gx,[1.,kernel_center,sigmap])
	cfl = convol(fl,kernel,/normalize,/edge_truncate)
	fsm = cfl
end


pro tram_projection2c, w, f, res, pixel_sampling, specimg, calw=calw, calf=calf, diag_out = diag_out, fiber_fractions = fiber_fractions, warray = warray, fiber_scale = fiber_scale, fiber_core_um = fiber_core_um, fiber_cladding_um = fiber_cladding_um, fiber_buffer_um = fiber_buffer_um, nfibers = nfibers, fiber_extra_sep_um = fiber_extra_sep_um, echellogram_file = echellogram_file, echellogram_wl_file = echellogram_wl_file, slitwidth_um = slitwidth_um

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
;;	mod_echel, xs, ys, ws, xs1, ys1, ws1
;;	xs_old = xs
;;	ys_old = ys
;;	ws_old = ws
;;	xs = xs1
;;	ys = ys1
;;	ws = ws1
	

	
	buffer = 100 ; the extra pixels (/2) on each side for the big arrays which get filled in and convolved
	upfactor = 2d
	
	ys = -1d * ys ;flip y around (stuart does this in his code)
	norders = (size(inp2))[2]
	
	warray=MAKE_ARRAY(2048, norders, /Double, Value=!Values.F_NAN)
	
	xs_pix = xs / pixel_size
	ys_pix = ys / pixel_size
	;xs,ys are the xy locations of the order sample points in mm
	;xs,ys_pix are the xy locations in pixels (with -1024 to 1023, not 0 to 2047)
	
	fs_base = dblarr((2048 + buffer)*upfactor,norders)
	fs_base2 = fs_base
	
	fs_det = dblarr(2048,2048)

	;read in the stellar spectrum
	input_wl = double(w * 1d3)
	input_fl = double(f)


	
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
	

	
	for i=0, norders-1 do begin
		;create highly sampled individual spectrum for integrating
		;have to upsample to get many pixels per output pixel
		;there is an option here to simply interpolate to the pixel values
		;below 1d5 was 1d6 (rct 4/2/2013)
		wl_order_highs = dindgen(1d6)/1d6 * (max(ws_base_right[*,i]) - min(ws_base_left[*,i]) + 10.) + min(ws_base_left[*,i]) - 5.
		midpts = (wl_order_highs[1:*] + wl_order_highs)/2.d
		xsis = midpts[1:*] - midpts
		xsis = [xsis[0],xsis,xsis[-1]]
		wl_order_highs_left = wl_order_highs - xsis/2d
		wl_order_highs_right = wl_order_highs + xsis/2d
		
;;		minl = min(wl_order_highs)
;;		maxl = max(wl_order_highs)
;;		res = 5d4
;;		pix = 3d
;;		pixwl = mean([minl,maxl]) / (res * pix)
;;		;pixwl = 1.15d / (res * pix)
;;		npn = round((maxl - minl)/pixwl) ;number of pixels in new array
;;		newwl = dindgen(npn)/double(npn)*(maxl - minl) + minl ;output wls

		;stop
		
		
		outlocs = where((max(input_wl) lt wl_order_highs) or (min(input_wl) gt wl_order_highs),noutlocs)
		
		fl_order_highs = interpol(input_fl,input_wl,wl_order_highs,/quadratic)
		;smooth_res_ryan,50000d,median(wl_order_highs),wl_order_highs,fl_order_highs,fl_order_highs2
		;stop
		;fl_order_highs1 = fl_order_highs
		;fl_order_highs = fl_order_highs2
		if noutlocs gt 0 then fl_order_highs[outlocs] = 0
		
		;print,minmax(fl_order_highs1)
		;print,minmax(fl_order_highs)
		
		;find where each output pixel falls in the upsampled array and integrate that part
		left_inds = value_locate(wl_order_highs_left,ws_base_left[*,i]) + 1
		right_inds = value_locate(wl_order_highs_right,ws_base_right[*,i])
		;for j=0, (size(fs_base))[1]-1 do fs_base[j,i] = tsum(wl_order_highs,fl_order_highs,left_inds[j]+1,right_inds[j]+1)
		for j=0, (size(fs_base))[1]-1 do begin
			totmid = total(fl_order_highs[left_inds[j]:right_inds[j]] * xsis[left_inds[j]:right_inds[j]],/double)
			li = left_inds[j]-1
			ri = right_inds[j]+1
			lf = (wl_order_highs_left[left_inds[j]] - ws_base_left[j,i]) / xsis[li]
			rf = (ws_base_right[j,i] - wl_order_highs_right[right_inds[j]]) / xsis[ri]
			ltot = lf * fl_order_highs[li] * xsis[li]
			rtot = rf * fl_order_highs[ri] * xsis[ri]
			tot = totmid + rtot + ltot
			fs_base[j,i] = tot
			;if j eq 200 then stop
			
;;			;fs_base[j,i] = int_tabulated(wl_order_highs[left_inds[j]+1:right_inds[j]+1],fl_order_highs[left_inds[j]+1:right_inds[j]+1],/double)
;;			fs_base[j,i] = total(wl_order_highs[left_inds[j]+1:right_inds[j]+1] * fl_order_highs[left_inds[j]+1:right_inds[j]],/double)
;;			left_extra_fraction = (wl_order_highs[left_inds[j]+1] - ws_base_left[j,i]) ;/ (xsis[left_inds[j]+1])
;;			left_extra = fl_order_highs[left_inds[j]+1] * left_extra_fraction
;;			right_extra_fraction = (wl_order_highs[right_inds[j]+1] - ws_base_right[j,i]); / (xsis[right_inds[j]+1])
;;			right_extra = fl_order_highs[right_inds[j]+1] * right_extra_fraction
;;			total_extra = left_extra + right_extra
;;			fs_base[j,i] = fs_base[j,i] + total_extra
			;if j eq 200 then stop
			
		endfor
		
;;		minl = min(ws_base[*,i])
;;		maxl = max(ws_base[*,i])
;;		nwl = reform(ws_base[*,i])
;;		rebin_ryan_debug,3d,50000d,wl_order_highs,fl_order_highs,wlr,flr,minl,maxl,new_wl_in = nwl
;;		
;;		fs_base2[*,i] = flr[0:2047+buffer]
;;
;;		stop
		;stop
		;;for j=0, (size(fs_base))[1]-1 do begin
			;fs_base[j,i] = tsum(wl_order_highs,fl_order_highs,left_inds[j],right_inds[j])
			
		
		;fill in the wavelength array
		warray[*,i] = (rebin(ws_base,2048+buffer,norders))[buffer/2 : 2047 + buffer/2,i]
		;warray[*,i] = wlr[buffer/2 : 2047+buffer/2]
		;fs_base[*,i] = flr[0:2047+buffer]
		
		
		
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
		;put flux only in the very center row (so the tri_surf knows to bring flux back down right away)
		fs_rect[*,100] = fs_base[*,i]
		
		
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
		
		for j=0, nfibers-1 do begin
			ycenter = (j + 1.) * kernel_sep_pix
			tk = psf_gaussian(npixel=[x_size_kernel,y_size_kernel],fwhm=[kernel_size_pix*upfactor,kernel_size_pix],centroid=[x_center_kernel,ycenter],/double,/normalize)/double(nfibers)
			if n_elements(slitwidth_um) ne 0 then begin
				slitwidth_pix = slitwidth_um * fiber_scale * upfactor
				slitkernel = psf_gaussian(npixel=x_size_kernel,fwhm=slitwidth_pix,centroid=x_center_kernel,/double,/normalize,ndimen=1)
				for k=0, y_size_kernel-1 do tk[*,k] = tk[*,k] * slitkernel
			endif
			kernel += tk
		endfor
		
	
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;convolve the kernel
;;		bigarr_lower = -1024 - buffer/2
;;		bigarr_upper = 1023 + buffer/2
;;		bigarr_lower_y = floor(min(ys_rect)) > bigarr_lower
;;		bigarr_upper_y = ceil(max(ys_rect)) < bigarr_upper
;;		
;;		det_lower = -1024
;;		det_upper = 1023
;;		ymin_ind = 0
;;		ymax_ind = bigarr_upper_y - bigarr_lower_y
;;		if det_lower gt bigarr_lower_y then begin
;;			extra = det_lower - bigarr_lower_y
;;			ymin_ind = extra
;;		endif
;;		if det_upper lt bigarr_upper_y then begin
;;			extra = bigarr_upper_y - det_upper
;;			ymax_ind = ymax_ind - extra
;;		endif
;;		if (det_upper lt bigarr_lower_y) or (det_lower gt bigarr_upper_y) then continue
;;		
;;		det_lower_y = det_lower > bigarr_lower_y
;;		det_upper_y = det_upper < bigarr_upper_y
;;
;;		det_lower_ind = det_lower + 1024
;;		det_upper_ind = det_upper + 1024
;;		det_lower_y_ind = det_lower_y + 1024
;;		det_upper_y_ind = det_upper_y + 1024
;;		
;;		
;;		bigarr_det_lower = buffer/2
;;		bigarr_det_upper = 2047 + buffer/2
;;		bigarr_det_lower_y = (buffer/2) 
;;		;tmp = tri_surf(fs_rect,xs_rect,ys_rect,gs=[1.,1.],bounds=[bigarr_lower,bigarr_lower_y,bigarr_upper,bigarr_upper_y]);,/linear)
;;		n_xs = bigarr_upper - bigarr_lower ;(size(fs_rect))[1]
;;		n_ys = bigarr_upper_y - bigarr_lower_y ;(size(fs_rect))[2]
;;		xsize = bigarr_upper - bigarr_lower + 1
;;		ysize = bigarr_upper_y - bigarr_lower_y + 1
;;		x0s = rebin(lindgen(n_xs),n_xs,n_ys) + bigarr_lower
;;		y0s = rebin(1#lindgen(n_ys),n_xs,n_ys) + bigarr_lower_y
;;		minx = min([x0s,xs_rect])
;;		miny = min([y0s,ys_rect])
;;		
;;		x0s_1 = x0s - bigarr_lower
;;		xs_rect_1 = xs_rect - bigarr_lower
;;		y0s_1 = y0s - bigarr_lower_y
;;		ys_rect_1 = ys_rect - bigarr_lower_y
;;		ys_rect_1 = reverse(ys_rect_1,2)
;;		y0s_1 = reverse(y0s_1,2)
;;		
;;		
;;	
;;		;tmp1 = warp_tri(x0s_1,y0s_1,xs_rect_1,ys_rect_1,fs_rect,/quintic);,output_size=[xsize,ysize])
;;		tmp2 = warp_tri(xs_rect_1,ys_rect_1,x0s_1,y0s_1,fs_rect,/quintic)
;;		
;;		if ysize lt 201 then tmp3 = tmp2[*,0:ysize-1] else tmp3 = tmp2
;;		
;;		
;;		tmp_conv = convol(tmp3,kernel,/center)
;;		;fs_det = fs_det + tmp_conv[bigarr_det_lower:bigarr_det_upper,bigarr_det_lower:bigarr_det_upper]
;;		t1 = fs_det[det_lower_ind:det_upper_ind,det_lower_y_ind:det_upper_y_ind]
;;		t2 = t1 + tmp_conv[bigarr_det_lower:bigarr_det_upper,ymin_ind:ymax_ind]
;;		fs_det[det_lower_ind:det_upper_ind,det_lower_y_ind:det_upper_y_ind] = t2
;;		print,'order ',i
		
		;oversized array (2048+buffer)^2 is called bigarr
		;this is where the resampling and convolution is done, to avoid edge effects
		;then take out the central region at the end
		
		;take out subarrays from bigarr defined by order coordinates with size 2148x(ys_rect_size)
		
		;first take the rectified image and shift the coordinates into _d system
		;where 0,0 is at lower corner of subarray d
		
		nx_rect = (size(fs_rect))[1]
		ny_rect = (size(fs_rect))[2]
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
		
		;resample onto the even grid
		;stop
		
		;tmp = convol(fs_rect,kernel,/center,/edge_truncate)
		;tmp_conv = warp_tri(xs_rect_d,ys_rect_d,x_grid_rect_d,y_grid_rect_d,tmp,output_size = [x_size_d,y_size_d])
		;if upfactor ne 1 then fs_rect_rebin = rebin(fs_rect,x_size_d/upfactor,ny_rect) else fs_rect_rebin = fs_rect
		
		tmp = warp_tri(xs_rect_d,ys_rect_d,x_grid_rect_d,y_grid_rect_d,fs_rect,output_size = [x_size_d,y_size_d],/quintic)
		;tmp = fs_rect
		
		tmp_conv = convol(tmp,kernel,/center,/normalize,/edge_truncate)
		;tmp_conv = tmp
		
		if upfactor ne 1 then tmp_conv_rebin = rebin(tmp_conv,x_size_d/upfactor,y_size_d) else tmp_conv_rebin = tmp_conv
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
		if y1_det lt 0 or y0_det gt 2047 then continue
		t1 = fs_det[xd0:xd1 , y0_det:y1_det]
		t2 = t1 + tmp_conv_rebin[x0:x1, y0_subarr:y1_subarr]
		fs_det[xd0:xd1 , y0_det:y1_det] = t2
		
		
		
		
		
		print,'order ',i

		
		
	endfor


	
	specimg = fs_det
	
end