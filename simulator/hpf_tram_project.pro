;+
; NAME:
;  hpf_tram_project
;
; PURPOSE:
;
;	Take a set of spectra and distribute flux onto an array of pixels. This is the meat 
;		of the simulator and has gone through many iterations. The current overall strategy:
;		- assign WAVELENGTH(X,Y) for each pixel according to an echellogram from Larry
;		- slice up the spectra according to these wavelength boundaries into 1-d spectra
;		- lay these out in "tramlines" onto a sub-detector array
;		- convolve these tramlines using a 2-d kernel
;		- warp these tramlines using tri_warp, into the correct curvature
;		- add this final order to the overall detector array
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;	hpf_tram_project, spec_params, optical_params, proj_params, det_params, out_file
;
; INPUTS:
;
;	spec_params: The spec_params structure as defined in hpf_initialize_spec_params
;
;	optical_params: The optical_params structure as defined in hpf_initialize_optical_params
;
;	proj_params: The proj_params structure as defined in hpf_process_optical_params
;
;	det_params: The det_params structure as defined in hpf_intialize_det_params (not needed)
;
;	out_file: The output fits location for the fluence array
;
; OUTPUTS:
;	
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-

pro hpf_tram_project, spec_params, optical_params, proj_params, det_params, out_file, generate_wlimg = generate_wlimg

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;PARAMETERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	buffer = 100 ; the extra pixels (/2) on each side for the big arrays which get filled in and convolved
	
	upfactor = optical_params.projection_upsample
	
	n_specs = n_elements(spec_params.spec_file)
	
	fs_det = dblarr(2048,2048) ; detector array
	
	norders = (size(proj_params.xs))[2]
	
	if keyword_set(generate_wlimg) then warray = dblarr(2048,norders)
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; FIBER and KERNEL PARAMETERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	fiber_core_um = optical_params.fiber_core_um
	fiber_buffer_um = optical_params.fiber_buffer_um
	fiber_scale = optical_params.fiber_scale
	fiber_extra_sep_um = optical_params.fiber_extra_sep_um
	slitwidth_um = optical_params.slitwidth_um

	;kernel size in pixels
	kernel_size_pix = fiber_core_um * fiber_scale
	fiber_extra_sep_pix = fiber_extra_sep_um * fiber_scale
	
	;separation between gaussians
	kernel_sep_pix = fiber_buffer_um * fiber_scale + fiber_extra_sep_pix

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; GET THE SPECTRA
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	input_wl = spec_params.shift_wl
	;for i=0, n_specs-1 do *(input_wl[i]) = double((1d3 * *(input_wl[i])))
	input_fl = spec_params.fl
	;for i=0, n_specs-1 do *(input_fl[i]) = double(*input_fl[i])
	;wl in nm
	;fl in photons / s / um

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;GET ECHELLOGRAM INFO
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	xs = proj_params.xs
	ys = proj_params.ys
	ws = proj_params.ws
	xs_pix = proj_params.xs_pix
	ys_pix = proj_params.ys_pix
	fs_base = proj_params.fs_base
	xs_base = proj_params.xs_base
	ws_base = proj_params.ws_base / 1d3
	xs_base_left = proj_params.xs_base_left
	xs_base_right = proj_params.xs_base_right
	ws_base_left = proj_params.ws_base_left / 1d3
	ws_base_right = proj_params.ws_base_right / 1d3
	ys_base = proj_params.ys_base
	xs_base_2d = proj_params.xs_base_2d
	
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;CREATE FLUX ARRAY
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	for i=0, norders-1 do begin
		print,'MEMORY IN TRAM PROJ, BEFORE ORDER',i
		print,memory()/1d6
		;create highly sampled individual spectrum for integrating
		;have to upsample to get many pixels per output pixel
		;there is an option here to simply interpolate to the pixel values
		;define high sample wavelengths/pixel edges
		;wl_order_highs = dindgen(1d6)/1d6 * (max(ws_base_right[*,i]) - min(ws_base_left[*,i]) + 10.) + min(ws_base_left[*,i]) - 5.
		
		
		;find where each output pixel falls in the upsampled array and integrate that part
		left_inds = ptrarr(n_specs,/allocate_heap)
		right_inds = ptrarr(n_specs,/allocate_heap)
		xsis = ptrarr(n_specs,/allocate_heap)
		input_wl_left = ptrarr(n_specs,/allocate_heap)
		input_wl_right = ptrarr(n_specs,/allocate_heap)
		for j=0, n_specs-1 do begin
			midpts = double(((*input_wl[j])[1:*]) + double(*input_wl[j]))/2d
			xsis1 = midpts[1:*] - midpts
			*xsis[j] = double([xsis1[0],xsis1,xsis1[-1]])
			*input_wl_left[j] = *input_wl[j] - *xsis[j]/2d
			*input_wl_right[j] = *input_wl[j] + *xsis[j]/2d
			*left_inds[j] = value_locate(*input_wl_right[j],ws_base_left[*,i]) + 1
			*right_inds[j] = value_locate(*input_wl_left[j],ws_base_right[*,i]) - 1
		endfor

		;stop
		;xsis /= 1d3 ; to get from nm to microns so fluxes are in /um
		;wl_order_highs_left /= 1d3
		;wl_order_highs_right /= 1d3 ;need to propagate constant units throughout code
		;stop
		
		for j=0, (size(fs_base))[1]-1 do begin
			for k=0, n_specs-1 do begin
				li_part = (*left_inds[k])[j]
				ri_part = (*right_inds[k])[j] + 1
				li = (*left_inds[k])[j]
				ri = (*right_inds[k])[j] + 1
				
				lf = (double((*input_wl_right[k])[li]) - double(ws_base_left[j,i])) / double((*xsis[k])[li_part])
				rf = (double(ws_base_right[j,i]) - double((*input_wl_left[k])[ri])) / double(((*xsis[k])[ri_part]))
				
				totmid = total((*input_fl[k])[li:ri] * (*xsis[k])[li:ri] ,/double)
				ltot = double(lf) * double((*input_fl[k])[li_part]) * double((*xsis[k])[li_part])
				rtot = double(rf) * double((*input_fl[k])[ri_part]) * double((*xsis[k])[ri_part])
				tot = totmid + rtot + ltot
				
				fs_base[j,i,k] = tot
				
				;if j eq 100 and k eq 2 then stop
				;if (k eq 2) and (ws_base[j,i] gt .82602) and (ltot + rtot gt 2.5) then stop
				;if j gt 1500 and k eq 2 then stop
				
			endfor

		endfor
		;stop
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;FILL IN WARRAY
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		if keyword_set(generate_wlimg) then warray[*,i] = (rebin(ws_base,2048+buffer,norders))[buffer/2 : 2047 + buffer/2,i]
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;CREATE/FILL IN RECTIFIED ARRAYS
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			
		xs_rect = dblarr((2048+buffer)*upfactor,200*1d)
		ys_rect = dblarr((2048+buffer)*upfactor,200*1d)
		fs_rect = dblarr((2048+buffer)*upfactor,200*1d)
		for j=0, 199 do begin
			xs_rect[*,j] = xs_base
			ypos = double(j) - 100d
			ys_rect[*,j] = ys_base[*,i] + ypos
		endfor
		
		;fill in the tramlines
		offset = 100 - (n_specs - 1d)/2d * kernel_sep_pix
		for j=0, n_specs-1 do begin
			ycenter = (j) * kernel_sep_pix + offset
			fs_rect[*,ycenter] = fs_base[*,i,j]
		endfor
		
		;create the kernel
		;y_size_kernel = fiber_buffer_um * fiber_scale + 2*fiber_extra_sep_pix
		x_size_kernel = (floor(3. * kernel_size_pix))*upfactor		
		if x_size_kernel mod 2 eq 0 then x_size_kernel += 1.
		y_size_kernel = floor(3. * kernel_size_pix)
		;if y_size_kernel mod 2 eq 0 then y_size_kernel += 1.
		;kernel = dblarr(x_size_kernel,y_size_kernel)
		x_center_kernel = floor(x_size_kernel / 2)
		y_center_kernel = floor(y_size_kernel / 2)
		
		case optical_params.kernel_type of
		'gaussian': kernel = psf_gaussian(npixel=[x_size_kernel,y_size_kernel],fwhm=[kernel_size_pix*upfactor,kernel_size_pix*1d],centroid=[x_center_kernel,y_center_kernel],/double,/normalize)
		
		'hrs': kernel = hrs_kernel_2d(x_size_kernel,y_size_kernel,kernel_size_pix)
		else:stop
		endcase

		
		if slitwidth_um gt 1 then begin
			slitwidth_pix = slitwidth_um * fiber_scale * upfactor 
			case optical_params.kernel_type of
			'gaussian': slitkernel = psf_gaussian(npixel=x_size_kernel,fwhm=slitwidth_pix,centroid=x_center_kernel,/double,/normalize,ndimen=1)
			'hrs': slitkernel = hrs_kernel_1d(x_size_kernel,slitwidth_pix)
			else:stop
			endcase
			
			for k=0, y_size_kernel-1 do kernel[*,k] = kernel[*,k] * slitkernel
			tot_kernel = total(kernel,/double)
			kernel /= tot_kernel
		endif	
			
		
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
		ys_rect_d = (ys_rect - ys_d_offset) * 1d
		
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
		
		gx = (max(xs_rect_d) - min(xs_rect_d))/double(x_size_d-1d)
		gy = (max(ys_rect_d) - min(ys_rect_d))/double(y_size_d-1d)
		tmp = tri_surf(fs_rect,xs_rect_d,ys_rect_d,gs=[gx,gy])
		print,'trisurf done'
		tmp_conv = convol(tmp,kernel,/center,/normalize,/edge_zero)
		;rebin back down if necessary
		st = size(tmp_conv,/dimen)
		if upfactor ne 1 then tmp_conv_rebin = rebin(tmp_conv,st[0]/upfactor,st[1]/1d) else tmp_conv_rebin = tmp_conv
		
		;stop
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;FILL IN TO DETECTOR ARRAY
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		;find the chunk of the detector that this corresponds to
		;y locations on detector
		y0_det = min(y_grid_rect) + 1024; + buffer/2
		y1_det = y0_det + y_size_d/1d - 1
		;y locations in the sub-array
		y0_subarr = 0
		y1_subarr = y_size_d/1d - 1
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
			y1_subarr = y_size_d/1d - 1 - extra/1d
		endif
		
		;if nothing is on the detector skip this part
		if y1_det lt 0 or y0_det gt 2047 then continue
		
		;fill in
		t1 = fs_det[xd0:xd1 , y0_det:y1_det]
		t2 = t1 + tmp_conv_rebin[x0:x1, y0_subarr:y1_subarr]
		fs_det[xd0:xd1 , y0_det:y1_det] = t2
		
		print,'order ',i
	
	endfor
		
	;sxaddpar,h
	spec_names = tag_names(spec_params)
	for i=0, n_elements(spec_params)-1 do begin
		for j=0, n_elements(spec_names)-1 do begin
			if typename(spec_params.(j)[i]) eq 'POINTER' then continue
			if strupcase(spec_names[j]) eq 'RV' then begin
				sxaddpar,h,'RV',spec_params.(j)[0]
				continue
			endif
			if n_elements(spec_params.(j)[i]) ne 1 then continue
			name = spec_names[j]+'_'+string(i,format='(I1)')
			sxaddpar,h,name,spec_params.(j)[i]
		endfor
	endfor
	optical_names = tag_names(optical_params)
	for i=0, n_elements(optical_params)-1 do begin
		for j=0, n_elements(optical_names)-1 do begin
			if typename(optical_params.(j)[i]) eq 'POINTER' then continue
			if n_elements(optical_params.(j)[i]) ne 1 then continue
			name = optical_names[j]
			sxaddpar,h,name,optical_params.(j)[i]
		endfor
	endfor
	

	mwrfits,fs_det,out_file,h,/create
	
	
	if keyword_set(generate_wlimg) then begin
		;wlimg_file = file_dirname(out_file,/mark_directory) + 'wlimg.fits'
		wlimg_file = optical_params.wlimg_file
		mwrfits,warray,wlimg_file,/create
		;optical_params.wlimg_file = wlimg_file
	endif

	
end