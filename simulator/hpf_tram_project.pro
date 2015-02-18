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
	convol_upfactor = optical_params.convol_upsample
	n_pixels = optical_params.n_pixels
	
	n_specs = n_elements(spec_params.spec_file)
	
	fs_det = dblarr(n_pixels,n_pixels) ; detector array
	
	norders = (size(proj_params.xs))[2]
	
	if keyword_set(generate_wlimg) then begin
		warray = dblarr(n_pixels,norders)
		bins_warray = dblarr(n_pixels,norders)
	endif
   
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
	slitwidth_pix = slitwidth_um * fiber_scale * convol_upfactor 
	
	;separation between gaussians
	kernel_sep_pix = fiber_buffer_um * fiber_scale + fiber_extra_sep_pix

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; GET THE SPECTRA
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	input_wl = spec_params.shift_wl
	input_fl = spec_params.fl
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
	
	ys_base_left = proj_params.ys_base_left
	ys_base_right = proj_params.ys_base_right
	ys_base_left_convol = proj_params.ys_base_left
	ys_base_right_convol = proj_params.ys_base_right
	
	ws_base_noup = proj_params.ws_base_noup /1d3
	ws_base_left_noup = proj_params.ws_base_left_noup /1d3
	ws_base_right_noup = proj_params.ws_base_right_noup /1d3
	
	fs_base_convol = proj_params.fs_base_convol
	xs_base_convol = proj_params.xs_base_convol
	xs_base_left_convol = proj_params.xs_base_left_convol
	xs_base_right_convol = proj_params.xs_base_right_convol
	ws_base_convol = proj_params.ws_base_convol / 1d3
	ws_base_left_convol = proj_params.ws_base_left_convol / 1d3
	ws_base_right_convol = proj_params.ws_base_right_convol / 1d3
	
	sfb = size(fs_base_convol,/dimen)
	upratio = double(upfactor) / double(convol_upfactor)
	if double(sfb[0]) * upratio mod 1 ne 0 then stop

	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;CREATE FLUX ARRAY
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	for i=0, norders-1 do begin
		print,'MEMORY IN TRAM PROJ, BEFORE ORDER',i
		print,memory()/1d6		
		
		;find where each output pixel falls in the upsampled array and integrate that part
		left_inds_i = ptrarr(n_specs,/allocate_heap)
		right_inds_i = ptrarr(n_specs,/allocate_heap)

    ;For all the spectra, map X/Y/WL
		for j=0, n_specs-1 do begin	
			wwii = dindgen(n_elements(*input_wl[j]))
			*left_inds_i[j] = interpol(wwii,*input_wl[j],ws_base_left_convol[*,i])
			*right_inds_i[j] = interpol(wwii,*input_wl[j],ws_base_right_convol[*,i])
		endfor
    
    ;For each pixel
		for j=0, (size(fs_base_convol))[1]-1 do begin
		  ;For each spectrum
			for k=0, n_specs-1 do begin
				tot2 = tsum(*input_wl[k],*input_fl[k],(*left_inds_i[k])[j],(*right_inds_i[k])[j])
				fs_base_convol[j,i,k] = tot2
			endfor
		endfor
		
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;FILL IN WARRAY
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		if keyword_set(generate_wlimg) then begin
			warray[*,i] = (rebin(ws_base_convol,n_pixels+buffer,norders))[buffer/2 : 2047 + buffer/2,i]
			wbins = abs(ws_base_right_noup - ws_base_left_noup)
			bins_warray[*,i] = wbins[buffer/2 : 2047 + buffer/2,i]
			;if i eq 10 then stop
		endif
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;CREATE/FILL IN RECTIFIED ARRAYS
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			
		xs_rect = dblarr((n_pixels+buffer)*upfactor,200*1d)
		ys_rect = dblarr((n_pixels+buffer)*upfactor,200*1d)
		fs_rect = dblarr((n_pixels+buffer)*convol_upfactor,200*1d)
		xs_rect_left = dblarr((n_pixels+buffer)*upfactor,200*upfactor)
		xs_rect_right = dblarr((n_pixels+buffer)*upfactor,200*upfactor)
		ys_rect_left = dblarr((n_pixels+buffer)*upfactor,200*upfactor)
		ys_rect_right = dblarr((n_pixels+buffer)*upfactor,200*upfactor)
		
		;for warping
		for j=0, 199 do begin
			xs_rect[*,j] = xs_base
			ypos = double(j) - 100d
			ys_rect[*,j] = ys_base[*,i] + ypos
		endfor
		
		;for polyclipping
		offsets = dindgen(200*upfactor)/(200*upfactor - 1.) * 199. - 99.5
		for j=0, 200.*upfactor-1 do begin
			xs_rect_left[*,j] = xs_base_left
			xs_rect_right[*,j] = xs_base_right
			ys_rect_left[*,j] = ys_base_left[*,i] + offsets[j]
			ys_rect_right[*,j] = ys_base_right[*,i] + offsets[j]
		endfor
		
		;fill in the tramlines
		offset = 100 - (n_specs - 1d)/2d * kernel_sep_pix
		for j=0, n_specs-1 do begin
			ycenter = round((j) * kernel_sep_pix + offset)
			fs_rect[*,ycenter] = fs_base_convol[*,i,j]
		endfor
		
		;create the kernel
		x_size_kernel = (floor(3. * kernel_size_pix))*convol_upfactor		
		if x_size_kernel mod 2 eq 0 then x_size_kernel += 1.
		y_size_kernel = floor(3. * kernel_size_pix)		
		if ((optical_params.kernel_type eq 'hrs') or (optical_params.kernel_type eq 'gaussian') or (optical_params.kernel_type eq 'step3')) and (convol_upfactor ne 1) then begin
			x_size_kernel = (floor(7. * kernel_size_pix))*convol_upfactor
			y_size_kernel = (floor(3. * kernel_size_pix));*convol_upfactor ;;;;;;;
			xyfact = round(x_size_kernel / y_size_kernel)
			y_size_kernel = x_size_kernel / xyfact
		endif

		;if y_size_kernel mod 2 eq 0 then y_size_kernel += 1.
		;kernel = dblarr(x_size_kernel,y_size_kernel)
		x_center_kernel = floor(x_size_kernel / 2)
		y_center_kernel = floor(y_size_kernel / 2)
		
		case optical_params.kernel_type of
		'gaussian': kernel = psf_gaussian(npixel=[x_size_kernel,y_size_kernel],fwhm=[kernel_size_pix*convol_upfactor,kernel_size_pix*1d],centroid=[x_center_kernel,y_center_kernel],/double,/normalize)
		'hrs': kernel = hrs_kernel_2d(x_size_kernel,y_size_kernel,kernel_size_pix)
		'circlehat': kernel = hpf_circlehat_kernel(x_size_kernel,y_size_kernel,kernel_size_pix * convol_upfactor)
		'step3': begin
			kernel_disp = hpf_disp_kernel(x_size_kernel,slitwidth_pix,ws_base_convol[*,i])
			kernel_xdisp = hpf_xdisp_kernel(y_size_kernel,kernel_size_pix)
			kernel_psf = hpf_camera_psf(x_size_kernel,y_size_kernel,convol_upfactor,1d)
		end
		'step3_tilt': begin
			kernel_disp = hpf_disp_kernel(x_size_kernel,slitwidth_pix,ws_base_convol[*,i])
			kernel_xdisp = hpf_xdisp_kernel(y_size_kernel,kernel_size_pix)
			kernel_psf = hpf_camera_psf(x_size_kernel,y_size_kernel,convol_upfactor,1d)
			kernel = hpf_combine_tilt_kernel(kernel_disp,kernel_xdisp,kernel_psf,optical_params.slit_angle, convol_upfactor=convol_upfactor)
			tot_kernel = total(kernel,/double)
			kernel /= tot_kernel
		end
		else:stop
		endcase

		
		if slitwidth_um gt 1 and (optical_params.kernel_type ne 'step3' and optical_params.kernel_type ne 'step3_tilt') then begin 
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
		
		nx_rect = (size(fs_rect))[1] / convol_upfactor
		ny_rect = (size(fs_rect))[2]
		x_grid_rect_d = rebin(lindgen(nx_rect),nx_rect,ny_rect)
		y_grid_rect_d = rebin(1#lindgen(ny_rect),nx_rect,ny_rect)
		
		;find offset between indices (0,2048) and coords (-1024,1023)
		xs_d_offset = min(xs_rect)
		ys_d_offset = min(ys_rect)
		
		;reposition the rect arrays to (0,2047)
		xs_rect_d = (xs_rect - xs_d_offset) * upfactor
		ys_rect_d = (ys_rect - ys_d_offset) * 1d
		
		xs_rect_d_left = (xs_rect_left - xs_d_offset)*upfactor
		xs_rect_d_right = (xs_rect_right - xs_d_offset)*upfactor
		ys_rect_d_left = (ys_rect_left - ys_d_offset)*upfactor
		ys_rect_d_right = (ys_rect_right - ys_d_offset)*upfactor
		
		;reposition the grid to (-1024,1023)
		x_grid_rect = x_grid_rect_d + xs_d_offset
		y_grid_rect = y_grid_rect_d + ys_d_offset
		
		;how big must the output be
		x_size_d = max(floor(xs_rect_d)) - min(floor(xs_rect_d))
		y_size_d = max(floor(ys_rect_d)) - min(floor(ys_rect_d))
		;if x_size_d mod 2 eq 1 then x_size_d += 1
		;if y_size_d mod 2 eq 1 then y_size_d += 1
		
		y_min = min(floor(ys_rect))
		y_max = max(floor(ys_rect))


		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;PERFORM THE WARPING / CONVOLUTION
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		gx = (max(xs_rect_d) - min(xs_rect_d))/double(x_size_d-1d)
		gy = (max(ys_rect_d) - min(ys_rect_d))/double(y_size_d-1d)
		
		;test = hpf_combine_tilt_kernel(kernel_disp,kernel_xdisp,kernel_psf)
		
		print,'begin convolve'
		if optical_params.kernel_type eq 'step3' then begin
			con1 = convol(fs_rect,kernel_disp,/center,/normalize,/edge_zero)
			con2 = convol(con1,transpose(kernel_xdisp),/center,/normalize,/edge_zero)
			tmp_conv_big = convol(con2,kernel_psf,/center,/normalize,/edge_zero)
			;tmp_conv_big = fs_rect
		endif else begin
			tmp_conv_big = convol(fs_rect,kernel,/center,/normalize,/edge_zero)
		endelse

		stc = size(tmp_conv_big,/dimen)
		if (upfactor/convol_upfactor) ne 1 then tmp_conv = frebin(temporary(tmp_conv_big),stc[0]*(upfactor/convol_upfactor),stc[1],/total) else tmp_conv = temporary(tmp_conv_big)
		tmp_conv_big = 0
		
		; interpolate in the y-direction if we need to upsample there?		
		tmp_conv1 = tmp_conv
		ss = stc[0]*upfactor/convol_upfactor
		tmp_conv = dblarr(ss,stc[1]*upfactor)
		xx1 = dindgen(stc[1])
		xx2 = dindgen(stc[1]*upfactor)/double(stc[1]*upfactor - 1d) * double(stc[1])
		for j=0, ss-1 do begin
			yy = tmp_conv1[j,*]
			oo = interpol(yy,xx1,xx2) / double(upfactor)
			tmp_conv[j,*] = oo
			;if j gt ss/2 then stop
		endfor

		print,'end convolve'
		;tmp_conv = tmp1
		;rebin back down if necessary
		print,'begin warp'
		
		if (optical_params.bypass_warp) and (optical_params.orders_shape eq 1 or optical_params.orders_shape eq 4) then begin
			tmp_warp = tmp_conv
		endif else begin
			case optical_params.warp_style of
			'tri_surf':tmp_warp = tri_surf(tmp_conv,xs_rect_d,ys_rect_d,gs=[gx,gy],linear=optical_params.linear_warp)
			'polyclip':begin
				print,'polyclipping begin'
				nnx = n_elements(xs_rect_d_left)
				nnx2 = nnx * 4.
				ppx = [[reform(xs_rect_d_left,nnx)], [reform(xs_rect_d_left,nnx)], [reform(xs_rect_d_right,nnx)], [reform(xs_rect_d_right,nnx)]] + 0.5d/double(upfactor)
				ppy = [[reform(ys_rect_d_left - 0.5,nnx)], [reform(ys_rect_d_left + 0.5,nnx)], [reform(ys_rect_d_right + 0.5,nnx)], [reform(ys_rect_d_right - 0.5,nnx)]]
				ppx = reform(transpose(temporary(ppx)),nnx2)
				ppy = reform(transpose(temporary(ppy)),nnx2)
				polyinds_in = dindgen(nnx2/4d + 1d)*4.
				polyclip_x_size = (size(xs_rect_d_left))[1]
				polyclip_y_size = floor(max(ys_rect_d_left) - min(ys_rect_d_left)) + 1.
				inds = polyfillaa(ppx,ppy,polyclip_x_size,polyclip_y_size,areas=areas,poly_indices = polyinds_in)
				tmp_warp = dblarr(polyclip_x_size,polyclip_y_size)
				;tot_areas_out = dblarr(polyclip_x_size,polyclip_y_size)s
				if nnx ne n_elements(tmp_conv) then message,'polyclip and fluence array differ nelements'
				for j=0, nnx-1 do begin
					is_edge = (j mod polyclip_x_size eq 0) or (floor(j / polyclip_x_size) eq 0) or (floor(j / polyclip_x_size) eq (polyclip_y_size - 1))
					ii1 = polyinds_in[j]
					ii2 = polyinds_in[j+1] - 1
					if ii1 gt ii2 then continue
					tt2 = double(areas[ii1:ii2]) * tmp_conv[j]
					test = total(areas[ii1:ii2],/double)
					if test ne 1 and ~is_edge then begin
						scale = 1d/test
						areas2 = areas[ii1:ii2] * scale
						tt2_scaled = areas2 * tmp_conv[j]
						tmp_warp[inds[ii1:ii2]] = temporary(tmp_warp[inds[ii1:ii2]]) + tt2_scaled
						;tot_areas_out[inds[ii1:ii2]] = total(areas2,/double)
					endif else begin
						tmp_warp[inds[ii1:ii2]] = temporary(tmp_warp[inds[ii1:ii2]]) + tt2
						;tot_areas_out[inds[ii1:ii2]] = total(areas[ii1:ii2],/double)
					endelse
				endfor
				print,'polyclipping end'
			end
			else:stop
			endcase
		endelse
		
		st = size(tmp_warp,/dimen)
		if upfactor ne 1 then tmp_conv_rebin = frebin(tmp_warp,st[0]/upfactor,st[1]/upfactor,/total) else tmp_conv_rebin = tmp_warp
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;FILL IN TO DETECTOR ARRAY
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		;find the chunk of the detector that this corresponds to
		;y locations on detector
		y0_det = min(y_grid_rect) + 1024; + buffer/2 ;4-30-14 trying 1023 instead of 1024
		y1_det = y0_det + y_size_d/1d - 1
		;y locations in the sub-array
		y0_subarr = 0
		y1_subarr = y_size_d/1d - 1
		;x locations in sub-array and detector always the same
		x0 = buffer/2
		x1 = x0 + 2047
		xd0 = 0
		xd1 = 2047
		
		;y0_det = round(y0_det)
		;y1_det = round(y1_det)
		
		;if there is any that spills off the end of the detector, adjust the two coordiantes
		if y0_det lt 0 then begin
			extra = abs(y0_det); try this + 1 4-30-14
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
		
		;y0_subarr = round(y0_subarr)
		;y1_subarr = round(y1_subarr)
		
		;fill in
		;t1 = fs_det[xd0:xd1 , y0_det:y1_det]
		t2 = temporary(fs_det[xd0:xd1 , y0_det:y1_det]) + temporary(tmp_conv_rebin[x0:x1, y0_subarr:y1_subarr])
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
	
	fs_det = 0
	
	if keyword_set(generate_wlimg) then begin
		;wlimg_file = file_dirname(out_file,/mark_directory) + 'wlimg.fits'
		wlimg_file = optical_params.wlimg_file
		out = {warray:warray, bins_warray:bins_warray}
		mwrfits,out,wlimg_file,/create
		;optical_params.wlimg_file = wlimg_file
	endif

	
end