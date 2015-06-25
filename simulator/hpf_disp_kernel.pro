function hpf_disp_kernel, npix, slitwidth_pix_in, wvl

	; make grating voigt
	c0     = 2.0056d
	c1     = 1.0593d
	fratio = 0.10d
	const  = 1.0d - c0*c1 + sqrt(fratio^2 + 2.0d*c1*fratio + c0*c0*c1*c1)
	mwl = median(wvl)
	dwl_det = median(abs(wvl[1:*] - wvl))
	wlrange = npix * dwl_det
	
	grating_res = 1.19d6
	
	dwl_grating = mwl / grating_res
	dwl_grating_pix = dwl_grating/3d
	wlrange_grating_kernel = (41d * dwl_grating_pix)
	wlrange_grating_kernel = wlrange
	npix_grating = wlrange / dwl_grating_pix
	npix_grating_kernel = 301d
	;npix_grating_kernel = npix_grating
	;if npix_grating mod npix_x ne 0 then begin
	fact = round(npix_grating / npix)
	;if fact mod 2 eq 0 then fact += 1
	npix_grating = fact * npix
	npix_grating_kernel = npix_grating
	;endif
	
	;upsample_xy = npix_x / npix_y

	slitwidth_pix = slitwidth_pix_in * fact
	;fiber_diameter_pix = fiber_diameter_pix_in * fact * upsample_xy
	
	grating_wvl = dindgen(npix_grating_kernel) * dwl_grating_pix + mwl - (wlrange_grating_kernel/2d)


	wis = dindgen(npix_grating_kernel) - floor(npix_grating_kernel/2)
	;ws1 = (dindgen(npix) - npix/2d)*dwl + mwl

	acoef = (npix_grating_kernel-1d) / (alog(max(grating_wvl)) - alog(min(grating_wvl)))
	ws = wis / acoef
	;res = mwl / (width * dwl)

	
	siglim = 12.0
	sigma  = 1.0 /(const*grating_res*2.0*sqrt(2.0*alog(2.0)))   ; sigma/c = dv/c / (2sqrt(2log2))
	;wvoigt = (findgen(nw) - float(nw)/2.0) / acoef
	;xv     = where(abs(wvoigt/sigma) le siglim)

	av     = sqrt(alog(2.0))*fratio
	;uv     = wvoigt[xv]/(sqrt(2.0)*sigma)
	uv = ws/(sqrt(2d) * sigma)
	vpsf   = voigt(av,uv)/(sigma*sqrt(2.0*!pi))

	;vint   = int_tabulated(wvoigt[xv],vpsf)
	;vint = int_tabulated(ws,vpsf)
	vint = total(vpsf)
	vpsf   = temporary(vpsf)/vint
	vsum   = total(vpsf)
	psf    = vpsf
	psfsum = vsum

	xs = abs(dindgen(npix_grating) - npix_grating/2 + 0.5)
	slit1 = dblarr(npix_grating)
	slit1[where(xs lt slitwidth_pix/2.)] = 1.
	
	conv = convol(slit1,psf,/center,/edge_truncate,/normalize)
	
	out1 = frebin(conv,npix,/total)
	
	if npix mod 2 eq 1 then out = out1 else begin
		final_x = dindgen(npix) - 0.5
		xx = dindgen(npix)
		out = interpol(out1,xx,final_x)
	endelse


	return,out
end