pro flux_checker
	exptime = 5d
	comparespec = '~/scratch/simulator/data/sp30/test1fl000_spec.fits'
	
	;;;FLUX LOAD
	c = 2.99792458d8	; m/s
	h = 6.626d-34 ;J s
	res = 50000.d
	if ~keyword_set(pixel_sampling) then pixel_sampling = 3d
	

	minl = .8d
	maxl = 1.4d

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;READ SPECTRUM
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	;spec = readfits('bt_34_extended.fits',head)
	;spec = transpose(spec)
	;cal = 0
	
	spec = mrdfits('~/work/spectra/ffp_30_300.fits')
	spec = transpose(spec)
	cal = 1

	wmod   = reform(spec[*,0])	; assume wavelength is in microns
	fmod   = reform(spec[*,1])	; assume flux is in Watts/m2/micron
	;dwmod  = 0.015d-4		; sampling of Hauschildt models in microns
	;sort
	w      = wmod[sort(wmod)]
	f      = fmod[sort(wmod)]
	
	;if flat then f[*] = mean(f)
	

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;CUT DOWN AND UPSAMPLE
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	g = where(w ge minl - .05d and w le maxl + .05d,ng)
	w = double(w[g])
	f = double(f[g])
	
	
	;upsample (ensures even/sufficient pixel sampling)
	if ~keyword_set(upsample_factor) then upsample_factor = -1
	if upsample_factor ne -1 then begin
		nel1 = n_elements(w)
		nel2 = upsample_factor*nel1 ;multiply this for upsampling
		neww = dindgen(nel2)/double(nel2) * (maxl - minl) + minl
		newf = interpol(f,w,neww)

		w = neww
		f = newf
		
		usf = string(nel2 / nel1,format='(I3)')
		
	endif
	ng = n_elements(w)
	
	w_orig = w
	f_orig = f
	

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;SCALE FLUX
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	jzp     = 3.31d-9	; Vega flux at J  (W/m2/micron)
	jvega   = 0.025d	; Vega mag at J
	deltaj  = 0.26d		; FWHM of J in microns
	isowj   = 1.215d	; Isophotal wavelength of J in microns
	rtel    = 10./2.0d	; radius of telescope in m
	telarea = double(!pi*rtel^2)    ; telescope area in m^2
	tput    = 0.025d		; total instrument throughput 

	;print, 'Scaling spectrum ...'
	
	jmag = 9.
	
	jmin   = isowj - deltaj/2.0d
	jmax   = isowj + deltaj/2.0d
	jndx   = where((w ge jmin) and (w le jmax))
	fmean  = int_tabulated(w[jndx],f[jndx])/deltaj
	jscale = (jzp/fmean)*10.0d^(-0.4d*(jmag - jvega))
	f      = temporary(f)*jscale		; scale spectrum to input J mag

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;CONVERT TO PHOTONS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;photons/s/m2/micron
	p = w * f / (h*c) * 1d-6

	pout = temporary(p)*telarea*tput ;photons/s/micron

	fold = f
	f = pout
	
	f = f * exptime
	
	if cal then begin
		w = w_orig
		f = f_orig
		p = w * f / (h*c) * 1d-6
		f = p * exptime
	endif



	;;;MODEL LOAD
	optical_model_file = 'support/model_ramsey_513.sav'
	restore,optical_model_file
	xs = model.xs
	ys = model.ys
	ws = model.ws
	lambdalow = model.lambdalow
	lambdahigh = model.lambdahigh
	
	li0=value_locate(w,lambdalow/1d3) ;what array index does lambdalow correspond to?
	li1=value_locate(w, lambdahigh/1d3) ;what array index does lambdahigh correspond to?


	
	norders = n_elements(model.orders)
	
	;make the warray, although this is not *required* for extraction now
	warray=MAKE_ARRAY(2048, norders, /Double, Value=!Values.F_NAN)
	
	;xs,ys are the xy locations of the order sample points in mm
	;xs,ys_pix are the xy locations in pixels (with -1024 to 1023, not 0 to 2047)
	pixel_size = 18d-3
	xs_pix = xs / pixel_size
	ys_pix = ys / pixel_size
	
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;PROCESS ECHELLOGRAM INFO
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;_base arrays are interpolated versions of the Barnes input to match what is required
	;for the array
	buffer = 0d
	upfactor = 1d
	
	fs_base = dblarr((2048 + buffer)*upfactor,norders)
	fs_base2 = fs_base
	
	xs_base = ((dindgen((2048+buffer)*upfactor) - (1024.+buffer/2d)*upfactor)/upfactor)
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
		;ws_base[*,i] = interpol(ws[*,i],xs_pix[*,i],xs_base,/spline)
		;same for left and right limits of each pixel
		;ws_base_left[*,i] = interpol(ws[*,i],xs_pix[*,i],xs_base_left,/spline)
		;ws_base_right[*,i] = interpol(ws[*,i],xs_pix[*,i],xs_base_right,/spline)
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
		
		;find the y for each x based on the x/y map alone
		;ys_base[*,i] = interpol(ys_pix[*,i],xs_pix[*,i],xs_base,/spline)
		;replicate the x base array for convolving and plotting
		xs_base_2d[*,i] = xs_base
	endfor
	
	
	
	;;; COMPARISON EXTRACTED SPECTRUM
	a = mrdfits(comparespec)
	cwl = reform(a[0,*])
	cfl = reform(a[2,*])
	totalfluxes_compare = dblarr(norders)
	iis = where(~finite(cfl))
	cfl2 = cfl
	cfl2[iis] = 0d
	
	for i=0, norders-1 do begin
		is = where( (cwl gt min(ws_base[*,i])/1d3) and (cwl lt max(ws_base[*,i])/1d3),nis)
		if nis eq 0 then stop
		intwl = cwl[is]
		intfl = cfl2[is]
		tot1 = int_tabulated(intwl,intfl,/double,/sort);total(intwl * intfl,/nan)
		;tot1 = total(intfl,/double)
		totalfluxes_compare[i] = tot1
	endfor
	totalflux_compare = total(totalfluxes_compare)

	;f2 = interpol(f,w,cwl)
	;iis = where(~finite(cfl))
	;f2[iis] = 0d
	
	totalfluxes = dblarr(norders)
	for i=0, norders-1 do begin
		print,'i ',i
		is = where( (w gt min(ws_base[*,i])/1d3) and (w lt max(ws_base[*,i])/1d3),nis)
		if nis eq 0 then stop
		intwl = w[is]
		intfl = f[is]
		tot1 = int_tabulated(intwl,intfl,/double,/sort);total(intwl * intfl,/nan)
		totalfluxes[i] = tot1
	endfor
	totalflux = total(totalfluxes)
	
	
	
	
	print,'total known flux: ',totalflux
	print,'total extracted flux: ',totalflux_compare
	
	
	
	stop
	
end
		
	
	
