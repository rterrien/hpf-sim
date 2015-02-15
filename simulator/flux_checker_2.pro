pro flux_checker_2

	; comparespec
	
	comparespec = '~/scratch/simulator/data/sp76/test1fl000_spec.fits'
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; OPTICAL PARAMETERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	optical_params = hpf_initialize_optical_params()
	print,'MEMORY AFTER LOAD SPEC'
	print,memory()/1d6
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;; PROJ PARAMETERS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	proj_params = hpf_process_optical_model(optical_params)
	print,'MEMORY AFTER PROC OPT'
	print,memory()/1d6

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;GET ECHELLOGRAM INFO
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	norders = (size(proj_params.xs))[2]

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
	

	;;; Put model into flux
	exptime = 5d
	;comparespec = '~/scratch/simulator/data/sp30/test1fl000_spec.fits'
	
	;;;FLUX LOAD
	c = 2.99792458d8	; m/s
	h = 6.626d-34 ;J s
	res = 50000.d
	pixel_sampling = 3d
	

	minl = .8d
	maxl = 1.4d

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;READ SPECTRUM
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	;spec = readfits('support/bt_34_extended.fits',head)
	;spec = transpose(spec)
	;cal = 0
	
	spec = mrdfits('~/work/spectra/ffp_30_80.fits')
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
		f = p * exptime * 1d-2 ;; filter amount must be changed manually
	endif
	
	;convolve? test
	dw = median(w[1:*] - w)
	fwhm = median(w) / 50000d
	fwhm_pix = fwhm/dw
	nkern = round(10d * fwhm_pix)
	if nkern mod 2 eq 0 then nkern += 1.
	kern = psf_gaussian(npixel = nkern, fwhm = fwhm_pix, /normalize, ndimen=1, /double)
	f1 = f
	f = convol(f1,kern,/normalize,/edge_truncate)
	
	;;; Extracted spec
	
	
	;;; COMPARISON EXTRACTED SPECTRUM
	a = mrdfits(comparespec)
	cwl = reform(a[0,*])
	cfl = reform(a[2,*])
	totalfluxes_compare = dblarr(norders)
	iis = where(~finite(cfl))
	cfl2 = cfl
	cfl2[iis] = 0d
	
	wwmins = dblarr(norders)
	wwmaxs = dblarr(norders)
	
	for i=0, norders-1 do begin
		is = where( (cwl gt min(ws_base[*,i])) and (cwl lt max(ws_base[*,i])) and finite(cfl),nis)
		wwmins[i] = min(cwl[is])
		wwmaxs[i] = max(cwl[is])
		if nis eq 0 then continue
		intwl1 = cwl[is]
		intfl1 = cfl2[is]
		ind = uniq(intwl1)
		intwl = intwl1[ind]
		intfl = intfl1[ind]
		if nis lt 10 then continue
		;tot1 = int_tabulated(intwl,intfl,/double,/sort);total(intwl * intfl,/nan)
		tot1 = tsum(intwl,intfl)
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
		;is = where( (w gt min(ws_base[*,i])) and (w lt max(ws_base[*,i])),nis)
		is = where( (w gt wwmins[i]) and (w lt wwmaxs[i]),nis)
		if nis eq 0 then continue
		intwl = w[is]
		intfl = f[is]
		;tot1 = int_tabulated(intwl,intfl,/double,/sort);total(intwl * intfl,/nan)
		tot1 = tsum(intwl,intfl)
		totalfluxes[i] = tot1
		;if i eq 10 then stop
	endfor
	totalflux = total(totalfluxes)
	
	
	
	
	print,'total known flux: ',totalflux
	print,'total extracted flux: ',totalflux_compare

	
	stop
	
end
