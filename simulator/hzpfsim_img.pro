;; hzpfsim_img
;; HPF DETECTOR SIMULATOR
;; This is the code to create the fluence array
;; 
;; CALLS
;; 
;; 
;; DIRECTLY MODIFIES
;;
;; CALLED BY
;; script_all_XX
;;
;; NOTES
;; 
;;
;; parameters:
;; vrad - RV in m/s
;; outfile - output file for fluence (.fits)
;; tellfile - telluric spectrum .sav to contaminate
;; specfile - fits file with input spectrum (if none provided it defaults to bt_34)

pro hzpfsim_img, vrad, outfile, tellfile=tellfile, specfile = specfile, diagfile = diagfile, calfile = calfile

	diag = diagfile ne !null
	if diag then begin
		openu,diaglun,diagfile,/get_lun,/append
		printf,diaglun,string(13B)+string(13B)+'############# BEGIN HZPFSIM_IMG.PRO ##############'
		printf,diaglun,'Run at: ',systime()
	endif


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;DEFINE PHYSICAL AND DETECTOR CONSTANTS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	c = 2.99792458d8	; m/s
	h = 6.626d-34 ;J s
	res = 50000.d
	
	gap=[0.0, 2.79, 2.67, 2.55, 2.45, 2.35, 2.25, 2.16, 2.08, 2.00, 1.93, 1.86, 1.79, 1.73, 1.67, 1.62, 1.56] * 1000 / 18. ;;pixels, center to center
	ordernum = LINDGEN(17)+46
	lambdalow=[13173, 12893, 12624, 12366, 12119, 11881, 11653, 11433, 11221, 11017, 10821, 10631, 10448, 10270, 10099, 9934, 9773]/1d4
	lambdahigh=[13390, 13105, 12832, 12570, 12319, 12078, 11845, 11622, 11407, 11199, 10999, 10806, 10620, 10440, 10266, 10098, 9935]/1d4
	
	minl = min(lambdalow)-.05d ;upper and lower limits with a little extra
	maxl = max(lambdahigh)+.05d

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;READ SPECTRUM
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	if keyword_set(specfile) then begin
		spec = mrdfits(specfile,0,head)
		spec = transpose(spec)
	endif else begin
		spec = readfits('bt_34.fits',head)
		spec = transpose(spec)
		specfile = 'bt_34.fits'
	endelse
	wmod   = reform(spec[*,0])	; assume wavelength is in microns
	fmod   = reform(spec[*,1])	; assume flux is in Watts/m2/micron
	;dwmod  = 0.015d-4		; sampling of Hauschildt models in microns
	;sort
	w      = wmod[sort(wmod)]
	f      = fmod[sort(wmod)]
	
	if diag then printf,diaglun,string(13B)+'Spec File: '+specfile
	

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;CUT DOWN AND UPSAMPLE
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	g = where(w ge minl - .05d and w le maxl + .05d,ng)
	w = double(w[g])
	f = double(f[g])

	;upsample by a factor of three (ensures even/sufficient pixel sampling)
	nel1 = n_elements(w)
	nel2 = 12.*nel1 ;multiply this for upsampling
	neww = dindgen(nel2)/double(nel2) * (maxl - minl) + minl
	newf = interpol(f,w,neww)
	
	w = neww
	f = newf
	ng = n_elements(w)
	
	usf = string(nel2 / nel1,format='(I3)')
	if diag then printf,diaglun,string(13B)+'Upsample Factor: '+usf


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;VELOCITY SHIFT
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	beta=double(vrad/c)
	w=double(w*SQRT( (1d + beta) / (1d - beta) ))
	
	if diag then printf,diaglun,string(13B)+'Vrad: '+string(vrad,format='(D+8.1)')

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;SCALE FLUX
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	jzp     = 3.31d-9	; Vega flux at J  (W/m2/micron)
	jvega   = 0.025d	; Vega mag at J
	deltaj  = 0.26d		; FWHM of J in microns
	isowj   = 1.215d	; Isophotal wavelength of J in microns
	rtel    = 7./2.0d	; radius of telescope in m
	telarea = double(!pi*rtel^2)    ; telescope area in m^2
	tput    = 0.04d		; total instrument throughput 

	;print, 'Scaling spectrum ...'
	
	jmag = 10.0
	
	jmin   = isowj - deltaj/2.0d
	jmax   = isowj + deltaj/2.0d
	jndx   = where((w ge jmin) and (w le jmax))
	fmean  = int_tabulated(w[jndx],f[jndx])/deltaj
	jscale = (jzp/fmean)*10.0d^(-0.4d*(jmag - jvega))
	f      = temporary(f)*jscale		; scale spectrum to input J mag
	if diag then begin
		printf,diaglun, 'Mean Model Flux in J band  = ', fmean, ' Watts/m2/micron'
		printf,diaglun, 'Mean Scaled Flux in J band = ', fmean*jscale, ' Watts/m2/micron'
	endif

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;TELnLURIC CONTAMINATE IF NECESSARY
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;do this while spectrum is still in energy
	telluric = 0
	if keyword_set(tellfile) then begin
		tellspec = mrdfits(tellfile)
		tellw = reform(tellspec[0,*])
		tellf = reform(tellspec[1,*])
		contam = interpol(tellf,tellw,w)
		f = f * contam
		telluric = 1
		if diag then printf,diaglun,'Telluric Contaminated File: '+tellfile
	endif else begin
		if diag then printf,diaglun,'Telluric Contaminated File: None'
	endelse
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;CONVERT TO PHOTONS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;photons/s/m2/micron
	p = w * f / (h*c) * 1d-6

	pout = temporary(p)*telarea*tput ;photons/s/micron

	fold = f
	f = pout
	
	if diag then printf,diaglun,'Avg Photons',mean(pout,/nan)
	
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;CONVERT TO PHOTONS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;Calibration Unit - added 8/2/2012, CFB
	;;Cal spectrum should be have wavelength in units of microns and flux
	;;in units of Watts/micron
	calflag=0
	IF SIZE(calfile,/Type) NE 0 THEN BEGIN
		cal = MRDFITS(calfile,0,calhead)
		calflag=1
		calwmod=REFORM(cal[0,*])
		calfmod=REFORM(cal[1,*])
		calw=calwmod[SORT(calwmod)]
		calf=calfmod[SORT(calwmod)]
		calg=WHERE(calw GE minl - 0.05d AND calw LE maxl + 0.05d,calng)
		IF calng EQ 0 THEN MESSAGE,'ERROR: Cal unit out of bounds!'
		calw=DOUBLE(calw[calg])
		calf=DOUBLE(calf[calg])
		calneww=DINDGEN(calng * 10L)/DOUBLE(calng * 10L) * (maxl - minl) + minl
		calnewf = INTERPOL(calf,calw,calneww)
		calw=calneww
		calf=calnewf
		calng=N_ELEMENTS(calw)
		;;Convert to photons / s / micron
		calp = calw * calf / (h*c) * 1d-6
		
		;;Attenuate (integrating sphere)
		calpout = calp * 1d-6
		
		;;Attenuate (spectrograph)
		calpout = calpout * tput
		
		;;Attenuate (ND filter)
		calpout = calpout * 1d-1
		
		;;Interpolate to wavelength vector of stellar spectrum
		
		;calf = INTERPOL(calpout,calw,w)
		calf = calpout
		;also calw
		
		if diag then printf,diaglun,'Cal File: ',calfile
	ENDIF

	if diag then printf,diaglun,'Cal Fiber Active: ',calflag
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;SMOOTH
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;smooth the entire array once for each order
	;this will produce n_orders copies of the spectrum, each smoothed to the resolution at the middle of order_n
	no = n_elements(lambdalow) ;number of orders
	fsms = dblarr(no,ng) ;smoothed array copies
	calfsms = DBLARR(no,n_elements(calf))   ;;added 8/2/2012, CFB
	;loop over orders
	for i=0, no-1 do begin
		if diag then printf,diaglun,'smoothing order ',i
		mwl = double(mean([lambdalow[i],lambdahigh[i]]))
		smooth_res_ryan,res,mwl,w,f,fstemp,diag_output = diag_output
		if diag then printf,diaglun,'Stellar Smoothing Output: '
		if diag then printf,diaglun,diag_output
		fsms[i,*] = double(fstemp)
		IF calflag THEN BEGIN  ;;added 8/2/2012, CFB
        	smooth_res_ryan,res,mwl,calw,calf,calfstemp,diag_output = diag_output
        	calfsms[i,*]=double(calfstemp)
        ENDIF
        if diag then printf,diaglun,'Cal Smoothing Output: '
        if diag then printf,diaglun,diag_output
	endfor

	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;PARAMETERS FOR MAKING THE SPEC IMAGE
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	nw      = n_elements(w)
	nypix = 22.
	nxarr = 2048L
	
	specimg=dblarr(nxarr, nxarr)
	wlimg = make_array(nxarr,nxarr,/double,value=!values.f_nan)
	y0pos= ROUND((nxarr - ROUND(TOTAL(gap)))/2.) ;what does this do?
	y0poscal=y0pos + nypix + 10 ;;origin of the cal fiber
	cgap=ROUND(TOTAL(gap, /Cumulative))
	norders=N_ELEMENTS(lambdalow)
	warray=MAKE_ARRAY(nxarr, norders, /Double, Value=!Values.F_NAN)
	li0=value_locate(w,lambdalow)
	li1=value_locate(w, lambdahigh)

	cal_li0 = value_locate(calw,lambdalow)
	cal_li1 = value_locate(calw,lambdahigh)

	
	;;;FOR TESTING
	;testout = dblarr(no,2,2048)

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;PROJECT INTO SPEC IMAGE (AND WLIMG)
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	for order=0, norders-1 do begin
		subimg = dblarr(nxarr,nypix) ;the entire order image
		subfls = (reform(fsms[order,*]))[li0[order]:li1[order]] ;fluxes for the order
		subwls = reform(w[li0[order]:li1[order]]) ;wls for the order
		
		;;lasers = [10830.,10640.,10530.]/1d4
;;		for k=0,2 do begin	
;;			diffs = abs(subwls - lasers[k])
;;			if min(diffs) gt 1d-5 then continue
;;			ind = where(diffs eq min(diffs))
;;			subfls[ind-5:ind+5] = 10.*max(subfls)
;;			print,'tripped ',k
;;			;stop
;;		endfor
		
		;get photons/s/m2 from photons/s/m2/micron
		;mp = (subwls + subwls[1:*])/2. ;midpoints
		;binsi = mp[1:*] - mp 
		;binsi = [binsi[0],binsi,binsi[-1]] ;bin sizes
		;subflsold = subfls
		;subfls = subfls; * binsi
		
		;rebin the sub-spectrum to a given pixel sampling
		rebin_ryan, 3.d, res, subwls, subfls, wr, fr, min(subwls,/nan), max(subwls,/nan)
		if diag then printf,diaglun,'n pixels for order ',n_elements(wr)

		;take 5 pixels off edges for safety
		wr = wr[4:-6]
		fr = fr[4:-6]
		
		;then just take first 2048 pixels to put into detector
		wsub = wr[0:2047]
		fsub = fr[0:2047]
		;;;;testing
		;testout[order,0,*] = wsub
		;testout[order,1,*] = fsub
		;;;;end testing
		farr = dblarr(2048,nypix)
		warr = dblarr(2048,nypix)
		for i=0, nypix-1 do begin
			farr[*,i] = fsub / double(nypix)
			warr[*,i] = wsub
		endfor
		specimg[0,-(y0pos + cgap[order] - nypix/2)] = farr
		wlimg[0,-(y0pos + cgap[order] - nypix/2)] = warr
		warray[0,order] = wsub
		
		IF calflag THEN BEGIN ;;added 8/2/2012, CFB
			subcalwls = reform(calw[cal_li0[order]:cal_li1[order]])
			subcalfls=(reform(calfsms[order,*]))[cal_li0[order]:cal_li1[order]]
			rebin_ryan, 3.d, res, subcalwls, subcalfls, calwr, calfr, min(subwls,/nan), max(subwls,/nan)
			calfr = calfr[4:-6]
			calfsub = calfr[0:2047]
			calfarr = dblarr(2048,nypix)
			for i=0, nypix-1 do calfarr[*,i] = calfsub / double(nypix)
			specimg[0,-(y0poscal + cgap[order] - nypix/2)] = calfarr
        ENDIF

	endfor
	

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;SAVE HEADER AND RESULT
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	h = []
	sxaddpar,h,'vrad',vrad
	sxaddpar,h,'tellcontam',telluric
	
	;writefits,outfile,specimg
	mwrfits,specimg,outfile,h,/create
	nlen = strlen(outfile)
	pref = strmid(outfile,0,nlen-5)
	writefits,pref+'_wlimg.fits',wlimg
	wfile = outfile+'_wl.sav'
	save,warray,filename=wfile
	
	if diag then begin
		printf,diaglun,string(13B)+string(13B)+'############# END HZPFSIM_IMG.PRO ##############'
		free_lun,diaglun
	endif

	
	;;;;testing
	;testout is order, 0/1 (w/f), pixel
;;	wl = reform(testout[-1,0,*])
;;	mp = (wl + wl[1:*])/2.
;;	xsi = mp[1:*] - mp
;;	xsi = [xsi[0],xsi,xsi[-1]]
;;	fl = reform(testout[-1,1,*])/xsi
;;	
;;	
;;	for i=15, 0, -1 do begin
;;		wl1 = reform(testout[i,0,*])
;;		mp1 = (wl1 + wl1[1:*])/2.
;;		xsi1 = mp1[1:*] - mp1
;;		xsi1 = [xsi1[0],xsi1,xsi1[-1]]
;;		fl1 = reform(testout[i,1,*])/xsi1
;;		a = where(wl1 gt wl[-1],na)
;;		if na lt 1 then stop
;;		wl = [wl, !values.f_nan,wl1[a]]
;;		fl = [fl, !values.f_nan,fl1[a]]
;;	endfor

;f1 = file_basename(outfile,'.fits')

;file = 'spt73/'+f1+'_spec.fits'

;mwrfits,[1#wl,1#fl],file,h,/create
;openw,1,file
;for i=0, n_elements(wl)-1 do begin
;	out = string(wl[i],format='(D15.8)')+' '+string(fl[i],format='(D25.4)')
;	printf,1,out
;endfor
;close,1

end

pro smooth_res_ryan, res, midwl, wl, fl, fsm, diag_output = diag_output
	;routine to do a simple convolution
	;with a gaussian of a fixed wavelength width
	;res = resolution
	;midwl = wavelength where resolution is converted to pixels
	;fsm = output flux
	;wl = wavelength array with even wavelength sampling
	;fl = flux array
	
	
	;figure out which pixel has midwl
	a = abs(wl - midwl)
	al = where(a eq min(a))
	
	;what wavelength interval is this pixel
	pixwl = abs(wl[al] - wl[al-1])
	
	;what should be the wavelength interval of the gaussian (FWHM)
	dl = midwl/res
	;in pixels
	dlpix = dl/pixwl
	;convert to sigma from FWHM
	sigmap = dlpix / 2.35
	
	;diag
	diag_output = 'midwl: '+string(midwl,format='(D7.5)')+string(13B)+'pixwl: '+string(pixwl,format='(G12.5)')+string(13B)+'dl: '+string(dl,format='(G12.5)')+string(13B)+'dlpix: '+string(dlpix,format='(G12.5)')

	
	;form kernel with 20X FWHM
	kernel_width = round(dlpix * 20)
	kernel_center = floor(kernel_width / 2)
	if kernel_width mod 2 eq 0 then kernel_width += 1
	;print, 'kernel width ',kernel_width
	gx = dindgen(kernel_width)
	kernel = gaussian(gx,[1.,kernel_center,sigmap])
	cfl = convol(fl,kernel,/normalize,/edge_truncate)
	fsm = cfl
	
	diag_output = diag_output+string(13B)+'Kernel Width: '+string(kernel_width,format='(I5)')
	
end

pro rebin_ryan, pix, res, wl, fl, wlr, flr, minl, maxl
	;rebin and conserve flux as found through riemann sum
	;pix = sampling
	;res = resolution
	;wlr/flr = output
	;minl/maxl = wavelength limits for output arrays
	
	;find pixel wavelength
	pixwl = mean([minl,maxl]) / (res * pix)
	;pixwl = 1.15d / (res * pix)
	npn = round((maxl - minl)/pixwl) ;number of pixels in new array
	newwl = dindgen(npn)/double(npn)*(maxl - minl) + minl ;output wls
	newfl = dindgen(npn) ;output fls

	np = n_elements(wl) ;number of pixels in old array
	
	;find midpoints for the input wavlength array
	mp = (wl + wl[1:*])/2. ;midpoints
	xsi = mp[1:*] - mp 
	xsi = [xsi[0],xsi,xsi[-1]] ;bin sizes
	rp = wl + xsi/2. ;right points
	lp = wl - xsi/2. ;left points
	
	;and for new array
	mpn = (newwl + newwl[1:*])/2. ;midpoints
	xsin = mpn[1:*] - mpn 
	xsin = [xsin[0],xsin,xsin[-1]] ;bin sizes
	rpn = newwl + xsin/2. ;right points
	lpn = newwl - xsin/2. ;left points
	
	;loop through all new pixels
	for i=0, npn-1 do begin
		ii = where(lp ge lpn[i] and rp le rpn[i],nii) ;all old pixels entirely within new one
		li = min(ii) - 1 ;leftmost old pixel
		ri = max(ii) + 1 ;rightmost old pixel
		tot = 0.d ;total for accumulating 
		if nii ne 0 then begin
			;add net flux (fl/wl * wl) from pixels completely within the limits
			tot += total(fl[ii]*xsi[ii])
			if li ge 0 then begin
				;add left fractional flux
				lf = (rp[li] - lpn[i])/xsi[li]
				tot += lf * fl[li] * xsi[li]
			endif
			if ri le np - 1 then begin
				;add right fractional flux
				rf = (rpn[i] - lp[ri])/xsi[ri]
				tot += rf * fl[ri] * xsi[ri]
			endif
		endif else begin
			;if new array pixles are smaller, or similar size
			;i haven't coded for this
			print,"NOT DEALT WITH YET"
			stop
		endelse
		newfl[i] = tot
	endfor
	wlr = newwl
	flr = newfl

end

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
