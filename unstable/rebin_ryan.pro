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
