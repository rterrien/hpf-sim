;; mask_cc3
;; HPF DETECTOR SIMULATOR
;; This code does the mask cross-correlation to find RV shifts
;; 
;; CALLS
;; 
;; 
;; DIRECTLY MODIFIES
;;
;; CALLED BY
;; test_mask_sb2
;;
;; NOTES
;; uses mpfitpeak with a gaussian+constant to fit the peaks
;;
;; parameters:
;; wl_in - input wavelengths
;; fl_in - input fluxes
;; rv_out - output in m/s
;; mask_lp_in - mask left points
;; mask_rp_in - mask right points
;; mask_we_in - mask weights
;; rv_err - estimate of the rv error (not used?)
;; nan - keyword to exclude mask points near a nan
;; tellfile - location of telluric .fits file, remove mask points in telluric features

pro mask_cc3, wl_in, fl_in, rv_out, mask_lp_in, mask_rp_in, mask_we_in, rv_err, nan=nan, tellfile = tellfile, diag_output = diag_output

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;DEFINE BASIC INPUTS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	diag_output = ''

	c=299792458.d ;m/s
	wl=double(wl_in)
	fl=double(fl_in)
	nel = n_elements(wl)
	
	mask_lp = mask_lp_in
	mask_rp = mask_rp_in
	mask_we = mask_we_in
	
	diag_output = diag_output + string(13B) + 'Total orig mask points: '+string(n_elements(mask_lp),format='(I6)')
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;REMOVE NAN MASK POINTS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;remove mask points that are near NAN groups so we don't have any which are partially included
	if keyword_set(nan) then begin
		;limit mask points - do not use any within 60 pixels of a nan
		mask_lp_final = []
		mask_rp_final = []
		;test all mask points
		nnan_tot = 0
		for i=0, n_elements(mask_lp)-1 do begin
			;find pixel centers
			lp = mask_lp[i]
			rp = mask_rp[i]
			cp = (lp + rp)/2.
			ds = abs(wl - cp)
			lm = where(ds eq min(ds))
			;look at all pixels within 60
			low = max([0,lm - 60])
			hi = min([n_elements(wl)-1,lm + 60])
			sub = fl[low:hi]
			;if any nearby pixels are NAN, do not include this
			nans = where(sub ne sub,nn)
			nnan_tot = nnan_tot + 1
			if nn eq 0 then begin
				mask_lp_final = [mask_lp_final,lp]
				mask_rp_final = [mask_rp_final,rp]
			endif
		endfor
		mask_lp_old = mask_lp
		mask_rp_old = mask_rp
		mask_lp = mask_lp_final
		mask_rp = mask_rp_final
	endif
	
	diag_output = diag_output + string(13B) + 'Total NaN-rejected Mask Points: '+string(nnan_tot,format='(I6)')

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;REMOVE MASK POINTS NEAR EDGES
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;do not use any points within 60 pixels of an edge either
	llim = wl[60]
	ulim = wl[-60]
	mcs = (mask_lp + mask_rp)/2.
	nori = n_elements(mcs)
	nedg = where(mcs gt llim and mcs lt ulim,nne)
	mask_lp = mask_lp[nedg]
	mask_rp = mask_rp[nedg]
	;print,'n orig ',nori
	;print,'n edge ',nne
	diag_output = diag_output + string(13B) + 'Edge-Rejected Mask Points: '+string(nne,format='(I6)')
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;REMOVE MASK POINTS NEAR TELLURICS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	if keyword_set(tellfile) then begin
		tellspec = mrdfits(tellfile)
		tellwl = reform(tellspec[0,*])
		tellfl = reform(tellspec[1,*])
		bad = where(tellfl lt .98, nbad)
		tellba = dblarr(n_elements(tellfl))
		tellba[*] = 0.d
		tellba[bad] = 1.d
		if nbad eq 0 then stop
		mask_lp_final = []
		mask_rp_final = []
		tot_tell = 0
		for i=0, n_elements(mask_lp)-1 do begin
			ll = where(tellwl ge mask_lp[i] and tellwl le mask_rp[i],nll)
			if nll eq 0 then stop
			low = min(ll)-30
			hi = max(ll)+30
			test = total(tellba[low > 0:hi < n_elements(tellba)-1])
			if test eq 0 then begin
				mask_lp_final = [mask_lp_final,mask_lp[i]]
				mask_rp_final = [mask_rp_final,masK_rp[i]]
				tot_tell += 1
			endif
		endfor
		n1 = n_elements(mask_lp)
		n2 = n_elements(mask_lp_final)
		;print,'n before tellurics ',n1
		;print,'n after tellurics ',n2
		mask_lp = mask_lp_final
		mask_rp = mask_rp_final
		
		diag_output = diag_output + string(13B) + string(tot_tell,format='(I4)')
		
	endif


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;FIND NECESSARY MIDPOINTS AND BINS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	nf = n_elements(mask_lp)
	
	mp = (wl + wl[1:*])/2.d ;midpoints
	xsi = double(mp[1:*] - mp)
	xsi = double([xsi[0],xsi,xsi[-1]]) ;bin sizes
	rp = double(wl + xsi/2.) ;right points
	lp = double(wl - xsi/2.) ;left points

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;ROUND 1 OF CCF FITTING
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;this will do a broad, shallow CCF (in terms of rv shifts) to find the central shift
	;to make a symmetric CCF for the second round of CCF fitting
	nvel = 501.d ;number of velocities to loop over
	velrange = 80000.d ;overall range
	vel_shifts = dindgen(nvel)/(nvel-1.)*velrange - (velrange/2.d)
	result = dblarr(nvel) ;holds CCF result
	result_weights = dblarr(nvel) ;
	
	diag_output = diag_output + string(13B) + '***** CCF ROUND 1 *****' + string(13B) + 'N velocities: ' + string(nvel,format='(I6)') + string(13B) + 'Velocity Range: ' + string(velrange,format='(I7)') + string(13B)

	;loop over velocity shifts, construct the CCF array
	for i=0, nvel-1 do begin
		;velocity shift the mask points
		beta = vel_shifts[i]/c
		lam_shift = double(sqrt( (1d + beta) / (1d - beta) ))
		mask_lp_shifted = reform(double(mask_lp*lam_shift))
		mask_rp_shifted = reform(double(mask_rp*lam_shift))
		
		;find which points are inside the mask points
		l_locs = value_locate(lp, mask_lp_shifted)
		r_locs = value_locate(rp, mask_rp_shifted)+1
		
		;check if there are any locations where there is not at least one full pixel
		;btwn left and right
		if min(r_locs - l_locs) lt 2 then begin
			print,'mask regions contain < 1 full pixel '
			stop
		endif
		
		;add everything together for this velocity shift
		sumarr = dblarr(nel)
		left_frac = (rp[l_locs] - mask_lp_shifted)/xsi[l_locs]
		right_frac = (mask_rp_shifted - lp[r_locs])/xsi[r_locs]
		for j=0, nf-1 do begin
			sumarr[l_locs[j]] += left_frac[j] * mask_we[j]
			sumarr[r_locs[j]] += right_frac[j] * mask_we[j]
			sumarr[l_locs[j]+1:r_locs[j]-1] += mask_we[j]
		endfor
		sum = total( sumarr * xsi * fl ,/nan,/double)
		result[i] = sum
	endfor

	;find the overall CCF peak
	ccf_fit1 = mpfitpeak(vel_shifts,result,a_ccf,nterms=4,sigma=a_ccf_errors,perror=a_ccf_errors);,estimates=[-50.,0.,5000.,200.])
	
	center1 = a_ccf[1]
	
	diag_output = diag_output + string(13B) + 'Center 1: ' + string(center1,format='(D+10.3)') + string(13B) + 'Peak Width: ' + string(a_ccf_errors[1],format='(D+10.3)')

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;ROUND 2 OF CCF FITTING
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;Now a finer scan of CCF
	velrange2 = 20000.d ; was 24000
	vel_shifts2 = (dindgen(501)/500.d - .5d)*velrange2 + center1
	nvel2 = n_elements(vel_shifts2)
	result2 = dblarr(nvel2)
	result_weights2 = dblarr(nvel2)
	
	diag_output = diag_output + string(13B) + '***** CCF ROUND 2 *****' + string(13B) + 'N velocities: ' + string(nvel2,format='(I6)') + string(13B) + 'Velocity Range: ' + string(velrange2,format='(I7)') + string(13B)


	for i=0, nvel2-1 do begin
		;velocity shift the mask points
		beta = vel_shifts2[i]/c
		lam_shift = double(sqrt( (1d + beta) / (1d - beta) ))
		
		;find pixels inside the mask points
		mask_lp_shifted = reform(double(mask_lp*lam_shift))
		mask_rp_shifted = reform(double(mask_rp*lam_shift))
		l_locs = value_locate(lp, mask_lp_shifted)
		r_locs = value_locate(rp, mask_rp_shifted)+1
		;check if there are any locations where there is not at least one full pixel
		;btwn left and right
		if min(r_locs - l_locs) lt 2 then begin
			print,'mask regions contain < 1 full pixel '
			stop
		endif
		
		;add everything together for this velocity shift
		sumarr = dblarr(nel)
		left_frac = (rp[l_locs] - mask_lp_shifted)/xsi[l_locs]
		right_frac = (mask_rp_shifted - lp[r_locs])/xsi[r_locs]
		for j=0, nf-1 do begin
			sumarr[l_locs[j]] += left_frac[j] * mask_we[j]
			sumarr[r_locs[j]] += right_frac[j] * mask_we[j]
			sumarr[l_locs[j]+1:r_locs[j]-1] += mask_we[j]
		endfor
		sum = total( sumarr * xsi * fl ,/nan,/double)
		result2[i] = sum
	endfor

	ccf_fit2 = mpfitpeak(vel_shifts2,result2,b_ccf,nterms=4,sigma=b_ccf_errors,perror=b_ccf_errors)
	
	center2 = b_ccf[1]
	
	diag_output = diag_output + string(13B) + 'Center 2 (RV Result): ' + string(center2,format='(D+10.3)') + string(13B) + 'Peak Width: ' + string(b_ccf_errors[1],format='(D+10.3)')

	
	rv_out = b_ccf[1]
	rv_err = b_ccf_errors[1]
	
end
	