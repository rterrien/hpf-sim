;+
; NAME:
;  hpf_mask_rv
;
; PURPOSE:
;
;  Perform the mask CCF RV analysis on a spectrum
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;	out = hpf_mask_rv(spec_params, mask_params, spec_struct)
;
; INPUTS:
;
;	spec_params: A structure as defined in hpf_initialize_spec_params
;
;	mask_params: A structure as defined in hpf_initialize_mask_params
;
;	spec_struct: A spectrum structure with {wl:dblarr, fl:dblarr}
;	
; OUTPUTS:
;	
;	A structure with the measured RVs and CCF peak fitting results
;
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-


function hpf_mask_rv, spec_params, mask_params, spec_struct

	out = {rv1_out:!values.f_nan,$
		rv1_err:!values.f_nan,$
		rv_out:!values.f_nan,$
		rv_err:!values.f_nan,$
		ccf1_res:dblarr(mask_params.ccf1_nterms),$
		ccf1_res_err:dblarr(mask_params.ccf1_nterms),$
		ccf2_res:dblarr(mask_params.ccf2_nterms),$
		ccf2_res_err:dblarr(mask_params.ccf2_nterms),$
		nan_reject:0d,$
		edge_reject:0d,$
		tell_reject:0d,$
		n_maskpts_orig:0d,$
		n_maskpts_final:0d}
		

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;DEFINE BASIC INPUTS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	diag_output = ''

	c=299792458.d ;m/s
	;nel = n_elements(spec_struct.wl)
	
	mask_lp = mask_params.lp
	mask_rp = mask_params.rp
	mask_we = mask_params.weights
	
	star_ind = where(spec_params.type eq 'STAR',nstar)
	if nstar ne 1 then message,'NO STAR or MORE THAN 1 STAR'
	
	wl1 = double(reform(spec_struct.wl))
	fl1 = double(reform(spec_struct.fl[star_ind[0],*]))
	
	wl = wl1
	fl = fl1 
	nel = n_elements(wl1)
	
	

	
;;	nel2 = 10.*nel1 ;multiply this for upsampling
;;	neww = dindgen(nel2)/double(nel2) * (maxl - minl) + minl
;;	;newf = interpol(fl1,wl1,neww)
;;	newf = hermite(wl1,fl1,neww)
;;	wl = neww
;;	fl = newf
;;	nel = nel2

	
	out.n_maskpts_orig = n_elements(mask_lp)
	

	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;REMOVE NAN MASK POINTS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;remove mask points that are near NAN groups so we don't have any which are partially included
	if mask_params.excl_nan then begin
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
	
	out.nan_reject = nnan_tot
	;diag_output = diag_output + string(13B) + 'Total NaN-rejected Mask Points: '+string(nnan_tot,format='(I6)')

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;REMOVE MASK POINTS NEAR EDGES
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;do not use any points within 60 pixels of an edge either
	llim = wl[150]
	ulim = wl[-150]
	mcs = (mask_lp + mask_rp)/2.
	nori = n_elements(mcs)
	nedg = where(mcs gt llim and mcs lt ulim,nne)
	mask_lp = mask_lp[nedg]
	mask_rp = mask_rp[nedg]
	;print,'n orig ',nori
	;print,'n edge ',nne
	;diag_output = diag_output + string(13B) + 'Edge-Rejected Mask Points: '+string(nne,format='(I6)')
	
	out.edge_reject = nne

	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;REMOVE MASK POINTS NEAR TELLURICS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;IMPLEMENT THIS LATER RCT 12-5-13 
	;NEED A PROCESSED TELLURIC SPECTRUM???
	if keyword_set(mask_params.excl_tell) then begin
		tellspec = mrdfits(mask_params.telluric_spec)
		tellwl = reform(tellspec[0,*])
		tellfl = reform(tellspec[1,*])
		
		bad = where(tellfl lt mask_params.tell_excl_level, nbad) ;originally .98
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
		
		;diag_output = diag_output + string(13B) + string(tot_tell,format='(I4)')
		
		out.tell_reject = tot_tell
		
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
	
;;	;try something more detailed
;;	xii = dindgen(n_elements(wl))
;;	xii_l = xii - 0.5d
;;	xii_r = xii + 0.5d
;;	lp2 = interpol(wl,xii,xii_l,/quadratic)
;;	rp2 = interpol(wl,xii,xii_r,/quadratic)
;;	xsi2 = rp2 - lp2
;;	xsi = xsi2
;;	rp = rp2
;;	lp = lp2
;;	;didnt' help

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;ROUND 1 OF CCF FITTING
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;this will do a broad, shallow CCF (in terms of rv shifts) to find the central shift
	;to make a symmetric CCF for the second round of CCF fitting
	nvel = mask_params.ccf1_nvel ;number of velocities to loop over
	velrange = mask_params.ccf1_velrange ;overall range
	vel_shifts = dindgen(nvel)/(nvel-1.)*velrange - (velrange/2.d)
	result = dblarr(nvel) ;holds CCF result
	result_weights = dblarr(nvel) ;
	
	;diag_output = diag_output + string(13B) + '***** CCF ROUND 1 *****' + string(13B) + 'N velocities: ' + string(nvel,format='(I6)') + string(13B) + 'Velocity Range: ' + string(velrange,format='(I7)') + string(13B)

	;loop over velocity shifts, construct the CCF array
	for i=0, nvel-1 do begin
		;velocity shift the mask points
		beta = vel_shifts[i]/c
		case strupcase(strcompress(spec_params.rv_type,/remove_all)) of
			'RELATIVISTIC':lam_shift = double(sqrt( (1d + beta) / (1d - beta) ))
			'NEWTONIAN':lam_shift = (1d + vel_shifts[i]/c)
		endcase
		
		mask_lp_shifted = reform(double(mask_lp*lam_shift))
		mask_rp_shifted = reform(double(mask_rp*lam_shift))
		
		;find which points are inside the mask points
		l_locs = value_locate(lp, mask_lp_shifted)
		r_locs = value_locate(rp, mask_rp_shifted)+1
		
		;check if there are any locations where there is not at least one full pixel
		;btwn left and right
	;;	if min(r_locs - l_locs) lt 2 then begin
;;			print,'mask regions contain < 1 full pixel '
;;			stop
;;		endif
		
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
	ccf_fit1 = mpfitpeak(vel_shifts,result,a_ccf,nterms=mask_params.ccf1_nterms,sigma=a_ccf_errors,perror=a_ccf_errors);,estimates=[-50.,0.,5000.,200.])
	
	center1 = a_ccf[1]
	
	out.ccf1_res = a_ccf
	out.ccf1_res_err = a_ccf_errors
	out.rv1_out = a_ccf[1]
	out.rv1_err = a_ccf_errors[1]
	
	diag_output = diag_output + string(13B) + 'Center 1: ' + string(center1,format='(D+10.3)') + string(13B) + 'Peak Width: ' + string(a_ccf_errors[1],format='(D+10.3)')

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;ROUND 2 OF CCF FITTING
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;Now a finer scan of CCF
	velrange2 = mask_params.ccf2_velrange ; was 24000
	nvel2 = mask_params.ccf2_nvel
	vel_shifts2 = (dindgen(nvel2)/(nvel2 - 1d) - .5d)*velrange2 + center1
	result2 = dblarr(nvel2)
	result_weights2 = dblarr(nvel2)
	
	;diag_output = diag_output + string(13B) + '***** CCF ROUND 2 *****' + string(13B) + 'N velocities: ' + string(nvel2,format='(I6)') + string(13B) + 'Velocity Range: ' + string(velrange2,format='(I7)') + string(13B)


	for i=0, nvel2-1 do begin
		;velocity shift the mask points
		beta = vel_shifts2[i]/c
		case strupcase(strcompress(spec_params.rv_type,/remove_all)) of
			'RELATIVISTIC':lam_shift = double(sqrt( (1d + beta) / (1d - beta) ))
			'NEWTONIAN':lam_shift = (1d + vel_shifts[i]/c)
		endcase
		
		;find pixels inside the mask points
		mask_lp_shifted = reform(double(mask_lp*lam_shift))
		mask_rp_shifted = reform(double(mask_rp*lam_shift))
		l_locs = value_locate(lp, mask_lp_shifted)
		r_locs = value_locate(rp, mask_rp_shifted)+1
		;check if there are any locations where there is not at least one full pixel
		;btwn left and right
		if min(r_locs - l_locs) lt 1 then begin
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

	ccf_fit2 = mpfitpeak(vel_shifts2,result2,b_ccf,nterms=mask_params.ccf2_nterms,sigma=b_ccf_errors,perror=b_ccf_errors)
	
	center2 = b_ccf[1]
	
	out.ccf2_res = b_ccf
	out.ccf2_res_err = b_ccf_errors
	out.rv_out = b_ccf[1]
	out.rv_err = b_ccf_errors[1]
	
	
	;diag_output = diag_output + string(13B) + 'Center 2 (RV Result): ' + string(center2,format='(D+10.3)') + string(13B) + 'Peak Width: ' + string(b_ccf_errors[1],format='(D+10.3)')

	return,out

end
	