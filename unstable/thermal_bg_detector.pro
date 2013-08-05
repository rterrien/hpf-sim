;; thermal_bg_detector
;; HPF DETECTOR SIMULATOR
;; This routine forms the thermal background array
;; 
;; CALLS
;; 
;; 
;; DIRECTLY MODIFIES
;;
;;
;; CALLED BY
;; expose.pro
;;
;; NOTES
;; 
;;
;; parameters:
;; det - detector structure
;; diag_ouput - string for diagnostics

function b_lambda_ergs, lambda_mic, temp
	;lambda in microns
	h = 6.626d-27 ;erg s
	k = 1.381d-16 ;erg / k
	c = 2.998d10  ;cm / s

	lambda_cm = lambda_mic * 1d-4
	
	b = 2d * h * c^2. / lambda_cm^5d / (exp(h*c / (lambda_cm * k * temp)) - 1d)
	
	;want output in / micron
	b = b * 1d / 1d4 ;cm / microns
	
	return,b
end


pro thermal_bg_detector, det, diag_output = diag_output

	diag_output = diag_output + string(13B) + '******************THERMAL_BG_DETECTOR.PRO OUTPUT **********************' + string(13B)
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;DEFINE PHYSICAL AND DETECTOR CONSTANTS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


	h = 6.626d-27 ;erg s
	k = 1.381d-16 ;erg / k
	c = 2.998d10  ;cm / s
	qe_wl = [.8,1.,1.23,1.5,2.0,3.5,4.4]
	wls = dindgen(1000)/1000d*(4.4 - .8) + .8
	pixel_area_cm = (18d-4)^2.
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;CHOOSE QE FUNCTION
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	if det.qe_flag ne 0b then begin
		qe_function_lookup,det,qearr
	endif else begin
		qearr = transpose([[0,2],[1.,1.]])
	endelse
	
	qes = interpol(qearr[1,*],qearr[0,*],wls)
	
	diag_output = diag_output + 'QE Function Choice: '+det.qe_loc + string(13B)


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;CONSTRUCT BLACKBODIES
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	bb_200_ergs = b_lambda_ergs(wls,det.bg_temp)
	
	diag_output = diag_output + 'Total BB ergs: '+string(int_tabulated(wls,bb_200_ergs,/double),format='(G20.10)') + string(13B)
	bb_filter_ergs = b_lambda_ergs(wls,80d)
	diag_output = diag_output + 'Filter is at 80K, total ergs from that: '+string(int_tabulated(wls,bb_filter_ergs,/double),format='(G20.10)') + string(13B)
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;ACCOUNT FOR FILTER LOSS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	if det.filter_flag then begin
		fes = interpol(det.filter_spec[1,*],det.filter_spec[0,*],wls)
		bb_200_ergs = bb_200_ergs * fes
		diag_output = diag_output + 'Filter loss mean: '+string(mean(fes),format='(D10.4)') + string(13B)
		diag_output = diag_output + 'Filter loss stddev: '+string(stddev(fes),format='(D10.4)') + string(13B)
	endif else begin
		diag_output = diag_output + 'No filter loss included ' + string(13B)
	endelse
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;ACCOUNT FOR QE
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	bb_200_ergs_qe = bb_200_ergs * qes
	
	bb_filter_ergs_qe = bb_filter_ergs * qes

	diag_output = diag_output + 'Total BB * QE ergs: ' + string(int_tabulated(wls,bb_200_ergs_qe,/double),format='(G20.10)') + string(13B)
	diag_output = diag_output + 'Total Filter * QE ergs: ' + string(int_tabulated(wls,bb_filter_ergs_qe,/double),format='(G20.10)') + string(13B)
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;CONVERT TO PHOTONS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	bb_200_photons_qe = bb_200_ergs_qe * wls * 1d-4 / (h * c)
	
	bb_filter_photons_qe = bb_filter_ergs_qe * wls * 1d-4 / (h * c)
	
	diag_output = diag_output + 'Total BB * QE photons: ' + string(int_tabulated(wls,bb_200_photons_qe,/double),format='(G20.10)') + string(13B)
	diag_output = diag_output + 'Total Filter * QE photons: ' + string(int_tabulated(wls,bb_filter_photons_qe,/double),format='(G20.10)') + string(13B)

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;SUM ALL PHOTONS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	bb_200_total_photons = int_tabulated(wls,bb_200_photons_qe,/double)
	
	bb_filter_total_photons = int_tabulated(wls,bb_filter_photons_qe,/double)
	
	;get total for window
	;3.05 steradians
	window_photons = bb_200_total_photons * 1.456 * pixel_area_cm * .1d ;emissivity
	
	;get total for grating
	grating_photons = bb_200_total_photons * .3 * pixel_area_cm * 1d ;emissivity
	
	;get total for filter
	filter_photons = bb_filter_total_photons * 4.48d * pixel_area_cm
	
	diag_output = diag_output + 'NOTE: using window at 1.456str, .1 emissivity' + string(13B)
	diag_output = diag_output + 'NOTE: using grating at .3str, 1 emissivity' + string(13B)
	diag_output = diag_output + 'NOTE: using filter at 4.48str, 1 emissivity' + string(13B)


	total_photons = window_photons + grating_photons + filter_photons
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;STORE RESULT
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	bgarr = dblarr(det.n,det.n)
	
	bgarr[*] = total_photons
	
	det.bgarr = bgarr
	
	diag_output = diag_output + 'Total photons: ' + string(total_photons,format='(G20.10)') + string(13B)
	diag_output = diag_output + '**********************************END THERMAL_BG_DETECTOR.PRO OUTPUT *******************************'
	
end
	
	
	
	
	
	
