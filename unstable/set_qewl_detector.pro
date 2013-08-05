;; set_qewl_detector
;; HPF DETECTOR SIMULATOR
;; This is the code to set the QE(wavelength) for the detector
;; 
;; CALLS
;; 
;; 
;; DIRECTLY MODIFIES
;;
;; CALLED BY
;; expose
;;
;; NOTES
;; 
;;
;; parameters:
;; det - detector structure
;; qewl - the quantum efficiency as a function of wavelength ([2,x] dblarr)

pro set_qewl_detector, det, qewl


	;find average qe to fill in the interorder regions
	avgqe = mean(qewl[1,*],/nan)
	
	;find indicies of wlarr with wavelength values and those with nans
	ifin = where(finite(det.wlim),nf,complement=inan,ncomplement=ni)
	
	;create qe array
	qearr = dblarr(det.n,det.n)
	
	;fill in order values
	qearr[ifin] = interpol(qewl[1,*],qewl[0,*],det.wlim[ifin])
	
	;fill in interorder values
	qearr[inan] = avgqe
	
	;store result
	det.qe = qearr

end








