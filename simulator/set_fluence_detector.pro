;; set_fluence_detector
;; HPF DETECTOR SIMULATOR
;; This is the code to set the fluence array of the detector
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
;; flux - a nxn array of fluxes in doubles

pro set_fluence_detector, det, flux

	det.fluence = flux

	
end
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	