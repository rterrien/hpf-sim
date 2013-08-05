;; reset_detector
;; HPF DETECTOR SIMULATOR
;; This code resets the detector array
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


pro reset_detector, det
	
	;Reset to value defined by Normal distribution with mean=0 and sigma = reset_noise
	common seeds, seed
	
	det.arr = double(randomu(seed, det.n, det.n, /normal) * det.reset_noise)
	det.signal_cts = 0.d
	
end
	







