;; set_wlim_detector
;; HPF DETECTOR SIMULATOR
;; This is the code to set the detector's wavelength image
;; 
;; CALLS
;; 
;; 
;; DIRECTLY MODIFIES
;;
;; CALLED BY
;; 
;;
;; NOTES
;; 
;;
;; parameters:
;; det - detector structure
;; wlim - nxn dblarr with the wavelengths for each pixel and NaNs in the interorder regions

pro set_wlim_detector, det, wlim
	
	det.wlim = wlim

end