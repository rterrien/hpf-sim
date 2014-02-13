;+
; NAME:
;  hpf_set_fluence_detector
;
; PURPOSE:
;
;	Set the fluence for the detector
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;	hpf_set_fluence_detector, det, flux
;
; INPUTS:
;
;	det: The detector structure
;
;	flux: The array of fluences
;
; OUTPUTS:
;	
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-

pro hpf_set_fluence_detector, det, flux

	det.fluence = flux

	
end
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	