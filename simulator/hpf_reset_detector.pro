;+
; NAME:
;  hpf_reset_detector
;
; PURPOSE:
;
;	Reset the detector array, with appropriate reset noise
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;	hpf_reset_detector, det
;
; INPUTS:
;
;	det: The detector structure
;
; OUTPUTS:
;	
;
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-

pro hpf_reset_detector, det
	
	;Reset to value defined by Normal distribution with mean=0 and sigma = reset_noise
	common seeds, seed
	
	det.arr = double(randomu(seed, det.n, det.n, /normal) * det.reset_noise)
	det.signal_cts = 0.d
	
end
	







