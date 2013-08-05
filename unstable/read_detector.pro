;; read_detector
;; HPF DETECTOR SIMULATOR
;; This is the code to readout the detector array
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
;; uses a double array output to enable use of NaNs
;; A/D conversion is a simple rounding at this point
;;
;; parameters:
;; det - detector structure
;; output - npix x npix array output in doubles


pro read_detector, det, output
	
	;common block for random draws
	common seeds, seed
	
	frame = det.arr
	
	;pull out array and persisted charge
	if det.persist_flag then begin
		frame += det.persist
	endif
	
	;form the read noise array
	noise = randomu(seed,det.n,det.n,/normal) * det.read_noise
	
	;apply read noise
	frame += noise
	
	if det.ipc_flag then begin
		;add also the "persisted" charge
		
		;add the ipc from top, right, bottom, left
		;origin = bottom left
		;column, row
		ipc = dblarr(det.n,det.n)
		;top
		ipc[*,0:-2] += double(det.ipc_coeff[*,0:-2,0] * frame[*,1:*])
		;right
		ipc[0:-2,*] += double(det.ipc_coeff[0:-2,*,1] * frame[1:*,*])
		;bottom
		ipc[*,1:*] += double(det.ipc_coeff[*,1:*,2] * frame[*,0:-2])
		;left
		ipc[1:*,*] += double(det.ipc_coeff[1:*,*,3] * frame[0:-2,*])
		
		;apply IPC
		frame += ipc
	endif
	
	;Done accumulating, now we can AD convert and go to double
	;skip A/D conversion if requested
	if det.ad_flag then frame = double(round(frame)) else frame = double(frame)
	
	if det.op_mask_flag then begin
		bad = where(det.op_mask eq 0b,nb)
		if nb gt 0 then frame[bad] = !values.f_nan
	endif
	
	output = frame
	
	
end




















