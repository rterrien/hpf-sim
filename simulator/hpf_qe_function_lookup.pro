;+
; NAME:
;  hpf_qe_function_lookup
;
; PURPOSE:
;
;  A table to translate QE function names into (wl,trans) arrays
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  hpf_qe_function_lookup
;
; INPUTS:
;
;	det_params: A det_params structure as defined in initialize_mask_params
;
; OUTPUTS:
;	
;	out: A 2-d array with (wavelength[microns],transmission)
;
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-


pro hpf_qe_function_lookup, det_params, out

	if det.qe_loc eq '' then begin
		out = -1
		return
	endif
	case det.qe_loc of
		'17': begin
			qe_wl = [.8,1.,1.23,1.5,2.0,3.5,4.4]
			qe_qe = [.5,.5,.7,.7,0.,0.,0.]
		end
		'25': begin
			qe_wl = [.8,1.,1.23,1.5,2.0,3.5,4.4]
			qe_qe = [.7,.7,.7,.7,.7,0.,0.]
		end
		'25adv': begin
			qe_wl = [.8,1.,1.23,1.5,2.0,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.5,4.4]
			qe_qe = [.7,.7,.7,.7,.7,.7,.7,.7,.53,.13,.02,.004,.006,.006,0.,0.]
			qearr=transpose([[qe_wl],[qe_qe]])
		end
		'17adv': begin
			qe_wl = [.81,1.,1.23,1.5,1.75,1.85,1.95,2.05,2.15,2.25,2.75,3.5,4.4]
			qe_qe = [.5,.5,.7,.7,.53,.13,.02,.004,.006,.006,0.,0.,0.]
			qearr=transpose([[qe_wl],[qe_qe]])
		end
		'h2rg_1.7qe_20micron': begin
			readcol,'support/h2rg_1.7qe_20micron',qe_wl,qe_qe,format='D,D'
			qearr=transpose([[qe_wl],[qe_qe]])
		end
		'h2rg_1.7qe_35micron': begin
			readcol,'support/h2rg_1.7qe_20micron',qe_wl,qe_qe,format='D,D'
			qearr=transpose([[qe_wl],[qe_qe]])
		end
		'h2rg_2.5qe_20micron': begin
			readcol,'support/h2rg_1.7qe_20micron',qe_wl,qe_qe,format='D,D'
			qearr=transpose([[qe_wl],[qe_qe]])
		end
		'h2rg_2.5qe_35micron': begin
			readcol,'support/h2rg_1.7qe_20micron',qe_wl,qe_qe,format='D,D'
			qearr=transpose([[qe_wl],[qe_qe]])
		end
		'h2rg_2.5qe_24micron': begin
			readcol,'support/h2rg_2.5qe_24micron',qe_wl,qe_qe,format='D,D'
			qearr=transpose([[qe_wl],[qe_qe]])
		end
		'h2rg_1.7qe_24micron': begin
			readcol,'support/h2rg_1.7qe_24micron',qe_wl,qe_qe,format='D,D'
			qearr=transpose([[qe_wl],[qe_qe]])
		end
	
		else: begin
			exists = file_test(det.qe_loc)
			if ~exists then begin
				print,'CANNOT FIND QE FILE: ',det.qe_loc
				stop
			endif
			testqe = file_search(det.qe_loc)
			readcol,det.qe_loc,qe_wl,qe_mask,qe_qe,qe_rqe_ir,format='D,D,D,D',comment='#'
			qearr=transpose([[qe_wl],[qe_qe/max(qe_qe,/nan)]])
		end
	endcase
	
	out = qearr

end