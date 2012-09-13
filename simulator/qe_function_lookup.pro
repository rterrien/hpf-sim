;; qe_function_lookup
;; HPF DETECTOR SIMULATOR
;; This is a table to translate qe function names (in det.qe_loc) into actual functions
;; 
;; CALLS
;; 
;; 
;; DIRECTLY MODIFIES
;; 
;;
;; CALLED BY
;; expose.pro
;; thermal_bg_detector.pro
;;
;; NOTES
;; 
;;
;; parameters:
;; det - detector structure
;; out - output qe-array [0,*] = wavelengths, [1,*] = qes


pro qe_function_lookup, det, out

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
			readcol,'h2rg_1.7qe_20micron',qe_wl,qe_qe,format='D,D'
			qearr=transpose([[qe_wl],[qe_qe]])
		end
		'h2rg_1.7qe_35micron': begin
			readcol,'h2rg_1.7qe_20micron',qe_wl,qe_qe,format='D,D'
			qearr=transpose([[qe_wl],[qe_qe]])
		end
		'h2rg_2.5qe_20micron': begin
			readcol,'h2rg_1.7qe_20micron',qe_wl,qe_qe,format='D,D'
			qearr=transpose([[qe_wl],[qe_qe]])
		end
		'h2rg_2.5qe_35micron': begin
			readcol,'h2rg_1.7qe_20micron',qe_wl,qe_qe,format='D,D'
			qearr=transpose([[qe_wl],[qe_qe]])
		end
		'h2rg_2.5qe_24micron': begin
			readcol,'h2rg_2.5qe_24micron',qe_wl,qe_qe,format='D,D'
			qearr=transpose([[qe_wl],[qe_qe]])
		end
		'h2rg_1.7qe_24micron': begin
			readcol,'h2rg_1.7qe_24micron',qe_wl,qe_qe,format='D,D'
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