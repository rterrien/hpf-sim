;; timestep_detector
;; HPF DETECTOR SIMULATOR
;; This routine progresses the detector structure through time
;; 
;; CALLS
;; 
;; 
;; DIRECTLY MODIFIES
;; det.signal_cts
;; det.arr
;; det.thermal_cts
;; det.persist
;; det.time
;;
;; CALLED BY
;; script_all_xx
;;
;; NOTES
;; 
;;
;; parameters:
;; det - detector structure
;; dt - timestep in seconds

pro timestep_detector, det, dt
	
	common seeds, seed
	
	

	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;MODIFY FLUENCE WITH QE IF NECESSARY
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	
	if det.qe_flag then $
		newflux = double(det.fluence) * double(det.qe) * double(dt) else $
		newflux = double(det.fluence) * double(dt)
		
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;MODIFY FLUENCE WITH FILTER IF NECESSARY
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	
	if det.filter_flag then $
		newflux = double(det.fluence) * double(det.filter_arr) * double(dt) 

		
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;BOOKKEEPING
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	det.signal_cts += mean(newflux[where(newflux gt 0)])
	dark_counts = double(det.dark_current) * double(dt)
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;UPDATE ACCUMULATOR ARRAY
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	if det.photon_noise_flag then begin
		photon_noise = randomn(seed,det.n,det.n,/double) * double(sqrt(newflux))
		det.arr += (double(newflux) + double(photon_noise)) / double(det.gain) + double(dark_counts)
	endif else begin
		det.arr += double(newflux) / double(det.gain) + double(dark_counts)
	endelse
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;BACKGROUND INCLUSION?
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;Note that randomu cannot draw poisson deviates with different means in an array
	;call- hence the nested for loop
	;the background noise and the above implementation of photon noise are inconsistent
	;is it necessary to account for the poisson nature of the process or can i just
	;use a normal distribution like above? (poisson limits to normal in large n)

	if det.bg_flag then begin
		thermcts = dblarr(det.n,det.n)
		for i=0, det.n-1 do begin
			for j=0, det.n-1 do begin
				b_noise = randomu(seed,1,poisson=det.bgarr[i,j],/double)
				thermcts[i,j] = b_noise
			endfor
		endfor
		det.arr += double(thermcts) * double(dt)
		det.thermal_cts += mean(thermcts) * dt
		
	endif
	
			
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;REPLACE SATURATED PIXELS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;check for pixels that exceed well depth, change these to well depth

	saturated = where(det.arr gt det.well_depth,nsat)
	if nsat gt 0 then det.arr[saturated] = det.well_depth
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;UPDATE PERSISTENCE ARRAY
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	if det.persist_flag then begin
		;update the persisted array from previous step
		;ie decay the previous amount and then add on the new charge
		det.persist = det.persist * exp(-1.d * dt / det.persist_time)
		det.persist += double(det.persist_percent * det.arr)
		
		;check for saturated pixels that exceed well depth, change these to well depth
		saturated = where(det.persist gt det.well_depth,nsat)
		if nsat gt 0 then det.persist[saturated] = det.well_depth
	endif
	
	;account timestep
	det.time += dt
end
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	