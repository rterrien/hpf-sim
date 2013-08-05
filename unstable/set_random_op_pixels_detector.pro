;; set_random_op_pixels_detector
;; HPF DETECTOR SIMULATOR
;; This is the code to set random bad pixels
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
;; percentage - percentage of bad pixels

pro set_random_op_pixels_detector, det, percentage

	;set a random selection of the pixels in detector det to inoperable
	
	n_total = det.n * det.n
	
	n_inop = round(n_total * percentage)
	
	all_inds = sort(randomu(seed,n_total))
	
	inop_inds = all_inds[0 : n_inop - 1]
	
	det.op_mask[*] = 1b
	det.op_mask[inop_inds] = 0b
	
end