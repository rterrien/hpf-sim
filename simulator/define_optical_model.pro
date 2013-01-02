;; define_optical_model
;; HPF DETECTOR SIMULATOR
;; This routine looks up the gaps, lambdalow/high, and number of orders for a given
;; optical model version and choice of orders
;; note that this has not been updated for the models from Dec 2012
;; 
;; CALLS
;; 
;; 
;; DIRECTLY MODIFIES

;;
;; CALLED BY
;; hzpfsim_img
;; initialize_detector
;;
;; NOTES
;; 
;;
;; parameters:
;; version - 1 or 2, the first or second optical model versions
;; orders_choice - blue or red, the bluest 22 or the reddest 17 orders
;; out_gap - the array of gaps between the orders in pixels
;; out_lambdalow/high - the lower and upper wavelength limits for each order
;; out_norders - the number of orders for the given model
;; out_ordernums - the order numbers (m)

pro define_optical_model, version, orders_choice=orders_choice, out_gap=out_gap, out_lambdalow=out_lambdalow, out_lambdahigh=out_lambdahigh, out_norders=out_norders, out_ordernums=out_ordernums

	case version of
	1:begin
		out_norders = 17
		out_gap = [0.0, 2.79, 2.67, 2.55, 2.45, 2.35, 2.25, 2.16, 2.08, 2.00, 1.93, 1.86, 1.79, 1.73, 1.67, 1.62, 1.56] * 1000 / 18. ;;pixels, center to center
		out_ordernums = LINDGEN(17)+46
		out_lambdalow = [13173, 12893, 12624, 12366, 12119, 11881, 11653, 11433, 11221, 11017, 10821, 10631, 10448, 10270, 10099, 9934, 9773]/1d4
		out_lambdahigh = [13390, 13105, 12832, 12570, 12319, 12078, 11845, 11622, 11407, 11199, 10999, 10806, 10620, 10440, 10266, 10098, 9935]/1d4
	end
	2:begin
		gap_all=[0.00, 2.78,2.66,2.55,2.44,2.34,2.25,2.16,2.08,2.00,1.93,1.86,1.79,1.73,1.67,1.62,1.56,1.51,1.47,1.42,1.38,1.33,1.30,1.26,1.22,1.19,1.15] * 1d3 / 18d
	lambdalow_all=[13223,12942,12672,12414,12166,11927,11698,11477,11264,11060,10862,10672,10488,10310,10138,9972,9811,9655,9504,9358,9216,9079,8945,8816,8690,8567,8448]/1d4
				lambdahigh_all=[13400,13114,12841,12579,12328,12086,11853,11630,11414,11207,11007,10814,10627,10447,10273,10105,9942,9784,9631,9483,9339,9200,9064,8933,8805,8681,8561]/1d4
		case orders_choice of
		'red':begin
			out_norders = 17
			ll = 0
			uu = 16
		end
		'blue':begin
			out_norders = 22
			ll = n_elements(gap_all)-22
			uu = n_elements(gap_all)-1
		end
		else:stop
		endcase
		out_gap = gap_all[ll:uu]
		out_lambdalow = lambdalow_all[ll:uu]
		out_lambdahigh = lambdahigh_all[ll:uu]
		;make sure the first gap is 0 (the code uses a cumulative sum to figure out the displacement)
		out_gap[0] = 0.
	end
	endcase

end