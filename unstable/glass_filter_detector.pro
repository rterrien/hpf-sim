;; glass_filter_detector
;; HPF DETECTOR SIMULATOR
;; This routine forms the filter array, which will be interpolated by the timestep
;; routine to reduce the number of recorded counts
;; 
;; CALLS
;; 
;; 
;; DIRECTLY MODIFIES
;;
;;
;; CALLED BY
;; expose
;;
;; NOTES
;; 
;;
;; parameters:
;; det - detector structure
;; out - the output wavelength,transmission array

function sellmeier, wl, arr

;;return sellmeier function to find index of refraction with a given
;;set of coefficients

;;arr = [c1,c2,c3,b1,b2,b3]

n =( 1. + arr[3] * wl^2./(wl^2.-arr[0]) + arr[4] * wl^2./(wl^2.-arr[1]) + arr[5] * wl^2./(wl^2.-arr[2]) )^(.5)

return, n

end

pro glass_filter_detector, det


	;prepare the arrays
	p50 = dblarr(2,3001)
	k50 = dblarr(2,3001)
	k12 = dblarr(2,3001)
	p12 = dblarr(2,3001)
	junk = strarr(86)
	
	;read in data
	openr,1,'kzfsn5_50_1nmstep.scan'
	readf,1,junk
	readf,1,k50
	close,1
	
	openr,1,'pk50_50_1nmstep.scan'
	readf,1,junk
	readf,1,p50
	close,1
	
	openr,1,'kzfsn5_12.5_1nmstep.scan'
	readf,1,junk
	readf,1,k12
	close,1
	
	openr,1,'pk50_12.5_1nmstep.scan'
	readf,1,junk
	readf,1,p12
	close,1
	
	;;smooth the scans by 25 pixels since they're quite noisy
	;;also cut them so they only go up to 3000nm (300:*), since they start
	;;from red side
	
	k50r = k50[*,300:*]
	p50r = p50[*,300:*]
	k12r = k12[*,300:*]
	p12r = p12[*,300:*]
	
	k50r[1,*] = smooth(10.^(-1.*k50r[1,*]),100)
	p50r[1,*] = smooth(10.^(-1.*p50r[1,*]),100)
	k12r[1,*] = smooth(10.^(-1.*k12r[1,*]),100)
	p12r[1,*] = smooth(10.^(-1.*p12r[1,*]),100)
	
	n = n_elements(p50r[0,*])
	
	;note that the pk50 50mm piece has a spike around 3microns, do not know if i trust this?
	;do not use unless very careful
	
	;;i know the sellmeier coeffs for kzfsn5 only at this point
	;;(these are for MICRONS!)
	ksell = [.00975488335,.0450495404,67.8786495,1.47727858,.191686941,.897333608]
	
	;pcoeffs = [2.2851698,-9.5978829d-3,9.9950907d-3,1.5661848d-4,-4.3866422d-6,3.1438773d-7]
	;kcoeffs = [2.6699840,-1.3941585d-2,2.2384056d-2,7.4780873d-4,-1.7341165d-5,3.4427318d-6]
	
	;;figure out the transmission for the interfaces
	n1 = 1.d ;;for air
	wlmic = k50r[0,*]/1000.
	n2 = sellmeier(wlmic,ksell)
	
	;;t = transmission coeff (as a function of wavelength)
	;; for one glass-air interface
	t = 4.*n1*n2 / (n1 + n2)^2.
	
	;;we will assume that this transmission coefficient is valid for both the kzfsn5
	;;and the pk50 since I can't find the sellmeier values for pk50
	;;(this is probably ok since it doesn't vary much across our wavelength range)
	
	;;add the interface losses back
	k50o = dblarr(2,n)
	p50o = dblarr(2,n)
	k12o = dblarr(2,n)
	p12o = dblarr(2,n)
	
	p50o[0,*] = p50r[0,*]
	p50o[1,*] = p50r[1,*]/t^2.
	k50o[0,*] = k50r[0,*]
	k50o[1,*] = k50r[1,*]/t^2.
	k12o[0,*] = k12r[0,*]
	k12o[1,*] = k12r[1,*]/t^2.
	p12o[0,*] = p12r[0,*]
	p12o[1,*] = p12r[1,*]/t^2.

	;;i do not trust the values the scanner reports around 1d-10 -
	;;these are small values but get scaled to non-negligible T's when the length
	;;of the glass is shortened
	;;so actually replace these with 0
	
	p50_rem10 = p50o[1,*]
	k50_rem10 = k50o[1,*]
	k12_rem10 = k12o[1,*]
	p12_rem10 = p12o[1,*]
	p50_rem10[where(p50_rem10 lt 2d-10)] = 2d-10
	k50_rem10[where(k50_rem10 lt 2d-10)] = 2d-10
	k12_rem10[where(k12_rem10 lt 2d-10)] = 2d-10
	p12_rem10[where(p12_rem10 lt 2d-10)] = 2d-10
	
	;;derive the absorption coeff per mm
	alphap = -1.*alog(p50_rem10)/50.8
	alphak = -1.*alog(k50_rem10)/50.8
	alphak12 = -1.*alog(k12_rem10)/12.5
	alphap12 = -1.*alog(p12_rem10)/12.5
	
	;;scale the 50mm T's to arbitrary thickness and re-remove transmission loss
	;;to find new T's
	;;use 50mm alphas?
	newp = p50_rem10 * exp(alphap12*(50.8 - det.filter_thick))*t^2.
	newk = k50_rem10 * exp(alphak*(50.8 - det.filter_thick))*t^2.
	newk12 = k12_rem10 * exp(alphak12 * (12.5 - det.filter_thick))*t^2.
	newp12 = p12_rem10 * exp(alphap12 * (12.5 - det.filter_thick))*t^2.
	
	;;if using antireflective coating, remove transmission loss
	if det.filter_ar_coating_flag then begin
		;scale the transmission coeffs to reach a maximum of 1 for each glass
		newp /= t^2.
		newp /= max(newp,/nan)
		newk /= t^2.
		newk /= max(newk,/nan)
		newk12 /= t^2.
		newk12 /= max(newk12,/nan)
		newp12 /= t^2.
		newp12 /= max(newp12,/nan)
	endif

	;;resample to form the new array
	wls = dindgen(1000)/1000d*(4.4 - .8) + .8
	tr = dblarr(1000)
	wave = reform(p50o[0,*])/1d3
	
	;;inside the range where there are scan values use those
	;;outside this range use 0
	inw = where(wls ge min(wave) and wls le max(wave),ninw,complement=outw,ncomplement=noutw)
	if ninw eq 0 or noutw eq 0 then stop
	
	;;choose the right glass
	case det.filter_glass of
		'kzfsn5': stop;trans = newk
		'pk50': stop;trans = newp
		'kzfsn5short': trans = newk12
		'pk50short': trans = newp12
		else: begin
			print,'bad glass name'
			stop
		end
	endcase
	
	;;interpolate where there are scan values
	tr[inw] = interpol(trans,wave,wls[inw])
	;tr[outw] = 0d
	tr[outw] = min(tr[inw])
	
	;;replace the values where I substituted 0 earlier (since I do not trust around 1d-10)
	tr[where(~finite(tr))]=0d
	
	;;project into an nxn array for multiplication by the timestep procedure
	res = dblarr(det.n,det.n)
	res[*] = !values.f_nan
	for i=0, det.n - 1 do begin
		if ~finite(det.wlim[0,i]) then continue
		resi = interpol(trans,wave,det.wlim[*,i])
		res[*,i] = resi
		
	endfor
	det.filter_arr = res
	
	;stop
	
	;;output the array to the detector structure
	det.filter_spec[0,*] = wls
	det.filter_spec[1,*] = tr
	
;;	plotsym,0,/fill
;;	openpps,'glassthickscale_125'
;;	plot,k12r[0,*]/1d3,k12r[1,*],/ylog,xtitle=textoidl('\lambda(\mum)'),ytitle='T',title='KZFSN5 Transmission Scaled from 12.5mm',chars=1.5
;;	oplot,wls,tr,co=fsc_color('red')
;;	legend,['12.5mm','6mm'],colors=[fsc_color('black'),fsc_color('red')],/bottom,/left,chars=1.2,psym=[8,8]
;;	closepps
	
end
	
	
		



