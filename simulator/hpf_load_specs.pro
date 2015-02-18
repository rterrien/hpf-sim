;+
; NAME:
;  hpf_load_specs
;
; PURPOSE:
;
;  Load and process the spectra into the spec_params structure
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  hpf_load_specs, spec_params, det_params
;
; INPUTS:
;
;	spec_params: A spec_params structure as defined in initialize_spec_params 	
;
;	det_params: A det_params structure as defined in initialize_det_params
;	
; OUTPUTS:
;	
;	The indicated spectra are loaded and processed (filtered, etc) into the spec_params structure
;
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-


pro hpf_load_specs, spec_params, det_params
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;DEFINE PHYSICAL AND DETECTOR CONSTANTS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	c = 2.99792458d8	; m/s
	h = 6.626d-34 ;J s
;	minl = .7d
;	maxl = 1.4d
	minl = hpf_spec_params.minl
	maxl = hpf_spec_params.maxl
	
	for i=0, n_elements(spec_params.spec_file)-1 do begin
		
		
		case strupcase(strcompress(spec_params.type[i],/remove_all)) of
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;; STAR
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		'STAR': begin
			extension = strupcase(strmid(spec_params.spec_file[i],3,/reverse_offset))
			if extension ne 'FITS' then message,'NOT FITS FILE: '+spec_params.spec_file[i]
			in = mrdfits(spec_params.spec_file[i])
			if (size(in))[2] gt 3 then in = transpose(in)
			if (size(in))[2] gt 3 then message,'CHECK ARRAY SIZE'
			w = reform(in[*,0])	; assume wavelength is in microns
			f = reform(in[*,1])	; assume flux is in Watts/m2/micron
			sw = sort(w)
			w = w[sw]
			f = f[sw]

			;if spec_params.flat[i] then f[*] = mean(f,/nan)

			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;;CUT DOWN AND UPSAMPLE
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

			g = where(w ge minl - .05d and w le maxl + .05d,ng)
			w = double(w[g])
			f = double(f[g])
	
			;upsample (ensures even/sufficient pixel sampling)
			upsample_factor = spec_params.upsample_factor[i]
			if ~finite(upsample_factor) then upsample_factor = 12.
			if upsample_factor ne -1 then begin
				nel1 = n_elements(w)
				nel2 = upsample_factor*nel1 ;multiply this for upsampling
				neww = dindgen(nel2)/double(nel2) * (maxl - minl) + minl
				newf = interpol(f,w,neww)
				;newf = hermite(w,f,neww)
				w = neww
				f = newf
			endif

			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;;SCALE FLUX
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
			jzp     = 3.31d-9	; Vega flux at J  (W/m2/micron)
			jvega   = 0.025d	; Vega mag at J
			deltaj  = 0.26d		; FWHM of J in microns
			isowj   = 1.215d	; Isophotal wavelength of J in microns
			rtel    = 10./2.0d	; radius of telescope in m
			telarea = double(!pi*rtel^2)    ; telescope area in m^2
			tput    = 0.025d		; total instrument throughput 
	
			jmag = spec_params.jmag[i]
	
			jmin   = isowj - deltaj/2.0d
			jmax   = isowj + deltaj/2.0d
			jndx   = where((w ge jmin) and (w le jmax))
			fmean  = int_tabulated(w[jndx],f[jndx])/deltaj
			jscale = (jzp/fmean)*10.0d^(-0.4d*(jmag - jvega))
			f      = temporary(f)*jscale		; scale spectrum to input J mag

			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;;TELLURIC CONTAMINATE IF NECESSARY
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;do this while spectrum is still in energy
			if spec_params.tellcontam_flag[i] then begin
				tellspec = mrdfits(spec_params.tellcontam_file[i]) ; expect wavelength in microns, y 0-1
				tellw = reform(tellspec[0,*])
				tellf = reform(tellspec[1,*])
				cl = where((w ge min(tellw)) and (w le max(tellw)),ncl)
				if ncl eq 0 then stop 
				contam = interpol(tellf,tellw,w[cl])
				f[cl] = f[cl] * contam
			endif

			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;; APPLY QE IF NECESSARY
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			if det_params.qe_flag then begin
				qe_wl = reform(*(det_params.qewl[0,*]))
				qe_tr = reform(*(det_params.qewl[1,*]))
				mmqewl = minmax(qe_wl)
				olow = where(w lt mmqewl[0],nlow)
				ohi = where(w gt mmqewl[1],nhi)
				qes = interpol(qe_tr,qe_wl,w)
				if nlow ne 0 then qes[olow] = qes[olow[-1]+1]
				if nhi ne 0 then qes[ohi] = qes[ohi[0]-1]
				f = f * qes
			endif
	
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;;CONVERT TO PHOTONS
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;photons/s/m2/micron

			p = w * f / (h*c) * 1d-6
			pout = temporary(p)*telarea*tput ;photons/s/micron
			fold = f
			f = pout
		end
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;; CAL
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		'CAL': begin
			extension = strupcase(strmid(spec_params.spec_file[i],3,/reverse_offset))
			if extension ne 'FITS' then message,'NOT FITS FILE: '+spec_params.spec_file[i]
			in = mrdfits(spec_params.spec_file[i])
			if (size(in))[2] gt 3 then in = transpose(in)
			if (size(in))[2] gt 3 then message,'CHECK ARRAY SIZE'

			w = reform(in[*,0])	; assume wavelength is in microns
			f = reform(in[*,1])	; assume flux is in Watts/m2/micron
			sw = sort(w)
			w = w[sw]
			f = f[sw]
			if spec_params.flat[i] then begin
				f[*] = mean(*spec_params.fl[0],/nan)
				if f[0] eq 0d then f[*] = 1d
			endif
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;;CUT DOWN AND UPSAMPLE
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

			g = where(w ge minl - .05d and w le maxl + .05d,ng)
			w = double(w[g])
			f = double(f[g])
	
			;upsample (ensures even/sufficient pixel sampling)
			upsample_factor = spec_params.upsample_factor[i]
			if ~finite(upsample_factor) then upsample_factor = 12.
			if upsample_factor ne -1 then begin
				nel1 = n_elements(w)
				nel2 = upsample_factor*nel1 ;multiply this for upsampling
				neww = dindgen(nel2)/double(nel2) * (maxl - minl) + minl
				newf = interpol(f,w,neww)
				w = neww
				f = newf
			endif
			

			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;; FILTER
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			
			if spec_params.filter[i] ne 1d then begin
				f = f * double(spec_params.filter[i])
			endif
			
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;; APPLY QE IF NECESSARY
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			if det_params.qe_flag then begin
				qe_wl = reform(*(det_params.qewl[0,*]))
				qe_tr = reform(*(det_params.qewl[1,*]))
				mmqewl = minmax(qe_wl)
				olow = where(w lt mmqewl[0],nlow)
				ohi = where(w gt mmqewl[1],nhi)
				qes = interpol(qe_tr,qe_wl,w)
				if nlow ne 0 then qes[olow] = qes[olow[-1]+1]
				if nhi ne 0 then qes[ohi] = qes[ohi[0]-1]
				f = f * qes
			endif
			
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;; Cut out dead regions of LFC
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			if spec_params.lfc_lims_flag eq 1 then begin
				gg = dblarr(n_elements(w))
				for j=0,(size(spec_params.lfc_lims,/dimen))[1]-1 do begin
					ggi = where( (w ge spec_params.lfc_lims[0,j]) and (w le spec_params.lfc_lims[1,j]),nggi)
					if nggi ne 0 then begin
						gg[ggi] = 1.
					endif
				endfor
				bbi = where(gg eq 0, nbbi)
				if nbbi ne 0 then f[bbi] = 1d-3 * median(f)
			endif
			
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;;CONVERT TO PHOTONS
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;photons/s/m2/micron

			p = w * f / (h*c) * 1d-6
			;pout = temporary(p)*telarea*tput ;photons/s/micron
			;fold = f
			f = p
		end
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;; SKY
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		'SKY': begin
			extension = strupcase(strmid(spec_params.spec_file[i],3,/reverse_offset))
			if extension ne 'FITS' then message,'NOT FITS FILE: '+spec_params.spec_file[i]
			in = mrdfits(spec_params.spec_file[i])
			if (size(in))[2] gt 3 then in = transpose(in)
			if (size(in))[2] gt 3 then message,'CHECK ARRAY SIZE'

			w = reform(in[*,0])	; assume wavelength is in microns
			f = reform(in[*,1])	; assume flux is in Watts/m2/micron
			sw = sort(w)
			w = w[sw]
			f = f[sw]
			if spec_params.flat[i] then begin
				f[*] = mean(*spec_params.fl[0],/nan)
				if f[0] eq 0d then f[*] = 1d
			endif
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;;CUT DOWN AND UPSAMPLE
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

			g = where(w ge minl - .05d and w le maxl + .05d,ng)
			w = double(w[g])
			f = double(f[g])
	
			;upsample (ensures even/sufficient pixel sampling)
			upsample_factor = spec_params.upsample_factor[i]
			if ~finite(upsample_factor) then upsample_factor = 12.
			if upsample_factor ne -1 then begin
				nel1 = n_elements(w)
				nel2 = upsample_factor*nel1 ;multiply this for upsampling
				neww = dindgen(nel2)/double(nel2) * (maxl - minl) + minl
				newf = interpol(f,w,neww)
				w = neww
				f = newf
			endif

			
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;; APPLY QE IF NECESSARY
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			if det_params.qe_flag then begin
				qe_wl = reform(*(det_params.qewl[0,*]))
				qe_tr = reform(*(det_params.qewl[1,*]))
				mmqewl = minmax(qe_wl)
				olow = where(w lt mmqewl[0],nlow)
				ohi = where(w gt mmqewl[1],nhi)
				qes = interpol(qe_tr,qe_wl,w)
				if nlow ne 0 then qes[olow] = qes[olow[-1]+1]
				if nhi ne 0 then qes[ohi] = qes[ohi[0]-1]
				f = f * qes
			endif


			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;;CONVERT TO PHOTONS
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;photons/s/m2/micron

			p = w * f / (h*c) * 1d-6
			pout = temporary(p)*telarea*tput ;photons/s/micron
			fold = f
			f = pout
		end
		'FLAT': begin
			w = *spec_params.wl[0]
			f = *spec_params.fl[0]
			f[*] = mean(*spec_params.fl[0],/nan)
			if f[0] eq 0d then f[*] = 1d
		
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;;TELLURIC CONTAMINATE IF NECESSARY
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;do this while spectrum is still in energy
			if spec_params.tellcontam_flag[i] then begin
			tellspec = mrdfits(spec_params.tellcontam_file[i]) ; expect wavelength in microns, y 0-1
			tellw = reform(tellspec[0,*])
			tellf = reform(tellspec[1,*])
			contam = interpol(tellf,tellw,w)
			f = f * contam
			endif

			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;; FILTER
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

			if spec_params.filter[i] ne 1d then begin
			f = f * double(spec_params.filter[i])
			endif


				;temporary make 1061 - 1065 illuminated
				;every .5nm alternate
			;;			mf = mean(*spec_params.fl[0],/nan)
			;;			test = (w*1d4) mod 10.
			;;			on = where(test gt 5.,complement = off)
			;;			f[on] = mf
			;;			f[off] = mf/2d
			;;			f[where(w lt 1.061)] = 0d
			;;			f[where(w gt 1.065)] = 0d
		end
		else: message,'unknown spectrum type'
		endcase
		
		*(spec_params.wl[i]) = w
		*(spec_params.fl[i]) = f
		*(spec_params.shift_wl[i]) = w
		
	endfor
	
end
	