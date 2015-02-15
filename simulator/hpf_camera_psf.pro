;+
; NAME:
;  hpf_camera_psf
;
; PURPOSE:
;
;  Return a two-dimensional kernel approximation of the HPF camera PSF (a gaussian 
;	with EE80 = 18microns
;
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  out = hpf_camera_psf
;
; INPUTS:
;	
; OUTPUTS:
;	
;	A 2-d double array corresponding to the 
;
; KEYWORD PARAMETERS:
;
;	init_params: A structure containing a subset of the parameters listed below, which will supersede the default parameters in this file
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 04-21-2014
;-


function hpf_camera_psf, npix_x, npix_y, npix_x_1, npix_y_1

	; ee80 to fwhm from experiment
	
	ef = [-0.27600846d,1.3293734]
	
	x_fwhm = poly(npix_x_1,ef)
	y_fwhm = poly(npix_y_1,ef)
	
	if npix_x mod 2 eq 0 then npix_x += 1.
	if npix_y mod 2 eq 0 then npix_y += 1.
	
	psf = psf_gaussian(npixel=[npix_x,npix_y],fwhm=[x_fwhm,y_fwhm],/double,/normalize)
	
	return,psf
	
end