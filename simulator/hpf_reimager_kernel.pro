function hpf_reimager_kernel, x_pixel_scale, y_pixel_scale
;;# Huygens PSF
;;# 1.0000 µm at 0.0000 mm.
;;# Data spacing is 1.000 µm.
;;# Data area is 64.000 by 64.000 µm.
;;# Strehl ratio: 0.994
;;# Pupil grid size: 128 by 128
;;# Image grid size: 64 by 64

; I will just hard code this for now, not sure if we need any degrees of freedom

; scales should be pix/um, expect 0.0286 / 4 / 100. for x and 0.0286 / 100. for y

	a = mrdfits('support/huygens_output_reimager_f3.3_f8.fits')
	len_f8 = 64d ;um
	scale = 3.3 / 8d
	len_f3 = len_f8 * scale
	
	lx_pix = round(len_f3 * x_pixel_scale)
	ly_pix = round(len_f3 * y_pixel_scale)
	
	out = congrid(a,lx_pix,ly_pix)
	stop
	return,out
	
end