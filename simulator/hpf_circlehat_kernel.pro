function hpf_circlehat_kernel, npix_x, npix_y, fiber_size_pix_x

	;if npix mod 2 eq 1 then message,'check to make sure xdisp will work'
	;if npix mod 2 eq 1 then offset = 0 else offset = 0.5
	;if npix mod 2 eq 0 then npix = npix + 1
	;offset = 0.5
	npix = npix_x
	xs = abs(dindgen(npix) - npix/2 + 0.5)
	
	fiber = dblarr(npix)
	fiber[where(xs lt fiber_size_pix_x/2.)] = 1.
	
	x = dindgen(npix)-npix/2 + 0.5
	;x = dindgen(npix-1)-npix/2 + 1.
	mu = 0d
	gamma = fiber_size_pix_x / 1d2
	
	napod = 10. * npix_x / npix_y
	;xapod = dindgen(napod)/double(napod-1) * !dpi / 2d
	;yapod = sin(xapod)
	xapod = dindgen(napod)
	ees = reverse(dindgen(napod)/(napod-1d) * 4.)
	yapod = 10d^(-1d * ees)
	
	apodsi = (npix - fiber_size_pix_x)/2
	xx = dindgen(napod)/double(napod-1) * apodsi ; (apodsi-1) + 1
	;xx = [0,xx]
	;yapod = [0,yapod]
	
	
	xxn = dindgen(apodsi)
	
	apod = dblarr(npix)
	apod[*] = 1d
	a1 = interpol(yapod,xx,xxn)
	apod[0:apodsi-1] = a1
	apod[-1*apodsi : *] = reverse(a1)
	
	
	;cauchy dist
	;ca1 = !dpi * gamma * (1 + (x - mu)^2D/gamma^2D)
	;ca = 1d / ca1
	ca_top = gamma / !dpi
	ca_bot = (x - mu)^2d + gamma^2.
	ca = ca_top / ca_bot
	
	conv = convol(fiber,ca,/center,/edge_truncate)
	;cgplot,conv,ps=1,/ylog
	out = conv * apod
	
	xs_2d = abs(dindgen(npix) - npix/2 + 0.5)

	xs2 = rebin(xs,npix_x,npix_x)
	ys2 = rebin(1#xs,npix_x,npix_x)
	
	dd = sqrt(xs2^2. + ys2^2.)
	yy = interpol(out,x,dd)
	
	; Calculate value for each distance
	zz = congrid(yy,npix_x,npix_y,/center,/interp)
	zztot = total(zz)
	zz = zz/zztot
	
	
	
	return,zz
	
end