function hpf_combine_tilt_kernel, kernel_disp, kernel_xdisp, kernel_psf, angle, convol_upfactor = convol_upfactor

	size_disp = size(kernel_disp,/dimen)
	size_xdisp = size(kernel_xdisp,/dimen)
	size_psf = size(kernel_psf,/dimen)
	
;	nx = 2*size_disp[0]
;	ny = 2*size_xdisp[0]
;  nx = 3*size_disp[0]
;  ny = 3*size_xdisp[0]
  nx = size_disp[0] ;* 2.
  ny = size_xdisp[0] ;* 2.

	
	mloc = round(size_xdisp)
	con1 = dblarr(nx,ny)
	;con1[size_disp/2 +1:size_disp/2 + size_disp-1 +1,mloc] = kernel_disp
	;con1[size_disp +1:size_disp + size_disp-1 +1,floor(1.5*mloc)] = kernel_disp
	
	con1[nx/2+1,ny/2] = 1d
	;con1[*,ny/2] = kernel_disp
	;con1[nx/2+1,*] = transpose(kernel_xdisp)
	;con1 = kernel_psf
	
	;con2 = convol(con1,transpose(kernel_xdisp),/normalize,/edge_zero)
	
	;con2 = convol(con1,kernel_disp,/center,/normalize,/edge_zero)
	
	con1b = convol(con1,kernel_disp,/center,/normalize,/edge_zero)
	
	con2 = convol(con1b,transpose(kernel_xdisp),/normalize,/edge_zero,/center)
	
	con3 = convol(con2,kernel_psf,/center,/normalize,/edge_zero)
	
	;con3 = convol(con2,transpose(kernel_xdisp),/normalize,/edge_zero,/center)	
	
	;angle = 30. * !dtor

	x_grid = rebin(dindgen(nx),nx,ny)
	y_grid = rebin(1#dindgen(ny),nx,ny)
	
	dx = y_grid * tan(angle*!dtor)
	
	if keyword_set(convol_upfactor) then dx *= convol_upfactor

	x_grid_warped = x_grid + dx
	y_grid_warped = y_grid


	;con3_warped = tri_surf(con3,x_grid_warped,y_grid_warped,gs=[1,1],/linear)
	
	nnx = n_elements(x_grid_warped)
	nnx2 = nnx * 4.
	ppx = [[reform(x_grid_warped - 0.5,nnx)],[reform(x_grid_warped - 0.5,nnx)],[reform(x_grid_warped + 0.5,nnx)],[reform(x_grid_warped + 0.5,nnx)]]
	ppy = [[reform(y_grid_warped - 0.5,nnx)],[reform(y_grid_warped + 0.5,nnx)],[reform(y_grid_warped + 0.5,nnx)],[reform(y_grid_warped - 0.5,nnx)]]
	ppx = reform(transpose(temporary(ppx)),nnx2) + 0.5
	ppy = reform(transpose(temporary(ppy)),nnx2) + 0.5
	polyinds_in = dindgen(nnx2/4d + 1d) * 4.
	polyclip_x_size = floor(max(x_grid_warped) - min(x_grid_warped)) + 1.
	polyclip_y_size = ny
	inds = polyfillaa(ppx,ppy,polyclip_x_size,polyclip_y_size,areas=areas,poly_indices = polyinds_in)
	
	tmp_warp = dblarr(polyclip_x_size,polyclip_y_size)
	for j=0, nnx-1 do begin
		ii1 = polyinds_in[j]
		ii2 = polyinds_in[j+1] - 1
		if ii1 gt ii2 then continue
		tt2 = double(areas[ii1:ii2]) * con3[j]
		tmp_warp[inds[ii1:ii2]] = temporary(tmp_warp[inds[ii1:ii2]]) + tt2
	endfor
  ;stop
	;return,tmp_warp
	return,con3
	
end
	

	