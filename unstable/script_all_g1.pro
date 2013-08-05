;; script_all_general
;; HPF DETECTOR SIMULATOR
;; Top-level script to run an HPF simulation
;; 
;; CALLS
;; hzpfsim_img
;; expose
;; extract
;; test_mask_sb2
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
;; 


pro script_all_g1, index, datdir = datdir, outdir = outdir, nvels = nvels

resolve_routine,'hzpfsim_img',/compile_full_file


ind = string(index,format='(I02)')

remove_outputs = 0

exposure_time = 900d



if ~keyword_set(nvels) then nvels = 2
vels=(randomu(seed,nvels,/uniform)-0.5) * 60d * 1000d

;restore,'/gpfs/home/rct151/scratch/simulator/out/maskrvs_66_0.sav'

;vels = rvs
;vels = rvs[0]

;vels = replicate(0d,nvels)
;restore,'result/vels28.sav'
;nvels = n_elements(vels)
;save,vels,filename='vels'+ind+'.sav'

;restore,'result/vels17.sav'

;restore,'vels14.sav'

num = indgen(n_elements(vels))
nnum = n_elements(num)

if ~(keyword_set(datdir) and keyword_set(outdir)) then begin
	computer = where_am_i()
	case computer of
	-1: begin
		print,'need datdir and outdir paths'
		stop
	end
	0: begin
		datdir = '/Volumes/RAID/HZPF/simulator/current_20130101/data/'
		outdir = '/Volumes/RAID/HZPF/simulator/current_20130101/result/'
	end
	1: begin
		datdir = '~/research/hzpf_sims/prelim/sbtest/test5/data/'
		outdir = '~/research/hzpf_sims/prelim/sbtest/test5/out/'
	end
	2: begin
		datdir = '/gpfs/scratch/rct151/simulator/data/'
		outdir = '/gpfs/scratch/rct151/simulator/out/'
	end
	endcase
endif


diag_output = outdir+'diag_output_'+ind+'.txt'

;define detector parameters

init_param={ $
	dark_current:0.00d, $ ;0.05
	qe_flag:0b, $
	qe_loc:'h2rg_2.5qe_24micron', $
	ad_flag:0b, $
	persist_flag:0b, $
	ipc_flag:0b, $
	ipc_mean:0.02, $
	ipc_sd:0.002, $
	read_noise:0L, $ ;18
	reset_noise:0L, $
	photon_noise_Flag:0b, $
	op_mask_flag:0b, $
	well_depth:1d30, $
	flat_flag:0b, $
	bg_flag:0b, $
	bg_temp:170.d, $
	filter_flag:0b, $
	filter_glass:'kzfsn5short',$
	filter_thick:6d,$
	filter_ar_coating_flag:1b,$
	optical_model_version:2,$
	optical_model_rangechoice:'red'}
	
;get the order params for hzpfsim_img
define_optical_model, init_param.optical_model_version, orders_choice=init_param.optical_model_rangechoice, out_gap=out_gap, out_lambdalow=out_lambdalow, out_lambdahigh=out_lambdahigh, out_norders=out_norders, out_ordernums=out_ordernums

;below is for testing
optical_model_file = 'model_ramsey_513.sav'
restore,optical_model_file
out_lambdahigh = model.lambdahigh/1d3
out_lambdalow = model.lambdalow/1d3
out_gap = replicate(0,n_elements(model.lambdalow))
out_norders = n_elements(model.lambdalow)
;above is for testing

;out_gap, out_lambdalow/high, out_norders, out_ordernums

;get the fiber params 
;;fiber_fractions = replicate(1/7d,7)
;;fiber_scale = 0.03d
;;fiber_core_um = 100d
;;fiber_cladding_um = 120d
;;fiber_buffer_um = 125d
;;nfibers = 7
;;fiber_extra_sep_um = 125d

fiber_fractions = [.5d,.5d]
fiber_scale = 0.03d
fiber_core_um = 300d
fiber_cladding_um = 0d
fiber_buffer_um = 500d
nfibers = 2
slitwidth_um = 100d
fiber_extra_sep_um = 200d


;other parameters
projection_type = 'tram1'; ;simple or 7_fibers or tram1
pixel_sampling=3.2
upsample_factor = 12
cal_upsample_factor=.1
velshift_style = 'relativistic' ;relativistic or newtonian

;file locations
specfile='bt_34_extended.fits'
maskfile='btmask2.sav'
weightsfile='btmask_weights2.sav'
contam_tellfile = '';'tellspec2.fits'
correct_tellfile = '';'tellspec2_detected.fits'
mask_tellfile = '';'tellspec2_detected.fits'
;calfile = 'cavity1_measured_watts_per_micron.fits'
;calfile = 'comb_120924_ll.fits'
;calfiles = file_search('/Volumes/RAID/HZPF/simulator/current_20120915/ffp_specs/*.fits',/fully_qualify_path)
;calfile = calfiles[3]
;echellogram_file = 'hpf demag=2-0x f8-5 2012dec15 v10-1-echelleogram.dat'
;echellogram_wl_file = 'hpf demag=2-0x f8-5 2012dec15 v10-1-wavelengths.dat'
optical_model = 'ramsey_513' ;choice of 'barnes_1212' or 'ramsey_513'

straight_orders = 1
upfactor = 2. ; this is for the warping/convolution array




openw,diaglun,diag_output,/get_lun
printf,diaglun,'################################# BEGIN DIAGNOSTIC RECORD ###############################'
printf,diaglun,'TIME RUN: ',systime()
free_lun,diaglun

fldir = datdir+'fl'+ind+'/'
imdir = datdir+'im'+ind+'/'
spdir = datdir+'sp'+ind+'/'
crdir = datdir+'cr'+ind+'/'

spawn,'mkdir ' + fldir
spawn,'mkdir ' + imdir
spawn,'mkdir ' + spdir
spawn,'mkdir ' + crdir

outs = fldir+'test1fl'+string(num,format='(I03)')+'.fits'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   MULTIPLE ROUNDS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if nnum gt 100 then begin
	if nnum mod 100 ne 0 then begin
		print,'use numbers divisible by 100'
		stop
	endif
	if nnum mod 100 gt 9 then stop
	for i=0, nnum/100 - 1 do begin
	
		remove_outputs = 1
	
		vs = vels[i*100 : i*100+99]
		
		for j=0, 99 do hzpfsim_img,vs[j],outs[i*n_elements(vs)+j],specfile=specfile,tellfile=contam_tellfile,diagfile = diag_output, calfile=calfile,projection_type = projection_type, pixel_sampling=pixel_sampling, upsample_factor = upsample_factor, cal_upsample_factor=cal_upsample_factor, velshift_style = velshift_style, orders_lambdahigh = out_lambdahigh, orders_lambdalow = out_lambdalow, orders_gaps = out_gap, fiber_fractions = fiber_fractions, fiber_scale = fiber_scale, fiber_core_um = fiber_core_um, fiber_cladding_um = fiber_cladding_um, fiber_buffer_um = fiber_buffer_um, slitwidth_um = slitwidth_um
		
		;calfile = 'cavity1_measured_watts_per_micron.fits'
			
		ind1 = string(i,format='(I1)')
		
		
		on = string(exposure_time,format='(I07)')
		
		expose,fldir,imdir,exptime=900d,iparam={dark_current:0.05d, qe_flag:1b, qe_loc:'25', persist_flag:0b, ipc_flag:1b, ipc_mean:0.02, ipc_sd:0.002, read_noise:18L, reset_noise:0L, photon_noise_Flag:1b, op_mask_flag:0b, well_depth:100000000, flat_flag:0b, bg_flag:0b, bg_temp:200.d, bg_str:1d},outnum=on,diagfile=diag_output
		
		expose,fldir,imdir,exptime=exposure_time,iparam=init_param,outnum=on,diagfile=diag_output
			
		extract,fldir,imdir,spdir,tellfile=correct_tellfile,diagfile=diag_output,projection_type = projection_type,orders_lambdahigh = out_lambdahigh, orders_lambdalow = out_lambdalow, orders_gaps = out_gap, fiber_scale = fiber_scale, fiber_core_um = fiber_core_um, fiber_cladding_um = fiber_cladding_um, fiber_buffer_um = fiber_buffer_um
		
		test_mask_sb2,spdir,'vels'+ind+'.sav','maskrvs_'+ind+'_'+ind1+'.sav',maskfile=maskfile,weightsfile=weightsfile,tellfile=mask_tellfile,diagfile=diag_output

		if remove_outputs then begin
			spawn,'rm -f '+imdir+'*'
			spawn,'rm -f '+fldir+'*'
			spawn,'mv '+spdir+'*.fits '+crdir
		endif
		
	endfor
	
endif else begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   SINGLE ROUND
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	
	vs = vels
		
	;for j=0, n_elements(vs)-1 do hzpfsim_img,vs[j],outs[j],specfile='bt_34_extended.fits',tellfile=contam_tellfile,diagfile = diag_output,calfile=calfile, projection_type = '7_fibers', model_version = 2, pixel_sampling=3.2, upsample_factor = 12, cal_upsample_factor = .1
	
	;calfiles = file_search('/Volumes/RAID/HZPF/simulator/current_20120915/ffp_specs/*.fits',/fully_qualify_path)
	;vs = replicate(0d,n_elements(calfiles))
	
	;for j=0, n_elements(vs)-1 do hzpfsim_img,vs[j],outs[j],specfile='bt_34_extended.fits',tellfile=contam_tellfile,diagfile = diag_output,calfile=calfiles[3], projection_type = '7_fibers', model_version = 2, pixel_sampling=3.2, upsample_factor = 12, cal_upsample_factor=.1, velshift_style = 'newtonian'
	
	for j=0, n_elements(vs)-1 do hzpfsim_img,vs[j],outs[j],specfile=specfile,tellfile=contam_tellfile,diagfile = diag_output,calfile=calfile, projection_type = projection_type, pixel_sampling=pixel_sampling, upsample_factor = upsample_factor, cal_upsample_factor=cal_upsample_factor, velshift_style = velshift_style, orders_lambdahigh = out_lambdahigh, orders_lambdalow = out_lambdalow, orders_gaps = out_gap, fiber_fractions = fiber_fractions, fiber_scale = fiber_scale, fiber_core_um = fiber_core_um, fiber_cladding_um = fiber_cladding_um, fiber_buffer_um = fiber_buffer_um, nfibers = nfibers, fiber_extra_sep_um = fiber_extra_sep_um, optical_model = optical_model, slitwidth_um = slitwidth_um, straight_orders = straight_orders, upfactor = upfactor

			
	ind1 = string(0,format='(I1)')
		
	on = string(exposure_time,format='(I07)')
	
	
	expose,fldir,imdir,exptime=exposure_time,iparam=init_param,outnum=on,diagfile=diag_output
			
	extract,fldir,imdir,spdir,tellfile=correct_tellfile,diagfile=diag_output,projection_type = projection_type,orders_lambdahigh = out_lambdahigh, orders_lambdalow = out_lambdalow, orders_gaps = out_gap, fiber_scale = fiber_scale, fiber_core_um = fiber_core_um, fiber_cladding_um = fiber_cladding_um, fiber_buffer_um = fiber_buffer_um, nfibers = nfibers, fiber_extra_sep_um = fiber_extra_sep_um, optical_model = optical_model, straight_orders = straight_orders
	
	
	test_mask_sb2,spdir,'vels'+ind+'.sav','maskrvs_'+ind+'_'+ind1+'.sav',maskfile=maskfile,weightsfile=weightsfile,tellfile=mask_tellfile,diagfile=diag_output
	
	if remove_outputs then begin
		spawn,'rm -f '+imdir+'*'
		spawn,'rm -f '+fldir+'*'
		spawn,'mv -f '+spdir+'*.fits '+crdir
	endif
	

endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   CLEANUP
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


spawn,'mv -f dv_maskrvs_*.sav '+outdir
spawn,'mv -f maskrvs_*.sav '+outdir
spawn,'mv -f sns_maskrvs_*.sav '+outdir
spawn,'rm -f *_tmp'
;spawn,'mv -f vels* '+outdir


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   SUMMARIZE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

dvfiles = file_search(outdir+'dv_maskrvs*'+ind+'*.sav',/fully_qualify_path)
rvfiles = file_search(outdir+'maskrvs*'+ind+'*.sav',/fully_qualify_path)

dv_all = []
rv_all = []
if n_elements(dvfiles) ne n_elements(rvfiles) then stop

for i=0, n_elements(dvfiles)-1 do begin
	restore,dvfiles[i]
	restore,rvfiles[i]
	dv_all = [dv_all,dv]
	rv_all = [rv_all,rvs]
endfor

openu,diaglun,diag_output,/get_lun,/append

printf,diaglun,string(13B)+'########### SUMMARY ##############'
printf,diaglun,'DV FILES: '
for i=0, n_elements(dvfiles)-1 do printf,diaglun,dvfiles[i]
printf,diaglun,'RV FILES: '
for i=0, n_elements(rvfiles)-1 do printf,diaglun,rvfiles[i]

printf,diaglun,string(13B) + 'DV MEAN: ', mean(dv_all)
printf,diaglun,'DV STDDEV: ', stddev(dv_all)
printf,diaglun,'RV MEAN: ' ,mean(rv_all)
printf,diaglun,'RV STDDEV: ', stddev(rv_all)
printf,diaglun,'EXP TIME: ',exposure_time

free_lun,diaglun

stop

end
