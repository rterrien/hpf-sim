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

contam_tellfile = '';'tellspec2.fits'
correct_tellfile = '';'tellspec2_detected.fits'
mask_tellfile = '';'tellspec2_detected.fits'

ind = string(index,format='(I02)')

remove_outputs = 1

exposure_time = 900d

calfile = 'cavity1_measured_watts_per_micron.fits'


if ~keyword_set(nvels) then nvels = 2
vels=(randomu(seed,nvels,/uniform)-0.5) * 60d * 1000d
save,vels,filename='vels'+ind+'.sav'

restore,'result/vels17.sav'

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
		datdir = '/Volumes/RAID/HZPF/simulator/current_20120908/data/'
		outdir = '/Volumes/RAID/HZPF/simulator/current_20120908/result/'
	end
	1: begin
		datdir = '~/research/hzpf_sims/prelim/sbtest/test5/data/'
		outdir = '~/research/hzpf_sims/prelim/sbtest/test5/out/'
	end
	endcase
endif


diag_output = outdir+'diag_output_'+ind+'.txt'


init_param={ $
	dark_current:0.05d, $ ;0.05
	qe_flag:1b, $
	qe_loc:'h2rg_2.5qe_24micron', $
	ad_flag:1b, $
	persist_flag:0b, $
	ipc_flag:0b, $
	ipc_mean:0.02, $
	ipc_sd:0.002, $
	read_noise:0L, $ ;18
	reset_noise:0L, $
	photon_noise_Flag:1b, $
	op_mask_flag:0b, $
	well_depth:100000000, $
	flat_flag:0b, $
	bg_flag:1b, $
	bg_temp:170.d, $
	filter_flag:1b, $
	filter_glass:'kzfsn5short',$
	filter_thick:6d,$
	filter_ar_coating_flag:1b}



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
		
		for j=0, 99 do hzpfsim_img,vs[j],outs[i*n_elements(vs)+j],specfile='bt_34.fits',tellfile=contam_tellfile,diagfile = diag_output, calfile=calfile
		
		;calfile = 'cavity1_measured_watts_per_micron.fits'
			
		ind1 = string(i,format='(I1)')
		
		
		on = string(exposure_time,format='(I07)')
		
		;expose,fldir,imdir,exptime=900d,iparam={dark_current:0.05d, qe_flag:1b, qe_loc:'25', persist_flag:0b, ipc_flag:1b, ipc_mean:0.02, ipc_sd:0.002, read_noise:18L, reset_noise:0L, photon_noise_Flag:1b, op_mask_flag:0b, well_depth:100000000, flat_flag:0b, bg_flag:0b, bg_temp:200.d, bg_str:1d},outnum=on,diagfile=diag_output
		
		expose,fldir,imdir,exptime=exposure_time,iparam=init_param,outnum=on,diagfile=diag_output
			
		extract,fldir,imdir,spdir,tellfile=correct_tellfile,diagfile=diag_output
		
		test_mask_sb2,spdir,'vels'+ind+'.sav','maskrvs_'+ind+'_'+ind1+'.sav',maskfile='btmask2.sav',weightsfile='btmask_weights2.sav',tellfile=mask_tellfile,diagfile=diag_output

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
		
	for j=0, n_elements(vs)-1 do hzpfsim_img,vs[j],outs[j],specfile='bt_34.fits',tellfile=contam_tellfile,diagfile = diag_output,calfile=calfile
			
	ind1 = string(0,format='(I1)')
		
	on = string(exposure_time,format='(I07)')
		
	expose,fldir,imdir,exptime=exposure_time,iparam=init_param,outnum=on,diagfile=diag_output
			
	extract,fldir,imdir,spdir,tellfile=correct_tellfile,diagfile=diag_output
	
	test_mask_sb2,spdir,'vels'+ind+'.sav','maskrvs_'+ind+'_'+ind1+'.sav',maskfile='btmask2.sav',weightsfile='btmask_weights2.sav',tellfile=mask_tellfile,diagfile=diag_output

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
spawn,'mv -f vels* '+outdir


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
