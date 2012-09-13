pro script_all_15

resolve_routine,'hzpfsim_img',/compile_full_file

;the top-level script to run an HPF simulation
;generate or set input velocitys shifts
;vels=(randomu(seed,20,/uniform)-0.5) * 60d * 1000d
;vels = [0d,0d]
;save,vels,filename='vels63.sav'

;restore,'vels1.sav'
restore,'vels14.sav'


;num = indgen(n_elements(vels))

;times = [1.,3.,10.,30.,100.,300.,1000.]; * 2.d ;,1500.] ;for hauschildt
;times = [1000.]*6.

;times = double([600,720,840,960,1080,1200])
times = [900d]

;num = n_elements(times)*n_elements(vels)
num = indgen(n_elements(times)*n_elements(vels))

;print, vels

ind = '15'
;datdir = '/Volumes/RAID/HZPF/simulator/test_20120621/'
datdir = '~/research/hzpf_sims/prelim/sbtest/test4/'

fldir = datdir+'fl'+ind+'/'
imdir = datdir+'im'+ind+'/'
spdir = datdir+'sp'+ind+'/'
crdir = datdir+'cr'+ind+'/'
spawn,'mkdir ' + fldir
spawn,'mkdir ' + imdir
spawn,'mkdir ' + spdir
spawn,'mkdir ' + crdir

outs = fldir+'test1fl'+string(num,format='(I03)')+'.fits'

;for i=0, 1 do begin
	
	;construct the simulated fluence arrays
	;vs = vels[i*100 : i*100+99]
	vs = vels
	;for j=0, 99 do hzpfsim_img,vs[j],outs[i*n_elements(vs)+j],specfile='bt_34.fits';,tellfile='tellspec2.fits'	
	;for j=0, 1 do hzpfsim_img,vs[j],outs[j],specfile='bt_34.fits'
	
	;low = i*100
	;hi = i*100+99
	i=0
	index = string(i,format='(I1)')
	;on = '0000000'
	on = string(times,format='(I07)')
	expose,fldir,imdir,exptime=900d,iparam={dark_current:0.05d, qe_flag:1b, qe_loc:'25', persist_flag:0b, ipc_flag:1b, ipc_mean:0.02, ipc_sd:0.002, read_noise:18L, reset_noise:0L, photon_noise_Flag:1b, op_mask_flag:0b, well_depth:100000000, flat_flag:0b, bg_flag:1b, bg_temp:200.d, bg_str:1d},outnum=on
		
	;run the extract script
	extract,fldir,imdir,spdir;,tellfile='tellspec2_detected.fits';,/crimg;,tellfile='tellspec_detected.fits'

	;run the mask ccf velocity extraction
	test_mask_sb2,spdir,'vels'+ind+'.sav','maskrvs_'+ind+'_'+index+'.sav',maskfile='btmask2.sav',weightsfile='btmask_weights2.sav';,tellfile='tellspec2_detected.fits'

	;spawn,'rm -f '+imdir+'*'
	;spawn,'rm -f '+fldir+'*'
	;spawn,'mv '+spdir+'*.fits '+crdir
;endfor

end