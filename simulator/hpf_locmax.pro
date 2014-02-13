;+
; NAME:
;  hpf_locmax
;
; PURPOSE:
;
;  Locate the maxima in a 1-d array
;
; CATEGORY:
;
;  HPF simulator
;
; CALLING SEQUENCE:
;
;  out = hpf_locmax(arr,width)
;
; INPUTS:
;
;	arr: The 1-d array to find the maxima in
;
;	width: The width to require a given max to have
;	
; OUTPUTS:
;	
;	Array same size as input, with original values at max and 0's elsewhere
;
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-03-2014
;-


function hpf_locmax, arr, width

;return an array of the same size as the input
;with original values where the mins are and 0's elsewhere

if n_elements(arr) lt width then message, 'arr is too small'
if (width mod 2) eq 0 then message, 'width must be odd'

ic = (width-1)/2
arrc = arr(ic:*)
b=bytarr(n_elements(arr)-width+1) + 1b
for i=1, ic do $
   b = b and (arrc ge arr(ic-i:*)) and (arrc ge arr(ic+i:*))
return, [arr(0:ic-1)*0, b*arrc, arr(0:ic-1)*0]

end 
