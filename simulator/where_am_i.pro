function where_am_i

	;0 for mlb
	;1 for ryan's laptop
	;2 for hammer
	;-1 for unknown
	
	a = getenv('USER')
	if a eq 'ryan' then return, 1
	
	b = getenv('HOST')
	if b eq 'mlb.astro.psu.edu' then return, 0
	
	c = getenv('HOSTNAME')
	if strmid(c,0,6) eq 'hammer' then return, 2
	if strmid(c,0,4) eq 'lion' then return, 3
	
	return, -1
	
end