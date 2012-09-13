function where_am_i

	;0 for mlb
	;1 for ryan's laptop
	;-1 for unknown
	
	a = getenv('USER')
	if a eq 'ryan' then return, 1
	
	b = getenv('HOST')
	if b eq 'mlb.astro.psu.edu' then return, 0
	
	return, -1
	
end