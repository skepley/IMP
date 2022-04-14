function arclen = arc_length(obj)
	darc_ds = ds(obj);
	rad = mtimes(darc_ds(1),darc_ds(1)) + mtimes(darc_ds(2),darc_ds(2)); 
	integrand = sqrt(rad);
	arclen = intds(integrand,[-1,1]);
end