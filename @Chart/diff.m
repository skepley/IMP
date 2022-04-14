function darc = diff(obj)
	darc_dt = [obj.var(1).diff;obj.var(2).diff];
end

function darc_ds = ds(obj)
	darc_ds = [obj.var(1).ds;obj.var(2).ds];
end