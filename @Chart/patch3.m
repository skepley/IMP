function patch3(obj,s,nticks,timespan,varargin)
	if nargin > 4
		patch_color = varargin{1};
	else
		patch_color = 'g';
	end
	t = linspace(0,1,nticks);
	tt = linspace(timespan(1),timespan(2),nticks);
	S = length(s); T = length(t);
	[Xs,Ys] = obj.mesh_eval(s,0);
	Zs = timespan(1)*ones(size(s));
	[Xs(S+1:S+T),Ys(S+1:S+T)] = obj.mesh_eval(s(end),t);
	Zs(S+1:S+T) = tt;
	[Xs(S+T+1:2*S+T),Ys(S+T+1:2*S+T)] = obj.mesh_eval(fliplr(s),1);
	Zs(S+T+1:2*S+T) = timespan(2)*ones(size(s));
	[Xs(2*S+T+1:2*(S+T)),Ys(2*S+T+1:2*(S+T))] = obj.mesh_eval(s(1),fliplr(t));
	Zs(2*S+T+1:2*(S+T)) = flip(tt);
	gcf;
	patch(Xs,Ys,Zs,patch_color,'EdgeColor',patch_color);
end