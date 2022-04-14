        
        function patch(obj,s,t,varargin)
            if nargin > 3
                patch_color = varargin{1};
            else
                patch_color = 'g';
            end
            S = length(s); T = length(t);
            [Xs,Ys] = obj.eval(s,t(1));
            [Xs(S+1:S+T),Ys(S+1:S+T)] = obj.mesh_eval(s(end),t);
            [Xs(S+T+1:2*S+T),Ys(S+T+1:2*S+T)] = obj.mesh_eval(fliplr(s),t(end));
            [Xs(2*S+T+1:2*(S+T)),Ys(2*S+T+1:2*(S+T))] = obj.mesh_eval(s(1),fliplr(t));
            gcf;
            patch(Xs,Ys,patch_color,'EdgeColor',patch_color);
        end