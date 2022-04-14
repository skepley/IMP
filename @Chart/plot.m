
        
        %% PLOTTING METHODS
        function plot(obj,s,t,varargin)
            [x,y] = obj.mesh_eval(s,t);
            x = x'; y = y';
            if nargin > 3
                plot_color = varargin{1};
                if isa(plot_color,'char')
                    plot(x,y,plot_color,'LineWidth',1.2);
                else
                    plot(x,y,'LineWidth',1.2,'Color',plot_color);
                end
            else
                plot(x,y,'LineWidth',1.2);
            end
        end