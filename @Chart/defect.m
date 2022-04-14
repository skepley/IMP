%% validation methods
function fout = defect(obj)
	XX = mtimes(obj.var(1),obj.var(1),Inf); %full X^2 including spillover terms
	XXY = mtimes(XX,obj.var(2),Inf); % full XXY including spillover terms
	XXY.coef(1:obj.modes(1),1:obj.modes(2)) = zeros(obj.modes);
	fout = XXY.norm; 
end
