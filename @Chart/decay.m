function objDecay = decay(obj)
	% - updated Mar18 2017
    error('This should not be called')
    if ~isprop(obj,'Coord') % all classes should use this going forward
        objDecay = obj.Coefficient.decay;
        return
    end

    % --------------------------- old syntax added so CR4BP will work  ----------------------------

    if isprop(obj,'Coord')
        phaseDecay = obj.Coord.decay; % phase space coordinates
    end

	if isprop(obj,'Variation') % variational equation coordinates
		variationDecay = obj.Variation.decay;
	else
		variationDecay = 0;
	end

	if isprop(obj,'AutoDiff') % automatical differentiation coordinates
		autoDiffDecay = obj.AutoDiff.decay;
	else
		autoDiffDecay = 0;
	end
    objDecay = max([max(phaseDecay(:)),max(variationDecay(:)),max(autoDiffDecay(:))]);
end
