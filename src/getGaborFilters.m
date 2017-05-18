function [gaborFilters] = getGaborFilters(HMAXparams, display)
    if nargin<1 || isempty(display)
        display = false;
    end
    
	gaborFilters = cell(1,HMAXparams.nscales);
	th = [0:HMAXparams.nth-1]*pi/HMAXparams.nth; % The orientations % TODO use parenthesis instead of brackets (does this work on Octave?)
	for scal = 1:HMAXparams.nscales
		nxy = HMAXparams.filter_sz(scal); % size of the filter
		xx = [-(nxy-1)/2:(nxy-1)/2];% TODO use parenthesis instead of brackets (does this work on Octave?)
		yy = xx;
		[x,y] = meshgrid(xx,yy);

		[filt] = gabor(x,y,th,scal,HMAXparams);
		mn = mean(mean(filt));
		for i=1:HMAXparams.nth
			filt(:,:,i) = filt(:,:,i) - mn(i); %centering
		end
		%filt = filt - mean(mean(filt)); % centering %Octave only
		if(is_octave)
			filt = filt ./ sqrt(sum(sumsq(filt))); % normalization (L2 norm)
		else
			filt = filt ./ sqrt(sum(sumsqr(filt))); % normalization (L2 norm)
		end

		% display the filters
        if display
			figure
            for i = 1:HMAXparams.nth
			  imas = filt(:,:,i);
			  nor = max(imas(:));
			  subplot(3,4,i)
			  imshow((imas/nor+1)/2)
            end
        end
		gaborFilters{scal} = filt;
	end
end
