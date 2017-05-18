function [gaborFilters] = getGaborFilters(display)
	if nargin<1 || isempty(display)
        display = false;
    end
	[sigma lambda gam nth nscales filter_sz pool_sz] = HMAXparameters();
	gaborFilters = cell(1,nscales);
	th = [0:nth-1]*pi/nth; % The orientations
	for scal = 1:nscales
		nxy = filter_sz(scal); % size of the filter
		xx = [-(nxy-1)/2:(nxy-1)/2];
		yy = xx;
		[x,y] = meshgrid(xx,yy);

		[filt] = gabor(x,y,th,scal);
		mn = mean(mean(filt));
		for i=1:nth
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
			for i = [1:nth]
			  imas = filt(:,:,i);
			  nor = max(imas(:));
			  subplot(3,4,i)
			  imshow((imas/nor+1)/2)
            end
		end

		gaborFilters{scal} = filt;
	end
end
