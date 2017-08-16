function [gaborFilters] = getGaborFilters(HMAXparams, display)
    %get gabor filters as output, given the parameters wanted
    if nargin==0 || isempty(HMAXparams)
        HMAXparams = HMAXparameters();
    end
    if nargin<2 || isempty(display)
        display = false;
    end

	gaborFilters = cell(1,HMAXparams.nscales);
	th = (0:HMAXparams.nth-1)*pi/HMAXparams.nth; % the orientations
	for scal = 1:HMAXparams.nscales
		nxy = HMAXparams.filter_sz(scal); % size of the filter
		xx = (-(nxy-1)/2:(nxy-1)/2);
		yy = xx;
		[x,y] = meshgrid(xx,yy);
		[filt] = gabor(x,y,th,scal,HMAXparams);
        if is_octave
            filt = filt - mean(mean(filt)); % centering, Octave only
        else
            mn = mean(mean(filt));
            for i=1:HMAXparams.nth
                filt(:,:,i) = filt(:,:,i) - mn(i); %centering
            end
        end
		
        % normalization (L2 norm)
        filt = filt ./ sqrt(sum(sum(filt.^2)));

		% display the filters
        if display
			figure
            for i = 1:HMAXparams.nth
			  imas = filt(:,:,i);
			  nor = max(imas(:));
			  subplot(HMAXparams.displayH,HMAXparams.displayW,i)
			  imshow((imas/nor+1)/2)
              title(sprintf('th = %d', i));
            end
            a = axes;
            t1 = title(sprintf('Gabor Filters | scale = %d', scal));
            set(a,'Visible','off');
            set(t1,'Visible','on');
        end
		gaborFilters{scal} = filt;
	end
end
