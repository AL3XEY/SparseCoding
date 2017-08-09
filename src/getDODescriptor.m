function [DO] = getDODescriptor(SO, HMAXparams, gaborFilters, display)
    if nargin<2 || isempty(HMAXparams)
        HMAXparams = HMAXparameters();
    end
    if nargin<3 || isempty(gaborFilters)
        gaborFilters = getGaborFilters(HMAXparams, false);
    end
    if nargin<4 || isempty(display)
        display = false;
    end
	%nth=HMAXparams.nth;
	nscales = HMAXparams.nscales;
	[h,w,nth,nchans] = size(SO{1});

	%gabor filters on each channel
	for scal = 1:nscales
	    for th = 1:nth
	        filtr = gaborFilters{scal}(:,:,th);
	        for chan=1:8
	            dochans{scal}(:,:,th,chan) = abs(filter2(filtr,SO{scal}(:,:,th,chan))); % filtered images
			end
		end
	end

	if display
		for scal=1:1%nscales
			for th=1:1%nth
				for chan=1:8
					figure
	                colormap gray
					imagesc(dochans{scal}(:,:,th,chan))
				end
			end
		end
	end

	%half-squaring
	for scal=1:nscales
		idx = find(dochans{scal}<0);
		dochans{scal}(idx) = 0;
		dochans{scal} = dochans{scal}.^2;
	end

	if display
		for scal=1:1%nscales
			for th=1:1%nth
				for chan=1:8
					figure
	                colormap gray
					imagesc(dochans{scal}(:,:,th,chan))
				end
			end
		end
	end

	%normalization
	DO=cell(1,nscales);
	for scal=1:nscales
		DO{scal}=zeros(h,w,nth,8);
	end
	k=1;
	sigma=0.225;
	sigma2=sigma^2;
	sm = sum(dochans{scal},3);
	for scal=1:nscales
		for th=1:nth
			DO{scal}(:,:,th,:) = sqrt((k*dochans{scal}(:,:,th,:))./(sigma2+sm));
		end
	end

	if display
		for scal=1:1%nscales
			for th=1:1%nth
				for chan=1:8
					figure
	                colormap gray
					imagesc(DO{scal}(:,:,th,chan))
				end
			end
		end
	end

end
