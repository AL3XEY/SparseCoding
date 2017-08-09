function [SO, ochans] = getSODescriptor(img, HMAXparams, gaborFilters, display) %img must be normalized between 0 and 255
    if nargin<2 || isempty(HMAXparams)
        HMAXparams = HMAXparameters();
    end
    if nargin<3 || isempty(gaborFilters)
        gaborFilters = getGaborFilters(HMAXparams, false);
    end
    if nargin<4 || isempty(display)
        display = false;
    end

    W = HMAXparams.W;
    Wexcit = (W>0).*W;
    Winhib = (W<0).*W;

	nchans = 2*HMAXparams.outChansHalf;
	img = double(img)./255;
	%imglms = rgb2lms(img);
	%img = imglms;
	[h,w,c] = size(img);

	nth=HMAXparams.nth;
	nscales = HMAXparams.nscales;
	ochans=cell(1,nscales);
	dochans=cell(1,nscales);
	excit=cell(1,nscales);
	inhib=cell(1,nscales);
	for scal=1:nscales
		ochans{scal} = zeros(h,w,nth,nchans);
	    dochans{scal} = zeros(h,w,nth,nchans);
		excit{scal} = zeros(h,w,nth,HMAXparams.inChans);
		inhib{scal} = zeros(h,w,nth,HMAXparams.inChans);
	end

	for scal = 1:nscales
	    for th = 1:nth
	        filterexcit = (gaborFilters{scal}(:,:,th)>0).*gaborFilters{scal}(:,:,th);
	        filterinhib = -(gaborFilters{scal}(:,:,th)<0).*gaborFilters{scal}(:,:,th);
	        for color=1:3
	            excit{scal}(:,:,th,color) = abs(filter2(filterexcit,img(:,:,color)));
	            inhib{scal}(:,:,th,color) = abs(filter2(filterinhib,img(:,:,color)));
	        end
	    end
	end

	for scal=1:nscales
		for th=1:nth
            for i=1:HMAXparams.outChansHalf
                for j=1:HMAXparams.inChans
                    ochans{scal}(:,:,th,i) = ochans{scal}(:,:,th,i) + Wexcit(j,i).*excit{scal}(:,:,th,j) + Winhib(j,i).*inhib{scal}(:,:,th,j);
                    ochans{scal}(:,:,th,i+HMAXparams.outChansHalf) = -(ochans{scal}(:,:,th,i) + Wexcit(j,i).*excit{scal}(:,:,th,j) + Winhib(j,i).*inhib{scal}(:,:,th,j));%-ochans{scal}(:,:,th,i);
                end
            end
            
            if display %TODO
				figure
				colormap gray
				imagesc(ochans{scal}(:,:,th,2))

				figure
				colormap gray
				imagesc(ochans{scal}(:,:,th,6))
            end
		end
	end

	%TODO display excit, inhib and then ochans/opponnencychannels

	%half-squaring
	for scal=1:nscales
		idx = find(ochans{scal}<0);
		ochans{scal}(idx) = 0;
		ochans{scal} = ochans{scal}.^2;
	end

	if display
		for scal=1:1%nscales
			for th=1:1%nth
	            for chan=1:8
	                figure
	                colormap gray
	                imagesc(ochans{scal}(:,:,th,chan))
	            end
			end
		end
	end

	%divisive normalization
	SO=cell(1,nscales);
	for scal=1:nscales
		SO{scal}=zeros(h,w,nth,8);
	end
	k=1;
	sigma=0.225;
	sigma2=sigma^2;
	sm = sum(ochans{scal},4);
    for scal=1:nscales
		for chan=1:8
			SO{scal}(:,:,:,chan) = sqrt((k*ochans{scal}(:,:,:,chan))./(sigma2+sm));
		end
    end

	if display
        %TODO show every orientation on the same figure and label it
		for scal=1:1%nscales
			for th=1:1%nth
				for chan=1:8
					figure
	                colormap gray
					imagesc(SO{scal}(:,:,th,chan))
				end
			end
		end
	end

end
