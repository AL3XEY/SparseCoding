function [ S1 ] = getS1( img, HMAXparams, gaborFilters, display )
    if nargin<2 || isempty(HMAXparams)
        HMAXparams = HMAXparameters();
    end
    if nargin<3 || isempty(gaborFilters)
        gaborFilters = getGaborFilters(HMAXparams, false);
    end
    if nargin<4 || isempty(display)
        display = false;
    end

    S1 = cell(1,HMAXparams.nscales);
	for scal = 1:HMAXparams.nscales
        for i = 1:HMAXparams.nth
		  filtr = gaborFilters{scal}(:,:,i);
		  S1{scal}(:,:,i) = abs(filter2(filtr,img)); % filtered images %TODO why abs?
        end

		% display filtered images
        if display
			figure
            for i = 1:HMAXparams.nth
			  imaf = S1{scal}(:,:,i);
			  vis = max(imaf(:));
			  subplot(3,4,i) %TODO only if 12 orientations %TODO in HMAXparams ?
			  imshow(uint8(255*(imaf/vis + 0.3)))
            end
         end
	end
end
