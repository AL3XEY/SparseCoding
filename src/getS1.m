function [ S1 ] = getS1( img, gaborFilters, HMAXparams, display )
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
			  subplot(3,4,i)
			  imshow(uint8(255*(imaf/vis + 0.3)))
            end
         end
	end
end
