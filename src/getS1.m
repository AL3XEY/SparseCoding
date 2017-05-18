function [ S1 ] = getS1( img, gaborFilters, display )
    if nargin<3 || isempty(display)
        display = false;
    end

    [sigma,lambda,gam,nth,nscales,filter_sz,pool_sz] = HMAXparameters();

    S1 = cell(1,nscales);
	for scal = 1:nscales
        for i = 1:nth
		  filtr = gaborFilters{scal}(:,:,i);
		  S1{scal}(:,:,i) = abs(filter2(filtr,img)); % filtered images
        end

		% display filtered images
        if display
			figure
            for i = 1:nth
			  imaf = S1{scal}(:,:,i);
			  vis = max(imaf(:));
			  subplot(3,4,i)
			  imshow(uint8(255*(imaf/vis + 0.3)))
            end
         end
	end
end
