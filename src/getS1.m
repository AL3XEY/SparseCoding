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
            %title(sprintf('scale = %d', scal));
            
            for i = 1:HMAXparams.nth
			  imaf = S1{scal}(:,:,i);
			  vis = max(imaf(:));
			  subplot(HMAXparams.displayH,HMAXparams.displayW,i) %TODO only if 12 orientations %TODO in HMAXparams ?
			  imshow(uint8(255*(imaf/vis + 0.3)))
              title(sprintf('th = %d', i));
            end
            a = axes;
            t1 = title(sprintf('S1 | scale = %d', scal));
            set(a,'Visible','off');
            set(t1,'Visible','on');
         end
	end
end
