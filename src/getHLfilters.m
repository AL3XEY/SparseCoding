function [HLfilters] = getHLfilters(imgs, nHL, display)
    if nargin<2 || isempty(nHL)
        nHL = 10;
        %nHL = 2; %the number of prototypes to take per image
        %nHL = 100;
        %nHL = 1000;
    end
    if nargin<3 || isempty(display)
        display = false;
    end

    if is_octave
		pkg load image;
    end
    
	HMAXparams = HMAXparameters();

	nimg = size(imgs,2);

    HLfilters = cell(1,HMAXparams.nscales-1);
	for scal=1:HMAXparams.nscales-1
		sz = HMAXparams.filter_sz(scal);
		HLfilters{scal} = zeros(sz, sz, HMAXparams.nth, nHL*nimg);
	end

	gaborFilters = getGaborFilters(HMAXparams, display); % build Gabor filters

	for imgcpt=1:nimg
		img = imgs{imgcpt};
		if size(img,3)==3
			img = double(rgb2gray(img));%/255; % convert it to grayscale
		end
		[dx,dy] = size(img);
        if display
            figure
            imshow(uint8(255*img)) % show original image
        end
        
		S1 = getS1(img, gaborFilters, HMAXparams, display);

		C1 = getC1(S1, dx, dy, HMAXparams, display);
		
        %Take random patches from C1
        rdx = 1+floor(rand(nHL,1).*(dx)); %TODO unique ?
        rdy = 1+floor(rand(nHL,1).*(dy));
        for scal=1:HMAXparams.nscales-1
			sz = HMAXparams.filter_sz(scal);
			[dx2,dy2,~]=size(C1{scal});
            maxX = dx2-sz; %TODO if maxX<0, stop
            maxY = dy2-sz; %TODO if maxY<0, stop
            ratio = dx2/dx;
			for cpt=1:nHL
                x = round(rdx(cpt)*ratio); %TODO calculate it as vector, before this loop ?
                y = round(rdy(cpt)*ratio);
                if x>maxX
                    x = maxX;
                end
                if y>maxY
                    y = maxY;
                end
				HLfilters{scal}(1:sz,1:sz,:,(imgcpt-1)*nHL+cpt) = C1{scal}(x:x+sz-1,y:y+sz-1,:);
			end
		end
	end
end
