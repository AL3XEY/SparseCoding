function [HLfilters] = getHLfilters(imgs, nHL, display)
    if nargin<2 || isempty(nHL)
        nHL = 1000; %the number of prototypes to take (total) across all images
    end
    if nargin<3 || isempty(display)
        display = false;
    end

    if is_octave
		pkg load image;
    end

    %Adapt the number of images to the number of HL filters wanted
    %TODO : this was easy to code but not optimal because if nHL >> nimg,
    %we store a big number of images (we could store nimg images and access
    %them when multiple times instead)
    %+ the random selection might have non-unique values (it's a very small
    %risk, but still, we generate new x and y every new image)
    %tl;dr : this solution is fast but costs a lot of memory
    nimg = size(imgs,2);
    images = cell(1,nHL);
    if nHL < nimg
		perms = randperm(nimg,nHL);
		for i=1:size(perms,2)
			images{i}=imgs{perms(i)};
		end
		imgs = images;
    end
    if nHL > nimg
		iterations = round(nHL/nimg);
		additions = mod(nHL,nimg);
        for i=1:nimg*iterations
            images{i} = imgs{mod(i-1,nimg)+1};
        end
        if additions
            perms = randperm(nimg,additions);
            for i=1:additions
                images{nimg*iterations+i}=imgs{perms(i)};
            end
        end
		imgs=images;
    end
    nimg = size(imgs,2);
    clear images
    
    HMAXparams = HMAXparameters();
    HLfilters = cell(1,HMAXparams.nscales-1);
	for scal=1:HMAXparams.nscales-1
		sz = HMAXparams.filter_sz(scal);
		HLfilters{scal} = zeros(sz, sz, HMAXparams.nth, nimg);
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
        rdx = 1+floor(rand().*(dx)); %TODO unicity ?
        rdy = 1+floor(rand().*(dy));
        for scal=1:HMAXparams.nscales-1
			sz = HMAXparams.filter_sz(scal);
			[dx2,dy2,~]=size(C1{scal});
            maxX = dx2-sz; %TODO if maxX<0, stop
            maxY = dy2-sz; %TODO if maxY<0, stop
            ratio = dx2/dx;
            x = 1+round(rdx*ratio);
            y = 1+round(rdy*ratio);
            if x>maxX
                x = maxX;
            end
            if y>maxY
                y = maxY;
            end
            HLfilters{scal}(1:sz,1:sz,:,imgcpt) = C1{scal}(x:x+sz-1,y:y+sz-1,:);
        end
	end
end
