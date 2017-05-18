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
			img = double(rgb2gray(img))/255; % convert it to grayscale
		end
		[dx,dy] = size(img);
		figure
		imshow(uint8(255*img)) % show original image

		S1 = getS1(img, gaborFilters, HMAXparams, display);

		C1 = getC1(S1, dx, dy, HMAXparams, display);

		for scal=1:HMAXparams.nscales-1
			%Take random patches from C1
			sz = HMAXparams.filter_sz(scal);
			%pxm = floor(dx/sz); pym = floor(dy/sz);
			[dx2,dy2,~]=size(C1{scal});
			rdx = 1+floor(rand(nHL,1).*(dx2-sz)); %TODO unique ?
			rdy = 1+floor(rand(nHL,1).*(dy2-sz)); %FIXME caution! if input image is too short, pym-sz or pxm-sz can be too small
			for cpt=1:nHL
				HLfilters{scal}(1:sz,1:sz,:,(imgcpt-1)*nHL+cpt) = C1{scal}(rdx(cpt):rdx(cpt)+sz-1,rdy(cpt):rdy(cpt)+sz-1,:);
			end
		end
		%clear S1,C1; %TODO necessary?
	end
end
