function [C2] = HMAXfunction(HLfilters, img, display)
	% Build the Gabor filters used by HMAX first layer (S1). For each of the 8 scales
	% the program can build the nth = 12 (here) filters corresponding to the different orientations.
    
    if nargin<3 || isempty(display)
        display = false;
    end
    
	[sigma lambda gam nth nscales filter_sz pool_sz] = HMAXparameters();
    
    gaborFilters = getGaborFilters(display); % build Gabor filters

	if size(img,3)==3
		img = double(rgb2gray(img))/255; % convert it to grayscale
	end
	[dx,dy] = size(img);
	figure
	imshow(uint8(255*img)) % show original image

	%%%%%%%%%%%%
	%%%  S1  %%%
	%%%%%%%%%%%%
	%% Build the filters for the different orientations %%
    S1 = getS1(img, gaborFilters, display);

	%%%%%%%%%%%%
	%%%  C1  %%%
	%%%%%%%%%%%%
    C1 = getC1(S1, dx, dy, display);
	for scal=1:nscales-1
		%%%%%%%%%%%%
		%%%  S2  %%%
		%%%%%%%%%%%%
		%%%%% S2 layer - compute the response of the HLfilters (prototypes randomly taken from C1 layers of a large dataset) to every C1 patch %%%%%
		%nHL is the number of prototypes
		nHL = size(HLfilters{scal},4);
        sz = filter_sz(scal);
		for cpt=1:nHL
			%calculate the response of each prototype over each patch of the C1 layer (see Mutch & Lowe 2008)
			% R = ||X - P||^2 / 2*sigma^2*alpha
			sigma2 = 1;
			alpha = (sz/4)^2;
            [dx2,dy2,~]=size(C1{scal});
			for x=1:dx2-sz
				for y=1:dy2-sz
					X = HLfilters{scal}(:,:,:,cpt);
					P = C1{scal}(x:x+sz-1,y:y+sz-1,:);
					R = X - P;
					R = R.^2;
					R = sum(sum(sum(R)));
					R = sqrt(R);
					R = R/(2*sigma2*alpha);
					S2{scal}(x,y,cpt) = R;
					%S2{scal}(x,y,cpt) = sqrt(sum(sum(sum((HLfilters{scal}(:,:,:,cpt) - C1{scal}(x:x+sz-1,y:y+sz-1,:)).^2))))/(2*sigma2*alpha);
				end
			end
		end
	end

	%Display the S2 layer
	if display
		for scal=1:nscales-2
			for cpt=1:nHL
				figure
				imgS2 = S2{scal}(:,:,cpt);
				vis = max(max(max(imgS2(:,:,:))));
				imshow(uint8(255*(imgS2/vis + 0.3)))
			end
		end
	end

	%%%%%%%%%%%%
	%%%  C2  %%%
	%%%%%%%%%%%%
	%%%%% C2 layer - max response from the S2 layer %%%%%
	C2 = zeros(nHL, 1);
	for scal=1:nscales-1
		S2b(scal,1:nHL) = max(max(S2{scal}(:,:,1:nHL))); %TODO or min ?
	end
	C2(1:nHL) = max(S2b(:,1:nHL));
	C2 = C2';
end

%close all;clear all;clc;display = false;img = imread('../res/Californie_m.JPG');load('HLfilters_2');tic;C2 = HMAXfunction(HLfilters, img, display);toc;
