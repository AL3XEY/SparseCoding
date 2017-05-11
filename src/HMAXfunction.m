function [C2] = HMAXfunction(HLfilters, img, nscales=8, norientations=12, display = false)
	%nscales=8;norientations=12;display=false;
	% Build the Gabor filters used by HMAX first layer (S1). For each of the 8 scales
	% the program can build the nth = 12 (here) filters corresponding to the different orientations.
	global filter_sz = [7,11,15,19,23,27,31,35];
	global sigma = [2.8,4.5,6.7,8.2,10.2,12.3,14.6,17.0];
	global lambda = [3.5,5.6,7.9,10.3,12.7,15.5,18.2,21.2];
	pool_sz = [8,10,12,14,16,18,22,24];
	global gam = 0.3;
	global nth = norientations;

	if size(img,3)==3
		img = rot90(rgb2hsv(img)(:,:,3),0); % convert it to grayscale
	end
	[dx,dy] = size(img);
	figure
	imshow(uint8(255*img)) % show original image

	%%%%%%%%%%%%
	%%%  S1  %%%
	%%%%%%%%%%%%
	%% Build the filters for the different orientations %%
	th = [0:nth-1]*pi/nth; % The orientations
	for scal = 1:nscales
		nxy = filter_sz(scal); % size of the filter
		xx = [-(nxy-1)/2:(nxy-1)/2];
		yy = xx;
		[x,y] = meshgrid(xx,yy);

		filt = gabor(x,y,th,scal);
		filt -= mean(mean(filt)); % centering
		filt ./= sqrt(sum(sumsq(filt))); % normalization (L2 norm)

		% display the filters
		if display
			figure
			for i = [1:nth]
			  imas = filt(:,:,i);
			  nor = max(imas(:));
			  subplot(3,4,i)
			  imshow((imas/nor+1)/2)
			endfor
		end

		for i = [1:nth]
		  filtr = filt(:,:,i);
		  S1{scal}(:,:,i) = abs(filter2(filtr,img)); % filtered images
		endfor

		% display filtered images
		if display
			figure
			for i = [1:nth]
			  imaf = S1{scal}(:,:,i);
			  vis = max(imaf(:));
			  subplot(3,4,i)
			  imshow(uint8(255*(imaf/vis + 0.3)))
			endfor
		end
	end

	%%%%%%%%%%%%
	%%%  C1  %%%
	%%%%%%%%%%%%
	for scal=1:nscales-1
		sz = pool_sz(scal);
		pxm = floor(dx/sz); pym = floor(dy/sz);

		for j = [1:nth] % pooling for every orientation
		  for px = [0:pxm-1]
		    for py = [0:pym-1]
		      C1{scal}(px+1,py+1,j) = max(max([S1{scal}(px*sz+1:(px+1)*sz,py*sz+1:(py+1)*sz,j) S1{scal+1}(px*sz+1:(px+1)*sz,py*sz+1:(py+1)*sz,j)]));
		    endfor
		  endfor
		endfor

		% display filtered and pooled images
		if display
			figure
			for i = [1:nth]
			  imaf = C1{scal}(:,:,i);
			  vis = max(imaf(:));
			  subplot(3,4,i)
			  imshow(uint8(255*(imaf/vis+0.3)))
			endfor
		end

		%%%%%%%%%%%%
		%%%  S2  %%%
		%%%%%%%%%%%%
		%%%%% S2 layer - compute the response of the HLfilters (prototypes randomly taken from C1 layers of a large dataset) to every C1 patch %%%%%
		%nHL is the number of prototypes
		nHL = size(HLfilters{scal},4);
		rdx = 1+floor(rand(nHL,1).*(pxm-sz)); %TODO unique ?
		rdy = 1+floor(rand(nHL,1).*(pym-sz)); %FIXME caution! if input image is too short, pym-sz or pxm-sz can be too small
		for cpt=1:nHL
			%calculate the response of each prototype over each patch of the C1 layer (see Mutch & Lowe 2008)
			% R = ||X - P||^2 / 2*sigma^2*alpha
			sigma2 = 1;
			alpha = (sz/4)^2;
			for x=1:pxm-sz
				for y=1:pym-sz
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

	%Display a small part of the S2 layer
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
	%%%  C1  %%%
	%%%%%%%%%%%%
	%%%%% C2 layer - max response from the S2 layer %%%%%
	C2 = zeros(nHL, 1);
	for scal=1:nscales-1
		S2b(scal,1:nHL) = max(max(S2{scal}(:,:,1:nHL))); %TODO or min ?
	end
	C2(1:nHL) = max(S2b(:,1:nHL));
	C2 = C2';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = gabor(x,y,thet,scale) %
	global sigma, global lambda; global gam; global nth;
	lamb = lambda(scale); % we chose the scale and the filter uses the parameters defined above
	sig = sigma(scale);
	[nx,ny] = size(x);
	x0 = reshape(permute(x,[3,1,2])' * cos(thet) + permute(y,[3,1,2])' * sin(thet),nx,ny,nth);
	y0 = reshape(permute(y,[3,1,2])' * cos(thet) - permute(x,[3,1,2])' * sin(thet),nx,ny,nth);
	out = exp( -0.5 * (x0.^2 + gam * y0.^2) / sig^2) .* cos(2 * pi * x0 / lamb);
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
nscales = 8;
norientations = 12;
display = false;
img = imread('../res/Californie_m.JPG');
load('HLfilters_2');
tic
C2 = HMAXfunction(HLfilters, img, nscales, norientations, display)
toc % print execution time

%close all;clear all;clc;nscales = 8;norientations = 12;display = false;img = imread('../res/Californie_m.JPG');load('HLfilters_2');tic;C2 = HMAXfunction(HLfilters, img, nscales, norientations, display);toc;
