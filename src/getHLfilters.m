function [HLfilters] = getHLfilters(imgs, nscales=8, norientations=12)
	%nscales=8; norientations=12;
	if is_octave
		pkg load image;
	end
	global filter_sz = [7,11,15,19,23,27,31,35]; %FIXME remove Matlab errors
	global sigma = [2.8,4.5,6.7,8.2,10.2,12.3,14.6,17.0];
	global lambda = [3.5,5.6,7.9,10.3,12.7,15.5,18.2,21.2];
	global gam = 0.3;
	global nth = norientations;

	display = false; % do we display everything or not

	pool_sz = [8,10,12,14,16,18,22,24];

	%nHL = 2; %the number of prototypes to take per image
	nHL = 100;
	%nHL = 1000;

	nimg = size(imgs,2);

	for scal=1:nscales-1
		sz = pool_sz(scal);
		HLfilters{scal} = zeros(sz, sz, nth, nHL*nimg);
	end

	for scal = 1:nscales %% Build the filters for the different orientations %%
		th = [0:nth-1]*pi/nth; % The orientations
		nxy = filter_sz(scal); % size of the filter
		xx = [-(nxy-1)/2:(nxy-1)/2];
		yy = xx;
		[x,y] = meshgrid(xx,yy);

		filtre = gabor(x,y,th,scal);
		filtre -= mean(mean(filtre)); % centering
		filtre ./= sqrt(sum(sumsq(filtre))); % normalization (L2 norm)
		filt{scal} = filtre;

		% display the filters
		if display
			figure
			for i = [1:nth]
			  imas = filt{scal}(:,:,i);
			  nor = max(imas(:));
			  subplot(3,4,i)
			  imshow((imas/nor+1)/2)
			endfor
		end
	end

	%imgs{1} = imread('../res/Californie_m.JPG');
	%imgs{2} = imread('../res/lena.ppm');
	%imgs{3} = imresize(imread('../res/ZDB/DSCN1572.JPG'),0.5);
	for imgcpt=1:nimg
		img = imgs{imgcpt};
		if size(img,3)==3
			img = rot90(rgb2hsv(img)(:,:,3),0); % convert it to grayscale
		end
		[dx,dy] = size(img);
		figure
		imshow(uint8(255*img)) % show original image

		for scal = 1:nscales
			for i = [1:nth]
			  filtr = filt{scal}(:,:,i);
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

			% taille des pools
		for scal=1:nscales-1
			sz = pool_sz(scal);
			pxm = floor(dx/sz); pym = floor(dy/sz);

			for j = [1:nth] % pooling for every orientation
			  for px = [0:pxm-1]
			    for py = [0:pym-1]
			      C1{scal}(px+1,py+1,j) = max(max([S1{scal}(px*sz+1:(px+1)*sz,py*sz+1:(py+1)*sz,j); S1{scal+1}(px*sz+1:(px+1)*sz,py*sz+1:(py+1)*sz,j)]));
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

			%Take random patches from C1
			rdx = 1+floor(rand(nHL,1).*(pxm-sz)); %TODO unique ?
			rdy = 1+floor(rand(nHL,1).*(pym-sz)); %FIXME caution! if input image is too short, pym-sz or pxm-sz can be too small
			%pxm
			%pym
			%sz
			%rdx
			%rdy
			for cpt=1:nHL
				HLfilters{scal}(1:sz,1:sz,:,(imgcpt-1)*nHL+cpt) = C1{scal}(rdx(cpt):rdx(cpt)+sz-1,rdy(cpt):rdy(cpt)+sz-1,:);
			end
		end
		clear S1,C1;
		if display
			clear imaf,vis;
		end
	end
end

function out = gabor(x,y,thet,scale) %
	global sigma, global lambda; global gam; global nth;
	lamb = lambda(scale); % we chose the scale and the filter uses the parameters defined above
	sig = sigma(scale);
	[nx,ny] = size(x);
	x0 = reshape(permute(x,[3,1,2])' * cos(thet) + permute(y,[3,1,2])' * sin(thet),nx,ny,nth);
	y0 = reshape(permute(y,[3,1,2])' * cos(thet) - permute(x,[3,1,2])' * sin(thet),nx,ny,nth);
	out = exp( -0.5 * (x0.^2 + gam * y0.^2) / sig^2) .* cos(2 * pi * x0 / lamb);
endfunction
