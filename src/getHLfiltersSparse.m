function [dicts] = getHLfilters(imgs, nHL, display)
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

    addpath('../fast-sc-2/code/');
    addpath('../fast-sc-2/code/sc2/');
    addpath('../fast-sc-2/code/sc2/nrf/');
        
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

		C1tmp = getC1(S1, dx, dy, HMAXparams, display);
        for scal=1:HMAXparams.nscales-1
            C1{scal}(:,:,(imgcpt-1)*(HMAXparams.nth)+1:(imgcpt)*(HMAXparams.nth)) = C1tmp{scal}(:,:,1:HMAXparams.nth);
        end
    end
    
    for scal=1:HMAXparams.nscales-1
		%%%%% Apply Sparse Coding algorithm %%%%%
        %C1 normalization
        mn = min(min(min(C1{scal})));
        mx = max(max(max(C1{scal})));
        C1{scal} = ((C1{scal}-mn)/(mx-mn))-0.5;

        num_bases = 128;
        beta = 0.1; %0.01
        batch_size = 1000;
        num_iters = 2;
        fname_save = [];
        Binit = [];
        sparsity_func= 'L1';
        epsilon = [];
        pars.display_images = true;%display;%true;%false;
        pars.display_every = 1;%int32(display);%1;%0;
        pars.save_every = 1;%int32(display);%1;%0;
        pars.save_basis_timestamps = false;%display;%false;%true;
        %winsize=8;
        %X = getdata_imagearray_all(C1{scal}, HMAXparams.filter_sz(scal)); %TODO or winsize ???
        X = getdata_imagearray(C1{scal}, HMAXparams.filter_sz(scal), 4096);
        [dicts{scal} S{scal} stat{scal}] = sparse_coding(X, num_bases, beta, sparsity_func, epsilon, num_iters, batch_size, fname_save, pars, Binit);
    end
end
