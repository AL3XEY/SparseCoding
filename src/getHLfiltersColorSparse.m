function [dicts] = getHLfiltersColorSparse(imgs, nHL, beta, iterations, oneDictPerScale, display)
    if nargin<2 || isempty(nHL)
        nHL = 10;
        %nHL = 1000; %the size of the dictionnary (number of atoms, or
        %prototypes)
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
    dicts = cell(HMAXparams.nscales-1,1);
    S = dicts;
    stat = dicts;
    nchans = 8; %FIXME
    C1 = cell(nimg*(HMAXparams.nscales-1)*nchans);
    nimg*(HMAXparams.nscales-1)*nchans

    %gaborFilters = getGaborFilters(HMAXparams, display); % build Gabor filters

    for imgcpt=1:nimg
		img = imgs{imgcpt};
        %if size(img,3)==3
		%	img = double(rgb2gray(img));%/255; % convert it to grayscale
        %end %TODO else, normalize between 0 and 1
        
		[dx,dy,c] = size(img);
		figure
		imshow(uint8(255*img)) % show original image

        %OpponencyCells
        SO = getSODescriptor(img, HMAXparams, display);
        DO = getDODescriptor(SO, HMAXparams, display);
        %TODO do we sum the complementary channels ?
        nchans = size(DO{1},4);
        for chan=1:nchans
            S1 = cell(1,HMAXparams.nscales);
            for scal=1:HMAXparams.nscales
               S1{scal}(:,:,:) = DO{scal}(:,:,:,chan);
            end
            C1tmp = getC1(S1, dx, dy, HMAXparams, display);
            %for scal=1:HMAXparams.nscales-1
            %    C1{scal}(:,:,(imgcpt-1)*(HMAXparams.nth)+1:(imgcpt)*(HMAXparams.nth)) = C1tmp{scal}(:,:,1:HMAXparams.nth);
            %end
            for scal=1:HMAXparams.nscales-1
                C1{(imgcpt-1)*(HMAXparams.nscales-1+nchans)+(chan-1)*(HMAXparams.nscales-1)+scal}(:,:,:) = C1tmp{scal}(:,:,:);
            end
        end
    end

    %%%%% Apply Sparse Coding algorithm %%%%%

    num_bases = nHL;
    %beta = 0.1; %0.01
    %TODO if beta empty, default value
    batch_size = 1000;
    num_iters = iterations;%2;
    %TODO if iterations empty, default value
    fname_save = [];
    Binit = [];
    sparsity_func= 'L1';
    epsilon = [];
    pars.display_images = true;%display;%true;%false;
    pars.display_every = 1;%int32(display);%1;%0;
    pars.save_every = 1;%int32(display);%1;%0;
    pars.save_basis_timestamps = false;%display;%false;%true;
    %winsize=8;
    if oneDictPerScale
        %FIXME
        for scal=1:HMAXparams.nscales-1 %FIXME
            %FIXME
            X = getdata_imagearray_all(C1{scal}, HMAXparams.filter_sz(scal)); %TODO or winsize ???)
            %C1 normalization
            mn = min(min(min(X)));
            mx = max(max(max(X)));
            X = ((X-mn)/(mx-mn))-0.5;
            [dicts{scal},S{scal},stat{scal}] = sparse_coding(X, num_bases, beta, sparsity_func, epsilon, num_iters, batch_size, fname_save, pars, Binit);
        end
    else
        cpt=0;
        for scalnimg=1:nimg*(HMAXparams.nscales-1)
            fsz = HMAXparams.filter_sz(1);
            [hh,ww,oo]=size(C1{scalnimg});
            hh = hh-fsz+1;
            ww = ww-fsz+1;
            cpt=cpt+ww*hh*oo-1;
        end
        X = zeros(fsz^2,cpt);
        oldValue=1;
        for scalnimg=1:nimg*(HMAXparams.nscales-1)
            [hh,ww,oo]=size(C1{scalnimg});
            hh = hh-fsz+1;
            ww = ww-fsz+1;
            newValue = oldValue+ww*hh*oo-1;
            X(:,oldValue:newValue) = getdata_imagearray_all(C1{scalnimg}, fsz);
            oldValue=newValue;
        end
        %C1 normalization
        mn = min(min(min(X)));
        mx = max(max(max(X)));
        X = ((X-mn)/(mx-mn))-0.5;
        [dicts,S,stat] = sparse_coding(X, num_bases, beta, sparsity_func, epsilon, num_iters, batch_size, fname_save, pars, Binit);
    end
end
