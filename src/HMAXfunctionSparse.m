function [C2] = HMAXfunctionSparse(dicts, imgs, display)
	% Build the Gabor filters used by HMAX first layer (S1). For each of the 8 scales
	% the program can build the nth = 12 (here) filters corresponding to the different orientations.

    if nargin<3 || isempty(display)
        display = false;
    end

    [HMAXparams] = HMAXparameters();
    
    addpath('../fast-sc-2/code/');
    addpath('../fast-sc-2/code/sc2/');
    addpath('../fast-sc-2/code/sc2/nrf/');

    gaborFilters = getGaborFilters(HMAXparams, display); % build Gabor filters
    nHL = size(dicts{1},2);
    nimg = size(imgs,2);
    C2 = zeros(nimg,nHL);
    for imgcpt=1:nimg
        img = imgs{1};
        if size(img,3)==3
            img = double(rgb2gray(img));%/255; % convert it to grayscale
        end
        [dx,dy] = size(img);
        figure
        imshow(uint8(255*img)) % show original image

        %%%%%%%%%%%%
        %%%  S1  %%%
        %%%%%%%%%%%%
        %% Build the filters for the different orientations %%
        S1 = getS1(img, gaborFilters, HMAXparams, display);

        %%%%%%%%%%%%
        %%%  C1  %%%
        %%%%%%%%%%%%
        C1 = getC1(S1, dx, dy, HMAXparams, display);

        %%%%%%%%%%%%
        %%%  S2  %%%
        %%%%%%%%%%%%
        %%%%% S2 layer - get sparse coefficients matrices %%%%%
        %gamma=0.01; %TODO
        %gamma=0.1; %TODO
        gamma=0.4;
        Sout = cell(HMAXparams.nscales,1);
        for scal=1:HMAXparams.nscales-1
            %C1 normalization
            mn = min(min(min(C1{scal})));
            mx = max(max(max(C1{scal})));
            C1{scal} = ((C1{scal}-mn)/(mx-mn))-0.5;
            for j=1:HMAXparams.nth
                X = getdata_imagearray_all(C1{scal}(:,:,j), HMAXparams.filter_sz(scal));
                Sout{scal}(:,:,j) = l1ls_featuresign (dicts{scal}, X, gamma);
                %for
                    %Sout{scal}(:,:,j)Xb(:, cpt) = reshape(I(i:i+winsize-1, j:j+winsize-1),winsize^2,1);
                %Sout{scal}(:,:,j) = l1ls_featuresign (dicts{scal}, C1{scal}(:,:,j), gamma); % TODO same dict for every orientation ?
            end
        end

        %%%%%%%%%%%%
        %%%  C2  %%%
        %%%%%%%%%%%%
        %%%%% C2 layer - max response from the S2 layer %%%%%
        %TODO we don't need the C2 layer for reconstruction
        for scal=1:HMAXparams.nscales-1
            Sout{scal} = permute(Sout{scal},[3,2,1]);
        end
        C2tmp = getC2(Sout, nHL, HMAXparams);
        C2(imgcpt,:) = C2tmp';
    end
end

%close all;clear all;clc;display = false;img = imread('../res/Californie_m.JPG');load('HLfilters_2');tic;C2 = HMAXfunction(HLfilters, img, display);toc;
