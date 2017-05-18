function [C2] = HMAXfunction(HLfilters, imgs, display)
	% Build the Gabor filters used by HMAX first layer (S1). For each of the 8 scales
	% the program can build the nth = 12 (here) filters corresponding to the different orientations.
    
    if nargin<3 || isempty(display)
        display = false;
    end
    
    [HMAXparams] = HMAXparameters();
        
    gaborFilters = getGaborFilters(HMAXparams, display); % build Gabor filters
    
    nimg = size(imgs,2);
    nHL = size(HLfilters{1},4);
    C2 = zeros(nimg,nHL);
    for imgcpt=1:nimg
        img = imgs{1};
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
        S1 = getS1(img, gaborFilters, HMAXparams, display);

        %%%%%%%%%%%%
        %%%  C1  %%%
        %%%%%%%%%%%%
        C1 = getC1(S1, dx, dy, HMAXparams, display);

        %%%%%%%%%%%%
        %%%  S2  %%%
        %%%%%%%%%%%%
        %%%%% S2 layer - compute the response of the HLfilters (prototypes randomly taken from C1 layers of a large dataset) to every C1 patch %%%%%
        S2 = getS2(C1,HLfilters, HMAXparams, display);

        %%%%%%%%%%%%
        %%%  C2  %%%
        %%%%%%%%%%%%
        %%%%% C2 layer - max response from the S2 layer %%%%%
        C2tmp = getC2(S2, nHL, HMAXparams);
        C2(nimg,:) = C2tmp';
    end
end

%close all;clear all;clc;display = false;img = imread('../res/Californie_m.JPG');load('HLfilters_2');tic;C2 = HMAXfunction(HLfilters, img, display);toc;
