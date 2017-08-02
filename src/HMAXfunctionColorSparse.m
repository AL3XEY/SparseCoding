function [C2] = HMAXfunctionColorSparse(dicts, imgs, gamma, display)
    if nargin<3 || isempty(display)
        display = false;
    end

    [HMAXparams] = HMAXparameters();

    addpath('../fast-sc-2/code/');
    addpath('../fast-sc-2/code/sc2/');
    addpath('../fast-sc-2/code/sc2/nrf/');

    gaborFilters = getGaborFilters(HMAXparams, display); % build Gabor filters
    if iscell(dicts)
        ndicts = size(dicts,2);
        nHL = size(dicts{1},2);
    else
        ndicts = 1;
        nHL = size(dicts,2);
    end
    nimg = size(imgs,2);
    C2 = zeros(nimg,nHL);
    nchans = 8; %FIXME
    for imgcpt=1:nimg
        img = imgs{imgcpt};
        %if size(img,3)==3
        %    img = double(rgb2gray(img));%/255; % convert it to grayscale
        %end
        [dx,dy,c] = size(img);
        figure
        imshow(uint8(255*img)) % show original image

        %%%%%%%%%%%%
        %%%  S1  %%%
        %%%%%%%%%%%%
%        S1 = getS1(img, gaborFilters, HMAXparams, display);
        %OpponencyCells
        SO = getSODescriptor(img, HMAXparams, display);
        DO = getDODescriptor(SO, HMAXparams, display);
        %TODO do we sum the complementary channels ?
        %%%%%%%%%%%%
        %%%  C1  %%%
        %%%%%%%%%%%%
        %nchans = size(DO{1},4);
        C1 = cell(nimg*(HMAXparams.nscales-1)*nchans);
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
                C1{(chan-1)*(HMAXparams.nscales-1)+scal}(:,:,:) = C1tmp{scal}(:,:,:);
            end
        end
%        C1 = getC1(S1, dx, dy, HMAXparams, display);

        %%%%%%%%%%%%
        %%%  S2  %%%
        %%%%%%%%%%%%
        %%%%% S2 layer - get sparse coefficients matrices %%%%%
        %gamma=0.01; %TODO
        %gamma=0.1; %TODO
        %gamma=0.4; %TODO default value

        %Sout = cell(HMAXparams.nscales,1);
        %for scal=1:HMAXparams.nscales-1
        %    %C1 normalization
        %    mn = min(min(min(C1{scal})));
        %    mx = max(max(max(C1{scal})));
        %    C1{scal} = ((C1{scal}-mn)/(mx-mn))-0.5;
        %    for j=1:HMAXparams.nth
        %        if iscell(dicts)
        %            dictionnary = dicts{scal};
        %        else
        %            dictionnary = dicts;
        %        end
        %        X = getdata_imagearray_all(C1{scal}(:,:,j), sqrt(size(dictionnary,1)));
        %        Sout{scal}(:,:,j) = l1ls_featuresign(dictionnary, X, gamma);
        %        %for
        %            %Sout{scal}(:,:,j)Xb(:, cpt) = reshape(I(i:i+winsize-1, j:j+winsize-1),winsize^2,1);
        %        %Sout{scal}(:,:,j) = l1ls_featuresign (dicts{scal}, C1{scal}(:,:,j), gamma); % TODO same dict for every orientation ?
        %    end
        %end

        if ndicts==1
            cpt=0;
            for scal=1:HMAXparams.nscales-1
                for chan=1:nchans
                    fsz = HMAXparams.filter_sz(1);
                    [hh,ww,oo]=size(C1{(chan-1)*(HMAXparams.nscales-1)+scal});
                    hh = hh-fsz+1;
                    ww = ww-fsz+1;
                    cpt=cpt+ww*hh*oo-1;
                end
            end
            X = zeros(fsz^2,cpt);
            oldValue=1;
            for scal=1:HMAXparams.nscales-1
                for chan=1:nchans
                    [hh,ww,oo]=size(C1{(chan-1)*(HMAXparams.nscales-1)+scal});
                    hh = hh-fsz+1;
                    ww = ww-fsz+1;
                    newValue = oldValue+ww*hh*oo-1;
                    X(:,oldValue:newValue) = getdata_imagearray_all(C1{(chan-1)*(HMAXparams.nscales-1)+scal}, fsz);
                    oldValue=newValue;
                end
            end
            %C1 normalization
            mn = min(min(min(X)));
            mx = max(max(max(X)));
            X = ((X-mn)/(mx-mn))-0.5;
            Sout = l1ls_featuresign (dicts, X, gamma);
        end %TODO other case

        %%%%%%%%%%%%
        %%%  C2  %%%
        %%%%%%%%%%%%
        %%%%% C2 layer - highest activation of every S2 atom %%%%%
        %for scal=1:HMAXparams.nscales-1
        %    Sout{scal} = permute(Sout{scal},[3,2,1]);
        %end
        %C2tmp = getC2(Sout, nHL, HMAXparams);
        %C2(imgcpt,:) = C2tmp';
        %%save('Sout.mat','Sout');
        if ndicts==1
            for i=1:size(dicts,2)
                C2(imgcpt,i) = max(abs(Sout(i,:)));
            end
        end
    end
end

%close all;clear all;clc;display = false;img = imread('../res/Californie_m.JPG');load('HLfilters_2');tic;C2 = HMAXfunction(HLfilters, img, display);toc;
