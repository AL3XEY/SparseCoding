function [C2] = HMAXfunction(imgs, HLfilters, HMAXparams, gaborFilters, display)
    if nargin<3 || isempty(HMAXparams)
        HMAXparams = HMAXparameters();
    end
    if nargin<4 || isempty(gaborFilters)
        gaborFilters = getGaborFilters(HMAXparams, false);
    end
    if nargin<5 || isempty(display)
        display = false;
    end
    if display
        shapeInserter = vision.ShapeInserter; %to draw rectangles around high-contribution areas on image
    end

    nimg = size(imgs,2);
    nHL = size(HLfilters{1},4);
    C2 = zeros(nimg,nHL);
    for imgcpt=1:nimg
        imgcpt
        img = imgs{imgcpt};
        if size(img,3)==3
            img = double(rgb2gray(img));%/255; % convert it to grayscale
        end
        [dx,dy] = size(img);
        if display
            figure
            imshow(uint8(255*img)) % show original image
        end
        %%%%%%%%%%%%
        %%%  S1  %%%
        %%%%%%%%%%%%
        %% Build the filters for the different orientations %%
        S1 = getS1(img, HMAXparams, gaborFilters, display);

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
        [C2tmp,coords] = getC2(S2, HMAXparams);
        C2(imgcpt,:) = C2tmp;

        if display
            %Contribution mapping
            %For each C2 value (high response to filters in S2), we highlight the region of the original image that later gave this result.
            for hl=1:nHL
                h = size(img,1);
                truc =coords(hl,:);
                maxX = truc(1);
                maxY = truc(2);
                scal = truc(3);
                h2 = size(S2{scal},1);
                ratio = h/h2;
                XStart = (maxX-1)*ratio;
                XEnd = (maxX+1)*ratio;
                YStart = (maxY-1)*ratio;
                YEnd = (maxY+1)*ratio;
                PTS = [[XStart YStart] [XEnd YEnd]];
                J = step(shapeInserter,img,PTS);
                figure;
                imshow(J);
                %TODO Octave :
                %hold ("on");
                %imshow (I);
                %plot (col, row, "ro");
                %hold ("off");
            end
        end
    end
end

%close all;clear all;clc;display = false;img = imread('../res/Californie_m.JPG');load('HLfilters_2');tic;C2 = HMAXfunction(HLfilters, img, display);toc;
