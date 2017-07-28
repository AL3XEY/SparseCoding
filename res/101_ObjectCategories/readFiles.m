close all
clear
clc
tic

%if is_octave
%    pkg load image;
%end

rootFolder = fullfile('./');
d = dir(rootFolder);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
imds = imageDatastore(fullfile(rootFolder, nameFolds), 'LabelSource', 'foldernames');
%tbl = countEachLabel(imds);
label=imds.Labels(1);
labelc=1;
nimgs = 100;
for jeej = 36:91
    images = cell(1,nimgs);
    grayImages = cell(1,nimgs);
    labels = zeros(1,nimgs);
    labelsString = zeros(1,nimgs);
    for i=1:100%nimgs
        j = (jeej-1)*100+i;
        file = imds.Files(j);
        img = imread(file{1});
        scaled = double(img)./255;
        [h,w,c]=size(scaled);
        [value,idx] = min([h,w]);
        %nscales = 4; %TODO get HMAX parameters
        %minFilterSz = 7;
        %minPoolSz = 8;
        %if (value/(2^nscales)) - minPoolSz < minFilterSz
        %    scaled = imresize(scaled, (value + minPoolSz*nscales)/(minFilterSz.*2^nscales)); %TODO BAD FORMULA
        %    [h,w,c]=size(scaled);
        %end
        if value < 512
            scaled = imresize(scaled, 512/value);
            [h,w,c]=size(scaled);
        end
        if c==3
            images{i} = scaled;
            %grayImages{i} = rgb2gray(scaled);
        else
            images{i} = zeros(h,w,3);
            images{i}(:,:,1) = scaled;
            images{i}(:,:,2) = scaled;
            images{i}(:,:,3) = scaled;
            %grayImages{i} = scaled;
        end
        if imds.Labels(j)~=label
            label = imds.Labels(j);
            labelc = labelc+1;
        end
        labels(i)=labelc;
        labelsString(i) = imds.Labels(j);
    end
    save(sprintf('caltech101-%d.mat',jeej),'images','labels','labelsString', '-v7.3');
    %images=grayImages;
    %save(sprintf('caltech101gray-%d.mat',jeej),'images','labels','labelsString', '-v7.3');
    jeej
end
jeej = 92
nimgs = size(imds.Files,1)
%nimgs = 10;
images = cell(1,nimgs-9100);
grayImages = cell(1,nimgs-9100);
labels = zeros(1,nimgs-9100);
labelsString = zeros(1,nimgs-9100);
for i=1:nimgs-9100
    j = (jeej-1)*100+i;
    file = imds.Files(j);
    img = imread(file{1});
    scaled = double(img)./255;
    [h,w,c]=size(scaled);
    [value,idx] = min([h,w]);
    %nscales = 4; %TODO get HMAX parameters
    %minFilterSz = 7;
    %minPoolSz = 8;
    %if (value/(2^nscales)) - minPoolSz < minFilterSz
    %    scaled = imresize(scaled, (value + minPoolSz*nscales)/(minFilterSz.*2^nscales)); %TODO BAD FORMULA
    %    [h,w,c]=size(scaled);
    %end
    if value < 512
        scaled = imresize(scaled, 512/value);
        [h,w,c]=size(scaled);
    end
    if c==3
        images{i} = scaled;
        %grayImages{i} = rgb2gray(scaled);
    else
        images{i} = zeros(h,w,3);
        images{i}(:,:,1) = scaled;
        images{i}(:,:,2) = scaled;
        images{i}(:,:,3) = scaled;
        %grayImages{i} = scaled;
    end
    if imds.Labels(j)~=label
        label = imds.Labels(j);
        labelc = labelc+1;
    end
    labels(i)=labelc;
    labelsString(i) = imds.Labels(j);
end
save(sprintf('caltech101-%d.mat',jeej),'images','labels','labelsString', '-v7.3');
%images=grayImages;
%save(sprintf('caltech101gray-%d.mat',jeej),'images','labels','labelsString', '-v7.3');
toc