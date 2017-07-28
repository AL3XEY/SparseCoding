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
nimgs = size(imds.Files,1)
images = cell(1,nimgs);
grayImages = cell(1,nimgs);
labels = zeros(1,nimgs);
label=imds.Labels(1);
labelc=1;
%nimgs
for i=(jeej-1)*100+1:jeej*100%nimgs
    %i
    file = imds.Files(i);
    img = imread(file{1});
    scaled = double(img)./255;
    [h,w,c]=size(scaled);
    [value,idx] = min([h,w]);
    nscales = 4; %TODO get HMAX parameters
    minFilterSz = 7;
    minPoolSz = 8;
    if (value/(2^nscales)) - minPoolSz < minFilterSz
        scaled = imresize(scaled, (value + minPoolSz)/(minFilterSz.*2^nscales));
        [h,w,c]=size(scaled);
    end
    if c==3
        images{i} = scaled;
        grayImages{i} = rgb2gray(scaled);
    else
        images{i} = zeros(h,w,3);
        images{i}(:,:,1) = scaled;
        images{i}(:,:,2) = scaled;
        images{i}(:,:,3) = scaled;
        grayImages{i} = scaled;
    end
    if imds.Labels(i)~=label
        label = imds.Labels(i);
        labelc = labelc+1;
    end
    labels(i)=labelc;
end

labelsString = imds.Labels;
   
clear 
toc