%read files and store them by category

close all
clear
clc
%tic

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
nimgs = size(imds.Labels,1);
[a,b]=hist(imds.Labels,unique(imds.Labels));
cpt=1;
for i=1:3%1:size(a,2) %TODO 1 2 3 5 7
    images = cell(1,a(i));
    grayImages = cell(1,a(i));
    labels = ones(1,a(i)).*i;
    labelString=b(i);
    labelString=labelString{1};
    for j=1:a(i)
        tic
        i
        j
        file = imds.Files(cpt);
        img = imread(file{1});
        scaled = double(img)./255;
        [h,w,c]=size(scaled);
        [value,idx] = min([h,w]);
        if value < 512
            scaled = imresize(scaled, 512/value);
            [h,w,c]=size(scaled);
        end
        if c==3
            images{j} = scaled;
            grayImages{j} = rgb2gray(scaled);
        else
            images{j} = zeros(h,w,3);
            images{j}(:,:,1) = scaled;
            images{j}(:,:,2) = scaled;
            images{j}(:,:,3) = scaled;
            grayImages{j} = scaled;
        end
        cpt=cpt+1;
    end
    toc
    tic
    save(sprintf('caltech101Category-%d-%s.mat',i,labelString),'images','labels','labelString', '-v7.3');
    images=grayImages;
    save(sprintf('caltech101CategoryGray-%d-%s.mat',i,labelString),'images','labels','labelString', '-v7.3');  
    toc
end
%toc