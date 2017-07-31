%Take N images per category for learning and another N images for testing
%(or less if there are less images in the category)

clear

N=5; %15, 30
load('../res/101_ObjectCategories/caltechHLfiltersgray.mat');
HLfiltersFinal = cell(1,size(HLfilters,2)); %3
cpt=1;
step=9;
for i=1:step:size(HLfilters{1},4) %
    for scal=1:size(HLfilters,2) %3
        HLfiltersFinal{scal}(:,:,:,cpt) = HLfilters{scal}(:,:,:,i);
    end
    cpt=cpt+1;
end
C2Learning = zeros(N*102,size(HLfiltersFinal{1},4));
C2Test = zeros(N*102,size(HLfiltersFinal{1},4));
labelsLearning = zeros(N*102);
labelsTest = zeros(N*102);
cptLearning = 1;
cptTest = 1;
HLfilters = HLfiltersFinal;
clear HLfiltersFinal

%TODO some images might be used for HLfilters and learning/test features

rootFolder = fullfile('../res/101_ObjectCategories/');
d = dir(rootFolder);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
imds = imageDatastore(fullfile(rootFolder, nameFolds), 'LabelSource', 'foldernames');
label=imds.Labels(1);
labelc=1;
nimgs = size(imds.Labels,1);
[a,b]=hist(imds.Labels,unique(imds.Labels));
oldCpt=0;
for i=1:size(b,2)%102
    tic
    i
    labelString=b(i);
    labelString=labelString{1};
    load(sprintf('../res/101_ObjectCategories/caltech101CategoryGray-%d-%s.mat',i,labelString));
    for j=1:N
        imagesLearning{j} = images{j};
        imagesTest{j} = images{N+j};
    end
    newCpt = oldCpt + size(imagesLearning,2); %size = N
    'entering C2 Learning'
    C2tmpLearning = HMAXfunction(HLfilters,imagesLearning,false);
    'entering C2 Test'
    C2tmpTest = HMAXfunction(HLfilters,imagesTest,false);
    C2Learning(oldCpt+1:newCpt,:) = C2tmpLearning;
    C2Test(oldCpt+1:newCpt,:) = C2tmpTest;
    labelsLearning(oldCpt+1:newCpt) = i;%labels(1:2:size(images,2));
    labelsTest(oldCpt+1:newCpt) = i;%labels(2:2:size(images,2));
    oldCpt = newCpt;

%     for j=1:size(images,2)
%         i
%         j
%         %imgs{1}=images{j};
%         C2tmp = HMAXfunction(HLfilters,imgs,false);
%         if mod(j,3)==2 %Training
%             C2Learning(cptLearning,:) = C2tmp;
%             labelsLearning(cptLearning)=labels(j);
%             %labelsStringLearning(cptLearning)=labelsString(j);
%             cptLearning = cptLearning+1;
%         end
%         if mod(j,3)==0 %Test
%             C2Test(cptTest,:) = C2tmp;
%             labelsTest(cptTest)=labels(j);
%             %labelsStringTest(cptTest)=labelsString(j);
%             cptTest = cptTest+1;
%         end
%     end
    toc
end
labelsString=b;
save('../res/101_ObjectCategories/caltech101FeaturesLearning.mat', 'C2Learning', 'labelsLearning','labelsString');
save('../res/101_ObjectCategories/caltech101FeaturesTest.mat', 'C2Test', 'labelsTest','labelsString');