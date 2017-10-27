% Sample script

display = false;
nHL = 10;
sparseLearning = true;
sparseCoding = true;
color = true;

%%%%% Load images - 3 groups %%%%%

imgsHmaxTraining{1} = imread('../res/boat.png');
imgsHmaxTraining{2} = imread('../res/peppers.png');
imgsSVMTraining{1} = imread('../res/baboon.png');
imgsSVMTraining{2} = imread('../res/fruits.png');
classes(1) = 1;
classes(2) = 2;
imgsSVMTesting{1} = imread('../res/lena.ppm');
imgsSVMTesting{2} = imread('../res/kids.ppm');

%%%%% learn HL filters on group 1 %%%%%

HMAXparams = HMAXparameters();
gaborFilters = getGaborFilters(HMAXparams, display);
beta=0.1;iterations=2;oneDictPerScale=false;
%%%HLfilters = getHLfiltersFusion(imgsHmaxTraining, nHL, sparseLearning, color ,HMAXparams, gaborFilters, display);
HLfilters = getHLfilters(imgsHmaxTraining, nHL, HMAXparams, gaborFilters, display);
%HLfilters = getHLfiltersSparse(imgsHmaxTraining, nHL, beta, iterations, oneDictPerScale, HMAXparams, gaborFilters, display);
%%HLfilters = getHLfiltersColor(imgsHmaxTraining, nHL, display);
%%HLfilters = getHLfiltersColorSparse(imgsHmaxTraining, nHL, beta, iterations, oneDictPerScale, HMAXparams, gaborFilters, display);

%%%%% apply HMAX on group 2 & 3 %%%%%

gamma=0.1;
%%%trainingFeatures = HMAXfunctionFusion(imgsSVMTraining, HLfilters, sparseLearning, color, sparseCoding, HMAXparams, gaborFilters, display);
%%%testingFeatures = HMAXfunctionFusion(imgsSVMTesting, HLfilters, sparseLearning, color, sparseCoding, HMAXparams, gaborFilters, display);
trainingFeatures = HMAXfunction(imgsSVMTraining, HLfilters, HMAXparams, gaborFilters, display);
testingFeatures = HMAXfunction(imgsSVMTesting, HLfilters, HMAXparams, gaborFilters, display);
%trainingFeatures = HMAXfunctionSparse(imgsSVMTraining, HLfilters, gamma, HMAXparams, gaborFilters, display);
%testingFeatures = HMAXfunctionSparse(imgsSVMTesting, HLfilters, gamma, HMAXparams, gaborFilters, display);
%%trainingFeatures = HMAXfunctionColor(imgsSVMTraining, HLfilters, HMAXparams, gaborFilters, display);
%%testingFeatures = HMAXfunctionColor(imgsSVMTesting, HLfilters, HMAXparams, gaborFilters, display);
%trainingFeatures = HMAXfunctionColorSparse(imgsSVMTraining, HLfilters, gamma, HMAXparams, gaborFilters, display);
%testingFeatures = HMAXfunctionColorSparse(imgsSVMTesting, HLfilters, gamma, HMAXparams, gaborFilters, display);

%%%%% SVM Training on group 2 %%%%%

%load trainingFeatures.mat
SVMModel = fitcsvm(trainingFeatures, classes,'KernelFunction','rbf','BoxConstraint',1,'ClassNames',[0,1]); %TODO BoxConstraint shouldn't be set to 1 because we have multiple classes
%TODO cross-validation

%%%%% SVM testing on group 3 %%%%%

[label,score] = predict(SVMModel,testingFeatures)