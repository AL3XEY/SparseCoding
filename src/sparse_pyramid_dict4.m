close all;
clear all;
clc;
if is_octave
	pkg load image;
end
addpath('../sepspyr/');
addpath('../fast-sc-2/code/');
addpath('../fast-sc-2/code/sc2/');
addpath('../fast-sc-2/code/sc2/nrf/');

%%%%% Parameters %%%%%

n_levels = 4; % nombre de niveaux de la pyramide.

%load ../res/att_faces/facesScaled.mat
%images = faces(:,:,1:1)+0.5;
load('../res/ZDB/zdbSmallScaled.mat');
images = zdb(:,:,1:1)+0.5;
[h w n_imgs] = size(images);

%% sparse coding
load('../fast-sc-2/results/dict/4/sc_L1_b128_beta0.4_20170228T202450.mat');
fooB = B;
clear B;
for j=1:n_levels+1
	%num_bases{j} = 128;
	beta{j} = 0.01; %0.05
	betaReconstruct{j} = beta{j};
	%batch_size{j} = 1000;
	%num_iters{j} = 3;
	%fname_save{j} = [];
	%Binit{j} = [];
	B{j} = fooB;
end
sparsity_func= 'L1';
epsilon = [];
windowsize = 8;
pars.display_images = false;
pars.display_every = 0;
pars.save_every = 0;
pars.save_basis_timestamps = true;

%%%%% Apply Laplacian pyramid algorithm to every image %%%%%

hgf = [1,4,6,4,1]/16; % demi filtre binomial à 5 points
[gfx,gfy] = meshgrid(hgf,hgf);
Gausfilt = gfx.*gfy; % filtre binomial
%for i=1:n_imgs
%	[Lplac] = laplacianPyramidDecomposition(images(:,:,i), n_levels);
%	for j=1:n_levels+1
%		decompositions{j}(:,:,i) = Lplac{j}(:,:,1);
%	end
%end

%%%%% Run through each scale and learn a dictionnary %%%%%

%arr = zeros(h, w, n_imgs)
%for j=1:n_levels+1
%	[hh ww nn] = size(decompositions{j});
%	arr = zeros(hh, ww, nn);
%	for i=1:n_imgs
%		arr(:,:,i) = decompositions{j}(:,:,i);
%	end
%	%X{j} = getdata_imagearray_all(arr, windowsize);
%	X{j} = getdata_imagearray(arr, windowsize, 1024);
%	%X{j} = getdata_imagearray(arr, windowsize, floor(1024/j/j/2));
%	[B{j} S{j} stat{j}] = sparse_coding(X{j}, num_bases{j}, beta{j}, sparsity_func, epsilon, num_iters{j}, batch_size{j}, fname_save{j}, Binit{j});
%end

%%%%% Apply Laplacian pyramid algorithm to test image %%%%%

img = images(:,:,1);
[Lplac] = laplacianPyramidDecomposition(img, n_levels);

%%%%% For each scale, for each signal, reconstruct it using corresponding dictionnary %%%%%

for j=1:n_levels+1
	[I In Iout{j} Sout{j} entropyI{j} entropyInoised entropyIout{j} PSNR_In PSNR_Iout{j} fresidue{j} fsparsity{j} sparsity{j}] = reconstruction(Lplac{j}, B{j}, betaReconstruct{j}, [], []);
end

%%%%% Reconstruct test image using pyramid %%%%%

imgout = laplacianPyramidReconstruction(Iout, n_levels);

%%%%% Display %%%%%

img = img-0.5;
imgout = imgout-0.5;
close all;
figure;
colormap gray;
imagesc(img, [-0.5 0.5]);
figure;
colormap gray;
imagesc(imgout, [-0.5 0.5]);
entropyImg = entropy(img)
entropyImgout = entropy(imgout)
PSNRImgout = psnr(imgout, img)

meanPSNR = 0;
meanSparsity = 0;
avgnzero = 0;
for j=1:n_levels+1
	meanPSNR = meanPSNR + PSNR_Iout{j};
	meanSparsity = meanSparsity + sparsity{j};
	avgnzero = avgnzero + sparsity{j}*size(B{j},2);
end
meanPSNR = meanPSNR / j
meanSparsity = meanSparsity / j
avgnzero = avgnzero / j %not really pertinent if levels have different dict size
