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

load ../res/att_faces/facesScaled.mat
images = faces(:,:,1:1)+0.5;
[h w n_imgs] = size(images);

%% sparse coding
num_bases = 128;
beta = 0.1;
batch_size = 1000;
num_iters = 1;
sparsity_func= 'L1';
epsilon = [];
fname_save = [];
Binit = [];
windowsize = 8;

%%%%% Apply Laplacian pyramid algorithm to every image %%%%%

hgf = [1,4,6,4,1]/16; % demi filtre binomial à 5 points
[gfx,gfy] = meshgrid(hgf,hgf);
Gausfilt = gfx.*gfy; % filtre binomial
n_levels = 2; % nombre de niveaux de la pyramide.
for i=1:n_imgs
	[Lplac] = laplacianPyramidDecomposition(images(:,:,i), n_levels);
	for j=1:n_levels+1
		decompositions{j}(:,:,i) = Lplac{j}(:,:,1);
	end
end

%%%%% Run through each scale and learn a dictionnary %%%%%

arr = zeros(h, w, n_imgs)
for j=1:n_levels+1
	[hh ww nn] = size(decompositions{j});
	arr = zeros(hh, ww, nn);
	for i=1:n_imgs
		arr(:,:,i) = decompositions{j}(:,:,i);
	end
	%X{j} = getdata_imagearray_all(arr, windowsize);
	X{j} = getdata_imagearray(arr, windowsize, 1024);
	[B{j} S{j} stat{j}] = sparse_coding(X{j}, num_bases, beta, sparsity_func, epsilon, num_iters, batch_size, fname_save, Binit);
end

%%%%% Apply Laplacian pyramid algorithm to test image %%%%%

img = images(:,:,1);
[Lplac] = laplacianPyramidDecomposition(img, n_levels);

%%%%% For each scale, for each signal, reconstruct it using corresponding dictionnary %%%%%

for j=1:n_levels+1
	[I In Iout{j} Sout] = reconstruction(Lplac{j}, B{j}, beta, [], []);
end

%%%%% Reconstruct test image using pyramid %%%%%

imgout = laplacianPyramidReconstruction(Iout, n_levels);

%%%%% Display %%%%%

img = img-0.5;
imgout = imgout;%-0.5;
close all;
figure;
colormap gray;
imagesc(img, [-0.5 0.5]);
figure;
colormap gray;
imagesc(imgout, [-0.5 0.5]);
