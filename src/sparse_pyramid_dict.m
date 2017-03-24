function [img imgout Lplac Iout meanPSNR meanSparsity avgnzero] = sparse_pyramid_dict(img, n_levels, B, beta)
if is_octave
	pkg load image;
end
addpath('../sepspyr/');
addpath('../fast-sc-2/code/');
addpath('../fast-sc-2/code/sc2/');
addpath('../fast-sc-2/code/sc2/nrf/');

%%%%% Parameters %%%%%

if nargin<2 || isempty(n_levels)
	n_levels = 2;
end
img = img+0.5;
[h w] = size(img);

%% sparse coding
if nargin<3 || isempty(B)
	load('../fast-sc-2/results/dict/4/sc_L1_b128_beta0.4_20170228T202450.mat');
	fooB = B;
	clear B;
	for j=1:n_levels+1
		beta{j} = 0.01; %0.05
		B{j} = fooB;
	end
	%TODO B could be a matrix, a list or a filename
end
if nargin<4 || isempty(beta)
	for j=1:n_levels+1
		beta{j} = 0.01;
	end
end

hgf = [1,4,6,4,1]/16; % demi filtre binomial à 5 points
[gfx,gfy] = meshgrid(hgf,hgf);
Gausfilt = gfx.*gfy; % filtre binomial

%%%%% Apply Laplacian pyramid algorithm to test image %%%%%

[Lplac] = laplacianPyramidDecomposition(img, n_levels);

%%%%% For each scale, for each signal, reconstruct it using corresponding dictionnary %%%%%

for j=1:n_levels+1
	[I In Iout{j} Sout{j} entropyI{j} entropyInoised entropyIout{j} PSNR_In PSNR_Iout{j} fresidue{j} fsparsity{j} sparsity{j}] = reconstruction(Lplac{j}, B{j}, beta{j}, [], []);
end

%%%%% Reconstruct test image using pyramid %%%%%

imgout = laplacianPyramidReconstruction(Iout, n_levels);

%%%%% Display %%%%%

for j=1:n_levels+1
	figure;
	colormap gray;
	imagesc(Lplac{j}, [0 1]);
end

for j=1:n_levels+1
	figure;
	colormap gray;
	imagesc(Iout{j}, [0 1]);
end

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
