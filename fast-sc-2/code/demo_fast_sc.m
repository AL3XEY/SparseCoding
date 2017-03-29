function demo_fast_sc(opt_choice)
% opt_choice = 1: use L1 penalty
% opt_choice = 2: use epslion-L1 penalty

if ~exist('opt_choice', 'var')
    opt_choice = 1;
end

if is_octave
    pkg load image;
end

% image data
%load ../data/IMAGES.mat
%X = getdata_imagearray(IMAGES, 14, 10000);
%load ../data/IMAGES_RAW.mat
%X = getdata_imagearray(IMAGESr, 8, 4096);
%X = getdata_imagearray(IMAGESr, 8, 100000);
%X = getdata_imagearray_all(IMAGESr, 8);
%X = getdata_imagearray_all(IMAGESr(:,:,4), 8);
load ../../res/att_faces/facesScaled.mat
winsize = 8;
X = getdata_imagearray(faces, winsize, 4096);
%X = getdata_imagearray_all(faces(:,:,1:20), winsize);
%load ../../res/ZDB/zdbSmallScaled.mat
%X = getdata_imagearray(zdb, 8, 4096);

% sparse coding parameters
num_bases = 128; %number of atoms in the dictionnary
beta = 0.1; %Lagrange multiplier
batch_size = 1000;
num_iters = 10; %number of iterations of the learning algorithm
if opt_choice==1
    sparsity_func= 'L1';
    epsilon = [];
elseif opt_choice==2
    sparsity_func= 'epsL1';
    epsilon = 0.01;
end

Binit = []; %Binit = getdata_imagearray(faces, 8, num_bases);
%Bi = rand(winsize^2,num_bases)-0.5;
%Bi = Bi - repmat(mean(Bi,1), size(Bi,1),1);
%Bi = Bi*diag(1./sqrt(sum(Bi.*Bi)));
%%Bi(:,1) = ones(winsize^2,1)*(-0.126);
%load ../results/dict/19/19.mat
%Bi(:,1) = B(:,129);
%Binit = Bi;

fname_save = sprintf('../results/sc_%s_b%d_beta%g_%s', sparsity_func, num_bases, beta, datestr(now, 30));

% run fast sparse coding
[B S stat] = sparse_coding(X, num_bases, beta, sparsity_func, epsilon, num_iters, batch_size, fname_save, [], Binit);
