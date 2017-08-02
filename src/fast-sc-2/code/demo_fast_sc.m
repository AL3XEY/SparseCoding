function demo_fast_sc(opt_choice)
% opt_choice = 1: use L1 penalty
% opt_choice = 2: use epslion-L1 penalty

if ~exist('opt_choice', 'var')
    opt_choice = 1;
end

if is_octave
    pkg load image;
end

winsize = 8;

%imgs(:,:,:,1) = rgb2lab2mat(imread('../../../res/lena.ppm'));
imgs(:,:) = double(imread('../../../res/lena.pgm'))./255-0.5;

channels = size(imgs,3);
%X = getdata_imagearray_all(imgs, winsize)
X = getdata_imagearray(imgs, winsize, 4096);

% sparse coding parameters
num_bases = 128; %number of atoms in the dictionnary
beta = 0.1; %Lagrange multiplier
batch_size = 1000;
num_iters = 2; %number of iterations of the learning algorithm
if opt_choice==1
    sparsity_func= 'L1';
    epsilon = [];
elseif opt_choice==2
    sparsity_func= 'epsL1';
    epsilon = 0.01;
end
Binit = [];
pars = [];
fname_save = sprintf('../results/sc_%s_b%d_beta%g_%s', sparsity_func, num_bases, beta, datestr(now, 30));

% run fast sparse coding
[B,S,stat] = sparse_coding(X, num_bases, beta, sparsity_func, epsilon, num_iters, batch_size, fname_save, pars, Binit, channels);
