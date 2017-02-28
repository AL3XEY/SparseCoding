function [I Sout Iout] = reconstruction(img, datas, winsize)

%load('../results/sc_L1_b128_beta0.4_20170227T121443.mat')
%load('../results/sc_L1_b128_beta0.4_20170227T171851.mat')
load(datas);
load('../data/IMAGES_RAW.mat');
I = IMAGESr(:,:,img);
figure;
imshow(mat2gray(I));
%I = imnoise(I, 'gaussian');
sigma = 0.2;
I = I + sigma*randn(size(I));

foo = winsize^2;

figure;
imshow(mat2gray(I));
%X = getdata_imagearray(I, 8, 4096);
X = getdata_imagearray_all(I, 8);
Sout = l1ls_featuresign (B, X, 1);
Xout = B*Sout;

cpt = 1;
for i=1:foo
	for j=1:foo
		Iout((i-1)*winsize+1:i*winsize, (j-1)*winsize+1:j*winsize) = reshape(Xout(:,cpt),winsize,winsize);
		cpt = cpt+1;
	end
end
figure;
imshow(mat2gray(Iout))
