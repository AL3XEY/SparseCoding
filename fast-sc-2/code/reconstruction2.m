function [I Sout Iout] = reconstruction(img, datas, winsize)

%load('../results/sc_L1_b128_beta0.4_20170227T121443.mat')
%load('../results/sc_L1_b128_beta0.4_20170227T171851.mat')
load(datas);
load('../data/IMAGES_RAW.mat');
I = IMAGESr(:,:,img);
[h w] = size(I);
figure;
imshow(mat2gray(I));

%I = imnoise(I, 'gaussian');
sigma = 0.2;
I = I + sigma*randn(size(I));

foo = h - winsize + 1;

figure;
imshow(mat2gray(I));
X = getdata_imagearray_all2(I, 8);
Sout = l1ls_featuresign (B, X, 1);
Xout = B*Sout;
meanCoef = zeros(h,w);

cpt = 1;
for i=1:foo
	for j=1:foo
		Iout(i:i+winsize-1, j:j+winsize-1) = reshape(Xout(:,cpt),winsize,winsize);
		meanCoef(i:i+winsize-1, j:j+winsize-1) = meanCoef(i:i+winsize-1, j:j+winsize-1)+1;
		cpt = cpt+1;
	end
end

Iout = Iout ./ meanCoef;

figure;
imshow(mat2gray(Iout))
