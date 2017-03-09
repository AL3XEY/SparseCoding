function [I Sout Iout] = reconstruction(img, datas, options, param) %TODO multiple options

% I : image used for reconstruction (noisy or not)
% Sout : output sparse coefficients
% Iout : image recovered using sparse coding
% img : filename (string), index in the dataset (integer 1 - 10) or actual image (matrix)
% datas : filename (string) for dictionnary or actual dictionnary (matrix)
% TODO options
% TODO param

if is_octave
	pkg load image;
	type = typeinfo(img);
	if type == 'scalar'
		load('../data/IMAGES_RAW.mat');
		I = IMAGESr(:,:,img);
	elseif type == 'sq_string'
		I = imread(img);
	elseif type == 'diagonal matrix' %TODO useful?
		I = img;
	else
		error('img is not a filename nor an index number nor a matrix') %TODO useful?
	end

	type = typeinfo(datas);
	if type == 'sq_string'
		load(datas); % load dictionnary B
	elseif type == 'diagonal matrix' %TODO useful?
		B = datas;
	else
		error('img is not a filename nor a matrix') %TODO useful?
	end
else
	%if isinteger(img) && size(img) == 1
	%	load('../data/IMAGES_RAW.mat');
	%	I = IMAGESr(:,:,img);
	if ischar(img)
		I = imread(img);
	elseif isnumeric(img)
		I = img;
		if size(img) == 1
			load('../data/IMAGES_RAW.mat');
			I = IMAGESr(:,:,img);
		end
	else
		error('img is not a filename nor an index number nor a matrix') %TODO useful?
	end

	if ischar(datas)
		load(datas); % load dictionnary B
	elseif isnumeric(datas) %TODO useful?
		B = datas;
	else
		error('img is not a filename nor a matrix') %TODO useful?
	end
end
[h w] = size(I);

%%%%%%%%%%%%%%%%%
%%%%% NOISE %%%%%
%%%%%%%%%%%%%%%%%
In = I;
if exist(options) && ~empty(options) && strcmp(options, 'noise')
	%I = imnoise(I, 'gaussian');

	%sigma = 0.2;
	%I = I + sigma*randn(size(I));

	%randnoise = reshape(round(rand(512^2,1)),512,512);
	%I = I.*randnoise;

	%https://fr.mathworks.com/help/stats/binornd.html?requestedDomain=www.mathworks.com
end
%%%%%%%%%%%%%%%%%

[szH szW] = size(B);
winsize = sqrt(szH);
foo = h - winsize + 1;
bar = w - winsize + 1;
X = getdata_imagearray_all(In, 8);
[Xh Xw] = size(X);
Sout = l1ls_featuresign (B, X, 1);

if ~(nargin<3)%exist(options) && ~empty(options)
	if strcmp(options, 'remove_last')
		if (nargin<4)%~exist(param) || empty(param)
			param = 1;
		end
		for i=1:param
			for j=1:Xw
				nzeros = find(Sout(:,j));
				tmp = zeros(size(nzeros,1),2);
				tmp(:,1) = nzeros(:);
				nzeros = tmp;
				for k=1:size(nzeros)
					nzeros(k,2) = Sout(nzeros(k,1),j);
				end
				[value index] = min(abs(nzeros(:,2)));
				Sout(nzeros(index,1),j) = 0;
			end
		end
	end
%	if strcmp(options, 'add_random')
%		if nargin<4
%			param = 1;
%		end
%		for i=1:param
%			for j=1:szW
%
%			end
%		end
%	end
end

Xout = B*Sout;
Iout = zeros(h,w);
meanCoef = zeros(h,w);

cpt = 1;
for i=1:foo
	for j=1:bar
		Iout(i:i+winsize-1, j:j+winsize-1) = Iout(i:i+winsize-1, j:j+winsize-1) + reshape(Xout(:,cpt),winsize,winsize);
		meanCoef(i:i+winsize-1, j:j+winsize-1) = meanCoef(i:i+winsize-1, j:j+winsize-1)+1;
		cpt = cpt + 1;
	end
end

Iout = Iout ./ meanCoef;

%%%%% Show and save the images and coefficients %%%%%

figure;
imshow(mat2gray(I));
imwrite(I, '../results/I.png');
entropyI = entropy(I)

figure;
imshow(mat2gray(In));
imwrite(In, '../results/In.png');
entropyInoised = entropy(In)
if is_octave
	PSNR_In = psnr(In, I) %FIXME not available for Matlab version < R2014a
end

figure;
imshow(mat2gray(Iout))
imwrite(Iout, '../results/Iout.png');
entropyIout = entropy(Iout)
if is_octave
	PSNR_Iout = psnr(Iout, I) %FIXME not available for Matlab version < R2014a
end
sparsity = sum(Sout(:)~=0)/length(Sout(:));
fprintf('sparsity = %g\n', sparsity);
save('Sout.mat', 'Sout', 'sparsity');
