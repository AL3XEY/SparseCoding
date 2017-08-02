function [I,Iout,Sout,entropyI,entropyIout,PSNR,fresidue,fsparsity,sparsity] = reconstruction(img, datas, gamma, options, param) %TODO multiple options

% I : image used for reconstruction (noisy or not)
% Iout : image recovered using sparse coding
% Sout : output sparse coefficients
% img : filename (string), index in the dataset (integer 1 - 10) or actual image (matrix)
% datas : filename (string) for dictionnary or actual dictionnary (matrix)
% gamma : Lagrange multiplier used for sparseness optimization

%interpretation of the input parameters
if is_octave
	pkg load image;
	type = typeinfo(img);
	if strcmp(type, 'sq_string')
		I = imread(img);
	else
		I = img;
	end

	type = typeinfo(datas);
	if strcmp(type, 'sq_string')
		load(datas); % load dictionnary B
	else
		B = datas;
	end
else
	if ischar(img)
		I = imread(img);
	else
		I = img;
	end

	if ischar(datas)
		load(datas); % load dictionnary B
	else
		B = datas;
	end
end
[h,w,c] = size(I);
%TODO if c==3 and pars.color = 3 ?
if c==3
	I = rgb2lab2mat(I);
end

[szH,szW] = size(B);
winsize = sqrt(szH/c);
foo = h - winsize + 1;
bar = w - winsize + 1;
X = getdata_imagearray_all(I, winsize); %get all the patches in the image as vectors
[Xh,Xw] = size(X);
if (nargin<3) || isempty(gamma)
	if exist('pars') && ~isempty('pars')
    	gamma = pars.beta/pars.sigma*pars.noise_var; %get gamma from input data
	else
		gamma = 0.1; %default gamma %TODO or throw error?
	end
end
Sout = l1ls_featuresign (B, X, gamma); %get the sparse coefficients matrix. "too many coefficients activated" means that the dictionnary is bad (minimizes error sparsity but not sparsity) or the gamma parameter is too low.

Xout = B*Sout; % get output signals (patches) using sparse coefficients matrix
[XoutH,XoutW] = size(Xout);
Xb = zeros(XoutH, XoutW);
Iout = zeros(h,w,c);
meanCoef = zeros(h,w,c);

% get output image by reshaping the vectors into matrices + get mean value for each pixel
cpt = 1;
if c==3
	for i=1:foo %TODO rename!
		for j=1:bar %TODO rename!
			Iout(i:i+winsize-1, j:j+winsize-1,1) = Iout(i:i+winsize-1, j:j+winsize-1,1) + reshape(Xout(1:XoutH/3,cpt),winsize,winsize);
	        %Xb(:, cpt) = reshape(I(i:i+winsize-1, j:j+winsize-1),winsize^2,1);
	        Iout(i:i+winsize-1, j:j+winsize-1,2) = Iout(i:i+winsize-1, j:j+winsize-1,2) + reshape(Xout((XoutH/3)+1:(2*XoutH/3),cpt),winsize,winsize);
	        Iout(i:i+winsize-1, j:j+winsize-1,3) = Iout(i:i+winsize-1, j:j+winsize-1,3) + reshape(Xout((2*XoutH/3)+1:XoutH,cpt),winsize,winsize);

			meanCoef(i:i+winsize-1, j:j+winsize-1,:) = meanCoef(i:i+winsize-1, j:j+winsize-1,:)+1;
			cpt = cpt + 1;
		end
	end
else
	for i=1:foo %TODO rename!
		for j=1:bar %TODO rename!
			Iout(i:i+winsize-1, j:j+winsize-1) = Iout(i:i+winsize-1, j:j+winsize-1) + reshape(Xout(:,cpt),winsize,winsize);
	        Xb(:, cpt) = reshape(I(i:i+winsize-1, j:j+winsize-1),winsize^2,1);
			meanCoef(i:i+winsize-1, j:j+winsize-1) = meanCoef(i:i+winsize-1, j:j+winsize-1)+1;
			cpt = cpt + 1;
		end
	end
end

Iout = Iout ./ meanCoef;
if c==3
	min(min(min(Iout)))
	max(max(max(Iout)))
	Iout = mat2lab2rgb(Iout);
	min(min(min(Iout)))
	max(max(max(Iout)))
	Iout = uint8(255*Iout);
	%min(min(min(Iout)))
	%max(max(max(Iout)))
	%Iout = Iout-0.5;

	I=uint8(mat2lab2rgb(I)*255);

	'I'
	min(min(min(I)))
	max(max(max(I)))
	'Iout'
	min(min(min(Iout)))
	max(max(max(Iout)))
end

%%%%% Display and save the images, coefficients and stats %%%%%
figure;
if c==3
	imshow(I);
	imwrite(I, '../results/I.png');
	entropyI = entropy(I);
else
	colormap gray;
	imagesc(I, [-0.5 0.5]);
	imwrite(uint8((I+0.5)*255), '../results/I.png');
end

figure;
if c==3
	imshow(Iout);
	imwrite(Iout, '../results/Iout.png');
else
	colormap gray;
	imagesc(Iout, [-0.5 0.5]);
	imwrite(uint8((Iout+0.5)*255), '../results/Iout.png');
end
entropyIout = entropy(Iout);
if is_octave || ~verLessThan('matlab', '8.3') %if Matlab R2014a and above
	PSNR = psnr(Iout, I);
else
	PSNR = -1;
end

%Idiff = abs(I - Iout);
Idiff = I - Iout;
%figure;
%colormap(parula);
%colormap(jet);
%imagesc(Idiff, [-1 1]);
%saveas(gcf, '../results/Idiff.png');

if exist('Iedges')
	figure;
	colormap gray;
	imagesc(Iedges);
end

sparsity = sum(Sout(:)~=0)/length(Sout(:));
avgnzero = nnz(Sout)/size(Sout,2);
nzero = nnz(Sout);

nzeromin = size(Sout,1);
nzerominidx = -1;
nzeromax = 0;
nzeromaxidx = -1;
for j=1:size(Sout,2) %Xw
	tempnzero = nnz(Sout(:,j));
	if tempnzero <= nzeromin
		nzeromin = tempnzero;
		nzerominidx = j;
	end
	if tempnzero >= nzeromax
		nzeromax = tempnzero;
		nzeromaxidx = j;
	end
end

addpath('sc2');
if exist('pars') % && ~isempty(pars)
	[fobj, fresidue, fsparsity] = getObjective2(B, Sout, Xb, pars.sparsity_func, pars.noise_var, pars.beta, pars.sigma, pars.epsilon)
else
	[fobj, fresidue, fsparsity] = getObjective2(B, Sout, Xb, 'L1', '1', gamma, '1', [])
end
fprintf('entropyI = %g\n', entropyI);
fprintf('entropyIout = %g\n', entropyIout);
fprintf('PSNR = %g\n', PSNR);
fprintf('fresidue = %g\n', fresidue);
fprintf('fsparsity = %g\n', fsparsity);
fprintf('sparsity = %g\n', sparsity);
fprintf('avgnzero = %g\n', avgnzero);
fprintf('nzero = %g\n', nzero);
fprintf('nzeromin = %g\n', nzeromin);
fprintf('nzeromax = %g\n', nzeromax);
if(is_octave) %TODO or always save with -v7 or -v7.3 ?
	save('../results/images.mat', 'I', 'Iout', 'Idiff', 'Sout', 'gamma', '-v7');
	save('../results/stats.mat', 'entropyI', 'entropyIout', 'sparsity', 'fresidue', 'fsparsity', 'PSNR', 'nzero', 'nzeromin', 'nzeromax', '-v7');
else
	save('../results/images.mat', 'I', 'Iout', 'Idiff', 'Sout', 'gamma');
	save('../results/stats.mat', 'entropyI', 'entropyIout', 'sparsity', 'fresidue', 'fsparsity', 'PSNR', 'nzero', 'nzeromin', 'nzeromax');
end
