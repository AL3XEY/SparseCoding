function [I,In,Iout,Sout,entropyI,entropyInoised,entropyIout,PSNR_In,PSNR_Iout,fresidue,fsparsity,sparsity] = reconstruction(img, datas, gamma, options, param) %TODO multiple options

% I : image used for reconstruction (noisy or not)
% Iout : image recovered using sparse coding
% Sout : output sparse coefficients
% img : filename (string), index in the dataset (integer 1 - 10) or actual image (matrix)
% datas : filename (string) for dictionnary or actual dictionnary (matrix)
% gamma : Lagrange multiplier used for sparseness optimization
% options : string, either 'remove_first', 'remove_last', 'add_random', 'noise' or 'edges' (see in the code below)
% param : parameter associated with the option selected

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
I = rgb2lab2mat(I);

%%%%% Noise addition option %%%%%
In = I;
if strcmp(options, 'noise') && ~(nargin<4) %exist(options) && ~isempty(options)
	%Type 1
	In = imnoise(I+0.5, 'speckle', 0.001)-0.5;

	%In = I + (randn(h, w)-0.5)/100; %random noise generation

	%randnoise = reshape(round(rand(512^2,1)),512,512);
	%In = I.*randnoise;

	%https://fr.mathworks.com/help/stats/binornd.html?requestedDomain=www.mathworks.com
end

[szH,szW] = size(B);
winsize = sqrt(szH/3);
foo = h - winsize + 1;
bar = w - winsize + 1;
X = getdata_imagearray_all(In, winsize); %get all the patches in the image as vectors
[Xh,Xw] = size(X);
if (nargin<3) || isempty(gamma)
	if exist('pars') && ~isempty('pars')
    	gamma = pars.beta/pars.sigma*pars.noise_var; %get gamma from input data
	else
		gamma = 0.1; %default gamma %TODO or throw error?
	end
end
Sout = l1ls_featuresign (B, X, gamma); %get the sparse coefficients matrix. "too many coefficients activated" means that the dictionnary is bad (minimizes error sparsity but not sparsity) or the gamma parameter is too low.

%%%%% Other options : remove_first, remove_last, add_random and edges %%%%%
if ~(nargin<4) && ~isempty(options)%exist(options) && ~isempty(options)
	if strcmp(options, 'remove_last') %remove the N coefficients with the smallest values (where N is the value of input parameter param, default = 1)
		if (nargin<5) && ~isempty(param)%~exist(param) || isempty(param)
			param = 1;
		end
		for i=1:param
			for j=1:Xw
				nzeros = find(Sout(:,j));
                if ~isempty(nzeros)
                    tmp = zeros(size(nzeros,1),2);
                    tmp(:,1) = nzeros(:);
                    nzeros = tmp;
                    for k=1:size(nzeros)
                        nzeros(k,2) = Sout(nzeros(k,1),j);
                    end
                    [value,index] = min(abs(nzeros(:,2)));
                    Sout(nzeros(index,1),j) = 0;
                end
			end
		end
	end
	if strcmp(options, 'remove_first') %remove the N coefficients with the highest values (where N is the value of input parameter param, default = 1)
		if (nargin<5) && ~isempty(param)%~exist(param) || isempty(param)
			param = 1;
		end
		for i=1:param
			for j=1:Xw
				nzeros = find(Sout(:,j));
                if ~isempty(nzeros)
                    tmp = zeros(size(nzeros,1),2);
                    tmp(:,1) = nzeros(:);
                    nzeros = tmp;
                    for k=1:size(nzeros)
                        nzeros(k,2) = Sout(nzeros(k,1),j);
                    end
                    [value,index] = max(abs(nzeros(:,2)));
                    Sout(nzeros(index,1),j) = 0;
                end
			end
		end
	end
	if strcmp(options, 'add_random') %add N random coefficients with the random values (lesser than the actual max coefficient value) (where N is the value of input parameter param, default = 1)
		if nargin<5 && ~isempty(param)
			param = 1;
		end
		for i=1:param
			for j=1:Xw
                zero = find(Sout(:,j)==0);
                if ~isempty(zeros)
                    nzeros = find(Sout(:,j));
                    if ~isempty(nzeros)
                        tmp = zeros(size(nzeros,1),2);
                        tmp(:,1) = nzeros(:);
                        nzeros = tmp;
                        for k=1:size(nzeros)
                            nzeros(k,2) = Sout(nzeros(k,1),j);
                        end
                        maxv = max(max(max(abs(nzeros(:,2))))); %TODO /2 ?
                    else
                        maxv = 1;
                    end
                    spl = datasample(zero, 1);
                    Sout(spl,j) = (rand()-rand())*maxv;
					%Sout(spl,j) = (rand()-rand())*maxv/2;
                end
			end
		end
	end
	if strcmp(options, 'edges') %TODO WIP
		%if nargin<5 && ~isempty(param)
		%	param = 1;
		%end
		Iedges = edge(I, 'sobel');
		Xedges = getdata_imagearray_all(Iedges, winsize);
		method = 1;
		for j=1:Xw
			%if(~isempty(imgEdges(imgEdges==1)))
			%if(~isempty(find(Xedges(:,j)~=0)))
			if(sum(sum(Xedges(:,j))) > 2)
				%disp(j)

				%%%
				if method==1
					nzeros = find(Sout(:,j));
	                if ~isempty(nzeros)
						Btmp = B;
		            	tmp = zeros(size(nzeros,1),2);
		            	tmp(:,1) = nzeros(:);
		            	nzeros = tmp;
		            	for k=1:size(nzeros)
		            		nzeros(k,2) = Sout(nzeros(k,1),j);
		            	end
		                [value index] = max(abs(nzeros(:,2)));
						Btmp(:,nzeros(index,1)) = 0;
						%Souttmp = l1ls_featuresign (Btmp, X(:,j), gamma);
						%Souttmp = l1ls_featuresign (Btmp, X(:,j), gamma*10);
						Souttmp = l1ls_featuresign (Btmp, X(:,j), gamma/10);
						Sout(:,j) = Souttmp(:, 1);
	                end
				end
				%%%

				%%%
				if method == 2 || method == 3 %remove first / last, but only on the edges
					for i=1:param
						nzeros = find(Sout(:,j));
		                if ~isempty(nzeros)
			                tmp = zeros(size(nzeros,1),2);
			                tmp(:,1) = nzeros(:);
			                nzeros = tmp;
			                for k=1:size(nzeros)
			                    nzeros(k,2) = Sout(nzeros(k,1),j);
			                end
							if method == 2 %remove first
			                	[value,index] = max(abs(nzeros(:,2)));
							else %remove last
								[value,index] = min(abs(nzeros(:,2)));
							end
							Sout(nzeros(index,1),j) = 0;
		                end
					end
				end
				%%%

				%Sout(:,j) =
            end
		end
	end
end

Xout = B*Sout; % get output signals (patches) using sparse coefficients matrix
[XoutH,XoutW] = size(Xout);
Xb = zeros(XoutH, XoutW);
Iout = zeros(h,w,c);
meanCoef = zeros(h,w,c);

% get output image by reshaping the vectors into matrices + get mean value for each pixel
cpt = 1;
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

Iout = Iout ./ meanCoef;
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

%%%%% Display and save the images, coefficients and stats %%%%%
figure;
%imshow(mat2lab2rgb(I));
imshow(I);
%imshow(mat2gray(I)), [-0.5 0.5]));
%colormap gray;
%imagesc(I, [-0.5 0.5]);
%imwrite(uint8((mat2lab2rgb(I))*255), '../results/I.png');
imwrite(I, '../results/I.png');
entropyI = entropy(I);

figure;
%imshow(mat2lab2rgb(In));
In=I;
imshow(In);
%imshow(mat2gray(In, [-0.5 0.5]));
%colormap gray;
%imagesc(In, [-0.5 0.5]);
%imwrite(uint8((mat2lab2rgb(In))*255), '../results/In.png');
imwrite(In, '../results/In.png');
entropyIn = entropy(In);
if is_octave || ~verLessThan('matlab', '8.3') %if Matlab R2014a and above
	PSNR_In = psnr(In, I);
else
	PSNR_In = -1;
end

figure;
%imshow(mat2lab2rgb(Iout));
imshow(Iout);
%imshow(mat2gray(Iout, [-0.5 0.5]));
%colormap gray;
%imagesc(Iout, [-0.5 0.5]);
%imwrite(uint8((mat2lab2rgb(Iout))*255), '../results/Iout.png');
imwrite(Iout, '../results/Iout.png');
entropyIout = entropy(Iout);
if is_octave || ~verLessThan('matlab', '8.3') %if Matlab R2014a and above
	PSNR_Iout = psnr(Iout, I);
else
	PSNR_Iout = -1;
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
fprintf('entropyIn = %g\n', entropyIn);
fprintf('PSNR_In = %g\n', PSNR_In);
fprintf('entropyIout = %g\n', entropyIout);
fprintf('PSNR_Iout = %g\n', PSNR_Iout);
fprintf('fresidue = %g\n', fresidue);
fprintf('fsparsity = %g\n', fsparsity);
fprintf('sparsity = %g\n', sparsity);
fprintf('avgnzero = %g\n', avgnzero);
fprintf('nzero = %g\n', nzero);
fprintf('nzeromin = %g\n', nzeromin);
fprintf('nzeromax = %g\n', nzeromax);
if(is_octave)
	save('../results/images.mat', 'I', 'In', 'Iout', 'Idiff', 'Sout', 'gamma', '-v7');
	save('../results/stats.mat', 'entropyI', 'entropyIn', 'entropyIout', 'sparsity', 'fresidue', 'fsparsity', 'PSNR_In', 'PSNR_Iout', 'nzero', 'nzeromin', 'nzeromax', '-v7');
else
	save('../results/images.mat', 'I', 'In', 'Iout', 'Idiff', 'Sout', 'gamma');
	save('../results/stats.mat', 'entropyI', 'entropyIn', 'entropyIout', 'sparsity', 'fresidue', 'fsparsity', 'PSNR_In', 'PSNR_Iout', 'nzero', 'nzeromin', 'nzeromax');
end
