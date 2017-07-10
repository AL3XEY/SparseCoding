close all
clear
clc
tic

W=ones(3,4);
W(2,1)=-1;
W(3,1)=0;
W(1,2)=2;
W(2,2)=-1;
W(3,2)=-1;
W(3,3)=-2;
W(:,1)=W(:,1)/sqrt(2);
W(:,2)=W(:,2)/sqrt(6);
W(:,3)=W(:,3)/sqrt(6);
W(:,4)=W(:,4)/sqrt(3);
WR = W(1,:);
WG = W(2,:);
WB = W(3,:);

%function [res] = f(x,y,lambda)

%end

nchans = 8;

%img = imread('../res/lena.ppm');
img = imread('../res/rgb-circles.bmp');
img = double(img)./255;
%
imgrgb = img;
%imglms = rgb2lms(img);
%img = imglms;
[h,w,c] = size(img);

%build gabor filters
HMAXparams = HMAXparameters();
nth=HMAXparams.nth;
nscales = HMAXparams.nscales;
gaborFilters = getGaborFilters(HMAXparams, false);
ochans=cell(1,nscales);
dochans=cell(1,nscales);
excit=cell(1,nscales);
inhib=cell(1,nscales);
for scal=1:nscales
	ochans{scal} = zeros(h,w,nth,8);
    dochans{scal} = zeros(h,w,nth,8);
	excit{scal} = zeros(h,w,nth,3);
	inhib{scal} = zeros(h,w,nth,3);
end

for scal = 1:nscales
    for th = 1:nth
        filtr = gaborFilters{scal}(:,:,th);
        filterexcit = (gaborFilters{scal}(:,:,th)>0).*gaborFilters{scal}(:,:,th);
        filterinhib = -(gaborFilters{scal}(:,:,th)<0).*gaborFilters{scal}(:,:,th);
        for color=1:3
            excit{scal}(:,:,th,color) = abs(filter2(filterexcit,img(:,:,color)));
            inhib{scal}(:,:,th,color) = abs(filter2(filterinhib,img(:,:,color)));
        end
    end
end

display=false;
for scal=1:nscales
	for th=1:nth
		ochans{scal}(:,:,th,1) = WR(1).*excit{scal}(:,:,th,1) + WG(1).*inhib{scal}(:,:,th,2);% + WB(1).*img(:,:,3);
		ochans{scal}(:,:,th,2) = WR(2).*excit{scal}(:,:,th,1) + WG(2).*inhib{scal}(:,:,th,2) + WB(2).*inhib{scal}(:,:,th,3);
		ochans{scal}(:,:,th,3) = WR(3).*excit{scal}(:,:,th,1) + WG(3).*excit{scal}(:,:,th,2) + WB(3).*inhib{scal}(:,:,th,3);
		ochans{scal}(:,:,th,4) = WR(4).*excit{scal}(:,:,th,1) + WG(4).*excit{scal}(:,:,th,2) + WB(4).*excit{scal}(:,:,th,3);
		%ochans{scal}(:,:,th,5:8) = -ochans{scal}(:,:,th,1:4);
		ochans{scal}(:,:,th,5) = -(WR(1).*inhib{scal}(:,:,th,1) + WG(1).*excit{scal}(:,:,th,2));% + WB(1).*img(:,:,3);
		ochans{scal}(:,:,th,6) = -(WR(2).*inhib{scal}(:,:,th,1) + WG(2).*excit{scal}(:,:,th,2) + WB(2).*excit{scal}(:,:,th,3));
		ochans{scal}(:,:,th,7) = -(WR(3).*inhib{scal}(:,:,th,1) + WG(3).*inhib{scal}(:,:,th,2) + WB(3).*excit{scal}(:,:,th,3));
		ochans{scal}(:,:,th,8) = -(WR(4).*inhib{scal}(:,:,th,1) + WG(4).*inhib{scal}(:,:,th,2) + WB(4).*inhib{scal}(:,:,th,3));
		if display
			figure
			colormap gray
			imagesc(ochans{scal}(:,:,th,2))

			figure
			colormap gray
			imagesc(ochans{scal}(:,:,th,6))
			%figure
			%colormap gray
			%imagesc(ochans{scal}(:,:,th,2))
			%figure
			%colormap gray
			%imagesc(ochans{scal}(:,:,th,3))
			%figure
			%colormap gray
			%imagesc(ochans{scal}(:,:,th,4))
		end
	end
end

%TODO display excit, inhib and then ochans/opponnencychannels

%half-squaring
for scal=1:nscales
	idx = find(ochans{scal}<0);
	ochans{scal}(idx) = 0;
	ochans{scal} = ochans{scal}.^2;
end



display = false;
if display
	for scal=1:1%nscales
		for th=1:1%nth
            for chan=1:8
                figure
                colormap gray
                imagesc(-ochans{scal}(:,:,th,chan))

                %figure
                %imagesc(ochans{scal}(:,:,th,1))

                %figure
                %imagesc(ochans{scal}(:,:,th,5))
            end
		end
	end
end

%divisive normalization
normchans=cell(1,nscales);
for scal=1:nscales
	normchans{scal}=zeros(h,w,nth,8);
end
k=1;
sigma=0.225;
sigma2=sigma^2;
%TODO see which loop order is the best
sm = sum(ochans{scal},4);
for scal=1:nscales
	for chan=1:8
		%ch = [1:8];
		%ch = ch(ch~=chan);
		%sm = sum(ochans{scal}(:,:,:,ch),4);
		normchans{scal}(:,:,:,chan) = sqrt((k*ochans{scal}(:,:,:,chan))./(sigma2+sm));
	end
end

display = false;
if display
	for scal=1:1%nscales
		for th=1:1%nth
			for chan=1:8
				figure
                colormap gray
				imagesc(-normchans{scal}(:,:,th,chan))

				%figure
				%colormap gray
				%imagesc((-normchans{scal}(:,:,th,1))

				%figure
				%colormap gray
				%imagesc(-normchans{scal}(:,:,th,5))
			end
		end
	end
end

%DO cells

%gabor filters on each channel

for scal = 1:nscales
    for th = 1:nth
        filtr = gaborFilters{scal}(:,:,th);
        for chan=1:8
            dochans{scal}(:,:,th,chan) = abs(filter2(filtr,normchans{scal}(:,:,th,chan))); % filtered images
		end
	end
end

display = false;
if display
	for scal=1:1%nscales
		%for th=1:1%nth
			for chan=1:8
				figure
                colormap gray
				imagesc(dochans{scal}(:,:,th,chan))
			end
		%end
	end
end

%half-squaring
for scal=1:nscales
	idx = find(dochans{scal}<0);
	dochans{scal}(idx) = 0;
	dochans{scal} = dochans{scal}.^2;
end

display = false;
if display
	for scal=1:1%nscales
		for th=1:1%nth
			for chan=1:8
				figure
                colormap gray
				imagesc(-dochans{scal}(:,:,th,chan))
			end
		end
	end
end

%normalization
donormchans=cell(1,nscales);
dosumchans=cell(1,nscales);
for scal=1:nscales
	donormchans{scal}=zeros(h,w,nth,8);
end
for scal=1:nscales
	dosumchans{scal}=zeros(h,w,nth,4);
end
k=1;
sigma=0.225;
sigma2=sigma^2;
sm = sum(dochans{scal},3);
for scal=1:nscales
	for th=1:nth
		donormchans{scal}(:,:,th,:) = sqrt((k*dochans{scal}(:,:,th,:))./(sigma2+sm));
	end
end

display = false;
if display
	for scal=1:1%nscales
		for th=1:1%nth
			for chan=1:8
				figure
                colormap gray
				imagesc(-donormchans{scal}(:,:,th,chan))
			end
		end
	end
end
display = false;
if display
	for scal=1:1%nscales
		for th=1:1%nth
			for chan=1:4
				figure
                colormap gray
				imagesc(-(donormchans{scal}(:,:,th,chan) + donormchans{scal}(:,:,th,chan+4)))
			end
		end
	end
end

%max pooling
pooled = cell(1,nscales);
for scal=1:nscales
	pooled{scal} = zeros(h,w,chan);
	pooled{scal} = max(donormchans{scal}, [], 3);
end

display=false;
if display
	for scal=1:1%nscales
		for chan=1:nchans
			figure
            colormap gray
			imagesc(-pooled{scal}(:,:,chan))
		end
	end
end
display=true;
if display
	for scal=1:1%nscales
		for chan=1:4
			figure
            colormap gray
			DOcells{scal}(:,:,chan) = pooled{scal}(:,:,chan) + pooled{scal}(:,:,chan+4);
			imagesc(-DOcells{scal}(:,:,chan))
		end
	end
end

toc
