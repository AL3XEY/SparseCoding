function [SO] = getSODescriptor(img, HMAXparams, display) %TODO ability to pass weight matrix as input

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

	nchans = 8;
	img = double(img)./255;
	%imglms = rgb2lms(img);
	%img = imglms;
	[h,w,c] = size(img);

	%build gabor filters
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
			if display %TODO
				figure
				colormap gray
				imagesc(ochans{scal}(:,:,th,2))

				figure
				colormap gray
				imagesc(ochans{scal}(:,:,th,6))
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

	if display
		for scal=1:1%nscales
			for th=1:1%nth
	            for chan=1:8
	                figure
	                colormap gray
	                imagesc(ochans{scal}(:,:,th,chan))
	            end
			end
		end
	end

	%divisive normalization
	SO=cell(1,nscales);
	for scal=1:nscales
		SO{scal}=zeros(h,w,nth,8);
	end
	k=1;
	sigma=0.225;
	sigma2=sigma^2;
	%TODO see which loop order is the best
	sm = sum(ochans{scal},4);
	for scal=1:nscales
		for chan=1:8
			SO{scal}(:,:,:,chan) = sqrt((k*ochans{scal}(:,:,:,chan))./(sigma2+sm));
		end
	end

	if display
		for scal=1:1%nscales
			for th=1:1%nth
				for chan=1:8
					figure
	                colormap gray
					imagesc(SO{scal}(:,:,th,chan))
				end
			end
		end
	end

end
