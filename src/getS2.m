function [ S2 ] = getS2( C1, HLfilters, HMAXparams, display )
    if nargin<3 || isempty(HMAXparams)
        HMAXparams = HMAXparameters();
    end
    if nargin<4 || isempty(display)
        display = false;
    end

    S2 = cell(1,HMAXparams.nscales-1);
    for scal=1:HMAXparams.nscales-1
		%nHL is the number of prototypes
		nHL = size(HLfilters{scal},4);
        sz = HMAXparams.filter_sz(scal);
        for cpt=1:nHL
			%calculate the response of each prototype over each patch of the C1 layer (see Mutch & Lowe 2008)
			% R = exp(||X - P||^2 / 2*sigma^2*alpha)
			sigma2 = 1;
			alpha = (sz/4)^2;
            [dx2,dy2,~]=size(C1{scal});
            for y=1:dy2-sz
                for x=1:dx2-sz
					X = HLfilters{scal}(:,:,:,cpt);
					P = C1{scal}(x:x+sz-1,y:y+sz-1,:);
					R = X - P;
					R = R.^2;
					R = sum(sum(sum(R)));
					R = sqrt(R);
					R = - R/(2*sigma2*alpha);
                    R = exp(R);
					S2{scal}(x,y,cpt) = R;
					%S2{scal}(x,y,cpt) = exp(-sqrt(sum(sum(sum((HLfilters{scal}(:,:,:,cpt) - C1{scal}(x:x+sz-1,y:y+sz-1,:)).^2))))/(2*sigma2*alpha));
                end
            end
        end
    end
    %Display the S2 layer
    if display
        for scal=1:HMAXparams.nscales-1
            for cpt=1:nHL
                figure
                imgS2 = S2{scal}(:,:,cpt);
                vis = max(max(max(imgS2(:,:,:))));
                imshow(uint8(255*(imgS2/vis + 0.3)))
            end
        end
    end
end
