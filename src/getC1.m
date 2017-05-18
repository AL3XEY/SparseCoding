function [ C1 ] = getC1( S1, dx, dy, display )
    if nargin<4 || isempty(display)
        display = false;
    end

    [sigma,lambda,gam,nth,nscales,filter_sz,pool_sz] = HMAXparameters();
    
    C1 = cell(1,nscales-1);
    for scal=1:nscales-1
		sz = pool_sz(scal);
		pxm = floor(dx/sz); pym = floor(dy/sz);

        for j = 1:nth % pooling for every orientation
            for px = 0:pxm-1
                for py = 0:pym-1
                    C1{scal}(px+1,py+1,j) = max(max([S1{scal}(px*sz+1:(px+1)*sz,py*sz+1:(py+1)*sz,j) S1{scal+1}(px*sz+1:(px+1)*sz,py*sz+1:(py+1)*sz,j)]));
                end
            end
        end

		% display filtered and pooled images
        if display
			figure
            for i = 1:nth
			  imaf = C1{scal}(:,:,i);
			  vis = max(imaf(:));
			  subplot(3,4,i)
			  imshow(uint8(255*(imaf/vis+0.3)))
            end
        end
    end
end

