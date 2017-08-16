function [HLfilters] = getHLfiltersFusion(imgs, nHL, sparseLearning, color ,HMAXparams, gaborFilters, display)
    if nargin<2 || isempty(nHL)
        nHL = 1000; %the number of prototypes to take (total) across all images
    end
    if nargin<5 || isempty(HMAXparams)
        HMAXparams = HMAXparameters();
    end
    if nargin<6 || isempty(gaborFilters)
        gaborFilters = getGaborFilters(HMAXparams, false);
    end
    if nargin<7 || isempty(display)
        display = false;
    end

    if is_octave
		pkg load image;
    end
    
    if sparseLearning
        addpath('../fast-sc-2/code/');
        addpath('../fast-sc-2/code/sc2/');
        addpath('../fast-sc-2/code/sc2/nrf/');
        nimg = size(imgs,2);
        HLfilters = cell(HMAXparams.nscales-1,1);
        S =HLfilters;
        stat = HLfilters;
        if color
            nchans = 2*HMAXparams.outChansHalf;
            C1 = cell(nimg*(HMAXparams.nscales-1)*nchans);
            nimg*(HMAXparams.nscales-1)*nchans
        else
            C1 = cell(nimg*(HMAXparams.nscales-1));
        end
    else
        %Adapt the number of images to the number of HL filters wanted
        %TODO : this was easy to code but not optimal because if nHL >> nimg,
        %we store a big number of images (we could store nimg images and access
        %them when multiple times instead)
        %+ the random selection might have non-unique values (it's a very small
        %risk, but still, we generate new x and y every new image)
        %tl;dr : this solution is fast but needs a lot of memory
        nimg = size(imgs,2);
        images = cell(1,nHL);
        if nHL < nimg
            perms = randperm(nimg,nHL);
            for i=1:size(perms,2)
                images{i}=imgs{perms(i)};
            end
            imgs = images;
        end
        if nHL > nimg
            iterations = round(nHL/nimg);
            additions = mod(nHL,nimg);
            for i=1:nimg*iterations
                images{i} = imgs{mod(i-1,nimg)+1};
            end
            if additions
                perms = randperm(nimg,additions);
                for i=1:additions
                    images{nimg*iterations+i}=imgs{perms(i)};
                end
            end
            imgs=images;
        end
        nimg = size(imgs,2);
        clear images 
    end
    HLfilters = cell(1,HMAXparams.nscales-1);
	for scal=1:HMAXparams.nscales-1
		sz = HMAXparams.filter_sz(scal);
		HLfilters{scal} = zeros(sz, sz, HMAXparams.nth, nimg);
    end
    
    for imgcpt=1:nimg
		img = imgs{imgcpt};
		if size(img,3)==3
			img = double(rgb2gray(img));%/255; % convert it to grayscale
		end
		[dx,dy,c] = size(img);
        if display
            figure
            imshow(uint8(255*img)) % show original image
        end

        if sparseLearning
            if color
                %OpponencyCells
                SO = getSODescriptor(img, HMAXparams, gaborFilters, display);
                DO = getDODescriptor(SO, HMAXparams, gaborFilters, display);
                %TODO do we sum the complementary channels ?
                nchans = size(DO{1},4);
                for chan=1:nchans
                    S1 = cell(1,HMAXparams.nscales);
                    for scal=1:HMAXparams.nscales
                       S1{scal}(:,:,:) = DO{scal}(:,:,:,chan);
                    end
                    C1tmp = getC1(S1, dx, dy, HMAXparams, display);
                    for scal=1:HMAXparams.nscales-1
                        C1{(imgcpt-1)*(HMAXparams.nscales-1+nchans)+(chan-1)*(HMAXparams.nscales-1)+scal}(:,:,:) = C1tmp{scal}(:,:,:);
                    end
                end
            else
                S1 = getS1(img, HMAXparams, gaborFilters, display);

                C1tmp = getC1(S1, dx, dy, HMAXparams, display);
                for scal=1:HMAXparams.nscales-1
                    C1{(imgcpt-1)*(HMAXparams.nscales-1)+scal}(:,:,:) = C1tmp{scal}(:,:,:);
                end
            end
        else
            S1 = getS1(img, HMAXparams, gaborFilters, display);

            C1 = getC1(S1, dx, dy, HMAXparams, display);

            %Take random patches from C1
            rdx = 1+floor(rand().*(dx)); %TODO unicity ?
            rdy = 1+floor(rand().*(dy));
            for scal=1:HMAXparams.nscales-1
                sz = HMAXparams.filter_sz(scal);
                [dx2,dy2,~]=size(C1{scal});
                maxX = dx2-sz; %TODO if maxX<0, stop
                maxY = dy2-sz; %TODO if maxY<0, stop
                ratio = dx2/dx;
                x = 1+round(rdx*ratio);
                y = 1+round(rdy*ratio);
                if x>maxX
                    x = maxX;
                end
                if y>maxY
                    y = maxY;
                end
                HLfilters{scal}(1:sz,1:sz,:,imgcpt) = C1{scal}(x:x+sz-1,y:y+sz-1,:);
            end            
        end
    end
    
    if sparseLearning
        %%%%% Apply Sparse Coding algorithm %%%%%
        num_bases = nHL;
        %beta = 0.1; %0.01
        %TODO if beta empty, default value
        batch_size = 1000;
        num_iters = iterations;%2;
        %TODO if iterations empty, default value
        fname_save = [];
        Binit = [];
        sparsity_func= 'L1';
        epsilon = [];
        pars.display_images = true;%display;%true;%false;
        pars.display_every = 1;%int32(display);%1;%0;
        pars.save_every = 1;%int32(display);%1;%0;
        pars.save_basis_timestamps = false;%display;%false;%true;
        %winsize=8;
        if oneDictPerScale
            %FIXME
%             for scal=1:HMAXparams.nscales-1 %FIXME
%                 %FIXME
%                 X = getdata_imagearray_all(C1{scal}, HMAXparams.filter_sz(scal)); %TODO or winsize ???)
%                 %C1 normalization
%                 mn = min(min(min(X)));
%                 mx = max(max(max(X)));
%                 X = ((X-mn)/(mx-mn))-0.5;
%                 [dicts{scal},S{scal},stat{scal}] = sparse_coding(X, num_bases, beta, sparsity_func, epsilon, num_iters, batch_size, fname_save, pars, Binit);
%             end
        else
            cpt=0;
            for scalnimg=1:nimg*(HMAXparams.nscales-1)
                fsz = HMAXparams.filter_sz(1);
                [hh,ww,oo]=size(C1{scalnimg});
                hh = hh-fsz+1;
                ww = ww-fsz+1;
                cpt=cpt+ww*hh*oo-1;
            end
            X = zeros(fsz^2,cpt);
            oldValue=1;
            for scalnimg=1:nimg*(HMAXparams.nscales-1)
                [hh,ww,oo]=size(C1{scalnimg});
                hh = hh-fsz+1;
                ww = ww-fsz+1;
                newValue = oldValue+ww*hh*oo-1;
                X(:,oldValue:newValue) = getdata_imagearray_all(C1{scalnimg}, fsz);
                oldValue=newValue;
            end
            %C1 normalization
            mn = min(min(min(X)));
            mx = max(max(max(X)));
            X = ((X-mn)/(mx-mn))-0.5;
            [dicts,S,stat] = sparse_coding(X, num_bases, beta, sparsity_func, epsilon, num_iters, batch_size, fname_save, pars, Binit);
        end
    end
end
