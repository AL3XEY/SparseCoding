function X = getdata_imagearray(IMAGES, winsize, type, option)
    [h,w,channels,num_images]=size(IMAGES);
    BUFF=4;
    if strcmp(type,'rand')
        num_patches = option;
    else
        if strcmp(type,'all')
            foo = h - winsize + 1;
            bar = w - winsize + 1;
        end
        if strcmp(type,'some')
            foo = floor(h / winsize);
            bar = floor(w / winsize);
        end
        patches_per_image = foo*bar;
        num_patches = num_images * patches_per_image;
    end
    cpt=1;
    % extract subimages at random from this image to make data vector X
    % Step through the images
    X= zeros((winsize^2)*channels, num_patches);
    for i=1:num_images,
        % Display progress
        fprintf('[%d/%d]',i,num_images);
        this_image=IMAGES(:,:,:,i);
        if strcmp(type,'rand')
           % Determine how many patches to take
            getsample = floor(num_patches/num_images);
            if i==num_images
                getsample = num_patches-cpt;
            end
            % Extract patches at random from this image to make data vector X
            for j=1:getsample
                r=BUFF+ceil((h-winsize-2*BUFF)*rand);
                c=BUFF+ceil((w-winsize-2*BUFF)*rand);
                cpt = cpt + 1;
                for chan=1:channels
                    temp = reshape(this_image(r:r+winsize-1,c:c+winsize-1,chan),winsize^2,1);
                    X((chan-1)*winsize^2+1:chan*winsize^2,cpt) = temp - mean(temp);
                end
            end 
        else
            for j=1:foo
                for k=1:bar
                    for chan=1:channels
                        X((chan-1)*winsize^2+1:chan*winsize^2,cpt) = reshape(this_image(j:j+winsize-1, k:k+winsize-1, chan),winsize^2,1)';
                    end
                    cpt = cpt+1;
                end
            end
        end
    end
    fprintf('\n');
end