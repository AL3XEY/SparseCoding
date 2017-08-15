function X = getdata_imagearray(IMAGES, winsize, type, option)
    if nargin<3 || isempty(type) || (~strcmp(type,'some') && ~strcmp(type,'rand'))
        %this is equivalent to previous 'all' type
        %selects every patch of the image(s)
        type = 'some';
        option = 1;
    end
    [h,w,channels,num_images]=size(IMAGES);
    if strcmp(type,'rand')
        BUFF=4;
        if nargin<4 || isempty(option)
            option = 4096;
        end
        num_patches = option;
    else
        if isempty(option)
            if nargin < 4
                option = winsize;
            else
                option = 1;
            end
        end
        if option==winsize
            foo = floor(h / option);
            bar = floor(w / option);
        end
        if option<winsize
            foo = floor((h - winsize + 1)/ option);
            bar = floor((w - winsize + 1)/ option);
        end
        if option>winsize
            foo = floor((h - winsize + 1)/ option)+1;
            bar = floor((w - winsize + 1)/ option)+1;
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
                        X((chan-1)*winsize^2+1:chan*winsize^2,cpt) = reshape(this_image((j-1)*option+1:(j-1)*option+winsize, (k-1)*option+1:(k-1)*option+winsize, chan),winsize^2,1)';
                    end
                    cpt = cpt+1;
                end
            end
        end
    end
    fprintf('\n');
end