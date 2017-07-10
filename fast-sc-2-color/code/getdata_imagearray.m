function X = getdata_imagearray(IMAGES, winsize, num_patches)

[h,w,c,num_images]=size(IMAGES);
sz= winsize;
BUFF=4;

totalsamples = 0;
% extract subimages at random from this image to make data vector X
% Step through the images
X= zeros((sz^2)*3, num_patches);
for i=1:num_images,

    % Display progress
    fprintf('[%d/%d]',i,num_images);

    this_image=IMAGES(:,:,:,i);

    % Determine how many patches to take
    getsample = floor(num_patches/num_images);
    if i==num_images, getsample = num_patches-totalsamples; end

    % Extract patches at random from this image to make data vector X
    for j=1:getsample
        r=BUFF+ceil((h-sz-2*BUFF)*rand);
        c=BUFF+ceil((w-sz-2*BUFF)*rand);
        totalsamples = totalsamples + 1;
        % X(:,totalsamples)=reshape(this_image(r:r+sz-1,c:c+sz-1),sz^2,1);
        temp =reshape(this_image(r:r+sz-1,c:c+sz-1,1),sz^2,1);
        X(1:sz^2,totalsamples) = temp - mean(temp);
        temp =reshape(this_image(r:r+sz-1,c:c+sz-1,2),sz^2,1);
        X((sz^2)+1:(sz^2)*2,totalsamples) = temp - mean(temp);
        temp =reshape(this_image(r:r+sz-1,c:c+sz-1,3),sz^2,1);
        X((sz^2)*2+1:(sz^2)*3,totalsamples) = temp - mean(temp);
    end
end
fprintf('\n');
