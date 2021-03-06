function X = getdata_imagearray_some(IMAGES, winsize)

num_images=size(IMAGES,3);
image_size=size(IMAGES,1);
BUFF=4;
foo = image_size / winsize
patches_per_image = foo^2
num_patches = num_images * patches_per_image
cpt = 1;
% extract subimages at random from this image to make data vector X
% Step through the images
X= zeros(winsize^2, num_patches);
for i=1:num_images,

    % Display progress
    fprintf('[%d/%d]',i,num_images);

    this_image=IMAGES(:,:,i);

%    % Determine how many patches to take
%    getsample = floor(num_patches/num_images);
%    if i==num_images, getsample = num_patches-totalsamples; end
%
%    % Extract patches at random from this image to make data vector X
%    for j=1:getsample
%        r=BUFF+ceil((image_size-sz-2*BUFF)*rand);
%        c=BUFF+ceil((image_size-sz-2*BUFF)*rand);
%        totalsamples = totalsamples + 1;
%        % X(:,totalsamples)=reshape(this_image(r:r+sz-1,c:c+sz-1),sz^2,1);
%        temp =reshape(this_image(r:r+sz-1,c:c+sz-1),sz^2,1);
%        X(:,totalsamples) = temp - mean(temp);
%    end
	for j=1:foo
		for k=1:foo
			X(:,cpt) = reshape(this_image((j-1)*winsize+1:j*winsize, (k-1)*winsize+1:k*winsize),winsize^2,1)';
			cpt = cpt+1;
		end
	end
end
fprintf('\n');
