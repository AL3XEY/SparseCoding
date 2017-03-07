function X = getdata_imagearray_all(IMAGES, winsize)

num_images=size(IMAGES,3);
image_size=size(IMAGES,1);
BUFF=4;
foo = image_size - winsize + 1
patches_per_image = foo^2
num_patches = num_images * patches_per_image
cpt = 1;
% extract subimages to make data vector X
% Step through the images
X= zeros(winsize^2, num_patches);
for i=1:num_images,

    % Display progress
    fprintf('[%d/%d]',i,num_images);

    this_image=IMAGES(:,:,i);

	for j=1:foo
		for k=1:foo
			X(:,cpt) = reshape(this_image(j:j+winsize-1, k:k+winsize-1),winsize^2,1)';
			cpt = cpt+1;
		end
	end
end
fprintf('\n');
