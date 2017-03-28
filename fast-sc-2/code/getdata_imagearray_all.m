function X = getdata_imagearray_all(IMAGES, winsize)

[h w num_images]=size(IMAGES);
BUFF=4;
foo = h - winsize + 1;
bar = w - winsize + 1;
patches_per_image = foo*bar;
num_patches = num_images * patches_per_image;
cpt = 1;
% extract subimages to make data vector X
% Step through the images
X= zeros(winsize^2, num_patches);
for i=1:num_images,

    % Display progress
    fprintf('[%d/%d]',i,num_images);

    this_image=IMAGES(:,:,i);

	for j=1:foo
		for k=1:bar
			X(:,cpt) = reshape(this_image(j:j+winsize-1, k:k+winsize-1),winsize^2,1)';
			cpt = cpt+1;
		end
	end
end
fprintf('\n');
