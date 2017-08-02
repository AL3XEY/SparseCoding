function X = getdata_imagearray_all(IMAGES, winsize)

[h,w,channels,num_images]=size(IMAGES);
foo = h - winsize + 1;
bar = w - winsize + 1;
patches_per_image = foo*bar;
num_patches = num_images * patches_per_image;
cpt = 1;
% extract subimages to make data vector X
% Step through the images
X= zeros((winsize^2)*channels, num_patches);
for i=1:num_images,

    % Display progress
    fprintf('[%d/%d]',i,num_images);

    this_image=IMAGES(:,:,:,i);

	for j=1:foo
		for k=1:bar
			for chan=1:channels
				X((chan-1)*winsize^2+1:chan*winsize^2,cpt) = reshape(this_image(j:j+winsize-1, k:k+winsize-1, chan),winsize^2,1)';
			end
	%		X(1:winsize^2,cpt) = reshape(this_image(j:j+winsize-1, k:k+winsize-1, 1),winsize^2,1)';
    %        X((winsize^2)+1:(winsize^2)*2,cpt) = reshape(this_image(j:j+winsize-1, k:k+winsize-1, 2),winsize^2,1)';
    %        X((winsize^2)*2+1:(winsize^2)*3,cpt) = reshape(this_image(j:j+winsize-1, k:k+winsize-1, 3),winsize^2,1)';
			cpt = cpt+1;
		end
	end
end
fprintf('\n');
