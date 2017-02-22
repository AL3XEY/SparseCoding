function [ pyr ] = pyramid( img, iterations=3, type='Gaussian', debug=false )
%pyramid Computes the Gaussian or Laplacian pyramid of an image
%   

pyr={};
gaussian = {};
gaussian{1} = img;

switch type
case 'Gaussian'
	for i=1:iterations
		foo = imsmooth(gaussian{i}, 'Gaussian');
		gaussian{i+1} = imresize(foo, 0.5);
	end
	pyr = gaussian;
case 'Laplacian'
	laplacian = {};
	laplacian{1} = img;
	for i=1:iterations
		foo = imsmooth(gaussian{i}, 'Gaussian');
		gaussian{i+1} = imresize(foo, 0.5);
		laplacian{i+1} = double(foo) - double(gaussian{i});
	end
	pyr = laplacian;
otherwise
	disp('not valid');
end

if debug
	for i=1:iterations+1
		figure;
		imshow(mat2gray(pyr{i}));
	end
end
