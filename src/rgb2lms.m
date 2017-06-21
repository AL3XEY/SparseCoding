function [lms] = rgb2lms(rgb)
	[h,w,c] = size(rgb);
	rgb = double(rgb);
	R = rgb(:,:,1);
	G = rgb(:,:,2);
 	B = rgb(:,:,3);
	L = 17.8824*R + 43.5161*G + 4.1193*B;
	M = 3.4557*R + 27.1554*G + 3.8671*B;
	S = 0.02996*R + 0.18431*G + 1.467*B;
	lms = zeros(h,w,c);
	lms(:,:,1) = L;
	lms(:,:,2) = M;
	lms(:,:,3) = S;
end
