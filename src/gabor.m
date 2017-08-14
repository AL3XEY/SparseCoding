function [out] = gabor(x,y,thet,scale,HMAXparams)
	%compute Gabor filters for orientation and scale/size given
	lamb = HMAXparams.lambda(scale); % we chose the scale and the filter uses the parameters defined above
	sig = HMAXparams.sigma(scale);
	[nx,ny] = size(x);
    x0 = zeros(nx,ny,HMAXparams.nth);
    y0 = x0;
	for i=1:HMAXparams.nth
        x0(:,:,i) = x .* cos(thet(i)) + y .* sin(thet(i));
		y0(:,:,i) = y .* cos(thet(i)) - x .* sin(thet(i));
        %XXX cos(pi/2) = cos(1.5708) = 6.123e-17 != 0
	end
	%x0 = reshape(permute(x,[3,1,2])' * cos(thet) + permute(y,[3,1,2])' *
	%sin(thet),nx,ny,nth); %Octave only
	%y0 = reshape(permute(y,[3,1,2])' * cos(thet) - permute(x,[3,1,2])' *
	%sin(thet),nx,ny,nth); %Octave only

	out = exp( -0.5 * (x0.^2 + HMAXparams.gam * y0.^2) / sig^2) .* cos(2 * pi * x0 / lamb);
end
