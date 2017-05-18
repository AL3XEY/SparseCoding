function [out] = gabor(x,y,thet,scale) %
	[sigma lambda gam nth] = HMAXparameters();
	lamb = lambda(scale); % we chose the scale and the filter uses the parameters defined above
	sig = sigma(scale);
	[nx,ny] = size(x);
	for i=1:nth
        x0(:,:,i) = x .* cos(thet(i)) + y .* sin(thet(i));
		y0(:,:,i) = y .* cos(thet(i)) - x .* sin(thet(i));
        %if(thet(i)==pi/2) % FIXME only way I found to counter the fact that cos(pi/2) = cos(1.5708) = 6.123e-17 != 0
        %    x0(:,:,i) = y;
        %    y0(:,:,i) = -x;
        %end
		%%x0(:,:,i) = x .* cosd(thet(i)) + y .* sind(thet(i));
		%%y0(:,:,i) = y .* cosd(thet(i)) - x .* sind(thet(i));
	end
	%x0 = reshape(permute(x,[3,1,2])' * cos(thet) + permute(y,[3,1,2])' *
	%sin(thet),nx,ny,nth); %Octave only
	%y0 = reshape(permute(y,[3,1,2])' * cos(thet) - permute(x,[3,1,2])' *
	%sin(thet),nx,ny,nth); %Octave only

	%gabA=2*pi*x0/lamb;
	%gabB=cos(gabA);
	%gabC=x0.^2;
	%gabD=gam*y0.^2;
	%gabE=(gabC+gabD)/sig^2;
	%gabF=exp(-0.5*gabE);
	%out = gabF.*gabB;
	out = exp( -0.5 * (x0.^2 + gam * y0.^2) / sig^2) .* cos(2 * pi * x0 / lamb);
end
