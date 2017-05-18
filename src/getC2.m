function [ C2 ] = getC2( S2, nHL )
    [sigma,lambda,gam,nth,nscales] = HMAXparameters();
    
    C2 = zeros(nHL, 1);
    S2b = zeros(nscales-1,nHL);
	for scal=1:nscales-1
		S2b(scal,1:nHL) = max(max(S2{scal}(:,:,1:nHL))); %TODO or min ?
	end
	C2(1:nHL) = max(S2b(:,1:nHL));
end

