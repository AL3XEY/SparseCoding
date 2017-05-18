function [ C2 ] = getC2(S2, nHL, HMAXparams)    
    C2 = zeros(nHL, 1);
    S2b = zeros(HMAXparams.nscales-1,nHL);
	for scal=1:HMAXparams.nscales-1
		S2b(scal,1:nHL) = max(max(S2{scal}(:,:,1:nHL))); %TODO or min ?
	end
	C2(1:nHL) = max(S2b(:,1:nHL));
end

