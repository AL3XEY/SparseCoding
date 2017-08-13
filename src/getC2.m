function [ C2,coords ] = getC2(S2, HMAXparams)
    if nargin<2 || isempty(HMAXparams)
        HMAXparams = HMAXparameters();
    end
    nHL = size(S2{1},3);
    C2tmp = zeros(nHL, 4);
    S2b = zeros(HMAXparams.nscales-1,nHL,3);
    for scal=1:HMAXparams.nscales-1
        for hl=1:nHL
            m=S2{scal}(:,:,hl);
            [value,idx] = max(m(:)); %TODO or min ?
            [maxX,maxY] = ind2sub(size(m),idx);
            %m(I1,I2)
            %Caution : max is always 1, I guess, if one of the HL filters is taken from current image
            S2b(scal,hl,:) = [value,maxX,maxY];
        end
    end
    for hl=1:nHL
        [val,scal] = max(S2b(:,hl,1),[],1);
        C2tmp(hl,1)=val;
        C2tmp(hl,2:3) = S2b(scal,hl,2:3);
        C2tmp(hl,4)=scal;
    end
    C2 = C2tmp(:,1)';
    coords = C2tmp(:,2:4);
end
