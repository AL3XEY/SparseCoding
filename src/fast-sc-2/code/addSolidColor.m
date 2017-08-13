function [B,pars] = addSolidColor(B,pars)
    [h,w]=size(B);
    Btmp = zeros(h,w+pars.channels);
    Btmp(1:h,1:w) = B;
    for chan=1:pars.channels
        Btmp(:,w+chan) = 0.25.*ones(h,1);
    end
    B = Btmp;
    pars.num_bases = size(B,2);
end
