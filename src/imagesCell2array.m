[arr] = function imagesCell2array(cel)
%take a cell and transform it to a matrix (for getdata_imagearray)
    s = size(cel,2);
    [h,w,c] = size(cel{1});
    arr = zeros(h,w,c,s);
    for i=1:s
        arr(:,:,:,i) = cel{i}(:,:,:);
    end
end
