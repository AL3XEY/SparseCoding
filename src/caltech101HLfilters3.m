%HLfilters=cell(1,3);
for i=1:92
    i
    tic
    load(sprintf('../res/101_ObjectCategories/caltech101gray-%d.mat',i));
    HLfilterstmp = getHLfilters(images,size(images,2),false);
    for j=1:size(HLfilterstmp{1},4)
        %HLfilters{((i-1)*100)+j} = HLfilterstmp{j};
        for scal=1:size(HLfilterstmp,2)
            HLfilters{scal}(:,:,:,((i-1)*100)+j) = HLfilterstmp{scal}(:,:,:,j);
        end
    end
    toc
end
save('../res/101_ObjectCategories/caltechHLfiltersgray.mat','HLfilters');
%save('../res/101_ObjectCategories/caltechHLfiltersgray1.mat','HLfilters');