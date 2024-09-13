function img0 = preprocess_imdilate(img)
%     min_val = min(img(:));
%     max_val = max(img(:));
%     img_ = (img - min_val)./(max_val - min_val);
len =size(img,3);
img0=img;
for i=1:len
    img0(:,:,i) = double(img(:,:,i));

    %% Image denoise
    se=strel('disk',4);%se=strel('ball',2,2,2);
    img0(:,:,i)=imdilate(  img0(:,:,i), se);

end

