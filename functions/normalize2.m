function input = normalize2(input)     
input = double(input);
[n1, n2, n3] = size(input);
for i =1:n3
    img = input(:,:,i);
    maxv = max(img(:));
    img_nozero = double(img>0);
    img = img.*(img_nozero);
    input(:,:,i) = img./(maxv).*255;
end
end