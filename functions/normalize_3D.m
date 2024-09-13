function Y =normalize_3D(X)
[n1 n2 n3]=size(X);
for i=1:n3
  b=X(:,:,i);
  b=b(:);
  b=normalize(b,'range');
  Y(:,:,i)=reshape(b,[n1 n2]);
end
end