function s = logsubexp(a, dim)

if nargin < 2
  dim = 1;
end

% subtract the largest in each column
[y, i] = max(real(a),[],dim);
dims = ones(1,ndims(a));
dims(dim) = size(a,dim);
a = a - repmat(y, dims);
a=exp(a);
if dim==1
    s = y + log(a(1,:)-a(2,:));
else
    s = y + log(a(:,1)-a(:,2));
end
i = find(~isfinite(y));
if ~isempty(i)
  s(i) = y(i);
end


