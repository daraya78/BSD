function samples = sample_discrete(vec,num)
    % vec is a column vector to be sampled like a discrete distribution
    % does not need to be normalized
    % vec is interpreted as p(1) = vec(1), p(2) = vec(2), etc.
    if ~exist('num','var'), num=1; end
    n=size(vec,2);
    if n==1
        vec=vec';
        n=size(vec,2);
    end
    cumvals = cumsum(vec);
    %sample = sum(rand() * cumvals(end) > cumvals) + 1;
    samples = sum(repmat(rand(num,1).*cumvals(end),1,n)>repmat(cumvals,num,1),2)+1;
end