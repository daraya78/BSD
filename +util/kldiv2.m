function KL = kldiv2(pVect1,pVect2,varargin)
%KLDIV Kullback-Leibler or Jensen-Shannon divergence between two distributions.
%   KLDIV(P1,P2) returns the Kullback-Leibler divergence between two
%   distributions specified over M variable values.  P1 is a
%   length-M vector of probabilities representing distribution 1, and P2 is a
%   length-M vector of probabilities representing distribution 2.
%   The Kullback-Leibler divergence is given by:
%
%       KL(P1(x),P2(x)) = sum[P1(x).log(P1(x)/P2(x))]
%
%   The elements of probability vectors P1 and P2 must 
%   each sum to 1 +/- .00001. Otherwise, the program will re-normalise both
%   P1 and P2
%
%   A "log of zero" warning will be thrown for zero-valued probabilities.
%   Handle this however you wish.  Adding 'eps' or some other small value 
%   to all probabilities seems reasonable.  (Renormalize if necessary.)
%
%   KLDIV(X,P1,P2,'sym') returns a symmetric variant of the Kullback-Leibler
%   divergence, given by [KL(P1,P2)+KL(P2,P1)]/2.  See Johnson and Sinanovic
%   (2001). Zero probability values can be handled as above.
%
%   KLDIV(X,P1,P2,'js') returns the Jensen-Shannon divergence, given by
%   [KL(P1,Q)+KL(P2,Q)]/2, where Q = (P1+P2)/2.  See the Wikipedia article
%   for "Kullback–Leibler divergence".  This is equal to 1/2 the so-called
%   "Jeffrey divergence."  See Rubner et al. (2000). In this case when both
%   P1 and P2 are zero it will produce NaNs. The program handles the
%   problem by dissregarding those values
%
%   EXAMPLE:  Let the event set and probability sets be as follow:
%                P1 = ones(5,1)/5;
%                P2 = [0 0 .5 .2 .3]' + eps;
%  
%                 KL = kldiv(X,P1,P2); 
%                 KL =  
%                      19.4899
%
%             Note that we avoided the log-of-zero warning by adding 'eps'
%             to all probability values in P2.  We didn't need to
%             renormalize because we're still within the sum-to-one
%             tolerance.
%  
%   REFERENCES:
%   1) Cover, T.M. and J.A. Thomas. "Elements of Information Theory," Wiley, 
%      1991.
%   2) Johnson, D.H. and S. Sinanovic. "Symmetrizing the Kullback-Leibler 
%      distance." IEEE Transactions on Information Theory (Submitted).
%   3) Rubner, Y., Tomasi, C., and Guibas, L. J., 2000. "The Earth Mover's 
%      distance as a metric for image retrieval." International Journal of 
%      Computer Vision, 40(2): 99-121.
%   4) <a href="matlab:web('http://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence','-browser')">Kullback–Leibler divergence</a>. Wikipedia, The Free Encyclopedia.
%
%   See also: MUTUALINFO, ENTROPY

% Note:
%   Modified by Nelson to account for zero probability values in either P1
%   or P2 (for the J-S divergence only, for the other cases just sum eps to
%   the two histograms) and to ensure that probabilities sum to one, i.e.
%   to ensure that P1 and P2 are normalised.

%if ~isequal(unique(varValue),sort(varValue)),
%    warning('KLDIV:duplicates','X contains duplicate values. Treated as distinct values.')
%end
%if ~isequal(size(varValue),size(pVect1)) || ~isequal(size(varValue),size(pVect2)),
%    error('All inputs must have same dimension.')
%end
% Check probabilities sum to 1:
if (abs(sum(pVect1) - 1) > .00001) || (abs(sum(pVect2) - 1) > .00001),
    %error('Probablities don''t sum to 1.')
    pVect1 = pVect1./sum(pVect1);
    pVect2 = pVect2./sum(pVect2);
end

if ~isempty(varargin),
    switch varargin{1},
        case 'js',
            
            % logQvect = log2((pVect2+pVect1)/2);
            Qvect = (pVect2+pVect1)/2;
            
            temppVect1 = pVect1.*(log2(pVect1./Qvect));
            temppVect1(isnan(temppVect1)) = 0;
            
            temppVect2 = pVect2.*(log2(pVect2./Qvect));
            temppVect2(isnan(temppVect2)) = 0;
            
%             KL = .5 * (sum(pVect1.*(log2(pVect1)-logQvect)) + ...
%                 sum(pVect2.*(log2(pVect2)-logQvect)));

            KL = .5 * (sum(temppVect1) + ...
                sum(temppVect2));

        case 'sym',
            
            KL1 = sum(pVect1 .* (log2(pVect1)-log2(pVect2)));
            KL2 = sum(pVect2 .* (log2(pVect2)-log2(pVect1)));

            KL = (KL1+KL2)/2;
            
        otherwise
            error(['Last argument' ' "' varargin{1} '" ' 'not recognized.'])
    end
else
    
    KL = sum(pVect1 .* (log2(pVect1)-log2(pVect2)));

end








