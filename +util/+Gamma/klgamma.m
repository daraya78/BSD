% KLGAMMA   Kullback-Leibler divergence for Gamma distributions.
%           Calculates KL(P||Q) where P and Q are Gamma distributions with
%           parameters {pa,pb} and {qa,qb}.
%
%           KL(P||Q) = \int d\pi P(\pi) ln { P(\pi) / Q(\pi) }.
%
% Usage:
%   >>  kl = KLGAMMA(pa,pb,qa,qb);
%
% Inputs:
%   pa    - [numeric] parameter a for distribution P.
%   pb    - [numeric] parameter b for distribution P.
%   qa    - [numeric] parameter a for distribution Q.
%   qb    - [numeric] parameter b for distribution Q.
%
% Additional notes:
%   This routine handles factorised P distributions, if their parameters
%   are specified multiply in either 'pa' or 'pb', as elements of a row
%   vector.
%
% Copyright (C) 2001 Matthew J. Beal
function [kl] = klgamma(pa,pb,qa,qb)
    n = max([size(pb,2) size(pa,2)]);

    if size(pa,2) == 1, pa = pa*ones(1,n); end
    if size(pb,2) == 1, pb = pb*ones(1,n); end
    qa = qa*ones(1,n); qb = qb*ones(1,n);

    kl = sum( pa.*log(pb)-gammaln(pa) ...
             -qa.*log(qb)+gammaln(qa) ...
         +(pa-qa).*(psi(pa)-log(pb)) ...
         -(pb-qb).*pa./pb );
end