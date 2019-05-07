function D=KLDirich(alpha1,alpha2)
    D = gammaln(sum(alpha1))-gammaln(sum(alpha2))-sum(gammaln(alpha1))+ ...
sum(gammaln(alpha2))+(alpha1-alpha2)*(psi(alpha1)-psi(sum(alpha1)))'; 
end

