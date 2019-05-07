function D=KLGamma2(shape1,scale1,shape2,scale2)
    D=(shape1-shape2)*digamma(shape1)-gammaln(shape1)+gammaln(shape2)+shape2*(log(scale2)-log(scale1))...
        +shape1*(scale1-scale2)/scale2; 
end

