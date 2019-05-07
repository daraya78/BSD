function re=logsubexpmat(A,B)
    nx=size(A,2);
    ny=size(A,1);
    Ar=reshape(A,[nx*ny 1]);
    Br=reshape(B,[nx*ny 1]);
    Cr=util.logsubexp([Ar Br],2);
    re=reshape(Cr,[ny nx]);
end
    
    
