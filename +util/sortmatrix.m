function sa = sortmatrix(mask,A)
for fil=1:size(mask,1)
    for col=1:size(mask,1)
        sa(mask(fil),mask(col))=A(fil,col);
    end
end


end

