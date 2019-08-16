function [out]=projection_cvxbody(in,threshold,model_PDE)

    
I=reshape(1e3./sqrt(in),model_PDE.n);
[m n]=size(I);
Ifill=I;
for i=1:m
    ind1=I(i,:)>threshold;
    indleft=min(find(ind1));
    indright=max(find(ind1));
    Ifill(i,indleft:indright)=threshold; 
end

out=max(Ifill,I);
out=1e6./(out(:).^2);
end