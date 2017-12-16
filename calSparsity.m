function [result]=calSparsity(abundance)
endmembNum=size(abundance,2);
imageSize=size(abundance,1);
total=0;
for l=1:endmembNum
    normDiv=norm(abundance(:,l),1)/norm(abundance(:,l),2);
    total=total+(sqrt(imageSize)-normDiv)/(sqrt(imageSize)-1);
end
result =total/sqrt(endmembNum);
end