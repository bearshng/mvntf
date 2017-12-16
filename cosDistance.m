function [meanDistance,Distance,sor]=cosDistance(res,original)
[row,col]=size(res);
meanDistance=Inf;
p=perms(1:col);
num=size(p,1);
sadMat=zeros(col,col);
for i=1:col
    temp=res(:,i);
    temp=repmat(temp,1,col);
%     size(temp)
%     size(diag(temp'*temp).*diag(original'*original))
%     dot(temp,original)
    sadMat(i,:)=acos(dot(temp,original)./(sqrt(diag(temp'*temp)).*sqrt(diag(original'*original)))');   
end
sadMat=sadMat';
meanDistance=Inf;
temp=[];
for i=1:num
    for j=1:col
        temp(j)=sadMat(j,p(i,j));
    end
     if mean(temp)<meanDistance
         sor=p(i,:);
         meanDistance=mean(temp);
         Distance=temp;
     end
end
  
%     temp=CalSAD(res,original,p(i,:));
%     if temp<sumDistance
%         sor=p(i,:);
%         sumDistance=temp;
%     end
% sumDistance=sumDistance;
end