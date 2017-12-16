function  [result]=HyperRmse(abundance,realAbundance,sor)
%description: this function calculates the root mean square error of
%two abudance maps
imageSize=size(realAbundance,1);
result=mean(sqrt(sum((realAbundance-abundance(:,sor)).^2,1)./imageSize));



% abundance=abundance'*permMatrix;
% for i = 1:6
%     RMSE_Lh_tmp(i) = rmse(abundance(:,i),realAbundance(:,i));
% end
% % AllRMSE_Lh(iter,:) = RMSE_Lh_tmp;
% result =  mean(RMSE_Lh_tmp,2);
  
end

% 
% function [re]=rmse(data,estimate)
% I = ~isnan(data) & ~isnan(estimate);
% data = data(I); estimate = estimate(I);
% 
% re=sqrt(sum((data(:)-estimate(:)).^2)/numel(data));
% end


