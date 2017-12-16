function  abundanceRes=normAbundance(btdRes,massNum,L)
rankVector=ones(1,massNum)*L;
for i=1:massNum
    W=btdRes{1}(:,sumRank(rankVector,i-1)+1:sumRank(rankVector,i));
    H=btdRes{2}(:,sumRank(rankVector,i-1)+1:sumRank(rankVector,i));
    abundanceRes(:,:,i)=W*H';
end
[width,len,height]=size(abundanceRes);
abundanceRes=reshape(abundanceRes,[width*len,massNum]);
temp=width*len;


%  for i=1:temp
%      abundanceRes(i,:)=abundanceRes(i,:)/sum(abundanceRes(i,:));
% end

end
