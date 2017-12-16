function Mat = myKr2(A,B,partA,partB)
temp=[];
partionNum=size(partB,2);

for i=1:partionNum
    t=A(:,sumRank(partA,i-1)+1:sumRank(partA,i))*B(:,sumRank(partB,i-1)+1:sumRank(partB,i))';
    temp=[temp t(:)];
end
Mat=temp;
end