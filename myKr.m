function Mat = myKr(A,B,partA,partB)
[J M]=size(A);
[K N]=size(B);
partionNum=size(partB,2);
indA=[0 cumsum(partA)];
indB=[0 cumsum(partB)];
indMat=[0 cumsum(partA.*partB)];
Mat=zeros(J*K,sum(partA.*partB));
for i=1:partionNum
     Mat(:,indMat(i)+1:indMat(i+1))= fast_kron(A(:,indA(i)+1:indA(i+1)) , B(:,indB(i)+1:indB(i+1)));
end
% Mat=temp;
end