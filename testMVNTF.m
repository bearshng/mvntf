load xiong_25dB.mat;
load syntheticDataIter1rank19.mat;
realEndmember=A;
realAbundance=s;
tensorData=reshape(x_n',64,64,224);
% btdInit=temp.btdInit;
btdInit{1}=max(btdInit{1},1e-4);
btdInit{2}=max(btdInit{2},1e-4);
btdInit{3}=max(btdInit{3},1e-4);
EndNum=6;
rankNumber=19;
options.derta=0.4;
options.convergeNum=5;
allSad=zeros(1,10);
allRmse=zeros(1,10);
t=mvntf(tensorData,EndNum,rankNumber,btdInit{1},btdInit{2},btdInit{3},options);
[sad,Distance,sor]=cosDistance(t{1},realEndmember);
[rmse]=HyperRmse(t{2},realAbundance,sor);
