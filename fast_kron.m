%*******************************************************************************
function C = fast_kron (A,B)
%FAST_KRON fast kronecker product
[I,L1]=size(A);
[J,L2]=size(B);

if (L1==1) && (L2==1)
    C=reshape(B*A.',I*J,1);
elseif (L1==1) && (L2>1)
    Bt=B.';
    C=reshape(Bt(:)*A.',L2,I*J).';
elseif (L2==1) && (L1>1)
    C=reshape(B*A(:).',I*J,L1);
else
    C=reshape(permute(reshape(B(:)*A(:).',[J,L2,I,L1]),[1 3 2 4]),[I*J,L1*L2]);
end
end
