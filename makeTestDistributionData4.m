function dataset=makeTestDistribution4(d,k,n,r)
%
% cluster centers are located randomly at R=1 positions.
% -> cluster centers are located randomly at R<1 positions.
% cluster members are located randomly from cluster center and the distances are less than r.
% -> using Gaussian distribution (normrnd(mu,sigma,m,n))
%
% %%%  how to make random distributeion %%%
%
% cluster centers are distributed around r=1\pm sigma 
% cluster members are distributed around each cluster center by Gaussian
%
%
% d: dimension
% k: number of clusters
% n: number of cluster members
% r: radius of distribution
%

counter=0;

% define cluster center vector (d \times k matrix)
C=makeNormalizedRand(d,k);
C=C*spdiags(normrnd(1,0.25,k,1),0,k,k);

X=[];
for i=1:k
X0=C(:,i)*ones(1,n)+2*r*makeNormalizedRand(d,n)*spdiags(normrnd(1,1,n,1),0,n,n);

X=horzcat(X,X0);
end

dataset.X=X;
dataset.C=C;
dataset.d=d;
dataset.k=k;
dataset.n=n;
dataset.r=r;
dataset.E=kron(speye(k,k),sparse(ones(1,n)));
return

function C=makeNormalizedRand(d,k)

A=rand(d,k)-1/2;
S=dot(A,A);
C=A*spdiags(sqrt(S').^-1,0,k,k);


return
