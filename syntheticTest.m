%Demo for synthetic data
% d: dimension, k:number of clusters, cn:number of members for each cluster
% r: radius, a: amplitude of Gaussian noise
d=2;
k=300;
cn=50;
r=.01;
a=.0;

skelname=sprintf('syntheticresult_d%05d_k_%05d_n_%d_r01',d,k,cn);
skeldata=sprintf('syntheticdata_d%05d_k_%05d_n_%d_r01',d,k,cn);

%if 1==0
if exist(sprintf('%s.mat',skeldata))==0

dataset=makeTestDistributionData4(d,k,cn,r);
X=dataset.X;
figure(3)
plot(X(1,:),X(2,:),'.') 
%keyboard
noiseX=makeNormalizedRand(d,k*cn)*a;
X=X+noiseX;
mX=mean(X,2);
mY=std(X,[],2); Z=diag(1./mY)*(X-mX*ones(1,size(X,2))); 
Xcorr=Z'*Z;

dataset.Xcorr=Xcorr;
save(skeldata,'dataset');

else

load(skeldata);

end
%load(skeldata);

Xcorr=dataset.Xcorr;
X=dataset.X;
C=dataset.C;
E=dataset.E;
n=size(Xcorr,1);
k=size(E,1);

%calc actual KMDLs
hv=[1e-2,1e-1,1,10,100];
%hv=[1e-8,1e-4,1,10,100];

[acKMDL1,acKMDL2]=computeKMCRonly(dataset,hv);

n=size(X,2);
kv=[ 2.^[2:floor(log2(n)) floor(log2(n))] n ];
computeKMCR(Xcorr,kv,hv,skelname);

%end %1==0


load(skeldata);
Xcorr=dataset.Xcorr;
X=dataset.X;
C=dataset.C;
E=dataset.E;
n=size(Xcorr,1);
k=size(E,1);

%calc actual KMDLs
hv=[1e-2,1e-1,1,10,100];
%hv=[1e-16,1e-12,1e-8,1e-4,1];

[acKMDL1,acKMDL2]=computeKMCRonly(dataset,hv);


%plotKMDL2tmp3(skelname,12,k,acKMDL1,acKMDL2)
plotKMCRresult(skelname,12,k,acKMDL1,acKMDL2)

