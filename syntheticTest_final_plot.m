d=10000;
k=300;
cn=50;
r=.01;
a=.0;

skelname=sprintf('testlooper0214_final2_d%05d_k_%05d_n_%d_r01',d,k,cn);
%skelname=sprintf('test');
skeldata=sprintf('datalooper0213_final1_d%05d_k_%05d_n_%d_r01',d,k,cn);


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
mY=std(X,[],2); 
Z=diag(1./mY)*(X-mX*ones(1,size(X,2))); 
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
%hv=[1e-16,1e-12,1e-8,1e-4,1];

[acKMDL1,acKMDL2]=calcKMDL(dataset,hv);

n=size(X,2);
if 1==0
kv=[ 2.^[2:floor(log2(n)) floor(log2(n))] n ];
%kv=[4024:100:8048]
%kmeans_arxiv_newtmp_h2(Xcorr,kv,hv,skelname);
%kmeans_arxiv_newtmp_h(Xcorr,kv,hv,skelname);
calc_KMCR_final(Xcorr,kv,hv,skelname);

%end %1==0

end %1==0
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

[acKMDL1,acKMDL2]=calcKMDL(dataset,hv);


plotKMDL2tmp3(skelname,12,k,acKMDL1,acKMDL2)

