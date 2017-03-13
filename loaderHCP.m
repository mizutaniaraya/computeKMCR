
fname=sprintf('HCP_anat_decimby%d.mat',decim);
a=load(fname);
% a.A, a.B, a.maskindv is loaded

mA=mean(a.A);
dA=std(a.A);
%devA=a.A-ones(size(a.A,1),1)*mA;
n=size(a.A,2);

zA=(a.A-ones(size(a.A,1),1)*mA)*spdiags(1./dA',0,n,n);

X=zA' * zA;

if 1==0

Nrep=1;
k=200;

for ii=1:Nrep
ii
n = size(X,2);
if 1==0
indl=find(a.B(1,:)<=0);
indr=find(a.B(1,:)>0);
initlabel=ones(1,n);
initlabel(indl)=ceil(k/2*rand(1,length(indl)));
initlabel(indr)=ceil(k/2*rand(1,length(indr))+k/2);
end
indlr=find(a.B(1,:)<=0 & a.B(2,:)<=0);
indlf=find(a.B(1,:)<=0 & a.B(2,:)>0);
indrr=find(a.B(1,:)>0 & a.B(2,:)<=0);
indrf=find(a.B(1,:)>0 & a.B(2,:)>0);
initlabel=ones(1,n);
llr=length(indlr);
llf=length(indlf);
lrr=length(indrr);
lrf=length(indrf);
lenall=llr+llf+lrr+lrf;

initlabel(indlr)=ceil(k/lenall*llr*rand(1,length(indlr)));
initlabel(indlf)=ceil(k/lenall*llf*rand(1,length(indlf))+k/lenall*llr);
initlabel(indrr)=ceil(k/lenall*lrr*rand(1,length(indrr))+k/lenall*(llr+llf));
initlabel(indrf)=ceil(k/lenall*lrf*rand(1,length(indrf))+k/lenall*(llr+llf+lrr));

%initlabel = ceil(k*rand(1,n));  % random initialization

[label,iter,mcoordmat]=litekmeans_rev2(X,k,a.B,initlabel);

[~,sidx]=sort(label);


k=length(unique(label));

E = sparse(1:n,label,1,n,k,n);  % transform label into indicator matrix
%mcoord = a.B*spdiags(mA',0,n,n)*(E*spdiags(1./sum(E,1)',0,k,k));    % compute m of each cluster
% simple average position 
mcoord = a.B*(E*spdiags(1./sum(E,1)',0,k,k));    % compute m of each cluster

figure(1)
plot(mcoord(2,:),mcoord(3,:),'.')
if ii==1
   hold on
end
drawnow

end

[~,sortmcoordidx]=sort(mcoord(2,:));
% sortmcorrdidx = 1...k

[~,sidx]=sort(sortmcoordidx(label));

figure(2)
imagesc(X(sidx,sidx))

end
