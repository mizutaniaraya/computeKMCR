function [MDL1,MDL2]=calcKMDL(dataset,amp)
% calcKMDL

Xcorr=dataset.Xcorr;
X=dataset.X;
C=dataset.C;
E=dataset.E;
n=size(Xcorr,1);
k=size(E,1);

% X=C*E+R;

R=X-C*E;

% Xcorr=CC*E+RR

Ninv=diag(sum(E,2))^-1;

CC=Xcorr*E'*Ninv;
RR=Xcorr-CC*E;

sumd=sum(dot(RR,RR));
a0  =sum(dot(Xcorr,Xcorr));

MDL1=a0/n*k+sumd;
MDL1=MDL1/a0;
h=(a0/n/n)./amp;
h=amp;

MDL2=k*n*log2((a0/n/n)./h+1)+n*n*log2((sumd/n/n)./h+1);
MDL2=MDL2./(n*n*log2((a0/n/n)./h+1));
MDL2=k*1*log2((a0/n/1)./h+1)+n*1*log2((sumd/n/1)./h+1);
MDL2=MDL2./(n*1*log2((a0/n/1)./h+1));



return
