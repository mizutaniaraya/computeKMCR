function computeKMCR_final(Xcorr,kv,hv,skelname)
%  computeKMCR_final(Xcorr,kv,hv,skelname)
%
%  Using litekmeans.m which is previously published in Matlab central
%  by Michael Chen
%
%  INPUT
%  Xcorr: NxP matrix (P column vectors are grouped by k)
%  kv   : vector of initial k (e.g. [ 2 3 4 5 ])
%  hv   : vector of quantization factor h
%  skelname: skelname.dat is the output file

%  OUTPUT
%  mink: 2x1 vector; ks that give minimum KMDL1/2 value
%  minkk: kks that give minimum KMDL1/2 value
%  minidx: 2xN matrix; index matrix for minimum KMDL1/2 value
%
%  OUTPUT FILE
%
%  k  kk  residual pMDL pMDL1 pMDL2
%
%  k: initial k
%  kk: resulting number of k
%  residual: residual
%  pMDLi: pseudo-MDL
%
if isempty(skelname)
    skelname='testMDL';
end

secondflag=1;

X2=dot(Xcorr,Xcorr);
a0=sum(X2);

for i=1:length(kv)
    
    k=kv(i);
    
    if k==0
        n=size(Xcorr,2);
        
        kk=0
        %    tic
        
        sumd=X2;
        ssum=sum(sumd);
        MDL=log(ssum);
        MDL1=    ssum;
        MDL2=    n*n*log2(ssum/n/n./h+1);
        MDL2new=MDL2./(n*n*log2(a0/n/n./h+1));
        
        yy=full([k kk a0 ssum MDL MDL1 MDL2new]);
        %save(sprintf('%s.dat',skelname),'yy','-append','-ascii')
        mink=0;minkk=0;minidx=0;
    else
        
        idx=litekmeans(Xcorr,k);
        n=length(idx)
        [u,~,label]=unique(idx);
        [~,sidx]=sort(idx);
        
        kk=length(u)
        %    tic
        E = sparse(1:n,label,1,n,kk,n);  % transform label into indicator matrix
        m = Xcorr*(E*spdiags(1./sum(E,1)',0,kk,kk));    % compute m of each cluster
        %    toc
        %    tic
        
        if kk~=1
            
            M=m(:,idx);
            m2=dot(m,m);
            M2=m2(idx);
            L2=bsxfun(@minus,dot(M,Xcorr),M2/2);
            sumd=X2-2*L2;
            ssum=sum(sumd);
            msum=sum(m2);
            MDL=log(ssum)+kk/n*log(a0);
            MDL1=a0/n*kk+    ssum;
            %    MDL1=msum/kk*kk+    ssum;
            MDL1=MDL1/a0; %normalization
            
            %h=a0/n/n/n*1;
            h=hv;
            %keyboard
            %%%% ( a0 is |X|^2_2 )
            
            MDL2=kk*n*log2(a0/n/n./h+1)+     n*n*log2(ssum/n/n./h+1);
            MDL2=kk*1*log2(a0/n./h+1)+     n*1*log2(ssum/n./h+1);
%            MDL2=kk*n*log2(a0/n/n./h)+     n*n*log2(ssum/n/n./h);
%            MDL2=kk*n*log2(msum/kk/n./h+1)+     n*n*log2(ssum/n/n./h+1);
%            MDL2new=MDL2./(n*log2(msum/kk/n./h));
            MDL2new=MDL2./(n*n*log2(a0/n/n./h+1));
            MDL2new=MDL2./(n*1*log2(a0/n/1./h+1));
%            normconst=(n*log2(msum/kk/n./h+1));
            normconst=(n*n*log2(a0/n/n./h+1));
            normconst=(n*1*log2(a0/n/1./h+1));
%    MDL2new=MDL2new*normconst;
%            normconst=1;
% normconst = kk*n*log2(a0/n/n/h+1)+     n*n*log2(ssum/n/n/h+1);
            
            um1=    n*n*log2(ssum/n/n./h+1)./normconst;
            
            %make new kv
            nmax=floor(log2(kk));
            %%%%    2nd stage
            if secondflag==1
                if nmax>1
                    disp('enter 2stage')
                    kkv=2.^[1:nmax];
                    mX2=dot(m,m);
                    %keyboard
                    for k2=kkv
                        iidx=litekmeans(m,k2);
                        nm=length(iidx)
                        [u,~,label]=unique(iidx);
                        [~,sidx]=sort(iidx);
                        
                        kkm=length(u)
                        if kkm~=1
                            E = sparse(1:nm,label,1,nm,kkm,nm);
                            mm = m*(E*spdiags(1./sum(E,1)',0,kkm,kkm));
                            M=mm(:,iidx);
                            m2=dot(mm,mm);
                            M2=m2(iidx);
                            L2=bsxfun(@minus,dot(M,m),M2/2);
                            sumd=mX2-2*L2;
                            mssum=sum(sumd);
                            ma0=sum(mX2);
                            mMDL1=a0/n*kkm + mssum;
                            mMDL2=kkm*n*log2(a0/n/n./h+1)+nm*n*log2(mssum/nm/n./h+1);
                            mMDL2=kkm*1*log2(a0/n/1./h+1)+nm*1*log2(mssum/nm/1./h+1);
                            %        um1=kkm*n*log2(a0/n/n/h+1);
                            um2=nm*n*log2(mssum/nm/n./h+1);
                            um3=n*n*log2(ssum/n/n./h+1);
                            
                            newMDL1=mMDL1+    ssum;
                            newMDL1=newMDL1/a0;
                            newMDL2=mMDL2+     n*n*log2(ssum/n/n./h+1);
                            newMDL2=mMDL2+     n*1*log2(ssum/n/1./h+1);
                            newMDL2new=newMDL2./normconst;
                            
                            yy=full([k kk a0 ssum MDL MDL1 MDL2new k2 kkm newMDL1 newMDL2new um1 um2 um3]);
                            save(sprintf('%s.dat',skelname),'yy','-append','-ascii')
                            
                        end %if kkm~=1
                        
                    end
                    
                    %      2nd stage
                end %1==0
                
            end %kk==1
            if secondflag~=1
            yy=full([k kk a0 ssum MDL MDL1 MDL2new]);
            save(sprintf('%s.dat',skelname),'yy','-append','-ascii')
            end 
            
        end
        %    yy=[k kk a0 ssum MDL MDL1 MDL2];
        
        if 1==0
            
            if i==1
                mink=[k;k];
                minkk=[kk;kk];
                minidx=[idx';idx'];
                lastMDL=[MDL1;MDL2'];
            else
                if MDL1<lastMDL(1)
                    mink(1)=k;
                    minkk(1)=kk;
                    minidx(1,:)=idx';
                    lastMDL(1)=MDL1;
                end
                if MDL2<lastMDL(2)
                    mink(2)=k;
                    minkk(2)=kk;
                    minidx(2,:)=idx';
                    lastMDL(2)=MDL2;
                end
            end
            
        end
        %save(sprintf('%s.dat',skelname),'yy','-append','-ascii')
        
        
    end %k==1
    
end %for


return


function label = litekmeans(X, k)
% Perform k-means clustering.
%   X: d x n data matrix
%   k: number of seeds
% Written by Michael Chen (sth4nth@gmail.com).
n = size(X,2);
last = 0;
label = ceil(k*rand(1,n));  % random initialization
%last = zeros(size(label));

while any(label ~= last)
    [u,~,label] = unique(label);   % remove empty clusters
    label=label(:)';
    k = length(u);
    E = sparse(1:n,label,1,n,k,n);  % transform label into indicator matrix
    m = X*(E*spdiags(1./sum(E,1)',0,k,k));    % compute m of each cluster
    last = label;
    [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1); % assign samples to the nearest centers
end
[~,~,label] = unique(label);

return

