function scores=IHRW(H)
%H: vertex-hyperedge incident matrix 
%scores: Predicted score. Specifically, the combination score of drug i,
%drug j, and drug k is scores(i,j,k).
%H=xlsread('lung_eff_dianbianguanlian.xlsx');
%H=xlsread('breast_eff_dianbianguanlian.xlsx');
%H=xlsread('colon_eff_dianbianguanlian.xlsx');

c=0.2;%Restart probability
[cc,dd]=size(H);

Dv=diag(sum(H,2));
De=diag(sum(H));
w=inv(Dv)*H*inv(De)*H';

for k=1:cc
    for i=1:cc
        p0=zeros(cc,1);
        p0(i)=1/2;
        p0(k)=1/2;
        pt(:,i,k)=(1-c)*inv(diag(ones(1,cc))-c*w')*p0;%%
    end
end
scores=zeros(cc,cc,cc);
for i=1:cc
    for j=1:cc
        for k=1:cc
            scores(i,j,k)=1/3*(pt(i,j,k)+pt(j,k,i)+pt(k,i,j));
        end
    end
end
for ii=1:cc
    scores(ii,:,ii)=0;
    scores(ii,ii,:)=0;
    scores(:,ii,ii)=0;
end
end
