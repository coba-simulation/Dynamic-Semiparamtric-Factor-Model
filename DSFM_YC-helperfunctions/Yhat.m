% ---------------------------------------------------------------------
% Quantlet:		Yhat
% ---------------------------------------------------------------------
% Description:		Function for DSFM_BS_YC
% ---------------------------------------------------------------------
% Author:		Piotr Majer, 20110714
% ---------------------------------------------------------------------

function y=Yhat(XY,Z,mcoefs,knotsmoney,kmoney,knotsmatur,kmatur,L,J);
tic

tt=unique(XY(:,4));
nominator=0;
tn=length(Z(:,1));
Yhat=0;
for t=1:length(tt)
    ind=find(XY(:,4)==tt(t));
    obsperday=length(ind);
    Xtmoney=XY(ind,1);
    Xtmatur=XY(ind,2);
    for k=1:obsperday
        model=0;
        SplineMon=spcol(knotsmoney,kmoney,Xtmoney(k));
        SplineMat=spcol(knotsmatur,kmatur,Xtmatur(k)).'; 
            for j=1:L+1
            model = model + Z(t,j).*SplineMon*mcoefs(:,:,j)*SplineMat;
            end
            Yhat=[Yhat; model];
            
    end
end
% denominator=sum((XY(:,3)-mean(XY(:,3))).^2);
 Yhat=Yhat(2:end);
 size(Yhat);
 y=reshape(Yhat,J,tn)';
