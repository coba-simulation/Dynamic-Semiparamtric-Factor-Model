% ---------------------------------------------------------------------
% Quantlet:		NCspDSFM
% ---------------------------------------------------------------------
% Description:	Function for DSFM_YC
% ---------------------------------------------------------------------
% Author:		Piotr Majer (based on codes of Szymon Borak), 20120331
% ---------------------------------------------------------------------


function sDSFMobj = NCspDSFM(IV,L,knotsmoney,kmoney,knotsmatur,kmatur,startingzeta,epsilon,maxiter,dim1,dim2)

%IV should contain 4 column moneyness, maturity, log iv, date

sDSFMobj = 0;

% number of splines
lkmat = length(knotsmatur)-kmatur;
lkmon = length(knotsmoney)-kmoney;

%coefitions of splines 
coefs = zeros(lkmon,lkmat,L+1);

%days which are in the sample
days = unique(IV(:,4));
T = length(days); %number of days

[ro,co]=size(startingzeta);
if  ro~=T
    error('Wrong number of days in statring betas')
end


% setting the parameters

minmatur = min(IV(:,2));
maxmatur = max(IV(:,2));
minmoney = min(IV(:,1));
maxmoney = max(IV(:,1));

%grid
 dim1 = 18; %dim2 = 40; %it seems to be enough 
gridmoney = linspace(minmoney,maxmoney,dim1);
% gridmatur = linspace(minmatur,maxmatur,dim2);
gridmatur=knotsmatur;
% 
z = startingzeta;

% matrices for the data on the splines
B2 = ones(lkmat,length(IV(:,1))); 
B1 = ones(lkmon,length(IV(:,1))); 

[IVU1 ,IIND1,IND1]=unique(IV(:,1));
[IVU2 ,IIND2,IND2]=unique(IV(:,2));

SMART1 = spcol(knotsmoney,kmoney,IVU1)' ;
SMART2 = spcol(knotsmatur,kmatur,IVU2)' ;

B1 = SMART1(:,IND1);  
B2 = SMART2(:,IND2);

c = ones(lkmon*lkmat*(L+1),1);



%Xbase is the very big matrix with spline results
Xbase = ones(length(IV(:,2)),lkmon*lkmat);
for ii = 1:lkmat
    for jj = 1:lkmon
        Xbase(:,jj + (ii-1)*lkmon) = (B1(jj,:).*B2(ii,:))';    
    end    
end

obsperday = ones(length(days),1);
for tdate = 1: length(days)
    ind = find(IV(:,4)==days(tdate));
    obsperday(tdate)=length(ind);    
end     

obsperday=[0;cumsum(obsperday)];


explainedvariance = 2*epsilon;
maxloop = maxiter;

description.explainedvariance = 100;
explainedvarianceold = explainedvariance+10*epsilon;


breakingrule = 2*epsilon;
riteration = IV(:,3);
r1iteration = 1000*epsilon*IV(:,3);

meanIV = mean(IV(:,3));
denominator =  sum((IV(:,3)-meanIV).^2);

description.breakingrule = breakingrule;
%%%%% 
%%%%commenent all from here
%


while (maxloop>0  &  breakingrule>epsilon )
    
    maxloop
    
    [Ynew, XXnew] = smartMatrix4sDSFM(Xbase,IV(:,3),z,obsperday);
    
    c = XXnewYnew;
    
    for i=1:(L+1)
        coefs(:,:,i)= xplreshape(c(((i-1)*lkmon*lkmat+1):((i)*lkmon*lkmat)),lkmon,lkmat);
    end
    
    
    % 
    %'regress on the time series'
    
    nominator=0;
%    denominator =0;

    modeliv=0;
    
    %riteration = 0;
    
    zz = z;
    
    for tdate = 1:T
        ind = find(IV(:,4)==days(tdate)); 
        Y = IV(ind,3);
        M0 = zeros(length(Y),1);
        M1L= zeros(length(Y),L);
        for j = 1:length(Y)
            MONMON = B1(:,ind(j))';
            MATMAT = B2(:,ind(j));
            %MONMON=spcol(knotsmoney,kmoney,IV(ind(j),1));
            %MATMAT=spcol(knotsmatur,kmatur,IV(ind(j),2)).';
%             M0(j)= MONMON*coefs(:,:,1)*MATMAT;  
            for ll = 1:(L+1)
                M1L(j,ll)=MONMON*coefs(:,:,ll)*MATMAT;
            end
        end
        ztemp = M1LY;
        
        z(tdate,1:(L+1))=ztemp';
        
        riteration(ind)=M1L*ztemp;
        
        nominator = nominator + sum((M1L*ztemp -Y).^2); 
        Yhat(tdate,:)=M1L*ztemp;
        
    end %for
    
    breakingrule = mean((riteration - r1iteration).^2);
    r1iteration = riteration;
    
    explainedvarianceold = explainedvariance;
    explainedvariance = nominator/denominator
    description.explainedvariance = [description.explainedvariance;explainedvariance]; 
    description.breakingrule=[description.breakingrule;breakingrule];
    
    
    maxloop=maxloop-1;  
    
end %while

%orthogonalization and norming part. 

mhat = ones(length(gridmoney),length(gridmatur),L+1);
for i = 1:(L+1)
    mhat(:,:,i) = spcol(knotsmoney,kmoney,gridmoney)*coefs(:,:,i)*spcol(knotsmatur,kmatur,gridmatur).'; 
end

%integral part
du = (gridmoney(2)-gridmoney(1))*(gridmatur(2)-gridmatur(1));

% create matrix for the integral purpose
% 122221
% 244442
% 244442
% 244442
% 122221

tempmatrix = 0*mhat(:,:,1)+1;
tempmatrix(2:(end-1),:)= 2*tempmatrix(2:(end-1),:);
tempmatrix(:,2:(end-1))= 2*tempmatrix(:,2:(end-1));

% matrices
GAMMA = ones(L+1);
% gamma = ones(L,1);

for i = 1:L+1
%     gamma(i)=sum(sum(tempmatrix.*mhat(:,:,1).*mhat(:,:,i+1)))*du/4;
    for j = 1:L+1
        GAMMA(i,j)=sum(sum(tempmatrix.*mhat(:,:,j).*mhat(:,:,i)))*du/4;
    end
end

%vectorize functions
MhatMatrix = xplvec(mhat(:,:,1));
for i=2:(L+1)
    MhatMatrix=[MhatMatrix,xplvec(mhat(:,:,i))]; 
end

%calculate znew and mhatnew
znew = z';
for i = 1:length(z(:,1))
    znew(1:end,i) = (GAMMA^0.5) * (z(i,1:end)' );
end
znew = znew';

MhatMatrixnew = MhatMatrix';
for i = 1:length(MhatMatrix(:,1))
%     MhatMatrixnew(1,i) = MhatMatrix(i,1)' -gamma'*GAMMA^-1*MhatMatrix(i,2:end)';   
    MhatMatrixnew(1:end,i) = (GAMMA^-0.5) * MhatMatrix(i,1:end)';
end
MhatMatrixnew=MhatMatrixnew';

% order wrt variance
B = znew(:,1:end)'*znew(:,1:end);
[V,D] = eigs(B);

zhat = znew';
for i = 1:length(znew(:,1))
    zhat(1:end,i) = V'*(znew(i,1:end))';
end
zhat = zhat';

MhatMatrix =  MhatMatrixnew';
for i = 1:length(MhatMatrixnew(:,1))
    MhatMatrix(1:end,i) = V' * MhatMatrixnew(i,1:end)';
end
MhatMatrix=MhatMatrix';

mhatort = mhat;
[r,c]=size(mhat(:,:,1));
for i = 1:(L+1)
    mhatort(:,:,i)=xplreshape(MhatMatrix(:,i),r,c) ; 
end

% make a spline representation of the functions
coefsotr = coefs;
for i = 1:(L+1)
    monmon=spap2(knotsmoney,kmoney,gridmoney,mhatort(:,:,i)');
    coefmonmon=monmon.coefs;
    matmon=spap2(knotsmatur,kmatur,gridmatur,coefmonmon.');
    coefsort(:,:,i) = matmon.coefs;
    values = spcol(knotsmoney,kmoney,gridmoney)*coefsort(:,:,i)*spcol(knotsmatur,kmatur,gridmatur).';
    max(max(values-mhatort(:,:,i)));
end

% give the results outside

description.eps = epsilon;
description.iter = maxiter-maxloop;

sDSFMobj.IV = IV;

sDSFMobj.z = zhat;
sDSFMobj.Yhat = Yhat;
sDSFMobj.mcoefs = coefsort;
sDSFMobj.knotsmoney=knotsmoney;
sDSFMobj.kmoney=kmoney;
sDSFMobj.knotsmatur=knotsmatur;
sDSFMobj.kmatur=kmatur;

sDSFMobj.mhat = mhatort;
sDSFMobj.gridmoney = gridmoney;
sDSFMobj.gridmatur = gridmatur;

sDSFMobj.description = description;