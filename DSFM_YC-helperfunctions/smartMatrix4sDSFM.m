% ---------------------------------------------------------------------
% Quantlet:		smartMatrix4sDSFM
% ---------------------------------------------------------------------
% Description:	Function for DSFM_BS_YC
% ---------------------------------------------------------------------
% Author:		Szymon Borak, 20080414
% ---------------------------------------------------------------------

function [Y, XX] = smartMatrix4sDSFM(Xbase,IV,z,obsperday)


[r,c] = size(z);

[rb,cb]=size(Xbase);

L = c-1;

if obsperday(end)~=length(IV)
    error('wrong number of observations in obsperday');
end

if r~=(length(obsperday)-1)
    error('wrong number of observations in z or obsperday ');
end


XX = zeros((L+1)*cb,(L+1)*cb);
Y = zeros((L+1)*cb,1);




for k = 2:length(obsperday)
    XX = XX + kron(z(k-1,:)'*z(k-1,:),Xbase((obsperday(k-1)+1):obsperday(k),:)'*Xbase((obsperday(k-1)+1):obsperday(k),:) );
    
    Y = Y+ kron( z(k-1,:)' , Xbase((obsperday(k-1)+1):obsperday(k),:)' * IV((obsperday(k-1)+1):obsperday(k)));
end


