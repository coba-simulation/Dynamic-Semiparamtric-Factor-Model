% ---------------------------------------------------------------------
% Quantlet:		xplreshape
% ---------------------------------------------------------------------
% Description:	Function for DSFM_BS_YC
% ---------------------------------------------------------------------
% Author:		Szymon Borak, 20080414
% ---------------------------------------------------------------------

function A = xplreshape(x,r,c)

if r*c ~= length(x)
error ('wrong numbers of rows and column')
end

A = ones(r,c);

for i = 1:c
    A(:,i)= x(((i-1)*r+1):(i)*r);
end
