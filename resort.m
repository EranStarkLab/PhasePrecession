% resort        indices to recover original order
%
% [ idx sx ] = resort( x )
%
% works on the columns of x
% gives the resorting indices, such that if
%
% sx = sort( x )
% 
% then
% 
% sx( idx ) == x

% 29-sep-13 ES

function [ idx, aa ] = resort( x )

[ aa, bb ] = sort( x );
[ ~, idx ] = sort( bb );

return

% EOF