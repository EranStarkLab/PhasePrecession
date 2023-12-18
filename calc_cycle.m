% calc_cycle       divide a phase vector to cycles 
%
% CALL             [cycPeriods, cycVec] = calc_cycle(phs)
%
%
% GETS              phs                  a vecotor of phases in radians 
%                   
%
% RETURNS           cycPeriods           the indices of the start and end of each cycle, phsFs
%                   cycVec               a vector in the size of phs with the cycle at each index
%
%                   NOTICE: if the vector starts with NaN values, that is
%                           considered cyc = 1
%
% written by        HES      23-Nov-21
% modified          HES      18-Dec-23

function [cyc_range, cyc]     = calc_cycle(phs, tol) 

%--------------------------------------------------------------------%
% check inputs
%--------------------------------------------------------------------%
nargs                          = nargin;

if nargs < 1 || isempty( phs ) 
    error( 'missing arguments' )
end

if nargin < 2 || isempty(tol)
    tol                         = 0.1;
end

phs                             = wrapTo2Pi( phs );
diffphs                         = diff( phs );
ch_cyc                          = find( abs( diffphs - min(diffphs) ) <= tol );
cyc_range1                      = [find(~isnan(phs),1); ch_cyc(:)+1];
cyc_range2                      = [cyc_range1(2:end)-1; find(~isnan(phs),1,'last')];
cyc_range                       = [cyc_range1, cyc_range2];
ncyc                            = size(cyc_range,1);

if nargout == 1
    return;
end

cyc                             = NaN( length(phs),1 );
for i = 1 : ncyc
    idx                         = cyc_range(i,1):cyc_range(i,2);
    cyc(idx)                    = i; 
end


return;

% EOF