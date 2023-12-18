% sortranges            to be a set of non-overlapping [ small large ] pairs
%
% CALL                  [ y, idx ] = sortranges( x, flag, dx )
%
% GETS                  x         a 2-column matrix of numbers corresponding to ranges
%                       flag      {1} combine ranges with consective borders  
%                                   0 does not combine (e.g. [ 1 2; 2 3 ] will remain)
%                       dx        bin size, defaults to 1 (for whole integers)
%
% RETUNRS               y         sorted version of x
%                       idx       indices of x correponding to y's rows. if all of x's rows are
%                                     non-overlapping, this is simply y = x( idx, : ) as in sort
%                                     if any rows were combined, idx( i ) corresponds to the 
%                                     last row in x contributing to y( i, : )
%
% example:
%   x = [ 10 20; 55 60; 50 58; 70 80 ];
%   [ y idx ] = sortranges( x )
%   x( idx, : )
%
% see also dilutesegments, getdatainranges, geteventsinranges, inranges, intersectranges, isoverlap, plotranges, setdiffranges, uniteranges

% linear (or log if unsorted) algorithm

% 03-mar-13 ES

% revisions
% 03-may-13 (1) non-integer ranges supported
%           (2) idxx output too
% 17-aug-19 cleaned up
% 10-sep-19 added help for flag and dx

function [ y, idxy, idxx ] = sortranges( x, flag, dx )

[ m, n ]                = size( x );
if isempty( x )
    y                   = [];
    idxy                = [];
    idxx                = [];
    return
end
if n ~= 2 || ~isreal( x )
    error( '%s: not a matrix of ranges\n', upper( mfilename ) )
end
if m == 1
    y                   = [ min( x ) max( x ) ];
    idxy                = 1;
    idxx                = 1;
    return
end

if sum( x( :, 1 ) <= x( :, 2 ) ) ~= m
    sidx                = x( :, 1 ) > x( :, 2 );
    t                   = x( sidx, : );
    x( sidx, : )        = t( :, [ 2 1 ] );
end

if sum( diff( x( :, 1 ) ) < 0 )
    [ x( :, 1 ), sidx ] = sort( x( :, 1 ) );
    x( :, 2 )           = x( sidx, 2 );
else
    sidx                = ( 1 : m )';
end

if exist( 'flag', 'var' ) && flag == 0
    y                   = x;
    idxy                = sidx;
    idxx                = sidx;
    return
end

y                       = zeros( m, 2 );
idxy                    = zeros( m, 1 );
idxx                    = zeros( m, 1 );
if all( x( : ) == round( x( : ) ) )
    dx                  = 1;
elseif ~exist( 'dx', 'var' ) || isempty( dx )
    dx                  = eps;
end
i                       = 1; 
j                       = 1;
y( j, : )               = x( i, : );
idxy( 1 )               = sidx( 1 );
idxx( 1 )               = 1;
for i                   = 1 : ( m - 1 )
    if x( i + 1, 1 ) > ( y( j, 2 ) + dx )
        j               = j + 1;
        y( j, : )       = x( i + 1, : );
        idxy( j )       = sidx( i + 1 );
    else
        y( j, 2 )       = max( y( j, 2 ), x( i + 1, 2 ) );
    end
    idxx( i + 1 )       = j;
end
if j <= m
    y( j + 1 : m, : )   = [];
    idxy( j + 1 : m )   = [];
end
idxx                    = idxx( sidx );

return

% EOF
