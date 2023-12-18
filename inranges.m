% inranges          determine which elements of a vector are in which range
%
% CALL              [ idx, seg, out ] = inranges( vec, mat, flag, forced )
%
% GETS              
%   vec:              vector (of indices). does not have to be sorted.
%   mat:              2-column matrix, each row is a set of two numbers
%                       indicating a range. does not have to be sorted
%   flag:             {0}; influences the 3rd output: if 1, will subtract
%                       lower range of the matching row of mat
%   forced:           {0}; if true, will force each vector element to a range.
%                       1: nearest range; 2: rounded up; -1: rounded down
%
% RETURNS 
% idx   indices of vec which are in one of the ranges of mat: 
%           vec( idx ) are the data in the ranges
% seg   the row number to which the element of idx corresponds s.t.
%           vec( idx( seg == SEG ) ) >= mat( SEG, 1 )
%           vec( idx( seg == SEG ) ) <= mat( SEG, 2 )
%           are all true
% out   the actual values in the ranges, given by out = vec( idx )
%           if flag is 1, then 
%               out( j ) = vec( idx( j ) ) - mat( i, 1 ) + 1
%           where i is the row of mat for which 
%               vec( idx( j ) ) >= mat( i, 1 ) 
%               vec( idx( j ) ) <= mat( i, 2 ) 
%           are both true
%
% NOTES 
% 1. if mat rows are not sorted (i.e. mat( i, 1 ) > mat( i  + 1, 1 ) for some i), 
%       then the seg will not be sorted either but the correspondence will be exact
% 2. if vec is not sorted, then the output will not be sorted either, but
%       the correspondence (with the sorted mat..) will be exact
% 
% USAGE: many... 
%   to count occurrences, check the length of the first output (not logical)
%   to generate a histogram, use the second argument
%   to create a PETH, call with flag on and use the third output
%
% CALLS         resort, sortranges
%
% see also dilutesegments, getdatainranges, geteventsinranges, intersectranges, isoverlap, plotranges, setdiffranges, sortranges, uniteranges

% linear algorithm O(n1+n2)

% 27-jan-12 ES

% revisions
% 04-sep-12 modified seg (bug)
% 26-dec-12 modified output in case mat = []; (should be empty)
% 03-mar-13 adapted for non-sorted data
% 11-jul-13 adapted seg also (inverse mapping)
% 26-dec-13 adapted for non-sorted input
% 09-oct-14 added a forced option (all vec elements are assigned to a range)
% 17-aug-19 cleaned up

function [ idx, seg, out ] = inranges( vec, mat, flag, forced )

% initialize output
idx                             = [];
seg                             = [];
out                             = [];

% arguments
nargs = nargin;
if isempty( vec ) || isempty( mat ) || nargs < 2
    return
end
if nargs < 3 || isempty( flag )
    flag                        = 0;
end
if nargs < 4 || isempty( forced )
    forced                      = 0;
end

vec                             = vec( : );
if size( mat, 2 ) ~= 2
    error( 'input size mismatch' )
end

% pre-process
n1                              = length( vec );
n2                              = size( mat, 1 );
if ~issorted( vec )
    mixed                       = 1;
    [ ridx, vec ]               = resort( vec );
else
    mixed                       = 0;
end
[ mat, sidx ]                   = sortranges( mat, 0 );

% actual algorithm
i1                              = 1;
i2                              = 1;
j                               = 0;
idx                             = zeros( n1, 1 );
seg                             = idx;
out                             = idx;
while i1 <= n1
    if vec( i1 ) > mat( i2, 2 )
        i2                      = i2 + 1;
    else
        if vec( i1 ) >= mat( i2, 1 )
            j                   = j + 1;
            idx( j )            = i1;
            seg( j )            = sidx( i2 );
            if flag
                out( j )        = vec( idx( j ) ) - mat( i2, 1 ) + 1;
            else
                out( j )        = vec( idx( j ) );
            end
        end
        i1                      = i1 + 1;
    end
    if i2 > n2
        break
    end
end

% post-process
if mixed
    seghat                      = zeros( n1, 1 );
    outhat                      = zeros( n1, 1 );
    idxhat                      = zeros( n1, 1 );
    seghat( idx( idx > 0 ) )    = seg( idx > 0 );
    outhat( idx( idx > 0 ) )    = out( idx > 0 );
    idxhat( idx( idx > 0 ) )    = idx( idx > 0 );
    idx                         = idxhat( ridx );
    seg                         = seghat( ridx );
    out                         = outhat( ridx );
    rmv                         = idx == 0;
    idx                         = find( ~rmv );
    seg( rmv )                  = [];
    out( rmv )                  = [];
else
    rmv                         = false( n1, 1 );
    rmv( ( j + 1 ) : n1 )       = 1;
    idx( rmv )                  = [];
    seg( rmv )                  = [];
    out( rmv )                  = [];
end

% force out of range values (inefficient code..):
if forced && length( idx ) ~= n1
    fidx                        = ( length( idx ) + 1 ) : n1;
    vec                         = vec( ridx );
    for i                       = 1 : length( fidx )
        dx                      = abs( vec( fidx( i ) ) - mat );
        if forced == 1              % closest
            row                 = find( sum( dx == min( dx( : ) ), 2 ) );
        elseif forced == 2          % above
            row                 = find( ( vec( fidx( i ) ) - mat( :, 2 ) ) < 0, 1, 'first' );
        elseif forced == -1         % below
            row                 = find( ( vec( fidx( i ) ) - mat( :, 2 ) ) < 0, 1, 'first' ) - 1;
        end
        if row < 1
            row                 = 1;
        elseif row > n2
            row                 = n2;
        end
        j                       = j + 1;
        idx( j, : )             = fidx( i ); 
        seg( j, : )             = row;
        if flag
            out( j, : )         = vec( fidx( i ) ) - mat( row, 1 );
        else
            out( j, : )         = vec( fidx( i ) );
        end
    end
end

return

% EOF

vec = ( 1 : 20 )' + 100;
mat = [ 1 5; 8 11; 16 18 ] + 100;
mat( 1, : ) = mat( 1, [ 2 1 ] );
mat( [ 2 3 ], : ) = mat( [ 3 2 ], : );
vechat = mixmat( vec );

[ idx seg ] = inranges( vechat, mat ); mat, for i = 1 : length( unique( seg( : ).' ) ), vechat( idx( seg == i ) ), end
[ idx seg ] = inranges( vec, mat ); mat, for i = 1 : length( unique( seg( : ).' ) ), vec( idx( seg == i ) ), end

mat = reshape( sort( round( rand( 4, 1 ) * 1e3 ) ), [ 2 2 ] )'; vec = round( rand( 100, 1 ) * 1e3 ); uvec = unique( vec );
