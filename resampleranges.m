% resampleranges          zero offset, keeping the total duration fixed
%
% out = resampleranges( in, p, q )
%
% in        matrix of ranges (or a vector of samples indices)
% p, q      up/downsampling integers (see resample.m)
%
% out       resampled ranges 
%
% example: to upsample a signal sampled at Fs by p, use either
% >> out = resampleranges( in, p * Fs, Fs );
% >> out = resampleranges( in, p, 1 );
% 
% it is guaranteed that 
% max( abs( in - resampleranges( resampleranges( in, p, 1 ), 1, p ) ) ) == 0
% max( abs( in - resampleranges( resampleranges( in, 1, p ), p, 1 ) ) ) < p
%
% see also          resample

% Algorithmic notes:
% the key is to realize that range upsampling is not a scalar function, i.e. for
%       out = resample( [ t0 t0 ], p, 1 ); a = out( 1 ); b = out( 2 )
% a != b

% 10-mar-13 ES

% revisions
% 19-may-13 previous version sorted the ranges, leading to potential
%               i/o mismatch. this version only sorts the rows
% 17-aug-19 cleaned up

function out = resampleranges( in, p, q )

% arguments
nargs           = nargin;
if nargs < 1 || isempty( in )
    out         = [];
    return
end 
siz             = size( in );
if length( siz ) > 1 && siz( 2 ) == 2
    in          = sort( in, 2 );
else
    in          = in( : );
end
if nargs < 2 || isempty( p )
    p           = 1;
end
if nargs < 3 || isempty( q )
    q           = 1;
end
if abs( round( p ) ) ~= p || p == 0
  error( '''P'' must be a positive integer' )
end
if abs( round( q ) ) ~= q || q == 0
  error( '''Q'' must be a positive integer' )
end

% actually resample
[ p, q ] = rat( p/q, 1e-12 );
up              = upsampleranges( in, p );
out             = round( downsampleranges( up, q ) );

if ~isequal( size( out ), siz )
    out         = reshape( out, siz );
end

return

%----------------------------------------------------------------------
% downsampleranges
%----------------------------------------------------------------------
function out = downsampleranges( in, fac )
if fac == 1
    out                 = in;
else
    out                 = ceil( in / fac );
end
return

%----------------------------------------------------------------------
% upsampleranges
%----------------------------------------------------------------------
function out = upsampleranges( in, fac )
if fac == 1
    out                 = in;
else
    sgn                 = sign( in( :, 1 ) );
    out( sgn ==  1, 1 ) = ( in( sgn ==  1, 1 ) - 1 ) * fac + 1;
    out( sgn == -1, 1 ) = ( in( sgn == -1, 1 ) + 0 ) * fac + 0;
    out( sgn ==  0, 1 ) = ( in( sgn ==  0, 1 ) + 0 ) * fac + 0;
    if size( in, 2 ) == 2
        out( :, 2 )     = in( :, 2 ) * fac;
    end
end

return

% EOF

