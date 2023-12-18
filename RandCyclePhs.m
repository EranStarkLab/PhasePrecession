% RandCyclePhs         randomize the phase of spikes within each theta cycle
%
% CALL                 [rand_phs, rand_cphs] = RandCyclePhs( spk, phs, periods)
%
%
% GETS                 spk                 time of spikes [phsFs]
%                      phs                 phase vector for the whole session [phsFs]
%                      periods             field times [phsFs]
% OPTIONAL  
%                      phsFs               {1250}
%                      model               {'uniform'} other options: LFP, spk, spk_rate
%                                               * spk_rate mode is not recommended
%                      nreps               {500} number of random repetitions 
%                      nbins               {50} number of bins for phase
%                      graphics            { 1 } 
%
% RETURNS
%                      rand_phs            phs of randomized spikes
%                                            * rows - each spike
%                                            * columns   - repetitions
%                      rand_cphs           cphs of randomized spikes
%                       
%
% CALLS                nothing
%
% written by           HES and ES 23-Nov-21
% modified             HS 18-Dec-23

function [rand_phs, rand_cphs, rand_cyc] = RandCyclePhs( spk, phs, periods, varargin )

%--------------------------------------------------------------------%
% check inputs
%--------------------------------------------------------------------%

[phsFs, model, nreps, nbins, graphics] = ParseArgPairs (...
    { 'phsFs', 'model', 'nreps', 'nbins', 'graphics'} ...
    ,{1250, 'uniform',500, 50, 0},varargin {:} );

if nargin < 3 || isempty(spk) || isempty(phs) || isempty(periods)
    error('Not enough input arguments');
end

nspk                            = length(spk);

%--------------------------------------------------------------------%
% Generate x axis
%--------------------------------------------------------------------%

bsize                           = (2*pi) / ( nbins -1 );
edges                           = (-bsize/2) : bsize : (2*pi + bsize/2);
binc                            = 0:bsize:(2*pi);


%--------------------------------------------------------------------%
% Generate distribution and draw phases
%--------------------------------------------------------------------%

% generate cdf
switch model
    
    case 'uniform'
        E                       = 0:bsize:(2*pi);
        B                       = E(1:end-1) + bsize/2;
        phs_pdf                 = pdf('Uniform',B,0,2*pi);
        nbins                   = nbins - 1;
        
    case 'normal'
                E                       = 0:bsize:(2*pi);
        B                       = E(1:end-1) + bsize/2;
        phs_pdf                 = pdf('normal',B, pi, pi/4);
        
    case 'LFP'
        field_phs               = phs(local_enumerate(periods));
        phs_pdf                 = histcounts(field_phs,edges,'Normalization','pdf');
        
    case 'spk'
        spk_phs                 = phs(spk);
        phs_pdf                 = histcounts(spk_phs,edges,'Normalization','pdf');
        
    case 'spk rate'
        field_phs               = phs(local_enumerate(periods));
        spk_phs                 = phs(spk);
        LFP_count               = histcounts(field_phs,edges);
        LFP_sec                 = LFP_count/phsFs;
        SPK_count               = histcounts(spk_phs,edges);
        SPK_rate                = SPK_count ./ LFP_sec;
        phs_pdf                 = SPK_rate / sum(SPK_rate);
    otherwise 
        error('mode not supported');
end

% calculate cdf
phs_pdf                         = phs_pdf( : ).';
phs_cdf                         = [0 cumsum( phs_pdf ) / sum( phs_pdf )] ;

% generate n random numbers, between 0 and 1
ndraw                           = nreps * nspk;
tidx                            = rand( ndraw, 1 );
CDF                             = phs_cdf( ones( ndraw, 1 ), : );

% sample the distribution (create the indices)
idx                             = tidx( :, ones( 1, nbins +1) ) > CDF;
rand_idx                        = sum( idx, 2 );

% get the phase correspondign to these indices
rand_phs                        = binc(rand_idx);
rand_phs                        = reshape(rand_phs,nspk,nreps);

%--------------------------------------------------------------------%
% Generate random cphases
%--------------------------------------------------------------------%

% compute phase cycle
[cyc_range, cyc]                = calc_cycle(phs);
spk_cyc                         = cyc(spk);
spk_phs                         = phs(spk);
spk_cyc_per                     = cyc_range(spk_cyc,:);

spk_DC                          = (spk_cyc - 1) * 2 * pi;
spk_cphs                        = spk_phs + spk_DC; 
rand_cphs                       = rand_phs + repmat(spk_DC,1,nreps);
rand_cyc                        = spk_cyc;


%--------------------------------------------------------------------%
% plot
%--------------------------------------------------------------------%

if ~graphics
    return;
end

cphs                            = phs + (cyc -1) * 2 * pi;
c                               = 100:113;

figure;
subplot(3,1,1);
plot(spk_cphs(c),1: length(c), '.r'); hold on;
plot(rand_cphs(c,1),1: length(c), '.b');
line( [cphs(spk_cyc_per(c,1)), cphs(spk_cyc_per(c,1))], ylim, 'color', 'k' )

ylabel('Spike number');
xlabel('Spike cphs');
set( gca, 'tickdir', 'out', 'box', 'off' )
legend({'Orig','rand','cycStart'});
title('Example spikes (one shuffled repetition)');


subplot(3,1,2);
plot(rand_cphs(c,:),1: length(c), '.b'); hold on;
plot(spk_cphs(c),1: length(c), '.r'); 
line( [cphs(spk_cyc_per(c,1)), cphs(spk_cyc_per(c,1))], ylim, 'color', 'k' )

ylabel('Spike number');
xlabel('Spike cphs');
set( gca, 'tickdir', 'out', 'box', 'off' )
legend({'Orig','rand','cycStart'});
title('Example spikes (all shuffled repetition)');

subplot(3,1,3);
plot(rand_cphs,1:nspk, '.b'); hold on;
plot(spk_cphs, 1:nspk, '.r'); 
line( [cphs(spk_cyc_per(c,1)), cphs(spk_cyc_per(c,1))], ylim, 'color', 'k' )
ylabel('Spike number');
xlabel('Spike cphs');
set( gca, 'tickdir', 'out', 'box', 'off' )
legend({'Orig','rand','cycStart'});
title('Example spikes (one shuffled repetition)');

figure;

edges                           = (-bsize/2) : bsize : (2*pi + bsize/2);
binc                            = 0:bsize:(2*pi);

% first row: before shuffle
subplot(2,4,1);
counts              = histcounts(spk_cyc,'BinWidth',1);
counts              = counts(counts~=0);
CYCedges            = 0.5 : 1 : 15.5;  
counts              = histcounts(counts,CYCedges);
histogram('binEdges', CYCedges, 'binCounts', counts);
xlim( [ min(CYCedges), max(CYCedges)] );

set( gca, 'tickdir', 'out', 'box', 'off' );
xlabel('Number of spikes');
ylabel('Number of cycles');
title('Spikes per cycle');
xlim( [ min(CYCedges), max(CYCedges)] );

subplot(2,4,2);
field_phs  = phs(local_enumerate(periods));
phi        = local_circ_mean( field_phs);
LFP_count  = histcounts(field_phs,edges);
LFP_sec    = LFP_count/phsFs;
ee         =  LFP_sec(1) + LFP_sec(end);
LFP_sec(1) = ee;
LFP_sec(end) = ee;
histogram('binEdges', edges, 'binCounts', LFP_sec);
xticks(0:(pi/2):(2*pi) );
xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
ylabel('Seconds spent');
set( gca, 'tickdir', 'out', 'box', 'off' );
line( [phi, phi], ylim, 'color', 'r','LineStyle','--' )

xlabel('Theta phase');
xlim( [ min(edges), max(edges)] );

title(sprintf('within-field LFP phase (phi=%.2f)',phi));

subplot(2,4,3);
SPK_count         = histcounts(spk_phs,edges);
ee                =  SPK_count(1) + SPK_count(end);
SPK_count(1) = ee;
SPK_count(end) = ee;
phi                 = local_circ_mean( spk_phs);
histogram('binEdges', edges, 'binCounts', SPK_count);
ylabel('Number of spikes');
set( gca, 'tickdir', 'out', 'box', 'off' )
xticks(0:(pi/2):(2*pi) );
xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
line( [phi, phi], ylim, 'color', 'r','LineStyle','--' )
xlabel('Theta phase');
xlim( [ min(edges), max(edges)] );
title(sprintf('Spike phase (phi=%.2f)',phi));

subplot(2,4,4);
SPK_rate = SPK_count ./ LFP_sec;
histogram('binEdges', edges, 'binCounts', SPK_rate);
phi  = local_circ_mean(binc(:), SPK_rate(:));
line( [phi, phi], ylim, 'color', 'r','LineStyle','--' )
xticks(0:(pi/2):(2*pi) );
xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel('Firing rate [spk/s]');
xlabel('Theta phase');
xlim( [ min(edges), max(edges)] );
title(sprintf('Spike rate (phi=%.2f)',phi));

% second row: after shuffle
% number of spikes per cycle
subplot(2,4,5);
rand_DC                         = rand_cphs - rand_phs;
rand_cyc                        = round( (  (rand_DC) / (2* pi) ) +1 );
counts                          = histcounts(rand_cyc(:),'BinWidth',1);
counts                          = counts(counts~=0);
counts                          = counts/nreps;
CYCedges                        = 0.5 : 1 : 15.5;  
counts                          = histcounts(counts,CYCedges);
histogram('binEdges', CYCedges, 'binCounts', counts);
xlim( [ min(CYCedges), max(CYCedges)] );
set( gca, 'tickdir', 'out', 'box', 'off' );
xlabel('Number of spikes');
ylabel('Number of cycles (per repetition)');
title('Spikes per cycle');

% LFP
subplot(2,4,6);
phi        = local_circ_mean( field_phs);
histogram('binEdges', edges, 'binCounts', LFP_sec);
xticks(0:(pi/2):(2*pi) );
xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
ylabel('Seconds spent');
set( gca, 'tickdir', 'out', 'box', 'off' );
line( [phi, phi], ylim, 'color', 'r','LineStyle','--' )
xlabel('Theta phase');
xlim( [ min(edges), max(edges)] );
title(sprintf('within-field LFP phase (phi=%.2f)',phi));

% spike 
subplot(2,4,7);
SPK_count         = histcounts(rand_phs(:),edges);
ee                =  SPK_count(1) + SPK_count(end);
SPK_count(1) = ee;
SPK_count(end) = ee;
phi                 = local_circ_mean( rand_phs(:));
histogram('binEdges', edges, 'binCounts', SPK_count);
ylabel('Number of spikes');
set( gca, 'tickdir', 'out', 'box', 'off' )
xticks(0:(pi/2):(2*pi) );
xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
line( [phi, phi], ylim, 'color', 'r','LineStyle','--' )
xlabel('Theta phase');
xlim( [ min(edges), max(edges)] );
title(sprintf('Spike phase (phi=%.2f)',phi));

subplot(2,4,8);
SPK_rate = SPK_count ./ LFP_sec;
histogram('binEdges', edges, 'binCounts', SPK_rate);
phi  = local_circ_mean(binc(:), SPK_rate(:));
line( [phi, phi], ylim, 'color', 'r','LineStyle','--' )
xticks(0:(pi/2):(2*pi) );
xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel('Firing rate [spk/s]');
xlabel('Theta phase');
xlim( [ min(edges), max(edges)] );
title(sprintf('Spike rate (phi=%.2f)',phi));

sgtitle(sprintf('Within-field phase rand: %s',model));

return;



%------------------------------------------------------------------------
% phi = local_circ_mean( t )
% compute mean direction
%------------------------------------------------------------------------
function phi = local_circ_mean( t, f )

% trigonometric functions
if nargin < 2 || isempty( f )
    nans                            = isnan( t );
    n                               = sum( ~nans, 1 );
    x                               = cos( t );
    y                               = sin( t );
else
    if size( t, 2 ) == 1
        t               = t * ones( 1, size( f, 2 ) );
    elseif size( f, 2 ) == 1
        f               = f * ones( 1, size( t, 2 ) );
    end
    if ~isequal( size( f ), size( t ) )
        error( 'input size mismatch' )
    end
    n                   = nansum( f );
    x                   = f .* cos( t );
    y                   = f .* sin( t );
end

% compute direction
sumx                            = nansum( x, 1 );
sumy                            = nansum( y, 1 );
C                               = sumx ./ n;
S                               = sumy ./ n;
phi                             = mod( atan2( S, C ), 2 * pi );

return % local_circ_mean

%------------------------------------------------------------------------
% vec = local_enumerate( mat )
% the inverse of parse
%------------------------------------------------------------------------


function vec = local_enumerate( mat )

% initialize output
vec                         = [];

% arguments
if isempty( mat )
    return
end

% code
mat                                 = sortranges( mat );
m                                   = size( mat, 1 );
durs                                = diff( mat, [], 2 ) + 1;
cdurs                               = cumsum( durs );
vec                                 = zeros( cdurs( m ), 1 );
idx                                 = [ 1; cdurs( 1 : m - 1 ) + 1 ];
for i                               = 1 : m
    vec( idx( i ) : cdurs( i ) )    = mat( i, 1 ) : mat( i, 2 );
end

return

