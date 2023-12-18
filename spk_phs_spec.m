% spk_phs_spec         Calculation of temporal precession by the phase spectrum method
%
% CALL                [ spkPhsPeak, phsSpec, phsSpec_x, nspk_output ] = spk_phs_spec (spk, phs, periods)
%
%
% GETS                 spk                 Time of spikes, in phsFs samples
% 		               phs 	               Theta phase in rad, sampled at phsFs
%                      periods             The start and end timepoint of each trial, in phsFs samples
% OPTIONAL
%                      phsFs               {1250}
%                      cphs_Fs             {150} number of bins within each phsFs
%                      M                   {50}   FFT length, nFFT= 2^floor( log2( M * Fs ) )
%                      dflag               {1} detrend flag for FFT
%                      fROI                {[ 7 13 ]} relevant frequencies
%                      dp_fROI             {[0.1 4]} phase ratio ROI
%                      rand_model          {'spk'} if 'shuffle' run shuff_spk_cyc, otherwise passed to RandCyclePhs
%                      subsample           {1} subsample the original spike train to the mean number of randomized spikes
%                      rand_nreps          {300} randomizatio repetitions
%                      spec_sig            {1} calculate significane
%                      nonsig_nreps        {10} if no sig test is required, rand_nreps is set to this number
%                      alpha               {0.05}
%                      graphics            {1}
%
% DOES
%                      Examines whether the spike of a specific units
%                      exhibit temporal precession in relation to LFP. The
%                      steps of the calculations are:
%                         (1) Extract the comulative phase of each spike
%                         (2) Calculate the spectrum of spike phase relative to LFP. to test significance, use
%                         randomization by LFP phase, which can detect either locking, precession, or procession
%                         (3) Generate randomized spike cphs trains. New spike phases are drawn from one the following
%                             distributions, as defined by rand_model:
%                               - 'uniform': A uniform distribution between 0 and 2pi
%                               - 'normal': A normal distribution centered around pi
%                               - 'LFP': The distribution of LFP oscillations recorded during this session
%                               - 'spk': The distribution of phases as a function of spike counts of the particular unit
%                               - 'spk rate': The distribution of phases as a function of spike firing rate of the particular unit
%                         (4) subsample the original spike train to  account for the effect of spike decimation on spectrum magnitude
%                         (5) Calculate the spectrum of each randomized train
%                         (6) Subtract the mean spiking phase randomized spectrum from the original spectrum
%                         (7) Find the peak frequency of the difference between the original and randomized spectra, dfp
%                         (8) Determine whether dfp is bigger than 1, and the closest local maxima is also bigger than 1
%                         (9) Compare the magnitude of the original spectrum at dfp, dfp_mag, to the distribution of magnitudes obtained from the randomized spectra, dfp_mag_rand
%
%                       A cell is determined to exhibit temporal precession if:
%                           (a) The spike phase spectrum exhibits a significant peak
%                           (b) The peak frequency of the difference spectrum, dfp, is larger than 1
%                           (c) The closest local maxima in the original spectrum is also bigger than 1
%                           (d) The probability of observing dfp_mag (or a higher magnitude) under empirical randomization was smaller than alpha
%
%
% RETURNS              spkPhsPeak          Spike phase spectrum peak. columns:
%                                            (1) dfp, peak frequency of spike phase relative to LFP
%                                            (2) dfp_localMax, the closest local maxima of the original spectrum to thea peak frequency
%                                            (3) dfp_mag, magnitude of the peak frequency
%                                            (4) dfp_mag_rand, the peak magnitude in randomizations
%                                            (5) peak_pval, p-value for the original spike spectrum peak
%                                            (6) rand_pval, randomization p-value
%                                            (7) teTPP_flag, occurance of temporal precession
%                      phsSpec             Full phase spectrum
%                      phsSpec_x           Phase spectrum x axis
%                      nspk_output         number of spikes: [median in rand, used for spectrum, original]
%
% CALLS                ParseArgPairs, calc_cycle, RandCyclePhs,
%                      calc_spectra, spike_spectra, myjet
%
% written by           HES and ES      24-Mar-22
% modified             HES             18-Dec-23
%
% based on a method by Mizuseki et al., 2009
% also used and developed by Geisler et al., 2007, Eliav et al., 2018, Qasim et al., 2021


function [ spkPhsPeak, phsSpec, phsSpec_x, nspk_output ] = spk_phs_spec(spk, phs, periods, varargin)


if nargin < 3 || isempty( spk ) || isempty( phs ) || isempty( periods )
    error( 'missing arguments' )
end


[cphs_Fs, M, dflag, dp_fROI,...
    rand_nreps, graphics, rand_model, subsample, ...
     spec_sig, nonsig_nreps, alpha] = ParseArgPairs (...
    { 'cphs_Fs',  'M','dflag', 'dp_fROI',...
    'rand_nreps', 'graphics', 'rand_model', 'subsample', ...
    'spec_sig','nonsig_nreps',  'alpha'} ...
    ,{200, 50, 1, [0.66 1.5], ...
    300, 1, 'spk', 1, ...
    1, 10,  0.05},varargin {:} );

if ~spec_sig
    rand_nreps                = nonsig_nreps;
end

%--------------------------------------------------------------------%
% Extract the comulative phase of each spike
%--------------------------------------------------------------------%
[~, cyc]                      = calc_cycle(phs);
cphsDC                        = (cyc-1) * 2 * pi;
cphs                          = phs + cphsDC;
spk                           = spk(:);
nspk                          = length(spk);
ospk                          = nspk;
spk_cphs                      = cphs(spk);
spk_cphs                      = double(spk_cphs);
T                             = 2 * pi;
half_bin                      = T / cphs_Fs / 2;
spk_bcphs                     = ceil( ( spk_cphs - half_bin ) * cphs_Fs / T ) + 1;

%--------------------------------------------------------------------%
% Generate randomized spike cphs trains
%--------------------------------------------------------------------%


[rand_phs_spk, rand_cphs_spk]   = RandCyclePhs( spk, phs, periods, ...
    'model',rand_model, 'nreps', rand_nreps, 'graphics',1);

rand_cphs_spk                   = double(rand_cphs_spk);
rand_bcphs_spk                  = ceil( ( rand_cphs_spk - half_bin ) * cphs_Fs / T ) + 1;

% subsample the original spike train to  account for the effect of spike decimation on spectrum magnitude
if subsample

    % find out the unique number of spikes in each repetition
    uspk                            = zeros(rand_nreps,1);
    for i = 1 : rand_nreps
        uspk(i)                     = length( unique( rand_bcphs_spk(:,i)  ) );
    end
    muspk                           = floor( mean(uspk) );

    % prune spikes to the number following randomization
    kidx                            = randperm( nspk, muspk);
    spk                             = spk(kidx);
    spk_bcphs                       = spk_bcphs(kidx);
    nspk                            = muspk;
end

%--------------------------------------------------------------------%
% Calculate the spectrum of spike phase relative to LFP
%--------------------------------------------------------------------%

% get phase periods
cperiods                      = cphs(periods); % periods in cphase
dcperiods                     = ceil( ( cperiods - half_bin ) * cphs_Fs / T ) + 1;

% calculate phase spectrum
[ phsSpec, phsSpec_x, o_peak ]       = spike_spectra( spk_bcphs, dcperiods, 'spkFs', ...
    cphs_Fs,'DSF',1,'M',M,'dflag',dflag, 'fROI', dp_fROI,'calcsig',spec_sig, 'graphics',1);
peak_pval                  = o_peak(5);

% look only in range
ROI_idx                         = ( phsSpec_x >= dp_fROI(1) ) & ( phsSpec_x <= dp_fROI(2) );
rPspecX                         = phsSpec_x(ROI_idx);
rphsSpec                        = phsSpec( ROI_idx );
norm_phsSpec                    = phsSpec / o_peak(4);

% find global maximum
[~, p_idx]                       = max(rphsSpec);
fp                               = rPspecX(p_idx);
hfp                              = o_peak(4);

% find local maximum
max_idx                          = local_find_local_max( rphsSpec );
local_max_fp                     = rPspecX(max_idx);

%--------------------------------------------------------------------%
% statistical test (randomization)
%--------------------------------------------------------------------%

% initialize
Rspec                           = NaN( length(phsSpec_x),rand_nreps );
rand_mag                        = NaN(rand_nreps,1);
rand_hfp                        = NaN(rand_nreps,1);
rand_fp                         = NaN(rand_nreps,1);
uspk                            = zeros(rand_nreps,1);

% generate an basline distribution
for i = 1 : rand_nreps
    [ Rspec(:,i),~, r_peak]= spike_spectra( rand_bcphs_spk(:,i), dcperiods, 'spkFs', ...
        cphs_Fs,'DSF',1,'M',M,'dflag',dflag, 'fROI', dp_fROI, 'calcsig', 0);
    rRspec                      = Rspec(ROI_idx,i);
    rand_mag(i)                 = rRspec(p_idx); % exact
    rand_hfp(i)                 = r_peak(4);
    uspk(i)                     = length( unique( rand_bcphs_spk(:,i) ) );
    [~, max_idx]                = max( rRspec );
    rand_fp(i)                  = rPspecX(max_idx);
end

% Summary
rspk                            = floor( mean(uspk) );
rRspec                          = Rspec(ROI_idx,:); %
RspecM                          = mean(Rspec,2);
RspecSD                         = std(Rspec,[], 2);
mrand_fp                           = nanmedian(rand_fp);


% test for peak in fp, fp found based on substracted data
dphsspec                           = phsSpec - RspecM;
rdspec                             = dphsspec(ROI_idx);
[~, dspec_pidx]                    = max( rdspec );
dfp                                = rPspecX(dspec_pidx);

Dmag                               = rphsSpec(dspec_pidx);
Drand_mag                          = rRspec(dspec_pidx,:);

[~, LC_idx]                        = min( abs(local_max_fp - dfp) );
dfp_localMax                       = local_max_fp(LC_idx);
if isempty(dfp_localMax)
    dfp_localMax                   = NaN;
end

% test for normalized peak in fp, fp found based on substracted data
dfp_mag                            = Dmag / hfp;
Nrand_mag_norm                     = Drand_mag' ./ rand_hfp;
dfp_mag_rand                       = nanmedian(Nrand_mag_norm);
rand_pval                          = ( sum( dfp_mag < Nrand_mag_norm ) + 1 )  / ( rand_nreps + 1 );

%--------------------------------------------------------------------%
% summary
%--------------------------------------------------------------------%

sigpeak                             = rand_pval <= alpha;
sigrand                             = peak_pval <= alpha;
onefD                               = dfp > 1; % fq is bigger than 1
alllocMax                           = dfp_localMax > 1; % the closest local maxima is also bigger than 1
teTPP_flag                          = sigpeak & sigrand & onefD & alllocMax;

spkPhsPeak                          = [ dfp, dfp_localMax, dfp_mag, dfp_mag_rand, peak_pval, rand_pval, teTPP_flag];
nspk_output                         = [ rspk, nspk, ospk ]; % random, raw used, original


%--------------------------------------------------------------------%
% graphics
%--------------------------------------------------------------------%
if ~graphics
    return;
end

figure;

% spike phase randomization
nbins = 50;
bsize                           = (2*pi) / ( nbins -1 );
edges                           = (-bsize/2) : bsize : (2*pi + bsize/2);
binc                            = 0:bsize:(2*pi);
spk_phs           = phs(spk);
SPK_count         = histcounts(spk_phs,edges);
ee                =  SPK_count(1) + SPK_count(end);
SPK_count(1)      = ee;
SPK_count(end)    = ee;
SPK_count_norm    = SPK_count / sum(SPK_count);
phi               = local_circ_mean( spk_phs);

rSPK_count         = histcounts(rand_phs_spk(:),edges);
ee                 =  rSPK_count(1) + rSPK_count(end);
rSPK_count(1)      = ee;
rSPK_count(end)    = ee;
rSPK_count_norm    = rSPK_count / sum(rSPK_count);
phiR               = local_circ_mean( rand_phs_spk(:));

subplot(2,3,1);
a(1)              = plot(binc, SPK_count_norm,'Color','b');
line( [phi, phi], ylim, 'color', [ 0 0 1 ] )
hold on;
a(2)              = plot(binc, rSPK_count_norm,'Color','r','LineStyle','--');
line( [phiR, phiR], ylim, 'color', 'red' ,'LineStyle','--')
ylabel('Number of spikes');

set( gca, 'tickdir', 'out', 'box', 'off' )
xticks(0:(pi/2):(2*pi) );
xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
xlabel('Theta phase');
xlim([min(edges), max(edges)]);
title('Spiking phase distribution ');
legend(a, { sprintf('orig phi=%.2f',phi), sprintf('rand phi=%.2f', phiR) } );

% phase spectrum
subplot(2,3,2)
PspecX2     = log2(phsSpec_x);
a(1) = plot( PspecX2, phsSpec, 'b');
hold on;
a(2) = plot( PspecX2, RspecM, '--r');
plot( PspecX2, RspecM+RspecSD, '--r');
plot( PspecX2, RspecM-RspecSD, '--r');
set( gca, 'tickdir', 'out', 'box', 'off');
xlabel( 'Frequency relative to LFP theta' );
ylabel( 'Power [AU]' );
mticks  =  [min(phsSpec_x(ROI_idx)), 1, max(phsSpec_x(ROI_idx))];
mticks = round( mticks, 2);
xticks(log2(mticks));
xticklabels(mticks);
xlim( [ log2(mticks(1)), log2( mticks(end) ) ] );
hold on; 
line( [log2(fp), log2(fp)], ylim, 'color', 'b' ,'LineStyle','--')
line( [log2(1), log2(1)], ylim, 'color', 'k' )
title( sprintf( 'Phase spectrum (pval=%.2g)', peak_pval ) );
legend(a, {sprintf('fp=%.3g',fp) , sprintf('rand fp=%.3g', mrand_fp) });

% all randomized spectra
subplot(2,3,3);
imagesc(Rspec(ROI_idx,:)' );
xx = find( phsSpec_x(ROI_idx) == 1);
colormap(myjet);
xticks2  =  [1, xx, numel(phsSpec_x(ROI_idx) )];
xlabels2 = round( rPspecX(xticks2), 2);
set( gca, 'xtick', xticks2);
set( gca, 'XTickLabel',xlabels2 );
[~, oneidx] = min( abs( rPspecX - 1 ) );
line( [oneidx, oneidx], ylim, 'color', [1 1 1],'LineStyle','--' )
% alines( oneidx,'x','LineStyle','--','Color',[1 1 1]);
line( [find(rPspecX == dfp), find(rPspecX == dfp)], ylim, 'color', 'r','LineStyle','--' )

title( 'Randomized spectra');
xlabel('Frequency relative to LFP theta');
ylabel('Randomization Repetitions');
set( gca, 'tickdir', 'out', 'box', 'off');
colorbar;


% Difference between raw and randomization
subplot(2,3,4)
plot( PspecX2, dphsspec, 'b');
xlabel( 'Frequency relative to LFP theta' );
ylabel( 'Power [AU]' );
mticks  =  [min(phsSpec_x(ROI_idx)), 1, max(phsSpec_x(ROI_idx))];
mticks = round( mticks, 2);
xticks(log2(mticks));
xticklabels(mticks);
xlim( [ log2(mticks(1)), log2( mticks(end) ) ] );
hold on; 
line( [log2(dfp), log2(dfp)], ylim, 'color', 'b','LineStyle','--' )
% alines(log2(dfp),'x','LineStyle','--','Color','b');
line( [log2(1), log2(1)], ylim, 'color', 'k' )
% alines(log2(1), 'x', 'Color','k');
title( sprintf( 'Spectra difference: dfp=%.2f, closest local max=%.2f', dfp, dfp_localMax ) );
set( gca, 'tickdir', 'out', 'box', 'off');


% histogram of hfp norm peak following normalization
subplot(2,3,5)
hist(Nrand_mag_norm,rand_nreps/3);
set( gca, 'tickdir', 'out', 'box', 'off');
xlabel('Peak magnitude at dfp');
ylabel('Number of randomized repetitions');
title( sprintf( 'Peak magnitude at dfp: ratio=%.2f, p=%.2g', dfp_mag/dfp_mag_rand, rand_pval) );
% alines(dfp_mag,'x','LineStyle','--','Color','r');
line( [dfp_mag, dfp_mag], ylim, 'color', 'r','LineStyle','--' )
% alines(dfp_mag_rand,'x','LineStyle','--','Color','k');
line( [dfp_mag_rand, dfp_mag_rand], ylim, 'color', 'k','LineStyle','--' )


sgtitle( sprintf('Temporal TPP=%d, model=%s, nreps = %d', teTPP_flag, rand_model, rand_nreps) );


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
% idx = local_find_local_max( x )
% detect all local maxima in a vector
%------------------------------------------------------------------------
function idx = local_find_local_max( x )

x                               = x( : );
d2                              = diff( sign( diff( x ) ) );
idx                             = find( d2 < -1 ) + 1;

return % local_find_local_max