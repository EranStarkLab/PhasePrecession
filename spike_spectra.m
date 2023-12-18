% spike_spectra        calculate the spectra of spikes
%
% CALL                 [ p1, fp, peak, pvals, sig_f ] = spike_spectra( spk, periods, varargin )
%
%
% GETS                 spk                 time of spikes
%
% OPTIONAL ARGUMENTS
%                      periods  {[]}         take spikes only within these periods [spkFs]
%                      spkFs    {20000}
%                      f0       {1}          minimal frequncy
%                      M        {10}         binSize factor
%                      DSF      {8}          downsampling factor for spikes
%                      dflag    {'constant'} normalization of spectra, passed to local_mtcsd1
%                      fROI     {[5 20]}     frequencies of interest
%                      nreps    {10}         for the calculation of pval
%                      alpha    {0.05}       significance level
%                      graphics {1}
%                      method   {'mt'}       fft method, can also be 'welch'
%
% RETURNS
%                      pow                   power at each frequency within ROI
%                      fp                    frequency vector
%                      peak                  peak of spectra [fpeak, apeak, gpeak, hfp]
%                                               - fpeak: frequncy at the peak
%                                               - apeak: amplitude at the peak
%                                               - gpeak: gain at the peak
%                                               - hfp: mean amplitude in high frequencies (used to calculate gpeak)
%                                               - pval: pvalue for the peak (bonferroni corrected)
%                      pvals_ROI             pvalue of spectra amplitude in each frequency according to fp
%                      sig_f                 frequencies with significant amplitude
%
% CALLS                resampleranges, inranges

% written by           ES and HES  04-May-21
% modified             HES       18-Dec-23


function [ pow, fp, peak, pvals, sig_f, nspk ] = spike_spectra( spk, periods, varargin )

% constants
NW                              = 3;                                        % MT: spike spectrum time-freq parameter
ZVAL                            = sqrt( eps );

%--------------------------------------------------------------------%
% (1) gather inputs
%--------------------------------------------------------------------%
[spkFs, f0, DSF, ...
    dflag, graphics, fROI, M, ...
    nreps, alpha, calcsig, method ] = ParseArgPairs (...
    {'spkFs', 'f0', 'DSF', 'dflag', 'graphics', 'fROI', 'M', ...
    'nreps', 'alpha', 'calcsig' , 'method'} ...
    , {20000, [], 8 , 1, 0, [ 5 20 ], 20, ...
    10, 0.05, 1, 'mt'} , varargin{ : } );


%--------------------------------------------------------------------%
% (2) preparations
%--------------------------------------------------------------------%

% downsample spk to save memory
Fs                              = round(spkFs / DSF);
spk                             = round( spk / DSF );
periods                         = resampleranges( periods, 1, DSF );


% compute some parameters
fHiROI                          = [ Fs / 2.5 Fs / 2 ];                      % in principle, may coincide with fROI
nFFTo                           = 2^floor( log2( M * Fs ) );
nFFT                            = nFFTo + (Fs - rem(nFFTo,Fs));
nWindow                         = nFFT;
fp0                             = ( 0 : Fs / nFFT : Fs / 2 )';
if isempty( f0 )
    f0                          = fp0( 2 );
end

%--------------------------------------------------------------------%
% (3) convert spike train to a continuous vector
%--------------------------------------------------------------------%

% initialize memory
npad                            = 1 / f0 * Fs;
dperiods                        = diff( periods, [], 2 ) + 1;
ninperiods                      = sum( dperiods );
nperiods                        = size( periods, 1 );
vec                             = zeros( ninperiods + npad * ( nperiods + 1 ), 1 );

% get the spikes in the periods
[ ~, seg, out ]                 = inranges( spk, periods );
cperiods                        = cumsum( [ 0; dperiods ] );

% initialize
pow0                            = NaN(length(fp0),1);
pvals0                          = NaN(length(fp0),1);
peak                            = NaN(5,1);
sig_f                           = NaN(length(fp0),1);
fidx                            = ( fp0 >= fROI(1) ) & ( fp0 <= fROI(2) );

leng                            = sum(fidx);
pow                             = NaN(leng,1);
fp                              = NaN(leng,1);
pvals                           = NaN(leng,1);

if isempty(seg)
    return;
end

% populate the output vector
for i                           = 1 : nperiods
    iidx                        = out( seg == i ) - periods( i, 1 ) + 1;
    sidx                        = iidx + npad * i + cperiods( i );
    imean                       = length( sidx ) / dperiods( i ); % trial-specific correction for mean
    vec( sidx )                 = 1 - imean;
end

%--------------------------------------------------------------------%
% (4) calculate spectra
%--------------------------------------------------------------------%

% compute the spectrum
switch method 
    case 'mt'
     cxx                        = local_mtcsd1( vec, nFFT, Fs, nWindow, nWindow / 2, NW, dflag );
     pow0                         = cxx( :, 1, 1 );
    case 'welch'
      pow0                        = my_spectrum( vec, nFFT, Fs, nFFT, nFFT / 2, 'none' );
end

% find hfp amplitude for normalization
hidx                           = (fp0 >= fHiROI(1) ) & (fp0 <= fHiROI(2) );
hfp                             = nanmean( pow0( hidx ) );
peak(4)                         = hfp;

% peak detection: find a global maximum in fROI
ffidx                           = find( fidx );
pow                             = pow0( ffidx );
fp                              = fp0( ffidx );
[apeak, gidx]                   = max(pow);
fpeakIdx                        = ffidx(gidx);
fpeak                           = fp0(fpeakIdx);

% make sure it is also a local maximum in the full spec
pidx                            = local_find_local_max( pow0 );
fpeak_local                     = fp0( pidx );

if ~any(fpeak_local == fpeak)
    fpeak                       = [];
end

if isempty(fpeak)
    return;
end
 
    
% scale the peak according to the HF spectrum
gpeak                           = apeak / hfp;
peak                            = [ fpeak apeak gpeak hfp ];

if any(isnan( hfp ) )
    keyboard;
end
   

%--------------------------------------------------------------------%
% (5) randomize for significance
%--------------------------------------------------------------------%

peak_pval                           = NaN;

if calcsig
    % (5.1) randomize the spike train N times (as in step 3)
    vec_hat                       	= zeros( ninperiods + npad * ( nperiods + 1 ), nreps );

    % populate the output vector
    for i                           = 1 : nperiods
        nspk                        = sum( seg == i );
        rng                         = [ 1 diff( periods( i, : ) ) + 1 ];
        for j                       = 1 : nreps
            iidx                 	= sort( randperm( rng( 2 ), nspk ) )';
            %iidx                        = out( seg == i ) - periods( i, 1 ) + 1;
            sidx                	= iidx + npad * i + cperiods( i );
            imean                	= length( sidx ) / dperiods( i ); % trial-specific correction for mean
            vec_hat( sidx, j )   	= 1 - imean;
        end
    end

    % (5.2) compute spectrum for each randomized train (as in step 4 )
    p1_hat                          = NaN( nFFT / 2 + 1, nreps );
    for j                           = 1 : nreps
        cxx_hat                  	= local_mtcsd1( vec_hat( :, j ), nFFT, Fs, nWindow, nWindow / 2, NW, dflag );
        p1_hat( :, j )          	= cxx_hat( :, 1, 1 );
    end

    % (3) compute point-wise p-values for each element in original train (based on a Gaussian assumption/empirical)
    myu                             = nanmean( p1_hat, 2 );
    sig                             = nanstd( p1_hat, [], 2 );
    if nreps < 1000
        pvals0                      = 1 - normcdf( pow0, myu, sig );
        zidx                        = pvals0 < ZVAL;
        pvals0( zidx )              = ZVAL;
    else
        pvals0                      = ( sum( repmat( pow0, [ 1 nreps ] ) <= p1_hat, 2 ) + 1 ) / ( nreps + 1 );
    end

    % (4) determine which frequency bins are low-probability
    TH                              = alpha / sum( fidx ); % bonferroni--corrected alfa level
    fsidx                           = pvals0 <= TH & fidx;
    sig_f                           = fp0( fsidx );

    peak_pval                       = pvals0(fpeakIdx) * sum(fidx);
end

pvals                               = pvals0(ffidx);
peak(1,5)                           = peak_pval;

%-------------------------------------------------------------------%
% plot
%--------------------------------------------------------------------%

if graphics
    
    figure
    plot( fp0, pow0, 'b')

    if calcsig
        hold on;    
        plot( fp0, myu, 'r', fp0, myu + 2 * sig, '--r', fp0, myu - 2 * sig, '--r' )
        hold on, plot( sig_f, pow0( fsidx ), '.r' )   
        plot( fpeak, apeak, 'or' )
    end    
    xlim( fROI )
    set( gca, 'tickdir', 'out', 'box', 'off', 'FontSize', 12 );
    xlabel( 'Frequency [Hz]' );
    ylabel( 'Power [AU]' );
    title(sprintf('fpeak:%.2f Hz; gpeak:%.2f, pval:%.3f',fpeak, gpeak, peak_pval ) );
    line( xlim, [hfp, hfp], 'color', [ 0 0 0 ],'LineStyle','--' )
    line( [fpeak, fpeak], ylim, 'color', [ 0 0 0 ],'LineStyle','--' )
    
end

return

% EOF

%------------------------------------------------------------------------
% [ y, f ] = local_mtcsd1( x, nFFT, Fs, WinLength, nOverlap, NW, Detrend )
% cross-spectra between one signal (first column) and all others
%------------------------------------------------------------------------
function [ y, f ] = local_mtcsd1( x, nFFT, Fs, WinLength, nOverlap, NW, Detrend )

nTapers                         = 2 * NW - 1;
winstep                         = WinLength - nOverlap;
nChannels                       = size( x, 2 );
nSamples                        = size( x, 1 );

% check for column vector input
if nSamples == 1
	x                           = x';
	nSamples                    = size( x, 1 );
	nChannels                   = 1;
end

% calculate number of FFTChunks per channel
nFFTChunks                      = 1 + floor( ( ( nSamples - WinLength ) / winstep ) );

% allocate memory
y                               = complex( zeros( nFFT, 2, nChannels ) );   % output array
Periodogram                     = complex( zeros( nFFT, nTapers, nChannels ) ); % intermediate FFTs
Temp1                           = complex( zeros( nFFT, nTapers ) );
Temp2                           = complex( zeros( nFFT, nTapers ) );
Temp3                           = complex( zeros( nFFT, nTapers ) );
eJ                              = complex( zeros( nFFT, 1 ) );

% calculate Slepian sequences
Tapers                          = dpss( WinLength, NW, nTapers, 'calc' );

% compute tapered periodogram with FFT 
TaperingArray                   = repmat( Tapers, [ 1 1 nChannels ] );
for j                           = 1 : nFFTChunks
	Segment                     = x( ( j - 1 ) * winstep + ( 1 : WinLength ), : );
	if ~isempty( Detrend )
		Segment                 = detrend( Segment, Detrend );
    end
	SegmentsArray               = permute(repmat(Segment, [ 1 1 nTapers ] ), [ 1 3 2 ] );
	TaperedSegments             = TaperingArray .* SegmentsArray;

	Periodogram( :, :, : )      = fft( TaperedSegments, nFFT );

	% products 
    Temp1                       = squeeze( Periodogram( :, :, 1 ) );        % Ch1
    for Ch2                     = 1 : nChannels 
        
        Temp2                   = squeeze( Periodogram( :, :, Ch2 ) );
        Temp2                   = conj( Temp2 );
        Temp3                   = Temp1 .* Temp2;
        eJ                      = sum( Temp3, 2 );
        y( :, 2, Ch2 )          = y( :, 2, Ch2 ) + eJ / ( nTapers * nFFTChunks ); % cross

        Temp3                   = squeeze( Periodogram( :, :, Ch2 ) ) .* Temp2;
        eJ                      = sum( Temp3, 2 );
        y( :, 1, Ch2 )          = y( :, 1, Ch2 ) + eJ / ( nTapers * nFFTChunks ); % auto
    end
    
end

% select subset of y, set up f array
if ~any( any( imag( x ) ) )    % x purely real
	if rem( nFFT, 2 )
		select                  = 1 : ( nFFT + 1 ) / 2;
	else
		select                  = 1 : nFFT / 2 + 1;
	end
	y                           = y( select, :, : );
else
	select                      = 1 : nFFT;
end
f                               = ( select - 1 )' * Fs / nFFT;

return % local_mtcsd1

%------------------------------------------------------------------------
% idx = local_find_local_max( x )
% detect all local maxima in a vector
%------------------------------------------------------------------------
function idx = local_find_local_max( x )

x                               = x( : );
d2                              = diff( sign( diff( x ) ) );
idx                             = find( d2 < -1 ) + 1;

return % local_find_local_max
