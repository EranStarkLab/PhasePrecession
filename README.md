# phasePrecession
Routines used to calculate temporal phase precession with the spike phase spectrum method, described in Sloin et al., 2023, bioRxiv.

**Overview
**The spike phase temporal precession algorithm described by Sloin et al was designed to statistically detect and describe temporal precession. The current method allows to reduce false positives created by phase lock and false negatives created by units exhibiting both phase lock.  The core algorithm is implemented by the MATLAB routine spk_phs_spec.m
All the routines necessary for spk_phs_spec.m are included in the functions folder. In addition, the MATLAB deep learning toolbox is necessary for running the routines.
Analysis
•	spk_phs_spec.m
o	computes spike phase spectrum and indicate the occurrence of temporal precession
•	calc_cycle.m
o	generates a vector of theta cycle from a theta phase vector
•	RandCyclePhs.m
o	randomize the phase of spikes 
•	spike_spectra.m
o	wrapper for calculating the spectrum of spike phase
Utilities
•	ParseArgPairs
o	flexible argument assigning
•	inranges
o	determine which elements of a vector are in which range
•	myjet
o	modified jet with extreme values in pure R,B
•	sortranges
o	to be a set of non-overlapping [ small large ] pairs
•	mixmat
o	mix matrix elements
•	my_spectrum
o	Welch spectrum for multiple signals. 
•	resampleranges
o	resample ranges from one Fs to another, while keeping the total duration fixed
•	resort
o	indices to recover original order
Data
•	1_precession_only.mat
o	spk: Spike time of an example unit exhibiting theta phase precession (same as Sloin et al., 2023, Fig. S9A). 
o	phs: LFP theta theta phase in radians, sampled at 1250 Hz.
o	periods: Start and end time of every crossing of the unit’s place field.
o	All three are sampled at 1250 Hz.
•	2_lock only.mat
o	Spike time, phs and periods of a unit exhibiting only theta phase lock (same as Sloin et al., 2023, Fig. S9B). 
•	3_precession_lock.mat
o	Spike time, phs and periods of a unit exhibiting both theta phase precession and theta phase lock (same as Sloin et al., 2023, Fig. S9C). 
Demo 
The wrapper spk_phs_spec_demo.m demonstrates the calculation of spike phase spectra by running it on three example units. One example unit exhibits phase locking, the second exhibits phase precession, and the third exhibits a combination of both. The examples are the same as the ones described in Sloin et al., 2023, bioRxiv, Fig. S2.
To demonstrate the algorithm for spike phase temporal precession, run spk_phs_spec_demo.m. To run spk_phs_spec_demo.m, you will need the following: 
To run the demo
•	Download all routines and data
•	In MATLAB, write spk_phs_spec_demo
