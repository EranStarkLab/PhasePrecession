**phasePrecession** 

This repository contains MATLAB routines used to calculate temporal phase precession by the spike phase spectrum method, described in Sloin et al., 2023, bioRxiv.
## **Overview**
The spike phase temporal precession algorithm described by Sloin et al was designed to statistically detect and describe temporal precession. The current method allows to reduce false positives created by phase lock and false negatives created by units exhibiting both phase lock.  The core algorithm is implemented by the MATLAB routine spk\_phs\_spec.m

All the routines necessary for spk\_phs\_spec.m are included in the functions folder. In addition, the MATLAB deep learning toolbox is necessary for running the routines.
### **Analysis**
- spk\_phs\_spec.m
  - computes spike phase spectrum and indicate the occurrence of temporal precession
- calc\_cycle.m
  - generates a vector of theta cycle from a theta phase vector
- RandCyclePhs.m
  - randomize the phase of spikes 
- spike\_spectra.m
  - wrapper for calculating the spectrum of spike phase
### **Utilities**
- ParseArgPairs
  - flexible argument assigning
- inranges
  - determine which elements of a vector are in which range
- myjet
  - modified jet with extreme values in pure R,B
- sortranges
  - to be a set of non-overlapping [ small large ] pairs
- mixmat
  - mix matrix elements
- my\_spectrum
  - Welch spectrum for multiple signals. 
- resampleranges
  - resample ranges from one Fs to another, while keeping the total duration fixed
- resort
  - indices to recover original order

**Data**

- 1\_precession\_only.mat
  - spk: Spike time of an example unit exhibiting theta phase precession (same as Sloin et al., 2023, Fig. S9A). 
  - phs: LFP theta theta phase in radians, sampled at 1250 Hz.
  - periods: Start and end time of every crossing of the unit’s place field.
  - All three are sampled at 1250 Hz.
- 2\_lock only.mat
  - Spike time, phs and periods of a unit exhibiting only theta phase lock (same as Sloin et al., 2023, Fig. S9B). 
- 3\_precession\_lock.mat
  - Spike time, phs and periods of a unit exhibiting both theta phase precession and theta phase lock (same as Sloin et al., 2023, Fig. S9C). 
## **Demo** 
The wrapper spk\_phs\_spec\_demo.m demonstrates the calculation of spike phase spectra by running it on three example units. One example unit exhibits phase locking, the second exhibits phase precession, and the third exhibits a combination of both. The examples are the same as the ones described in Sloin et al., 2023, bioRxiv, Fig. S2.

To demonstrate the algorithm for spike phase temporal precession, run spk\_phs\_spec\_demo.m. To run spk\_phs\_spec\_demo.m, you will need the following: 
### To run the demo
- Download all routines and data
- In MATLAB, write spk\_phs\_spec\_demo
### Demo results

