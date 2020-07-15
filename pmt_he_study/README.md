# supernemopmts
Helium Permeation Project

This repository is for the collation of the data logging and analysis code for this project.

Project Description:

The SuperNEMO experiment is design to search for the BSM event called Neutrinoless Double Beta Decay 0nBB. This is a "would be" rare event of the already rare Double Beta Deacy event 2nBB. SuperNEMO is unique in design as it can create a full event mapping topology. It is designed to also have excellent energy resolution down to the MeV level.

The experiment is made up of several parts, two of which are the Tracker and Calorimeter. The Tracker contains a gas medium neccessary for measuring the trajectories of the produced electrons from XnBB in the source foil (and from background elsewhere). This gas is primarily (95%) Helium. The tracker is then surrounded by hundreds of optical modules (OMs), made up of scintillating blocks and PhotoMultiplier Tubes (PMTs), the latter being the focus of this study.

Helium is known to permeate through glass and ruin the vacuum conditions of a PMT. This poses a threat to the effective lifetime of the PMT. As Helium in used in high concentations and in large volume within SuperNEMO, the study of the effect of afterpulsing must be further investigated.

The aims of this study are as follows:
- Quantatatively measure the resolution of a test Hamamatsu 8" PMT under controlled conditions as a function of its exposure to Helium gas
- Measure the rate of afterpulsing as a function of exposure

What should be contained in this repository?:
- The Bi Spectrum analysis code Area_Spectrum.cpp
- Afterpulse signal processing code sweep_shape_index.cpp
- Any template waveforms that have collected

- Perhaps collate some of the results on here or save on the UCL.hep cluster.
