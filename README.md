# PMT-ShapeAnalysis
Calorimeter Commissioning Data PMT waveform shape analysis
NOTE: This is for use on CRD files. This work will be integrated into the SNFEE software. I will note when the very last commit is made at the end of the this file

AverageWaveforms.py

    This file takes the raw CRD files and creates a series of Average Waveforms for each OM in the file. The OM is specified as follows:

    Slot X    - (0 - 29 for the main wall)
    Channel Y - (0 - 12 for the main wall)

    You can define X and Y in the topology part of the main() function.

    For the X-Wall and G-Veto OMs, specify this in the main() function and set the topology to the maximal Slot and Channel numbers. The output plots will reflect the actual OM layout

    The High Threshold HT has a value of 1 or 0 - 1 means triggers both LT and HT
    The Low Threshold Only LTO has 1 or 0 aswell - 1 means trigger low Threshold

    This code maps the event rate within the data files for both the LTO and HT individually and for the sum. The Average waveforms will loop over a selected sample (defined in the function getAverageWaveforms).

    Once the Average waveforms are collected, the function checkWaveforms() is called. This creates a template (defined when calling createTemplate function within checkWaveforms) OM average waveform take from one of the OMs. It then performs an inner product between the template and the rest of the average waveforms. These numbers are recorded in the shape 2D histogram. With a good choice of template any poor SHAPE INDEX (SI) values may indicate a poor PMT and the average waveforms should be checked. Further investigation of issues can be found by looking at some of the recorded sample waveforms.

The other files are less relevant to my analysis. "shape_index.cpp" is a C++ version of AverageWaveforms.py. Spectra.cpp collects the areas of each spectrum which can be converted to charge (and energy) and creates a spectra for each OM

Mapping is redundant but is included for completeness. It maps the events recorded in the data file (as done in AverageWaveforms)

Pulse_Modelling.cpp is in it's infancy but will be the code to run when modelling the PMT pulses.
