#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"

using namespace std;

void waveform_reader() {
	  
	ifstream resultsFile;
	resultsFile.open("../archive/March_1.4K_Ch1_waveforms.dat");
  
	ofstream currentOutFile;
  
	// Make an output file to store selected waveforms:
    TFile* histFile = new TFile("all_nonempty_waveforms.root","RECREATE");
  
	if(resultsFile.is_open())
    {
		string data;
		
		/* 
		#############################################################################
		
		This is a defining section
		
		#############################################################################
		*/
		
		// Number of waveforms to scan:
		int nWaveforms_Requested;
		// Number of Bins
		int nBins;
		int max_area = 2000;
		
		cout << " -- " << "The maximum area is set as "<< " -- " << max_area <<" -- " << endl;
		cout << " -- " << "Define the number of bins  "<< " -- ";
		cin >> nBins;
		cout << " -- " << "The bin Width is "<< " -- " << double(max_area/nBins) << " -- " << endl;
		cout << " -- -- -- -- -- " << endl;
		cout << " -- " << "Define how many waveforms to be analysed " << " -- ";
		cin >> nWaveforms_Requested;
		cout << " -- -- -- -- -- " << endl;
		
		// Find the full window of the waveform (at what point does the left column go back to 0?)
		// Counter to get number of total lines in file:
		int pointsCounter = 0;

		// Counter to fill trees and arrays
		int counter = 0;

		// Counter to count number of waveforms
		int waveform_counter = 0;

		// If the timebase doesn't start at 0 we want to find out what it does start at:
		Double_t timebase = 0.;
		
		// Define a counter for the number high energy and low energy samples 
		int high_energy_counter = 0;
		int low_energy_counter = 0;
		
		/*
		#############################################################################
		
		End of defining section
		
		#############################################################################
		*/

		// Open file to obtain the record length of the waveforms:       
		while (getline(resultsFile, data))
        {
     
			stringstream dataReadin;
			dataReadin << data;
      
			Double_t timeStamp = 0;
			Double_t ADC = 0;
      
			//dataReadin >> timeStamp >> ADC >> baseline >> trigger >> longGate >> shortGate >> zero;
			dataReadin >> timeStamp >> ADC;
			// If the first line starts at zero, then set that to be the start of the timebase. When the next timebase is reached (i.e. the start of another waveform) break out of the loop.
			if (pointsCounter == 0) {
				timebase = timeStamp;
			}
			else if (timebase == timeStamp) {
				break;
			}
      
			// Increment counter to find out how many points per waveform
			++pointsCounter;

			// currentOutFile << timeStamp*pow(10,-9) << " " << ADC*pow(10,-3) << std::endl;
      
		}
		
		/* 
		#############################################################################
		
		This is a defining section
		
		
		#############################################################################
		*/
		
		cout << "Number of points per waveform: " << pointsCounter << endl;
		
		// Define an integer to have the number of points per waveform to use to define arrays
		const unsigned int nLines = pointsCounter;
    
		// Go back to the beginning of the file to read the first waveform again
		resultsFile.clear();
		resultsFile.seekg(0, ios::beg);

		// Define graph/hist and its arrays
		double* timebase_value;
		timebase_value = new double[nLines];

		double* adc_value;
		adc_value = new double[nLines];

		// Define Histograms for areas:
		TH1D* total_area = new TH1D("total_area", "Area", nBins, 0,max_area);
		
		
		/*
		#############################################################################
		
		End of defining section
		
		#############################################################################
		*/
		
		// Re-read in file variables
	
		while (getline(resultsFile, data) && waveform_counter < nWaveforms_Requested) {
			
			stringstream dataReadin;
			dataReadin << data;
      
			Double_t timeStamp = 0;
			Double_t ADC = 0;
      
			dataReadin >> timeStamp >> ADC;
			
			// Fill arrays:
			timebase_value[counter] = timeStamp;
			adc_value[counter] = ADC;
			
			double pre_trigger = 590;
			
			
			if (counter == pointsCounter - 1)
            {
			
				++waveform_counter;

				double sampling_rate = 1.; //scope sampling rate in GHz
				
				// Define doubles for area calculations:
				double area = 0;
                double Total_Area_Actual = 0;
				int integration_window = 0;
				
				
				// Average the pre-trigger region to get the baseline:
				double baseline_average = 0.;
				double baseline_calculated = 0.;
				for (int j = 0; j < pre_trigger - 20.; j++) {
					baseline_average += adc_value[j];
				}
				baseline_calculated = baseline_average/(pre_trigger - 20.);
				//cout << "The baseline value is " << baseline_calculated << endl;
				
				for (int i = 0; i < nLines; i++) {
                    
					// hist->SetBinContent(timebase_value[i], adc_value[i]);
					// hist->SetBinError(timebase_value[i], 1.); // assume an error of +/- 1 mV on the amplitude
	
					// Calculate the area of the initial pulse
					if (timebase_value[i] >= pre_trigger && integration_window < 40){
						area = area + ((adc_value[i] - baseline_calculated)/50);
						integration_window++;
    
					}
				}
		
				//calculate the area of the inverted pulse
                Total_Area_Actual = abs(area);
                //cout << "Area is: " << Total_Area_Actual << endl;
                
				//fill the actual areas into a histogram
                if (Total_Area_Actual > 0)
                {
                
                    total_area->Fill(Total_Area_Actual);
                    
                    //set the title and label the axis
                    total_area->SetTitle("Full Bi-207 Spectrum");
                    total_area->GetXaxis()->SetTitle("charge/pC");
                    total_area->GetYaxis()->SetTitle("ADC-count/mV");
                
                }
				// Do I want to save this waveform for viewing ?
                //3000
                //500
				if (Total_Area_Actual > 10. && high_energy_counter < 10000) {

                    // Make a unique histogram name:
                    char histname[20];
                    sprintf(histname,"waveform_%d",high_energy_counter);

                    TH1F* hist = new TH1F(histname,histname, (nLines-1)/sampling_rate, 0, (nLines-1)/sampling_rate);
                    hist->SetTitle(histname);
                    hist->GetXaxis()->SetTitle("Timestamp/ns");
                    hist->GetYaxis()->SetTitle("ADC-count/mV");
                    
                    for (int i = 0; i < nLines; i++) {
                        hist->SetBinContent(timebase_value[i], adc_value[i]);
                        hist->SetBinError(timebase_value[i], 1.); // assume an error of +/- 1 mV on the amplitude
                    }
                    hist->Write();
					high_energy_counter++;
				}
//                else {
//                    delete hist;
//                }
//
				
				cout << " -- " << "waveform counter: "<<waveform_counter << " -- " << endl;
				//waveform_tree->Fill();
				counter = 0;
                
				
			}
			
			else
            {
				++counter;
			}
			
		} //close while loop to read in file
	
		
		cout << " -- " << endl;
		cout << " -- " << endl;
		
		cout << "-- Successful loop --" << endl;
		
		cout << " -- " << endl;
		cout << " -- " << endl;
		
		cout << "-- The number of waveforms in this file is : " << waveform_counter << " --"<< endl;
		
		total_area->Write();
		resultsFile.close();
		histFile->cd();
		histFile->Close();
		
	}//close if loop to see if file is open

	else
    {
    cout << " -- Error opening file -- ";
	}
	
	cout << " -- " << endl;
	cout << " -- " << endl;
	
	cout << "-- Code terminated --" << endl;
	
	cout << " -- " << endl;
	cout << " -- " << endl;
  
  
}

