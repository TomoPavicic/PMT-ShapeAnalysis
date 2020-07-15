#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <math.h>
#include <numeric>
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include <TRandom3.h>
#include <vector>
#include "TMultiGraph.h"


/**************************************************************
/
/    CHECK
/   ->APULSE REGION START
/   ->AREA RANGE FROM THE BI SPECTRUM
/   -> PRE-TRIGGER
/
 ***************************************************************/



using namespace std;

void template_reader()
{
    
    ifstream resultsFile;
    resultsFile.open("/Users/adilaislam/helium_pmt/1.4kV_2Ch_March/channel_1/archive/March_1.4K_Ch1_waveforms.dat");
    
    // Make an output file to store selected waveforms:
    TFile* histFile = new TFile("/Users/adilaislam/helium_pmt/1.4kV_2Ch_March/channel_1/ROOT_files/template_array_March1.4kV_CH1.root","RECREATE");
    
    // Make an output file to store array:
    ofstream fl("/Users/adilaislam/helium_pmt/1.4kV_2Ch_March/channel_1/template_files/template_array_March1.4kV_CH1.dat");
    
    if (!fl)
    {
        cout << "file could not be open for writing ! " <<endl;
        
    }
    
    if(resultsFile.is_open())
    {
        string data;
        
        /*
         #############################################################################
         
         This is a defining section
         
         #############################################################################
         */
        
        // Number of waveforms to scan:
        int nTemplate_Requested = 100;
        cout << "Number of waveforms templates to be averaged: " << nTemplate_Requested << endl;;
        
        // Afterpulsing region:
        int apulse_region_start = 1000;
        
        // Sweeping window width:
        const unsigned int pulse_width = 40.;
        
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
            
            dataReadin >> timeStamp >> ADC ;
            
            // If the first line starts at zero, then set that to be the start of the timebase. When the next timebase is reached (i.e. the start of another waveform) break out of the loop.
            if (pointsCounter == 0)
            {
                timebase = timeStamp;
            }
            else if (timebase == timeStamp)
            {
                break;
            }
            
            // Increment counter to find out how many points per waveform
            ++pointsCounter;
            
        }
        
        cout << " -- -- -- -- -- " << endl;
        cout << " -- First loop succesfully iterated -- " << endl;
        cout << " -- -- -- -- -- " << endl;
        cout << " -- Number of points per waveform: " << pointsCounter << endl;
        cout << " -- -- -- -- -- " << endl;
    
        /*
         #############################################################################
         
         This is a defining section
         
         
         #############################################################################
         */
        
        
        // Define an integer to have the number of points per waveform to use to define arrays
        const unsigned int nLines = pointsCounter;
        
        // Go back to the beginning of the file to read the first waveform again
        resultsFile.clear();
        resultsFile.seekg(0, ios::beg);
        
        // Define arrays for first read through
        double* timebase_value;
        timebase_value = new double[nLines];
        
        double* adc_value;
        adc_value = new double[nLines];
        
        double* vector_template;
        vector_template = new double[pulse_width];
        
        for (int n = 0; n < pulse_width; n++)
        {
            vector_template[n] = 0;
            //cout << "starting vector_template: " << vector_template[n] << endl;
        }
                
        double temp_num = 0;
        
        // Define Histograms for areas:
        TH1D* total_area = new TH1D("total_area", "Area", 1100, 0, 100);
        
        
        /*
         #############################################################################
         
         End of defining section
         
         #############################################################################
         */
        
        while (getline(resultsFile, data) && temp_num < nTemplate_Requested)
        {
            stringstream dataReadin;
            dataReadin << data;
            
            Double_t timeStamp = 0;
            Double_t ADC = 0;
            
            dataReadin >> timeStamp >> ADC ;
            
            timebase_value[counter] = timeStamp;
            adc_value[counter] = ADC;
            
            double charge_area = 0;
            double Total_Absolute_Area = 0;
            int integraton_window_iterator = 0;
            
            // Define the Pre-Trigger region:
            double pre_trigger = 590.0;

            
            if (counter == pointsCounter - 1)
            {
                // HERE WE HAVE A FILLED ARRAY. NOW ANALYSE THIS WAVEFORM.
                
                // Also find the pulse minimum within the first 700 samples
                double pulse_minimum = 1000000;
                int pulse_minimum_location = 0;
                for (int j = 0; j < 700; j++)
                {
                    if (adc_value[j] < pulse_minimum){
                        pulse_minimum = adc_value[j];
                        pulse_minimum_location = j;
                    }
                }
               // cout << "For this pulse, the minimum value is " << pulse_minimum << " at location " << pulse_minimum_location << endl;

                
                // NOW DEFINE WINDOW BASED ON WHERE THE MINIMUM IS
                // cout << pre_trigger << endl;
                pre_trigger = pre_trigger + (pulse_minimum_location-605.0);
                //cout << pre_trigger << endl;

                
                // Average the pre_trigger region to get the baseline.
                double baseline_average = 0.;
                double baseline_calculated = 0.;
                for (int j = 0; j < pre_trigger - 40.; j++)
                {
                    baseline_average += adc_value[j];
                }

                baseline_calculated = baseline_average/(pre_trigger - 40.);
                
                
                
                
                double charge_area = 0.;
                double Total_Absolute_Area = 0.;
                int integraton_window_iterator = 0.;
                
                
                // Make a unique histogram name:
                char histname[20];
                sprintf(histname,"waveform_%d",waveform_counter);
                
                ++waveform_counter;
                
                double sampling_rate = 1.; //scope sampling rate in GHz
                
                TH1F* hist = new TH1F(histname,histname, (nLines-1)/sampling_rate, 0, (nLines-1)/sampling_rate);
                hist->SetTitle(histname);
                hist->GetXaxis()->SetTitle("Timestamp/ns");
                hist->GetYaxis()->SetTitle("ADC-count/mV");
                
                
                for (int i = 0; i < nLines; i++)
                {
                    hist->SetBinContent(timebase_value[i], adc_value[i]);
                    hist->SetBinError(timebase_value[i], 1.); // assume an error of +/- 1 mV on the amplitude
                    
                    //calculate the area of the initial pulse
                    if (timebase_value[i] >= pre_trigger && integraton_window_iterator < 40)
                    {
                        charge_area = (charge_area + (adc_value[i] - baseline_calculated)/50.);
                        integraton_window_iterator++;
                    }
                }
                
                
                //calculate the area of the inverted pulse
                Total_Absolute_Area = abs(charge_area);
                //cout << "Total Absolute Area: " << Total_Absolute_Area << endl;
                
                //fill the actual areas into a histogram
                total_area->Fill(Total_Absolute_Area);
                
                //set the title and label the axis
                total_area->SetTitle("Area Histogram");
                total_area->GetXaxis()->SetTitle("charge/pC");
                total_area->GetYaxis()->SetTitle("ADC-count/mV");
                
                // Do I want to save this waveform for viewing ?
                
                if (Total_Absolute_Area > 2 && Total_Absolute_Area < 3 && high_energy_counter < 20)
                {
                    hist->Write();
                    high_energy_counter++;
                }

                else
                {
                    delete hist;
                }
                
                
                if (Total_Absolute_Area > 2 && Total_Absolute_Area < 3) //area range of the K(976keV) peak on the Bi spectrum
                {
                    for (int j = pre_trigger; j < (pre_trigger+pulse_width); j++)
                    {
                       int q = j - pre_trigger;
                        vector_template[q] = vector_template[q] + 10*(adc_value[j] - baseline_calculated)/(abs(charge_area));
                        //cout << " -- " << vector_template[q]<< endl;
                    }
                    
                    
                    temp_num++;
                    //cout << " -- -- -- -- -- " << endl;
                    cout << " -- Template Number " << temp_num << endl;
                }
                
                //cout << " -- " << "waveform counter: "<< waveform_counter << " -- " << endl;
                
                counter = 0;
                
            }
            
            else
            {
                ++counter;
                
            }
            
        }
        
        cout << " -- -- -- -- -- " << endl;
        cout << " -- Second loop successfully iterated" << endl;
        cout << " -- -- -- -- -- " << endl;
        
        total_area->Write();
        resultsFile.close();
        histFile->cd();
        histFile->Close();
        
        //create the array by averaging out the vector template over the number of waveform templates.
        for (int h = 0; h < pulse_width; h++)
        {
            vector_template[h] = vector_template[h]/(nTemplate_Requested);
            cout << " -- " << vector_template[h] << endl;
            fl<<vector_template[h]<<endl;
        }
        
        fl.close();
        
        cout << " -- -- -- -- -- " << endl;
        cout << " -- Array successfully created" << endl;
        cout << " -- -- -- -- -- " << endl;
        
    }//close if loop to see if file is open
    
    else
    {
        cout << "-- Error opening file --";
    }
    
    cout << " -- " << endl;
    cout << " -- " << endl;
    
    cout << "-- Code terminated --" << endl;
    
    cout << " -- " << endl;
    cout << " -- " << endl;
    
}
