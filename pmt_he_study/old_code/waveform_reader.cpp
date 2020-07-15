#include "TF1.h"
#include "TStyle.h"
#include "TObjString.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/* To compile:
 root [] gSystem->CompileMacro("waveform_reader.cpp")
 root [] waveform_reader()
 
 NOTE: It is much quicker to compile the code for one channel at a time for large files!!! (Otherwise, the code takes >1hr to create the files).
 */

void waveform_reader()
{
    
    TString pathFileName("/Users/adilaislam/helium_pmt/data/A1000_B1000_t1346_long.xml");
    //cout << "Reading waveform from file : " << pathFileName << endl;
    string LINE;
    ifstream infile;
    infile.open(pathFileName);
    
    
    // Make an output file to store the waveforms:
    /*ofstream fl_Ch0("../archive/22.2.19/long/channel_0/22_02_Ch0_waveforms.dat");
    
    if (!fl_Ch0)
    {
        cout << "file could not be open for writing ! " <<endl;
        
    }*/
    
    ofstream fl_Ch1("../archive/22.2.19/long/channel_1/22_02_Ch1_waveforms.dat");
    
    if (!fl_Ch1)
    {
        cout << "file could not be open for writing ! " <<endl;
        
    }

    int skip_header_lines = 45;
    
    for (int iSkip=0; iSkip<skip_header_lines; ++iSkip)
    {
        if (!infile.eof()) getline(infile,LINE);
        
    }
    
    int trace_counter;
    
    TString separator = " ";
    
    while(!infile.eof())
    {
        getline(infile,LINE);
        // std::cout << LINE << std::endl;
        
        TString LineString(LINE);
        TObjArray* words = LineString.Tokenize(separator.Data());
        
        cout << "Number of words = " << words->GetEntries() << endl;
        
        if (words->GetEntries() > 10)
        {
            // This looks like a line containing a trace
            
            // Extract channel number :
            TObjString channelWord(*((TObjString*)((*words)[1])));
            
            /*if (channelWord.GetString().Contains("=\"0\""))
            {
                cout << "CHANNEL 0 !!" << endl;
                
            }
            
           if (channelWord.GetString().Contains("=\"1\""))
           {
                cout << "CHANNEL 1 !!" << endl;
               
           }*/
            
            for (int iword=2; iword < words->GetEntries(); iword++)
            {
                TObjString thisWord(*((TObjString*)((*words)[iword])));
                //cout << iword-2 << " " << thisWord.String() << endl;
                
                //write into files
                /*if (channelWord.GetString().Contains("=\"0\""))
                {
                    fl_Ch0 << iword-2 <<" "<< thisWord.String() << endl;
                    
                }*/
                
                if (channelWord.GetString().Contains("=\"1\""))
                {
                    fl_Ch1 << iword-2 <<" "<< thisWord.String() << endl;
                    
                }
            }
            
            trace_counter++;
            

        }
        
    }
    
    
    //fl_Ch0.close();
    fl_Ch1.close();

    cout << " -- " << endl;
    cout << " -- " << endl;
    
    cout << "-- Files created --" << endl;
    
    cout << " -- " << endl;
    cout << " -- " << endl;
    
    cout << "-- Code Terminated --" << endl;

    cout << " -- " << endl;
    cout << " -- " << endl;
    

}



