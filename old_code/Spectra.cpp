#include "TF1.h"
#include "TStyle.h"
#include "TObjString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDatime.h"
#include <iostream>
#include <fstream>
#include <string>

void ROOT_file_Topology(std::vector<int> & topology, std::string root_filename);
void Read_waveforms(std::vector<int> & topology, std::string data_filename,std::string root_filename);
Double_t getTotalArea(TObjArray* ADC_values,int pre_PULSE_region,int integration_window);
double checkWidth(TObjArray* ADC_values, int pre_PULSE_region);
void getFilenames(std::string Filenames_file,std::vector<std::string> & filenames);

void Spectra_116()
{
    TDatime().Print();
    // Create vector to define the map of the wall
    std::vector<int> topology(2);
    topology[0] = 20; // Columns
    topology[1] = 13; // Rows

    std::vector<std::string> filenames;
    getFilenames("../data/filenames_116.txt",filenames);

    ROOT_file_Topology(topology,"../ROOT_files/Spectra_116.root");

    for (int filenum = 0; filenum < filenames.size(); filenum++)
    {
        std::cout << filenames.at(filenum) << std::endl;
        Read_waveforms(topology,filenames.at(filenum),"../ROOT_files/Spectra_116.root");
        TDatime().Print();
    } 
}

void getFilenames(std::string Filenames_file,std::vector<std::string> & filenames)
{

    TString pathFileName(Filenames_file.c_str());
    string LINE;
    ifstream infile;
    infile.open(pathFileName);
    TString separator = " ";
    if (!infile)
    {
        std::cout << "Unable to open filenames.txt" << std::endl;
    }
    while(!infile.eof())
    {
        getline(infile,LINE);
        //std::cout << LINE << std::endl;
        filenames.push_back("../data/" + LINE);
    }
    infile.close();
}

void Read_waveforms(std::vector<int> & topology, std::string data_filename,std::string root_filename)
{
    TFile* root_file = new TFile(root_filename.c_str(),"UPDATE");
    if (root_file == nullptr)
    {
        std::cout << "null pointer" << std::endl;
    }

    TString pathFileName(data_filename.c_str());
    string LINE;
    ifstream infile;
    infile.open(pathFileName);
    
    int skip_header_lines = 9;
    
    for (int iSkip=0; iSkip<skip_header_lines; ++iSkip)
    {
        if (!infile.eof()) getline(infile,LINE);
    }
    
    int line_number = 0;
    
    TString separator = " ";

    Int_t Slot;
    Int_t Channel;
    Int_t Signal_boolean;
    
    while(!infile.eof())
    {
        getline(infile,LINE);
        
        TString LineString(LINE);
        TObjArray* words = LineString.Tokenize(separator.Data());
        
        if (line_number ==1)
        {
            // Extract channel number :
            TObjString slotWord(*((TObjString*)((*words)[1])));
            TObjString channelWord(*((TObjString*)((*words)[3])));
            TObjString Signal_word(*((TObjString*)((*words)[7])));

            Slot = (int)slotWord.String().Atof();
            Channel = (int)channelWord.String().Atof();
            Signal_boolean = (int)Signal_word.String().Atof();

            if (Channel == 13){Signal_boolean=0;}

        }
        
        if (line_number == 2)
        {
            
            root_file->cd();
            line_number = -1;
            
            if (Signal_boolean != 0)
            {
                char Spectra_histname[15];
                sprintf(Spectra_histname,"Spectra_%d_%d",Slot,Channel);

                char Width_histname[15];
                sprintf(Width_histname,"Width_%d_%d",Slot,Channel);

                char Slot_Number[15];
                sprintf(Slot_Number,"Slot%d",Slot);

                char Channel_Number[15];
                sprintf(Channel_Number,"Channel%d",Channel);

                root_file->GetDirectory(Slot_Number)->GetDirectory(Channel_Number)->cd();

                TH1D* spectra_hist = (TH1D*)root_file->GetDirectory(Slot_Number)->GetDirectory(Channel_Number)->Get(Spectra_histname);
                TH1D* width_hist = (TH1D*)root_file->GetDirectory(Slot_Number)->GetDirectory(Channel_Number)->Get(Width_histname);
                
                Double_t Total_Integrated_Area = getTotalArea(words,160,400);

                double width = checkWidth(words,100);
                width_hist->Fill(width);
                width_hist->Write("",TObject::kOverwrite);
                delete width_hist;
                
                spectra_hist->Fill(Total_Integrated_Area);
                spectra_hist->Write("",TObject::kOverwrite);
                delete spectra_hist;
            }
           
        }
        delete words;
        line_number++;
    }
    root_file->Close();
}

void ROOT_file_Topology(std::vector<int> & topology, std::string root_filename)
{
    TFile* root_file = new TFile(root_filename.c_str(),"RECREATE");
    for (int slot_num = 0; slot_num < topology[0]; slot_num++)
    {
        root_file->cd();
        char Slot_Number[10];
        sprintf(Slot_Number,"Slot%d",slot_num);
        root_file->mkdir(Slot_Number);

        for (int channel_num = 0; channel_num < topology[1]; channel_num++)
        {
            char Channel_Number[10];
            sprintf(Channel_Number,"Channel%d",channel_num);
            root_file->GetDirectory(Slot_Number)->mkdir(Channel_Number);

            char Spectra_Number[15];
            sprintf(Spectra_Number,"Spectra_%d_%d",slot_num,channel_num);

            char Width_Number[15];
            sprintf(Width_Number,"Width_%d_%d",slot_num,channel_num);

            root_file->GetDirectory(Slot_Number)->GetDirectory(Channel_Number)->cd();
            TH1D* spectra_hist = new TH1D(Spectra_Number,Spectra_Number,100,0,100000);
            TH1D* width_hist = new TH1D(Width_Number,Width_Number,100,0,1);

            spectra_hist->Write();
            width_hist->Write();

            delete spectra_hist;
            delete width_hist;
        }
    }
    root_file->Close();
    std::cout << ">>> ROOT file created with the wall topology" << std::endl;
}

Double_t getTotalArea(TObjArray* ADC_values,int pre_PULSE_region,int integration_window)
{
    Double_t intermediate_baseline = 0.0;
    int pre_TRIGGER_region = pre_PULSE_region - 40;
    for (int i = 0; i < pre_TRIGGER_region; i++) 
    {
        TObjString thisWord(*((TObjString*)((*ADC_values)[i])));
        Double_t ADC_value = (Double_t)thisWord.String().Atof();
        intermediate_baseline += ADC_value;
    }
    Double_t ADC_baseline_calculated = intermediate_baseline/(Double_t)pre_TRIGGER_region;

    Double_t intermediate_integrated_area = 0.0;
    for (int i = pre_PULSE_region; i < pre_PULSE_region + integration_window; i++)
    {
        //std::cout << i << std::endl;
        TObjString thisWord(*((TObjString*)((*ADC_values)[i])));
        Double_t ADC_value = (Double_t)thisWord.String().Atof();
        intermediate_integrated_area += (Double_t)ADC_value - (Double_t)ADC_baseline_calculated;
    }
    Double_t Total_Integrated_Area = std::abs(intermediate_integrated_area);
    return Total_Integrated_Area;
}

double checkWidth(TObjArray* ADC_values, int pre_PULSE_region)
{
    Double_t intermediate_baseline = 0.0;
    int pre_TRIGGER_region = pre_PULSE_region - 40;
    
    int position;
    double amplitude = 2000.0;
    for (int i = 0; i < ADC_values->GetEntries()-1; i++)
    {
        TObjString thisWord(*((TObjString*)((*ADC_values)[i])));
        Double_t ADC_value = (Double_t)thisWord.String().Atof();
        if (i < pre_TRIGGER_region)
        {
            intermediate_baseline += (double)ADC_value;
        }
        if (amplitude > (double)ADC_value)
        {
            //std::cout << " amplitude " << amplitude << " - " << ADC_value <<  std::endl;
            amplitude = (double)ADC_value;
            position = i;
        }
    }
    Double_t ADC_baseline_calculated = intermediate_baseline/(Double_t)pre_TRIGGER_region;
    amplitude = amplitude - ADC_baseline_calculated;

    //std::cout << " position " << position << std::endl;

    int start,stop;

    for (int i = 0; i < ADC_values->GetEntries()-1; i++)
    {
        TObjString thisWord(*((TObjString*)((*ADC_values)[i])));
        Double_t ADC_value = (Double_t)thisWord.String().Atof();
        if (i < position && ADC_value - ADC_baseline_calculated > amplitude*0.5)
        {
            start = i;
            //std::cout << " start " << start << std::endl;
        }
        if (i > position && ADC_value - ADC_baseline_calculated < amplitude*0.5)
        {
            stop = i;
            //std::cout << " amplitude " << amplitude << " - " << ADC_value - ADC_baseline_calculated <<  std::endl;
            //std::cout << " stop " << stop << std::endl;
        }
    }
    double width = ((double)stop - (double)start)/(std::abs(amplitude));
    return width;
}

