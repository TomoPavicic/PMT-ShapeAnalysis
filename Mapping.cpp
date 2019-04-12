#include "TF1.h"
#include "TStyle.h"
#include "TObjString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2I.h"
#include "TDatime.h"
#include <iostream>
#include <fstream>
#include <string>

void Read_waveforms(std::vector<int> & topology, std::string data_filename,std::string root_filename);
Double_t getTotalArea(TObjArray* ADC_values,int pre_PULSE_region,int integration_window);
int getFilenames(std::string Filenames_file,std::vector<std::string> & filenames);
void ROOT_file_Topology(std::vector<int> & topology, std::string root_filename);

void Mapping_116()
{
    TDatime().Print();
    // Create vector to define the map of the wall
    std::vector<int> topology(2);
    topology[0] = 20; // Columns
    topology[1] = 13; // Rows

    std::vector<std::string> filenames;
    int file_check = getFilenames("../data/filenames_116.txt",filenames);

    if (file_check == 1)
    {
        ROOT_file_Topology(topology,"../ROOT_files/Mapping_116.root");

        for (int filenum = 0; filenum < filenames.size(); filenum++)
        {
            std::cout << filenames.at(filenum) << std::endl;
            Read_waveforms(topology,filenames.at(filenum),"../ROOT_files/Mapping_116.root");
            TDatime().Print();
        } 
    }else{
        return;
    }
}

int getFilenames(std::string Filenames_file,std::vector<std::string> & filenames)
{

    TString pathFileName(Filenames_file.c_str());
    string LINE;
    ifstream infile;
    infile.open(pathFileName);
    TString separator = " ";
    if (!infile)
    {
        std::cout << "Unable to open filenames.txt" << std::endl;
        return 0;
    }
    while(!infile.eof())
    {
        getline(infile,LINE);
        //std::cout << LINE << std::endl;
        filenames.push_back("../data/" + LINE);
    }
    infile.close();
    return 1;
}

void Read_waveforms(std::vector<int> & topology, std::string data_filename,std::string root_filename)
{
    TFile* root_file = new TFile(root_filename.c_str(),"UPDATE");
    if (root_file == nullptr)
    {
        std::cout << "null pointer" << std::endl;
    }

    TH2I* Mapping_hist = (TH2I*)root_file->Get("Mapping_hist");
    TH2I* Amp_hist = (TH2I*)root_file->Get("Amp_hist");

    TString pathFileName(data_filename.c_str());
    string LINE;
    ifstream infile;
    infile.open(pathFileName);

    std::vector<int> waveform_nums;
    std::vector<int> amp_nums;
    for (int veci = 0; veci < 280; veci++)
    {
        waveform_nums.push_back(0);
        amp_nums.push_back(0);
    }
    
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
    Int_t Amp;
    
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
            TObjString Amp_word(*((TObjString*)((*words)[23])));

            Slot = (int)slotWord.String().Atof();
            Channel = (int)channelWord.String().Atof();
            Signal_boolean = (int)Signal_word.String().Atof();
            Amp = std::abs((int)Amp_word.String().Atof());
            //std::cout << Amp << std::endl;

            if (Channel == 13){Signal_boolean=0;}
            

        }
        
        if (line_number == 2)
        {
            int OM_num = Slot + (Channel*20);
            
            root_file->cd();
            line_number = -1;
            
            if (Signal_boolean != 0)
            {
                waveform_nums.at(OM_num)++;
                Mapping_hist->Fill(Slot,Channel,1);
                amp_nums.at(OM_num) += Amp;
                
            }
           
        }
        delete words;
        line_number++;
    }
    for (int slot_num = 0; slot_num < topology.at(0);slot_num++)
    {
        for (int channel_num = 0; channel_num < topology.at(1);channel_num++)
        {
            int OM_NUM = slot_num + (channel_num*20);
            if (waveform_nums.at(OM_NUM) != 0)
            {
                Amp_hist->Fill(slot_num,channel_num,(double)amp_nums.at(OM_NUM)/(double)waveform_nums.at(OM_NUM));
            }

        }
    }
    Amp_hist->Write();
    Mapping_hist->Write("",TObject::kOverwrite);
    delete Amp_hist;
    delete Mapping_hist;
    root_file->Close();
}

void ROOT_file_Topology(std::vector<int> & topology, std::string root_filename)
{
    TFile* root_file = new TFile(root_filename.c_str(),"RECREATE");

    TH2I* Mapping_hist = new TH2I("Mapping_hist","Mapping_hist",topology.at(0),0,topology.at(0),topology.at(1),0,topology.at(1));
    TH2I* Amp_hist = new TH2I("Amp_hist","Amp_hist",topology.at(0),0,topology.at(0),topology.at(1),0,topology.at(1));

    Mapping_hist->Write();
    Amp_hist->Write();
    delete Mapping_hist;
    delete Amp_hist;
    
    root_file->Close();
    std::cout << ">>> ROOT file created with the wall topology" << std::endl;
}