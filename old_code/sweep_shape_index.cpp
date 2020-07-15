#include "TF1.h"
#include "TStyle.h"
#include "TObjString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDatime.h"
#include <iostream>
#include <fstream>
#include <string>

void Read_waveforms(std::vector<int> & topology, std::string data_filename,std::string root_filename,std::vector<double> & template_vector);
int getFilenames(std::string Filenames_file,std::vector<std::string> & filenames);
void ROOT_file_Topology(std::vector<int> & topology, std::string root_filename);
void Sweep(TObjArray* ADC_values,std::vector<double> & template_vector,Double_t calculated_baseline, TH1D* shape_hist, TH1D* amp_hist);
int getTemplate(int Slot,int Channel,std::string root_filename,std::vector<double> & template_vector);
void Normalise_vector(std::vector<double> & input_vector);
Double_t getInnerProduct(std::vector<double> & input_vector1, std::vector<double> & input_vector2);
Double_t getBaseline(TObjArray* ADC_values,int pre_PULSE_region);

void sweep_shape_index()
{
    TDatime().Print();
    // Create vector to define the map of the wall
    std::vector<int> topology(2);
    topology[0] = 20; // Columns
    topology[1] = 13; // Rows

    std::vector<std::string> filenames;
    int file_check = getFilenames("/Users/willquinn/Documents/PhD/SNEMO_ComData_Analysis/French_CaloWall/data/filenames_100.txt",filenames);

    std::vector<double> template_vector;
    int template_check = getTemplate(0,0,"/Users/willquinn/Documents/PhD/SNEMO_ComData_Analysis/French_CaloWall/ROOT_files/Average_Waveforms_run100.root",template_vector);
    if (template_check == 0)
    {
        std::cout << "Couldn't grab Template " << std::endl;
        return;
    }
    TH1D* template_hist = new TH1D("template_hist","template_hist",template_vector.size(),0,template_vector.size());
    for (int itemp = 0; itemp < template_vector.size(); itemp++)
    {
        template_hist->SetBinContent(itemp,template_vector.at(itemp));
    }
    
    //std::cout << "Got template " << std::endl;
    Normalise_vector(template_vector);
    //std::cout << "Normalised template " << std::endl;
    if (file_check != 0)
    {
        ROOT_file_Topology(topology,"/Users/willquinn/Documents/PhD/SNEMO_ComData_Analysis/French_CaloWall/ROOT_files/sweep_100.root");

        TFile* root_file = new TFile("/Users/willquinn/Documents/PhD/SNEMO_ComData_Analysis/French_CaloWall/ROOT_files/sweep_100.root","UPDATE");
        root_file->cd();

        template_hist->Write();
        delete template_hist;
        root_file->Close();

        for (int filenum = 0; filenum < filenames.size(); filenum++)
        {
            std::cout << filenames.at(filenum) << std::endl;
            Read_waveforms(topology,filenames.at(filenum),"/Users/willquinn/Documents/PhD/SNEMO_ComData_Analysis/French_CaloWall/ROOT_files/sweep_100.root",template_vector);
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
        filenames.push_back("/Users/willquinn/Documents/PhD/SNEMO_ComData_Analysis/French_CaloWall/data/" + LINE);
    }
    infile.close();
    return 1;
}

void Read_waveforms(std::vector<int> & topology, std::string data_filename,std::string root_filename,std::vector<double> & template_vector)
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

    std::vector<int> waveform_nums;
    for (int inum = 0; inum < 280; inum++)
    {
        waveform_nums.push_back(0);
    }
    
    int skip_header_lines = 9;
    
    for (int iSkip=0; iSkip<skip_header_lines; ++iSkip)
    {
        if (!infile.eof()) getline(infile,LINE);
    }
    
    int line_number = 0;
    int waveform_num = 0;
    
    TString separator = " ";

    Int_t Slot;
    Int_t Channel;
    Int_t Signal_boolean;
    
    while(!infile.eof() && waveform_num < 100)
    {
        getline(infile,LINE);

        //TString LineString(LINE);
        //std::cout << "Line " << LineString << std::endl;
        
        if (line_number == 1)
        {
            TString LineString(LINE);
            //std::cout << "Line " << LineString << std::endl;
            TObjArray* words = LineString.Tokenize(separator.Data());

            // Extract channel number :
            TObjString slotWord(*((TObjString*)((*words)[1])));
            TObjString channelWord(*((TObjString*)((*words)[3])));
            TObjString Signal_word(*((TObjString*)((*words)[7])));

            Slot = (int)slotWord.String().Atof();
            Channel = (int)channelWord.String().Atof();
            Signal_boolean = (int)Signal_word.String().Atof();

            if (Channel == 13){Signal_boolean=0;}

            delete words;

        }
        
        if (line_number == 2)
        {
            line_number = -1;
            
            if (Signal_boolean != 0)
            {
                int OM_num = Slot + (Channel*20);
                char Shape_Number[30];
                sprintf(Shape_Number,"Shape_%d_%d_%d",Slot,Channel,waveform_nums.at(OM_num));
                char Amplitude_Number[30];
                sprintf(Amplitude_Number,"Amplitude_%d_%d_%d",Slot,Channel,waveform_nums.at(OM_num));

                //std::cout << waveform_num << std::endl;
                TString LineString(LINE);
                TObjArray* words = LineString.Tokenize(separator.Data());

                TH1D* shape_hist = new TH1D(Shape_Number,Shape_Number,words->GetEntries() - template_vector.size(),0,words->GetEntries() - template_vector.size());
                TH1D* amp_hist = new TH1D(Amplitude_Number,Amplitude_Number,words->GetEntries() - template_vector.size(),0,words->GetEntries() - template_vector.size());

                Double_t calculated_baseline = getBaseline(words,100);
                Sweep(words,template_vector,calculated_baseline,shape_hist,amp_hist);

                shape_hist->Write();
                amp_hist->Write();
                delete words;
                delete shape_hist;
                delete amp_hist;
                waveform_num++;
                waveform_nums.at(OM_num)++;
            }
        }
        line_number++;
    }
    root_file->Close();
}

void ROOT_file_Topology(std::vector<int> & topology, std::string root_filename)
{
    TFile* root_file = new TFile(root_filename.c_str(),"RECREATE");
    root_file->Close();
    std::cout << ">>> ROOT file created with the wall topology" << std::endl;
}

void Sweep(TObjArray* ADC_values,std::vector<double> & template_vector,Double_t calculated_baseline, TH1D* shape_hist, TH1D* amp_hist)
{
    //std::cout << "In sweep " << std::endl;
    for (int isweep = 0; isweep < ADC_values->GetEntries() - template_vector.size(); isweep++)
    {   
        std::vector<double> test_vector;
        for (int ivec = 0; ivec < template_vector.size(); ivec++)
        {
            TObjString thisWord(*((TObjString*)((*ADC_values)[isweep + ivec])));
            Double_t ADC_value = (Double_t)thisWord.String().Atof();
            test_vector.push_back(ADC_value - calculated_baseline);
        }
        Double_t amp_index = getInnerProduct(template_vector,test_vector);
        Normalise_vector(test_vector);
        Double_t shape_index = getInnerProduct(template_vector,test_vector);
        //std::cout << "shape_index " << shape_index << std::endl;
        shape_hist->SetBinContent(isweep,shape_index);
        amp_hist->SetBinContent(isweep,amp_index);
    }
    
}

int getTemplate(int Slot,int Channel,std::string root_filename,std::vector<double> & template_vector)
{
    //std::cout << "getting Template" << std::endl;
    TFile* root_file = new TFile(root_filename.c_str(),"READ");
    if (root_file == nullptr)
    {
        return 0;
    }
    
    char Waveform_Number[30];
    sprintf(Waveform_Number,"Waveform_%d_%d_average",Slot,Channel);

    TH1D* template_hist = (TH1D*)root_file->Get(Waveform_Number);
    if (template_hist == nullptr)
    {
        std::cout << "No such template called " << Waveform_Number << " in file " << root_filename << std::endl;
        return 0;
    }
    //std::cout << "Got Histogram" << std::endl;

    for (int itemp = (int)template_hist->GetMinimumBin() - 60; itemp < ((int)template_hist->GetMinimumBin() + 150); itemp++)
    {
        template_vector.push_back(template_hist->GetBinContent(itemp));
        //std::cout << itemp << " " << template_hist->GetMinimumBin() << std::endl;
        if (template_vector.size() > 1000)
        {
            std::cout << "Error: Template exceding 50 bins" << std::endl;
            return 0;
        }
    }
    return 1;
}

void Normalise_vector(std::vector<double> & input_vector)
{
    double Normalisation_factor = 0.;
    for (int inorm = 0; inorm < input_vector.size(); inorm++)
    {
        Normalisation_factor += input_vector.at(inorm)*input_vector.at(inorm);
    }
    Normalisation_factor = sqrt(Normalisation_factor);

    for (int inorm = 0; inorm < input_vector.size(); inorm++)
    {
        input_vector.at(inorm) = input_vector.at(inorm)/Normalisation_factor;
    }
    return;
}

Double_t getInnerProduct(std::vector<double> & input_vector1, std::vector<double> & input_vector2)
{//Both vectors must be normalised
    if (input_vector1.size() != input_vector2.size())
    {
        std::cout << "Error: sizes of vectors do not match " << std::endl;
        return 0.0;
    }
    Double_t shape_index = 0.;
    for (int iprod = 0; iprod < input_vector1.size(); iprod++)
    {
        shape_index += (Double_t)(input_vector1.at(iprod)*input_vector2.at(iprod));
    }
    return shape_index;
}

Double_t getBaseline(TObjArray* ADC_values,int pre_PULSE_region)
{   
    Double_t calculated_baseline = 0.;
    int iword = 0;
    while (iword < pre_PULSE_region)
    {
        TObjString thisWord(*((TObjString*)((*ADC_values)[iword])));
        Double_t ADC_value = (Double_t)thisWord.String().Atof();
        calculated_baseline += ADC_value;
        iword++;
    }
    calculated_baseline = calculated_baseline/(Double_t)iword;
    return calculated_baseline;
}




