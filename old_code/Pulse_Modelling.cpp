#include "TF1.h"
#include "TStyle.h"
#include "TObjString.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include <string>

void ROOT_file_Topology(std::vector<int> & topology, std::string root_filename);
void Read_waveforms(std::vector<int> & topology, std::string data_filename,std::string root_filename);
Double_t getTotalArea(TObjArray* ADC_values,int pre_PULSE_region,int integration_window);
int GetParameter(TH1D* p, double* par);
double checkWidth(TObjArray* ADC_values, int pre_PULSE_region);
double mypulse(double *x, double *par);
TFitResultPtr pulse_fit(TH1D* data, double* initpar);

void Spectra()
{
    // Create vector to define the map of the wall
    std::vector<int> topology(2);
    topology[0] = 20; // Columns
    topology[1] = 13; // Rows

    ROOT_file_Topology(topology,"../ROOT_files/Waveforms_BackUP_Run104_Spectra.root");
    Read_waveforms(topology,"../data/run104_1.dat","../ROOT_files/Waveforms_BackUP_Run104_Spectra.root");
}

void Read_waveforms(std::vector<int> & topology, std::string data_filename,std::string root_filename)
{
    TFile* root_file = new TFile(root_filename.c_str(),"UPDATE");
    if (root_file == nullptr)
    {
        std::cout << "null pointer" << std::endl;
    }

    TH1D* Chi2_hist = new TH1D("Chi2_hist","Chi2_hist",50,0,4);
    Chi2_hist->SetXTitle("Chi2_reduced");

    TH2D* Chi2_hist2 = new TH2D("Chi2_hist2","Chi2_hist2",50,-2500,0,50,0,3);

    TString pathFileName(data_filename.c_str());
    //std::cout << "Reading waveform from file : " << pathFileName << std::endl;
    string LINE;
    ifstream infile;
    infile.open(pathFileName);
    
    int skip_header_lines = 9;
    
    for (int iSkip=0; iSkip<skip_header_lines; ++iSkip)
    {
        if (!infile.eof()) getline(infile,LINE);
    }
    
    int line_number = 0;
    //int line_counter = 0;
    
    TString separator = " ";

    Int_t Slot;
    Int_t Channel;
    Int_t Signal_boolean;

    int waveform_num = 0;
    std::vector<int> waveform_nums(280,0);
    
    while(!infile.eof() && waveform_num < 1000000)
    {
        getline(infile,LINE);
        
        TString LineString(LINE);
        TObjArray* words = LineString.Tokenize(separator.Data());
        //words->Print();
        
        if (line_number ==1)
        {
            // Extract channel number :
            TObjString slotWord(*((TObjString*)((*words)[1])));
            TObjString channelWord(*((TObjString*)((*words)[3])));
            TObjString Signal_word(*((TObjString*)((*words)[7])));

            Slot = (int)slotWord.String().Atof();
            Channel = (int)channelWord.String().Atof();
            Signal_boolean = (int)Signal_word.String().Atof();

            //std::cout << "Slot " << Slot << " Channel " << Channel << std::endl;
            //std::cout << "Signal Boolean " << Signal_boolean << std::endl;
            if (Channel == 13){Signal_boolean=0;}

        }
        
        if (line_number == 2)
        {
            
            root_file->cd();
            //std::cout << "checking line" << line_counter << std::endl;
            line_number = -1;
            //line_counter ++;
            //std::cout << Slot << " " << Channel << std::endl;
            //std::cout << Channel << " " << channel_num << std::endl;
            if (Signal_boolean != 0)
            {
                char Spectra_histname[15];
                sprintf(Spectra_histname,"Spectra_%d_%d",Slot,Channel);

                char Width_histname[15];
                sprintf(Width_histname,"Width_%d_%d",Slot,Channel);

                char waveform_histname[30];
                sprintf(waveform_histname,"Waveform_%d_%d_%d",Slot,Channel,waveform_nums.at(Slot + (Channel*20)));

                char Slot_Number[15];
                sprintf(Slot_Number,"Slot%d",Slot);

                char Channel_Number[15];
                sprintf(Channel_Number,"Channel%d",Channel);

                root_file->GetDirectory(Slot_Number)->GetDirectory(Channel_Number)->cd();

                //TH1D* spectra_hist = (TH1D*)root_file->GetDirectory(Slot_Number)->GetDirectory(Channel_Number)->Get(Spectra_histname);
                TH1D* width_hist = (TH1D*)root_file->GetDirectory(Slot_Number)->GetDirectory(Channel_Number)->Get(Width_histname);

                TH1D* waveform_hist = new TH1D(waveform_histname,waveform_histname,words->GetEntries()-2,0,words->GetEntries()-2);

                
                Double_t Total_Integrated_Area = getTotalArea(words,100,400);

                for (int iword=0; iword < words->GetEntries()-2; iword++)
                {
                    TObjString thisWord(*((TObjString*)((*words)[iword])));
                    Double_t ADC_value = (Double_t)thisWord.String().Atof();
                    waveform_hist->SetBinContent(iword,ADC_value);
                    //std::cout << iword << " - " << ADC_value << std::endl;
                }

                double width = checkWidth(words,100);
                width_hist->Fill(width);
                width_hist->Write("",TObject::kOverwrite);
                delete width_hist;


                //std::cout << width << std::endl;
                
                //spectra_hist->Fill(Total_Integrated_Area);
                
                //spectra_hist->Write("",TObject::kOverwrite);

                double* guess_parameters;
                guess_parameters = new double[5];

                int test = GetParameter(waveform_hist,guess_parameters);

                for (int i = 0; i < 5; i++)
                {
                   //std::cout << guess_parameters[i] << std::endl; 
                }
                
                TFitResultPtr fitted_parameters = (TFitResultPtr)pulse_fit(waveform_hist,guess_parameters);
                
                //std::cout << fitted_parameters << std::endl;;

                waveform_num++;
                //std::cout << waveform_num << std::endl;

                waveform_nums.at(Slot + (Channel*20))++;
                if (test == 0)
                {
                    waveform_hist->Write();
                    //std::cout << (double)fitted_parameters->Chi2()/(1023. - 5.) << std::endl;
                }
                Chi2_hist->Fill((double)fitted_parameters->Chi2()/(1023. - 5.));
                Chi2_hist2->Fill((double)fitted_parameters->Parameter(0),(double)fitted_parameters->Chi2()/(1023. - 5.));
                
                delete waveform_hist;
            }
           
        }
        delete words;
        line_number++;
    }
    root_file->cd();
    Chi2_hist->Write();
    Chi2_hist2->Write();
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
            TH1D* spectra_hist = new TH1D(Spectra_Number,Spectra_Number,50,0,100000);
            TH1D* width_hist = new TH1D(Width_Number,Width_Number,10,0,1);

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

TFitResultPtr pulse_fit(TH1D* data, double* initpar)
{
    double mypulse(double* x, double* par);
    int npoints = data->GetEntries();

    TF1* fitfunction = new TF1("func",mypulse,0.0,(double)npoints,5);
    fitfunction->SetParName(0, "Amplitude");
    fitfunction->SetParName(1, "Rise Time");
    fitfunction->SetParName(2, "Decay Time");
    fitfunction->SetParName(3, "Pulse Onset");
    fitfunction->SetParName(4, "Baseline");
    fitfunction->SetParameter(0, initpar[0]);
    fitfunction->SetParameter(1, initpar[1]);
    fitfunction->SetParameter(2, initpar[2]);
    fitfunction->SetParameter(3, initpar[3]);
    fitfunction->SetParameter(4, initpar[4]);
    fitfunction->SetParLimits(0,-10000,0);
    //fitfunction->SetParLimits(1,0,1000);
    //fitfunction->SetParLimits(2,10,10000);
    //fitfunction->SetParLimits(3,0,7);
    //fitfunction->SetParLimits(4,180,200);
    //fitfunction->SetParLimits(5,1500,2500);
    //fitfunction->SetParLimits(6,0,1);
    TFitResultPtr fr = data->Fit(fitfunction,"SNQ");

    delete fitfunction;
    return fr;
}

double mypulse(double *x, double *par)
{
    double xx = x[0];
    if (xx < par[3])
        {return (par[4]);}
    double f = par[0]*(TMath::Exp(-(xx - par[3])/par[2]) - TMath::Exp(-(xx - par[3])/par[1]));
    return (f + par[4]);
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

int GetParameter(TH1D* p, double* par)
{
    int Pulselength = p->GetNbinsX();
    int fiveperc = (int)(0.04*Pulselength);
    double baseline = p->Integral(0,fiveperc);
    baseline /= (double)fiveperc;
    double amplitude = baseline;
    int position = 0;
    for (int i = 0; i < Pulselength-1; i++)
    {
        if (amplitude > (double)p->GetBinContent(i))
        {
            amplitude = (double)p->GetBinContent(i);
            position = (int)p->GetBin(i);
        }
    }

    par[4] = baseline; // storage place in par array for baseline parameter

    // corrupt pulses I
    if (baseline<=5.0) 
    { // adjust as needed
        par[0] = par[1] = par[2] = par[3] = 0.0;
        //std::cout << "Bad Pulse I" << std::endl;
        return 0;  
    }

    par[0] = amplitude - baseline; // amplitude in [3]
    //par[4] = (double)position; // position in [4]

    // corrupt pulses II
    if (position<=(2*fiveperc)) 
    { // early pulse max - reject
        par[1] = par[2] = 0.0;
        //std::cout << "Bad Pulse II" << std::endl;
        return 0;  
    }

    int start=0;
    int stop=0;
    // Onset 30%
    // stop at 70% height
    for (int i=fiveperc;i<=position;i++) 
    {
        if (((double)(p->GetBinContent(i))-baseline) >= (0.3*(amplitude - baseline)))
        {
            start = i;
            //std::cout << start << std::endl;
        }
        if (((double)(p->GetBinContent(i))-baseline) <= (0.7*(amplitude - baseline)))
        {
            stop = i;
            //std::cout << stop << std::endl;
        }
    }
    par[1] = (double)(stop - start); 
    //std::cout << stop << " " << start << std::endl;// rise time parameter in [2]
    par[3] = (double)start;  // onset in [1]

    // corrupt pulses III
    if (stop<=start || start<=0.0) 
    {
        par[1] = par[3] = par[0] = par[2] = 0.0;
        //std::cout << "Bad Pulse III" << std::endl;
        return 0;  
    }

      // decay time
    for (int i=position;i<Pulselength;i++) 
        if (((double)(p->GetBinContent(i))-baseline) >= ((amplitude - baseline)/TMath::Exp(1.0)))
            {stop = i;}

    par[2] = (double)(stop - position); // in [5]
      
    return 1;
}
