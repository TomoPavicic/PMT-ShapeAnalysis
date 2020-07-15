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

void fit_function(std::string data_file_name_prefix,std::string date,int channel_classifier) 
{
	double shift = channel_classifier*4; 

	std::string root_filename = "../../ROOT_files/" + date + "_" + data_file_name_prefix + "_areas.root";
	TFile* file = new TFile(root_filename.c_str());
	char Spectra_histname[15];
	sprintf(Spectra_histname,"Spectra_Ch%d",channel_classifier);

	TH1D* Total_Area1 = (TH1D*)file->Get(Spectra_histname);
	//TH1F* Total_Area2 = (TH1F*)file->Get("total_area");

	TCanvas* c1 = new TCanvas();

	TF1* fit1 = new TF1("fit1","[0]*(7.08*TMath::Gaus(x,[1],[2]) + 1.84*TMath::Gaus(x,[1]*(1 + 72.144/975.651),[2]*1.036) + 0.44*TMath::Gaus(x,[1]*(1 + 84.154/975.651),[2]*1.042)) + [3]*(exp([4]*x)/(1 + exp((x - [5])/[6])))",29+shift,38+shift);
    
    //TF1* fit2 = new TF1("fit2","[0]*(1.537*TMath::Gaus(x,[1],[2]) + 0.111*TMath::Gaus(x,[1]*(1 + 83.1538/481.6935),[2]*1.084)) + [3]*(exp([4]*x)/(1 + exp((x - [5])/[6]))) + [7]*(exp([8]*x)/(1 + exp((x - [9])/[10])))",10,20); //the last exponent removed 

	//TF1* fit2 = new TF1("fit2","[0]*(1.537*TMath::Gaus(x,[1],[2]) + 0.442*TMath::Gaus(x,[1]*(1 + 72.144/481.6935),[2]*1.072) + 0.111*TMath::Gaus(x,[1]*(1 + 83.1538/481.6935),[2]*1.084)) + [3]*(exp([4]*x)/(1 + exp((x - [5])/[6]))) + [7]*(exp([8]*x)/(1 + exp((x - [9])/[10])))",10,20);
    

	fit1->SetParNames("A_1","mu_1","sigma_1","B1","C1","ce_1","slope_1");
	
	fit1->SetParLimits(0,0,400);
	fit1->SetParLimits(1,30+shift,32+shift);
	fit1->SetParLimits(2,0.8,2);
	fit1->SetParLimits(3,0,600);
	fit1->SetParLimits(4,0,0.09);
	fit1->SetParLimits(5,20+shift,30+shift);
	fit1->SetParLimits(6,0,10);
	//fit1->SetParLimits(7,0,1000);
	

	fit1->SetParameters(319,31+shift,1.09,10,0.001,26+shift,1);
	Total_Area1->Fit("fit1","","",28+shift,36+shift);

	//Total_Area1->GetXaxis()->SetRange(2,400);
	Total_Area1->SetXTitle("Charge /pC");
	Total_Area1->SetYTitle("Counts");
	Total_Area1->SetTitle("Bi Integrated Charge Spectrum");

	Total_Area1->Draw();
	c1->SetGrid();
	c1->Update();

/*
	TCanvas* c2 = new TCanvas();

	fit2->SetParNames("A_2","mu_2","sigma_2","B2","C2","ce_2","slope_2","B1","C1","ce_1","slope_1");
	
	fit2->SetParLimits(0,200,700);
	fit2->SetParLimits(1,11,13);
	fit2->SetParLimits(2,0.5,0.7);
	fit2->SetParLimits(3,0,1000);
	fit2->SetParLimits(4,0,1);
	fit2->SetParLimits(5,8,11);
	fit2->SetParLimits(6,0,2);
	fit2->SetParLimits(7,69,70);
	fit2->SetParLimits(8,0.1,0.12);
	fit2->SetParLimits(9,25,26);
	fit2->SetParLimits(10,0.5,0.6);
	

	fit2->SetParameters(300,11.5,0.6,100,0.001,9,4,69.8004,0.117317,25.6946,0.511693);
	Total_Area2->Fit("fit2","","",11,15);

	Total_Area2->GetXaxis()->SetRange(2,350);
	Total_Area2->SetXTitle("Charge /pC");
	Total_Area2->SetYTitle("Counts");
	Total_Area2->SetTitle("Bi Integrated Charge Spectrum");

	Total_Area2->Draw();
	c2->SetGrid();
	c2->Update();
	*/
	
}
