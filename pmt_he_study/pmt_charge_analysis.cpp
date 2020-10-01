// Standard
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

// Argument parser
#include <boost/program_options.hpp>

// ROOT
#include <TTree.h>
#include <TDatime.h>
#include <TFile.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TDatime.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>

namespace po = boost::program_options;


typedef struct {
    Int_t sweep_start, pre_trigger, trigger, trig_tolerance;
    Double_t shape_cut, amp_cut, charge_cut, resistance;
    std::vector<Double_t> integration;
} CONF;

typedef struct {
    std::string date, timestamp, prefix;
    Int_t voltage, tot_event_ch0, tot_event_ch1;
} DESC;

std::vector<Double_t> process_line( const std::string & s, char delimiter );
Double_t get_baseline( std::vector<Double_t> &vec, CONF &config );
Double_t get_charge( std::vector<Double_t> &vec, Double_t baseline, CONF &config, Int_t peak_cell );
Double_t get_amplitude( std::vector<Double_t> &vec, Double_t baseline );
Int_t get_peak_cell( std::vector<Double_t> &vec );
std::vector<std::string> split( const std::string& s, char delimiter );
CONF read_config( std::string filename );
DESC process_name( std::string &s );
void fit_function(TFile* rootfile,int channel_classifier, std::vector<Double_t> & results, DESC  &description);
Double_t bi_func_ch0(Double_t *x, Double_t *par);
Double_t bi_func_ch1(Double_t *x, Double_t *par);

Int_t main(Int_t argc, char* argv[])
{
    std::string input_file;
    std::string config_file;
    std::string output_file;

    po::options_description desc("Allowed options"); desc.add_options()
            ("help", "produce help message")
            ("i", po::value(&input_file), "Name of input file")
            ("c", po::value(&config_file), "Name of configuration file")
            ("o", po::value(&output_file), "Name of output file");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 0;
    }

    // Process the arguements

    const std::string tag_ch0 = "trace channel=\"0\"";
    const std::string tag_ch1 = "trace channel=\"1\"";

    // ==================================================================
    // Get the relevant strings
    DESC description = process_name( input_file );

    std::cout << ">>> File process      : " << input_file << std::endl;
    std::cout << ">>> Data file prefix  : " << description.prefix << std::endl;;
    std::cout << ">>> Output filename   : " << output_file << std::endl;;
    std::cout << ">>> Time-stamp        : " << description.timestamp << std::endl;;
    std::cout << ">>> Date              : " << description.date << std::endl;;
    std::cout << ">>> Voltage           : " << description.voltage << std::endl;;

    // Read the config file
    CONF config_object = read_config( config_file );
    std::cout << ">>> Settings          : " << std::endl;
    std::cout << ">>> Amp cut           : " << config_object.amp_cut << std::endl;
    std::cout << ">>> Shape cut         : " << config_object.shape_cut << std::endl;
    std::cout << ">>> Charge cut        : " << config_object.charge_cut << std::endl;
    std::cout << ">>> Pre-trigger       : " << config_object.pre_trigger << std::endl;
    std::cout << ">>> Trigger           : " << config_object.trigger << std::endl;
    std::cout << ">>> Trigger tolerance : " << config_object.trig_tolerance << std::endl;
    std::cout << ">>> Integration       : " << std::endl;
    std::cout << ">>> Low/High edge     : " << config_object.integration[0] << "/" << config_object.integration[1] << std::endl;
    std::cout << ">>> Sweep start       : " << config_object.sweep_start << std::endl;
    std::cout << ">>> Resistance        : " << config_object.resistance << std::endl;

    std::vector<std::string> results_files;
    results_files.push_back("/unix/nemo4/PMT_He_Study_nemo4/Results/" + std::to_string(description.voltage) + "/resolution_vs_date_Ch0.dat");
    results_files.push_back("/unix/nemo4/PMT_He_Study_nemo4/Results/" + std::to_string(description.voltage) + "/resolution_vs_date_Ch1.dat");

    TFile root_file( output_file.c_str(), "RECREATE" );
    // ==================================================================

    // Create a container to hold all the waveform information
    // this will be the construct to fill Int_to the TTree
    typedef struct {
        Int_t OM_ID;
        Double_t pulse_charge, pulse_amplitude, baseline;
    } EVENTN;
    EVENTN eventn;

    TTree tree("T","Tree containing waveform information");
    tree.Branch("pulse_charge",&eventn.pulse_charge);
    tree.Branch("pulse_amplitude",&eventn.pulse_amplitude);
    tree.Branch("pulse_baseline",&eventn.baseline);
    tree.Branch("OM_ID",&eventn.OM_ID);

    // Create the histograms to hold the information we want
    TH1D area_spectra_hist_0("Area_Spectra_Ch0","Area_Spectra_Ch0",200,0,60);
    TH1D area_spectra_hist_1("Area_Spectra_Ch1","Area_Spectra_Ch1",200,0,60);
    TH1D amp_spectra_hist_0("Amplitude_Spectra_Ch0","Amplitude_Spectra_Ch0",200,0,1000);
    TH1D amp_spectra_hist_1("Amplitude_Spectra_Ch1","Amplitude_Spectra_Ch1",200,0,1000);
    TH1D baseline_hist_0("Baseline_Spectra_Ch0","Baseline_Spectra_Ch0",200,900,1000);
    TH1D baseline_hist_1("Baseline_Spectra_Ch1","Baseline_Spectra_Ch1",200,900,1000);
    std::vector<TH1D> area_hists;
    area_hists.push_back(area_spectra_hist_0);
    area_hists.push_back(area_spectra_hist_1);
    std::vector<TH1D> amp_hists;
    amp_hists.push_back(amp_spectra_hist_0);
    amp_hists.push_back(amp_spectra_hist_1);
    std::vector<TH1D> baseline_hists;
    baseline_hists.push_back(baseline_hist_0);
    baseline_hists.push_back(baseline_hist_1);

    std::ifstream data_file( input_file );
    std::string data_line;

    if (!data_file.good())
    {
        std::cout << "matched_filter.cpp : ERROR opening input file : " << input_file << std::endl;
        std::cout << "EXIT" << std::endl;
        exit(1);
    }

    Int_t channel_indicator{0};

    std::vector<int> channel_event_num( 2, 0 );
    std::vector<int> thousand_counter( 2, 0 );
    std::vector<int> channel_waveform_num( 2, 0 );

    TDatime().Print();
    std::cout << ">>> Beginning data read..." << std::endl;

    // Read the data file
    while ( std::getline( data_file, data_line ) && !data_file.eof() )
    {
        // If the line is a data line we can filter it out to process
        if (data_line.find(tag_ch0) != std::string::npos) {
            channel_indicator = 0;
            //std::cout << "Channel: " << channel_indicator << std::endl;
        } else if (data_line.find(tag_ch1) != std::string::npos) {
            channel_indicator = 1;
            //std::cout << "Channel: " << channel_indicator << std::endl;
        } else { continue; } // All other lines are garbage

        // Process      =================================================
        std::vector <Double_t> data = process_line(data_line, ' ');
        //std::cout << "Data vector length: " << data.size() << std::endl;

        // Analysis     =================================================
        channel_event_num[channel_indicator] ++;

        Int_t peak_cell = get_peak_cell( data );

        if ( peak_cell > config_object.trigger + config_object.trig_tolerance || peak_cell < config_object.trigger - config_object.trig_tolerance )
        {
            // Waveform is likely empty so will not store
            continue;
        }

        Double_t baseline = get_baseline( data, config_object );
        Double_t pulse_charge = get_charge( data, baseline, config_object, peak_cell );

        // std::cout << "charge : " << pulse_charge << std::endl;

        if ( pulse_charge < config_object.charge_cut )
        {
            // Pulse is either too small or just noise
            continue;
        }

        channel_waveform_num[channel_indicator]++;

        Double_t pulse_amplitude = get_amplitude( data, baseline );

        // Output       =================================================

        amp_hists[channel_indicator].Fill(pulse_amplitude);
        area_hists[channel_indicator].Fill(pulse_charge);
        baseline_hists[channel_indicator].Fill(baseline);

        eventn.pulse_amplitude  = pulse_amplitude;
        eventn.pulse_charge     = pulse_charge;
        eventn.baseline         = baseline;
        eventn.OM_ID            = channel_indicator;
        tree.Fill();

        // Show the progress of the read file
        if ( channel_waveform_num[channel_indicator] % 1000 == 0)
        {
            thousand_counter[channel_indicator]++;
            std::cout << std::endl;
            TDatime().Print();
            std::cout << ">>> #1000 waveforms analysed Ch" << channel_indicator << " : " << thousand_counter[channel_indicator] << std::endl;
        }
    }

    // Write the TTree to the file
    std::cout << ">>> Number of events in file: " << channel_event_num[0] + channel_event_num[1] << std::endl;
    std::cout << ">>> Events per channel: " << channel_event_num[0] << " Ch0, " << channel_event_num [1] << " Ch1" << std::endl;
    std::cout << ">>> After cuts : " << std::endl;
    std::cout << ">>> Events in file : " << channel_waveform_num[0] + channel_waveform_num[1] << std::endl;
    std::cout << ">>> Events per channel: " << channel_waveform_num[0] << " Ch0, " << channel_waveform_num[1] << " Ch1" << std::endl;
    description.tot_event_ch0 = channel_event_num[0];
    description.tot_event_ch1 = channel_event_num[1];
    tree.Fill();

    // Output       =================================================
    std::vector<Double_t> results;
    for (int i_channel = 0; i_channel<2 ; i_channel++)
    {
        std::string canvas_name =
                "/unix/nemo4/PMT_He_Study_nemo4/ROOT_files/Area_Spectrum_Files/" + std::to_string(description.voltage) + "V/PDFs/" +
                description.date + "_Canvas_Spectra_Ch" + std::to_string(i_channel) + ".pdf";

        std::ofstream resolution_file;
        results.clear();
        fit_function(&root_file,i_channel, results, description);

        if (results.size() > 0) //filter the not used channels
        {
            resolution_file.open(results_files[i_channel], std::ios_base::app);
            std::cout << ">>> Writing to output file: " << results_files[i_channel];
            resolution_file << TDatime().Get()
                            << ","
                            << description.date
                            << ","
                            << description.timestamp
                            <<","
                            << results.at(0)
                            <<","
                            <<results.at(1)
                            <<","
                            <<results.at(2)
                            <<","
                            <<results.at(3)
                            <<","
                            << results.at(2)/results.at(0) * 100.0
                            <<","
                            << results.at(4) << std::endl;
        }
    }
    //              =================================================
    root_file.Write();
}


std::vector<Double_t> process_line( const std::string & s, char delimiter )
{
    // Split the string of data Int_to a vector of Double_ts
    std::vector<Double_t> vec;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        if ( token.find ( ">" ) != std::string::npos )
        {
            continue;
        }else if( token.find ( "<" ) != std::string::npos )
        {
            continue;
        }else{
            if ( (Int_t)token.length() < 4 && (Int_t)token.length() > 0 ){
                vec.push_back( (Double_t)std::atoi(token.c_str()) );
            }
        }
    }
    return vec;
}
std::vector<std::string> split( const std::string& s, char delimiter )
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}
DESC process_name( std::string &s )
{
    DESC description;
    std::string date, time_stamp, prefix;
    Int_t voltage;

    std::vector<std::string> dir = split( s, '/' );
    date = dir[dir.size() - 2];

    // Remove the .xml
    std::vector<std::string> temp = split(dir[dir.size() - 1], '.');
    // Isolate
    std::vector<std::string> temp1 = split(temp[0], '_');
    // Get the voltage
    voltage = std::stoi(split( temp1[0], 'A' )[1]);

    time_stamp = split( temp1[2], 't' )[1];

    prefix = temp[0];

    description.date = date;
    description.timestamp = time_stamp;
    description.prefix = prefix;
    description.voltage = voltage;

    return description;
}
CONF read_config( std::string filename )
{
    CONF config;

    std::ifstream file(filename);
    std::string line;

    if (!file.good())
    {
        std::cout << "matched_filter.cpp : ERROR opening configuration file : " << filename << std::endl;
        std::cout << "EXIT" << std::endl;
        exit(1);
    }

    while ( std::getline( file, line ) && !file.eof() )
    {
        if ( line.find ( "#" ) != std::string::npos )
        {
            // A comment so ignore
            continue;

        }
        else if  (line.empty())
        {
            // Empty line so ignore
            continue;
        }else{
            std::vector<std::string> settings = split( line, ':' );
            if ( settings[0] == "integration" ) { config.integration[0] = std::stod(settings[1]); config.integration[1] = std::stod(settings[2]); }
            else if ( settings[0] == "sweep_start" ) { config.sweep_start = std::stoi(settings[1]); }
            else if ( settings[0] == "pre_trigger" ) { config.pre_trigger = std::stoi(settings[1]); }
            else if ( settings[0] == "shape_cut" ) { config.shape_cut = std::stod(settings[1]); }
            else if ( settings[0] == "amp_cut" ) { config.amp_cut = std::stod(settings[1]); }
            else if ( settings[0] == "charge_cut" ) { config.charge_cut = std::stod(settings[1]); }
            else if ( settings[0] == "trigger" ) { config.trigger = std::stod(settings[1]); }
            else if ( settings[0] == "trig_tolerance" ) { config.trig_tolerance = std::stod(settings[1]); }
            else if ( settings[0] == "resistance" ) { config.resistance = std::stod(settings[1]); }
            else { continue; }
        }
    }

    return config;
}
Double_t get_amplitude( std::vector<Double_t> &vec, Double_t baseline )
{
    Double_t amplitude = vec[0];
    for ( Int_t i = 0 ; i < (Int_t)vec.size() ; i++ )
    {
        if ( vec[i] < amplitude )
        {
            amplitude = vec[i];
        }
    }
    return (-1)*(amplitude - baseline);
}
Int_t get_peak_cell( std::vector<Double_t> &vec )
{
    Int_t peak_cell = 0;
    Double_t temp = vec[0];
    for ( Int_t i = 0 ; i < (Int_t)vec.size() ; i++ )
    {
        if ( vec[i] < temp )
        {
            temp = vec[i];
            peak_cell = i;
        }
    }
    return peak_cell;
}
Double_t get_charge( std::vector<Double_t> &vec, Double_t baseline, CONF &config, Int_t peak_cell )
{
    Double_t charge{0.0};

    Int_t start,end;
    for ( Int_t i = (Int_t)config.pre_trigger; i < peak_cell ; i++ )
    {
        if ( (vec[i] - baseline) < (vec[peak_cell] - baseline)*config.integration[0] )
        {
            start = i;
            break;
        }else{ continue; }
    }
    for ( Int_t i = peak_cell ; i < (Int_t)vec.size() ; i++ )
    {
        if ( (vec[i] - baseline) > (vec[peak_cell] - baseline)*config.integration[1] )
        {
            end = i;
            break;
        }else{ continue; }
    }

    for ( Int_t i = start ; i < end ; i++ ) { charge += (vec[i] - baseline); }
    return (-1.0)*charge/( config.resistance );
}
Double_t get_baseline( std::vector<Double_t> &vec, CONF &config )
{
    Int_t pre_trigger = config.pre_trigger;
    Double_t baseline = 0;
    for ( Int_t i = 0 ; i < pre_trigger ; i++ )
    {
        baseline += vec[i];
    }
    return (Double_t)baseline/(Double_t)pre_trigger;
}
void fit_function(TFile* rootfile,int channel_classifier, std::vector<Double_t> & results, DESC &description)
{
    rootfile->cd();

    std::string spectra_name = "Area_Spectra_Ch" + std::to_string(channel_classifier);
    std::string canvas_name = description.date+"_Canvas_Spectra_Ch" + std::to_string(channel_classifier);
    std::string pdf_name = "/unix/nemo4/PMT_He_Study_nemo4/ROOT_files/Area_Spectrum_Files/" +
                                std::to_string(description.voltage)+"V/PDFs/" + description.date +
                                "_Canvas_Spectra_Ch" + std::to_string(channel_classifier) +".pdf";

    TH1D* charge_spectrum = (TH1D*)rootfile->Get(spectra_name.c_str());

    if (charge_spectrum->GetEntries() == 0)
    {
        delete charge_spectrum;
        return;
    }
    else{

        int bi207_1MeV_peak_position = charge_spectrum->GetMaximumBin()*charge_spectrum->GetBinWidth(2);
        TCanvas c1 = new TCanvas(canvas_name.c_str(),canvas_name.c_str());

        if (channel_classifier == 0)
        {

            TF1 fit("fit",bi_func_ch0,bi207_1MeV_peak_position - 3,bi207_1MeV_peak_position + 5);
            fit.SetParNames("A","mu","sigma");

            fit.SetParLimits(0,0,400);
            fit.SetParLimits(1,bi207_1MeV_peak_position - 2,bi207_1MeV_peak_position + 1);
            fit.SetParLimits(2,0.8,2);
            fit.SetParameters(319,bi207_1MeV_peak_position,1.09);

            charge_spectrum->Fit("fit","","",bi207_1MeV_peak_position - 3,bi207_1MeV_peak_position + 5);

            Double_t chi2       = fit.GetChisquare()/fit.GetNDF();
            Double_t mu         = fit.GetParameter(1);
            Double_t mu_err     = fit.GetParError(1);
            Double_t sigma      = fit.GetParameter(2);
            Double_t sigma_err  = fit.GetParError(2);

            results.push_back(mu);
            results.push_back(mu_err);
            results.push_back(sigma);
            results.push_back(sigma_err);
            results.push_back(chi2);
        }
        else if (channel_classifier == 1)
        {

            TF1 fit("fit",bi_func_ch1,bi207_1MeV_peak_position - 2,bi207_1MeV_peak_position + 5);
            fit.SetParNames("A","mu","sigma");

            fit.SetParLimits(0,0,400);
            fit.SetParLimits(1,bi207_1MeV_peak_position - 2,bi207_1MeV_peak_position + 1);
            fit.SetParLimits(2,0.8,2);

            fit.SetParameters(319,bi207_1MeV_peak_position,1.09);

            charge_spectrum->Fit("fit","","",bi207_1MeV_peak_position - 2,bi207_1MeV_peak_position + 5);

            Double_t chi2       = fit.GetChisquare()/fit.GetNDF();
            Double_t mu         = fit.GetParameter(1);
            Double_t mu_err     = fit.GetParError(1);
            Double_t sigma      = fit.GetParameter(2);
            Double_t sigma_err  = fit.GetParError(2);

            results.push_back(mu);
            results.push_back(mu_err);
            results.push_back(sigma);
            results.push_back(sigma_err);
            results.push_back(chi2);
        }

        charge_spectrum->SetXTitle("Charge /pC");
        charge_spectrum->SetYTitle("Counts");
        charge_spectrum->SetTitle("Bi Integrated Charge Spectrum");

        charge_spectrum->Draw();
        c1->SetGrid();
        c1->Update();
        c1->Draw();
        gStyle->SetOptFit(1);
        gStyle->SetStatY(0.9);
        gStyle->SetStatX(0.9);
        gStyle->SetStatW(0.8);
        gStyle->SetStatH(0.1);
        c1->SaveAs(pdf_name.c_str(),"pdf");
        delete charge_spectrum;
        delete c1;
    }
    return;
}
Double_t bi_func_ch1(Double_t *x, Double_t *par)
{
    Float_t xx =x[0];
    Double_t f = par[0]*(7.08*TMath::Gaus(xx,par[1],par[2]) +
                        1.84*TMath::Gaus(xx,par[1]*(1 + 72.144/975.651),par[2]*1.036) +
                        0.44*TMath::Gaus(xx,par[1]*(1 + 84.154/975.651),par[2]*1.042)) +
                                0.515*(exp(0.2199*xx)/(1 + exp((xx - 31.68)/2.48)));
    return f;
}
Double_t bi_func_ch0(Double_t *x, Double_t *par)
{
    Float_t xx =x[0];
    Double_t f = par[0]*(7.08*TMath::Gaus(xx,par[1],par[2]) +
                        1.84*TMath::Gaus(xx,par[1]*(1 + 72.144/975.651),par[2]*1.036) +
                        0.44*TMath::Gaus(xx,par[1]*(1 + 84.154/975.651),par[2]*1.042)) +
                                0.464*(exp(0.254*xx)/(1 + exp((xx - 28.43)/2.14)));
    return f;
}

