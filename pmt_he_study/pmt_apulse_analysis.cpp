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

namespace po = boost::program_options;


typedef struct {
    Int_t sweep_start, pre_trigger, trigger, trig_tolerance;
    Double_t shape_cut, amp_cut, charge_cut, resistance;
    std::vector<Double_t> integration;
} CONF;

typedef struct {
    std::string date, timestamp, prefix, template_file, input_file, output_file, config_file;
    Int_t voltage, tot_event_ch0, tot_event_ch1;
} DESC;

typedef struct {
    Int_t apulse_num;
    std::vector<Int_t> apulse_times;
    std::vector<Double_t> apulse_amplitudes;
    std::vector<Double_t> mf_shapes;
    std::vector<Double_t> mf_amps;
} MATCHFILTER;

std::vector<Double_t> process_line( const std::string & s, char delimiter );
Double_t get_baseline( std::vector<Double_t> &vec, CONF &config );
Double_t get_charge( std::vector<Double_t> &vec, Double_t baseline, CONF &config, Int_t peak_cell );
Double_t get_apulse_charge( std::vector<Double_t> &vec, Double_t baseline, CONF &config );
Double_t get_amplitude( std::vector<Double_t> &vec, Double_t baseline );
Int_t get_peak_cell( std::vector<Double_t> &vec );
std::vector<std::string> split( const std::string& s, char delimiter );
CONF read_config( std::string filename );
DESC process_name( std::string &s );
MATCHFILTER sweep( std::vector<Double_t> &vec, CONF &config, Int_t channel , Double_t baseline, std::vector<Double_t>& temp );
std::vector<std::vector<Double_t>> get_template_pulses( DESC &description );
Double_t get_inner_product( std::vector<Double_t> &vec1, std::vector<Double_t> &vec2 );

Int_t main(Int_t argc, char* argv[])
{
    std::string input_file;
    std::string config_file;
    std::string output_file;
    std::string template_file;

    po::options_description desc("Allowed options"); desc.add_options()
            ("help", "produce help message")
            ("i", po::value(&input_file), "Name of input file")
            ("c", po::value(&config_file), "Name of configuration file")
            ("o", po::value(&output_file), "Name of output file")
            ("t", po::value(&template_file), "Name of template file");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 0;
    }

    // Process the arguements
    std::cout << std::endl;

    const std::string tag_ch0 = "trace channel=\"0\"";
    const std::string tag_ch1 = "trace channel=\"1\"";

    // ==================================================================
    // Get the relevant strings
    DESC description            = process_name( input_file );
    description.input_file      = input_file;
    description.output_file     = output_file;
    description.config_file     = config_file;
    description.template_file   = template_file;

    std::cout << ">>> File process      : " << description.input_file << std::endl;
    std::cout << ">>> Data file prefix  : " << description.prefix << std::endl;
    std::cout << ">>> Output filename   : " << description.output_file << std::endl;
    std::cout << ">>> Config filename   : " << description.config_file << std::endl;
    std::cout << ">>> Template filename : " << description.output_file << std::endl;
    std::cout << ">>> Time-stamp        : " << description.timestamp << std::endl;
    std::cout << ">>> Date              : " << description.date << std::endl;
    std::cout << ">>> Voltage           : " << description.voltage << std::endl;

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
    results_files.push_back("/unix/nemo4/PMT_He_Study_nemo4/Results/" + std::to_string(description.voltage) + "/apulseNUM_vs_date_Ch0.dat");
    results_files.push_back("/unix/nemo4/PMT_He_Study_nemo4/Results/" + std::to_string(description.voltage) + "/apulseNUM_vs_date_Ch1.dat");

    TFile root_file( output_file.c_str(), "RECREATE" );

    std::vector<std::vector<Double_t>> template_vectors = get_template_pulses( description );
    // ==================================================================

    // Create a container to hold all the waveform information
    // this will be the construct to fill Int_to the TTree
    typedef struct {
        Int_t OM_ID;
        Double_t apulse_charge, pulse_charge ,pulse_amplitude, baseline;
    } EVENTN;
    EVENTN eventn;

    MATCHFILTER matchfilter;

    TTree tree("T","Tree containing waveform information");
    tree.Branch("apulse_times",&matchfilter.apulse_times);
    tree.Branch("apulse_amplitudes",&matchfilter.apulse_times);
    tree.Branch("mf_shapes",&matchfilter.mf_shapes);
    tree.Branch("mf_amps",&matchfilter.mf_amps);
    tree.Branch("apulse_num",&matchfilter.apulse_num);
    tree.Branch("apulse_charge",&eventn.apulse_charge);
    tree.Branch("pulse_charge",&eventn.pulse_charge);
    tree.Branch("pulse_baseline",&eventn.baseline);
    tree.Branch("pulse_amplitude",&eventn.pulse_amplitude);

    // Create the histograms to hold the information we want

    TH1D area_spectra_hist_0 = TH1D("Area_Spectra_Ch0","Area_Spectra_Ch0",200,0,60);
    TH1D area_spectra_hist_1 = TH1D("Area_Spectra_Ch1","Area_Spectra_Ch1",200,0,60);
    TH1D amp_spectra_hist_0 = TH1D("Amplitude_Spectra_Ch0","Amplitude_Spectra_Ch0",200,0,1000);
    TH1D amp_spectra_hist_1 = TH1D("Amplitude_Spectra_Ch1","Amplitude_Spectra_Ch1",200,0,1000);
    TH1D baseline_hist_0 = TH1D("Baseline_Spectra_Ch0","Baseline_Spectra_Ch0",200,900,1000);
    TH1D baseline_hist_1 = TH1D("Baseline_Spectra_Ch1","Baseline_Spectra_Ch1",200,900,1000);
    TH1D apulse_area_spectra_hist_0 = TH1D("Apulse_Area_Spectra_Ch0","Apulse_Area_Spectra_Ch0",200,0,60);
    TH1D apulse_area_spectra_hist_1 = TH1D("Apulse_Area_Spectra_Ch1","Apulse_Area_Spectra_Ch1",200,0,60);
    std::vector<TH1D> area_hists;
    area_hists.push_back(area_spectra_hist_0);
    area_hists.push_back(area_spectra_hist_1);
    std::vector<TH1D> amp_hists;
    amp_hists.push_back(amp_spectra_hist_0);
    amp_hists.push_back(amp_spectra_hist_1);
    std::vector<TH1D> baseline_hists;
    baseline_hists.push_back(baseline_hist_0);
    baseline_hists.push_back(baseline_hist_1);
    std::vector<TH1D> apulse_area_hists;
    area_hists.push_back(apulse_area_spectra_hist_0);
    area_hists.push_back(apulse_area_spectra_hist_1);

    int nbins{7000-config_object.sweep_start};

    std::vector<std::vector<TH1D>> hist_vector;

    TH1D h_apulse_num_0_hh("h_apulse_num_Ch0_hh","h_apulse_num_Ch0_hh",20,0,20);
    TH1D h_apulse_amplitudes_0_hh("h_apulse_amp_Ch0_hh","h_apulse_amp_Ch0_hh",200,10,100);
    TH1D h_apulse_times_0_hh("h_apulse_times_Ch0_hh","h_apulse_times_Ch0_hh",nbins,(Double_t)config_object.sweep_start,7000);
    std::vector<TH1D*> channel0_vector(h_apulse_num_0_hh, h_apulse_amplitudes_0_hh,h_apulse_times_0_hh);

    TH1D h_apulse_num_1_hh("h_apulse_num_Ch1_hh","h_apulse_num_Ch1_hh",20,0,20);
    TH1D h_apulse_amplitudes_1_hh("h_apulse_amp_Ch1_hh","h_apulse_amp_Ch1_hh",200,10,100);
    TH1D h_apulse_times_1_hh("h_apulse_times_Ch1_hh","h_apulse_times_Ch1_hh",nbins,(Double_t)config_object.sweep_start,7000);
    std::vector<TH1D*> channel1_vector(h_apulse_num_1_hh, h_apulse_amplitudes_1_hh,h_apulse_times_1_hh);

    hist_vector.push_back(channel0_vector);
    hist_vector.push_back(channel1_vector);

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
        } else if (data_line.find(tag_ch1) != std::string::npos) {
            channel_indicator = 1;
        } else { continue; } // All other lines are garbage

        // Process      =================================================
        std::vector <Double_t> data = process_line(data_line, ' ');

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

        if ( pulse_charge < config_object.charge_cut )
        {
            // Pulse is either too small or just noise
            continue;
        }

        channel_waveform_num[channel_indicator]++;

        Double_t pulse_amplitude        = get_amplitude( data, baseline );
        Double_t apulse_charge          = get_apulse_charge( data, baseline, config_object );
        MATCHFILTER temp_match_filter   = sweep( data, config_object, baseline, template_vectors[channel_indicator] );

        // Output       =================================================

        matchfilter.apulse_times        = temp_match_filter.apulse_times;
        matchfilter.apulse_amplitudes   = temp_match_filter.apulse_amplitudes;
        matchfilter.mf_amps             = temp_match_filter.mf_amps;
        matchfilter.mf_shapes           = temp_match_filter.mf_shapes;
        matchfilter.apulse_num          = temp_match_filter.apulse_num;

        eventn.pulse_amplitude  = pulse_amplitude;
        eventn.pulse_charge     = pulse_charge;
        eventn.baseline         = baseline;
        eventn.apulse_charge    = apulse_charge;
        eventn.OM_ID            = channel_indicator;

        for ( Int_t i = 0; i < matchfilter.apulse_amplitudes.size() ; i++ )
        {
            hist_vector[channel_indicator][1].Fill(matchfilter.apulse_amplitudes[i]);
            hist_vector[channel_indicator][1].Fill(matchfilter.apulse_times[i]);
        }
        hist_vector[channel_indicator][0].Fill(matchfilter.apulse_num);
        amp_hists[channel_indicator].Fill(pulse_amplitude);
        area_hists[channel_indicator].Fill(pulse_charge);
        baseline_hists[channel_indicator].Fill(baseline);
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
    // Write all the afterpulse information to the analysis file
    for (int ichannel =0;ichannel<2;ichannel++)
    {
        Int_t results;
        Double_t percentage_apulse;

        std::ofstream apulse_file;
        std::string output_filename = results_file_directory + "apulseNUM_vs_date_Ch" + std::to_string(ichannel) + ".dat";

        if (hist_vector[ichannel][0].GetEntries() > 0)
        {
            for (Int_t ibin = 2; ibin < hist_vector[ichannel][0].GetNbinsX(); ibin++)
            {
                results += (Int_t) hist_vector[ichannel][0].GetBinContent(ibin);
            }
            percentage_apulse = Double_t(results)/Double_t(hist_vector[ichannel][0].GetEntries())*100.0;
        }


        std::cout << "Opening output file... " << output_filename << std::endl;
        apulse_file.open(output_filename, std::ios_base::app);
        std::cout << ">>> Results to be appended to file: "<< std::endl;
        std::cout << ">>> Percentage of pulses that have afterpulses "<< percentage_apulse << std::endl;
        std::cout << ">>> Standard Deviation of the apulse_num plot "<< hist_vector[ichannel][0].GetStdDev() << std::endl;
        std::cout << ">>> " << std::endl;
        apulse_file << description.voltage
                    << ","
                    << TDatime().Get()
                    << ","
                    << description.date
                    << ","
                    << description.timestamp
                    << ","
                    << percentage_apulse
                    << ","
                    << hist_vector[ichannel][0].GetStdDev()
                    << std::endl;

        std::cout << "Closing output file..." << std::endl;

    }
    //              =================================================
    root_file.Write();
    root_file.Close();
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
    Double_t error = sqrt(baseline);
    Int_t gate = config.integration;

    Int_t start,end;
    for ( Int_t i = (Int_t)config.pre_trigger; i < peak_cell ; i++ )
    {
        if ( (vec[i] - baseline) < (vec[peak_cell] - baseline)*config.integration[0] )
        {
            start = i;
            break;
        }else{ continue; }
    }
    for ( Int_t i = peak_cell ; i < vec.size() ; i++ )
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
Double_t get_apulse_charge( std::vector<Double_t> &vec, Double_t baseline, CONF &config )
{
    Double_t apulse_charge = 0;
    for ( Int_t i = config.sweep_start ; i < vec.size() ; i++ )
    {
        apulse_charge += (vec[i] - baseline);
    }
    return (-1.0)*apulse_charge
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
MATCHFILTER sweep( std::vector<Double_t> &vec, CONF &config, Int_t channel , Double_t baseline, std::vector<Double_t>& temp )
{
    MATCHFILTER matchfilter;
    Int_t sweep_start = config.sweep_start;
    Double_t shape_cut = config.shape_cut;
    Double_t amp_cut = config.amp_cut;

    std::vector<Double_t> apulse_amp_vec;
    std::vector<Int_t> apulse_time_vec;
    Int_t current_apulse = 0;
    Int_t previous_apulse = sweep_start - (Int_t)temp[channel].size()*2;

    // Create containers for the shape and amplitude sweeps for storing
    std::vector<Double_t> shape_vec;
    std::vector<Double_t> amp_vec;

    for ( Int_t isweep = sweep_start; isweep < (Int_t)vec.size() - (Int_t)temp[channel].size(); isweep++ )
    {
        std::cout << "In sweep" << std::endl;
        std::vector<Double_t> test;
        for ( Int_t i_vec = 0; i_vec < (Int_t)temp[channel].size(); i_vec++ )
        {
            test.push_back( vec[i_vec] - baseline );
        }
        std::cout << "Got test vector" << std::endl;
        Double_t norm = get_inner_product( test, test );
        std::cout << "Normalised test vector" << std::endl;
        Double_t amplitude = get_inner_product( test, temp[channel]);
        std::cout << "Got amplitude" << std::endl;
        Double_t shape = amplitude/norm;

        shape_vec.push_back(shape);
        amp_vec.push_back(amplitude);

        if ( shape > shape_cut && amplitude > amp_cut) // We have an afterpulse
        {
            Int_t distance_to_nearest_afterpulse = isweep - previous_apulse;

            // Check whether still on same afterpulse by defining a 2*template length temporal distance
            if (distance_to_nearest_afterpulse > (Int_t)temp[channel].size()*2)
            {
                // This is a new afterpulse:
                apulse_amp_vec.push_back( amplitude );
                apulse_time_vec.push_back( isweep );
                previous_apulse = isweep;
                current_apulse = apulse_amp_vec.size()-1;
            } else{
                // We are still analysing the same afterpulse
                if (amplitude > apulse_amp_vec[current_apulse])
                {
                    apulse_amp_vec[current_apulse] = amplitude;
                    apulse_time_vec[current_apulse] = isweep;
                }
            }
        }
    }
    matchfilter.apulse_amplitudes = apulse_amp_vec;
    matchfilter.apulse_times = apulse_time_vec;
    matchfilter.mf_amps = amp_vec;
    matchfilter.mf_shapes = shape_vec;
    matchfilter.apulse_num = (Int_t)apulse_time_vec.size();
    return matchfilter;
}
std::vector<std::vector<Double_t>> get_template_pulses( DESC &description )
{
    std::vector<std::vector<Double_t>> template_pulses;
    TFile temp_root_file{description.template_file.c_str(), "READ"};
    for (Int_t itemp = 0; itemp < 2; itemp++)
    {
        std::vector<Double_t> temp_vector; // Define a temporary filling vector
        //Get the template histogram from the file
        std::string histname = "A" + std::to_string(description.voltage)
                                + "_B" + std::to_string(description.voltage)
                                + "_Ch" + std::to_string(itemp)
                                + "_Template";

        TH1D* template_hist = (TH1D*)temp_root_file.Get(hist_name.c_str());

        for (Int_t ihist = 0; ihist < template_hist->GetEntries(); ihist++)
        {
            if ( template_hist->GetBinContent(ihist) > 0.0 )
            {
                temp_vector.push_back(template_hist->GetBinContent(ihist));
            }
        }
        delete template_hist;
        template_pulses.push_back(temp_vector);
    }
    temp_root_file.Close();

    return template_pulses;
}
Double_t get_inner_product( std::vector<Double_t> &vec1, std::vector<Double_t> &vec2 )
{
    Double_t inner_product;
    for ( Int_t i_vec = 0 ; i_vec < (Int_t)vec1.size() ; i_vec++ )
    {
        inner_product += vec1[i_vec]*vec2[i_vec];
    }
    return inner_product;
}



