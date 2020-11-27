// Standard library:
#include <iostream>
#include <exception>
#include <cstdlib>


#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"


#include <snfee/snfee.h>
#include <snfee/io/multifile_data_reader.h>
#include <snfee/data/raw_trigger_data.h>


// This project:
#include <sncabling/sncabling.h>
#include <sncabling/om_id.h>
#include <sncabling/calo_hv_id.h>
#include <sncabling/calo_hv_cabling.h>
#include <sncabling/label.h>

#include <sncabling/service.h>
#include <sncabling/calo_signal_cabling.h>

typedef struct {
    Int_t OM_ID, row, col, wall, ID;
} EVENTN;

typedef struct {
    Int_t low_edge = 30;
    Int_t high_edge = 100;
    Int_t temp_length = 130;
} TEMP_INFO;

Double_t get_inner_product( std::vector<Double_t> &vec1, std::vector<Double_t> &vec2 );
std::vector<std::vector<Double_t>> get_template_pulses( std::string template_file , Int_t n_temp );
void update_temp_vector( std::vector<std::vector<Double_t>> &template_vectors, std::vector<Double_t> new_vector, Int_t OM_ID );
Int_t get_peak_cell( std::vector<Double_t> &vec );
void write_templates( std::vector<std::vector<Double_t>> &template_vectors );
Double_t get_baseline( std::vector<Double_t> &vec );
Int_t get_max_value( std::vector<Double_t> &vec );
void draw_waveform( std::vector<Double_t> &vec, Int_t n_samples, Double_t baseline, EVENTN &eventn, std::string output_directory);
void draw_pulse( std::vector<Double_t> &temp, std::vector<Double_t> &test, Int_t i, Double_t convo, Double_t sample_time, EVENTN &eventn);
void save_hist( std::vector<Double_t> &vec, std::string x_label, std::string y_label, std::string title, std::string file_name, Int_t n_bins, Double_t min_bin, Double_t max_bin, TFile* root_file);
Double_t get_pulse_time_mf(std::vector<Double_t> &vec);
std::vector<Double_t> read_energy_coef( std::string filename );
std::vector<std::string> split( const std::string& s, char delimiter );

bool debug = true;


void usage()
{
    std::clog<<std::endl;
    std::clog<<"+--------------------------------------------------+"<<std::endl;
    std::clog<<"| SuperNEMO calorimeter commissioning tutorial lv0 |"<<std::endl;
    std::clog<<"+--------------------------------------------------+"<<std::endl;

    std::clog<<">>> How to use: "<<std::endl;
    std::clog<<">>> -help "<<std::endl;
    std::clog<<">>> -i  -  std::string /input/file/path/.gz "<<std::endl;
    std::clog<<">>> -o  -  std::string /output/file/path/.root "<<std::endl;
    std::clog<<">>> -t  -  BOOL create template, def:false "<<std::endl;
    std::clog<<">>> -OM -  INT def: 1000 (no plots), chosen OM to plot examples "<<std::endl;
    std::clog<<">>> -W  -  BOOL analyse waveforms def:false "<<std::endl;
    std::clog<<std::endl;
}

// Main program

int main(int argc, char **argv)
{
    sncabling::initialize();
    int error_code = EXIT_SUCCESS;

    std::string input_file_name, output_file_name;
    bool do_template = false;
    int chosen_OM = 1000;
    bool do_waveforms = false;

    try {
    
        if (argc > 0)
        {
            for (int i{1}; i < argc;  i++)
	        {
                std::string s{argv[i]};

                if ( s == "-help" )
                {
                    usage();
                    return 0;
                }
                else if ( s == "-i" )
                {
                    input_file_name = std::string(argv[i+1]);
                }
                else if ( s == "-o" )
                {
                    output_file_name = std::string(argv[i+1]);
                }
                else if ( s == "-t" )
                {
                    if ( std::string(argv[i+1]) == "true" )
                    {
                        do_template = true;
                    }else{
                        do_template = false;
                    }
                }
                else if ( s == "-OM" )
                {
                    chosen_OM = std::stoi(argv[i+1]);
                }
                else if ( s == "-W" )
                {
                    if ( std::string(argv[i+1]) == "true" )
                    {
                        do_waveforms = true;
                    }else{
                        do_waveforms = false;
                    }
                }
	        }
        }
    
        if (input_file_name.length() < 1)
        {
	        std::clog<<"Invalid input file"<<std::endl;
	        return 0;
        }

        std::clog<<"Input file name : "<<input_file_name<<std::endl;

        std::vector<Double_t> energy_coefs = read_energy_coef("Resultat_mean_energie_charge_france_435.txt");

        TEMP_INFO template_info;
        std::vector<std::vector<Double_t>> template_vectors;
        if ( do_template )
        {
            for (int k = 0; k < 260; ++k)
            {
                std::vector<Double_t> temp(template_info.temp_length, 0.0);
                template_vectors.push_back(temp);
            }

            std::cout << "Initialise template vectors" << std::endl;
            for (int j = 0; j < (Int_t)template_vectors.size(); ++j)
            {
                std::cout << std::endl;
                std::cout << "Template " << j << std::endl;
                for (int i = 0; i < (Int_t)template_vectors[j].size(); ++i)
                {
                    std::cout << "( " << i << " , " << template_vectors[j][i] << " )" << std::endl;
                }
            }

        } else {
            template_vectors = get_template_pulses( "templates.root", 260 );
        }

        sncabling::service snCabling;
        snCabling.initialize_simple();
        // Access to the calorimeter signal readout cabling map:
        const sncabling::calo_signal_cabling & caloSignalCabling = snCabling.get_calo_signal_cabling();

        TFile* output_file = new TFile(output_file_name.c_str(), "RECREATE");
    
        Int_t event_num;
        Int_t row;
        Int_t column;
        Int_t OM_ID;
        Int_t charge;
        Int_t amplitude;
        Int_t baseline;
        Double_t fall_time;
        Double_t rise_time;
        Double_t peak_time;
        Int_t calo_hit_num;
        unsigned long long calo_tdc;
        Int_t run_num;
        Int_t wall_num;
        Int_t trig_id;
        Double_t calo_time;
        std::vector<uint16_t> waveform;
        Double_t energy;

        // Create a ROOT Tree
        TTree tree("T","Tree containing simulated vertex data");
        tree.Branch("event_num",&event_num);
        tree.Branch("row",&row);
        tree.Branch("column",&column);
        tree.Branch("OM_ID",&OM_ID);
        tree.Branch("charge",&charge);
        tree.Branch("baseline",&baseline);
        tree.Branch("amplitude",&amplitude);
        tree.Branch("rise_time",&rise_time);
        tree.Branch("fall_time",&fall_time);
        tree.Branch("peak_time",&peak_time);
        tree.Branch("calo_hit_num",&calo_hit_num);
        tree.Branch("calo_tdc",&calo_tdc);
        tree.Branch("run_num",&run_num);
        tree.Branch("wall_num",&wall_num);
        tree.Branch("trig_id",&trig_id);
        tree.Branch("calo_time",&calo_time);
        tree.Branch("energy",&energy);
        //tree.Branch("waveform",&waveform);

        bool cont = true;

        // Configuration for raw data reader
        snfee::io::multifile_data_reader::config_type reader_cfg;
        reader_cfg.filenames.push_back(input_file_name);

        // Instantiate a reader:
        snfee::io::multifile_data_reader rtd_source(reader_cfg);

        // Working RTD object --> Raw Trigger Data
        // 1 record per trigger composed by few CaloHit
        snfee::data::raw_trigger_data rtd;

        event_num = 0;

        EVENTN eventn;

        std::size_t rtd_counter = 0;
        while ( rtd_source.has_record_tag() )
        {
            rtd_counter++;
      
            // Load the next RTD object:
            rtd_source.load(rtd);
            // General informations:
            int32_t trigger_id = rtd.get_trigger_id();
            int32_t run_id     = rtd.get_run_id();
      
            if(rtd_counter %10000 == 0 )std::clog<<"In Run : "<<run_id<<" Trigger # "<<trigger_id <<std::endl;
      
            std::size_t calo_counter = 0;
            // Loop on calo hit records in the RTD data object:
            for (const auto & p_calo_hit : rtd.get_calo_hits())
            {
	            // Dereference the stored shared pointer oin the calo hit record:
	            const snfee::data::calo_hit_record & calo_hit = *p_calo_hit;
	            calo_counter++;
	            uint64_t tdc             = calo_hit.get_tdc();        // TDC timestamp (48 bits)
	            int32_t  crate_num       = calo_hit.get_crate_num();  // Crate number (0,1,2)
	            int32_t  board_num       = calo_hit.get_board_num();  // Board number (0-19)
	            //if (board_num >= 10){ board_num++; };                 // convert board_num  from [10-19] to [11-20]
	            int32_t  chip_num        = calo_hit.get_chip_num();   // Chip number (0-7)
	            auto     hit_num         = calo_hit.get_hit_num();

                /*
                if(rtd_counter < 100 ){
                  std::clog<<"   |-> tdc      : "<< tdc <<std::endl;
                  std::clog<<"   |-> calo data from CaloFEB : "<<crate_num<<"."<<board_num<<"."<<chip_num<<std::endl;
                }
                */

	            // Extract SAMLONG channels' data:
	            // 2 channels per SAMLONG
	            for (int ichannel = 0; ichannel < snfee::model::feb_constants::SAMLONG_NUMBER_OF_CHANNELS; ichannel++)
	            {
	                const snfee::data::calo_hit_record::channel_data_record & ch_data = calo_hit.get_channel_data(ichannel);
	                bool    ch_lt           {ch_data.is_lt()};            // Low threshold flag
	                bool    ch_ht           {ch_data.is_ht()};            // High threshold flag
	                int32_t ch_baseline     {ch_data.get_baseline()};     // Computed baseline       (LSB: ADC unit/16)
	                int32_t ch_peak         {ch_data.get_peak()};         // Computed peak amplitude (LSB: ADC unit/8)
	                int32_t ch_peak_cell    {ch_data.get_peak_cell()};    // Computed peak cell
	                int32_t ch_charge       {ch_data.get_charge()};       // Computed charge
	                int32_t ch_rising_cell  {ch_data.get_rising_cell()};  // Computed rising cell
	                int32_t ch_falling_cell {ch_data.get_falling_cell()}; // Computed falling cell

	                Double_t ch_rising_cell_  = Double_t(ch_rising_cell);
	                Double_t ch_falling_cell_ = Double_t(ch_falling_cell);
	                Double_t ch_peak_cell_    = Double_t(ch_peak_cell);

	                Double_t rising_actual    = (ch_rising_cell_*6.25)/256.0;
	                Double_t falling_actual   = (ch_falling_cell_*6.25)/256.0;
	                Double_t peak_actual      = ch_peak_cell_*6.25/8.0;

	                sncabling::calo_signal_id readout_id(sncabling::CALOSIGNAL_CHANNEL,
	                        crate_num, board_num,
	                        snfee::model::feb_constants::SAMLONG_NUMBER_OF_CHANNELS * chip_num + ichannel);

	  
	                if (caloSignalCabling.has_channel(readout_id))
	                {
	                    if (ch_ht)
	                    {
		                    const sncabling::om_id & calo_id = caloSignalCabling.get_om(readout_id);
                            row = calo_id.get_row();
                            column = calo_id.get_column();
                            OM_ID = row + column*13;
                            Double_t energy_t = -1.0 * (Double_t)ch_charge * energy_coefs[OM_ID];

                            if (energy_t > 0.7)
                            {
                                amplitude = ch_peak;
                                baseline  = ch_baseline;
                                charge    = ch_charge;
                                rise_time = rising_actual;
                                fall_time = falling_actual;
                                peak_time = peak_actual;
                                calo_hit_num = hit_num;
                                calo_tdc = tdc;
                                run_num = run_id;
                                wall_num = crate_num;
                                trig_id = trigger_id;
                                energy = energy_t;

                                eventn.OM_ID = OM_ID;
                                eventn.col = column;
                                eventn.row = row;
                                eventn.wall = crate_num;
                                eventn.ID = event_num;

                                calo_time = -1000.0;

                                if (do_waveforms)
                                {
                                    if ( ch_peak_cell > 1024 - 160){ continue; }
                                    if ( chosen_OM == 1000 ){} else if ( OM_ID != chosen_OM ){ continue; }

                                    uint16_t waveform_number_of_samples = calo_hit.get_waveform_number_of_samples();
                                    std::vector<Double_t> waveform_adc;
                                    for (uint16_t isample = 0; isample < waveform_number_of_samples; isample++)
                                    {
                                        uint16_t adc = calo_hit.get_waveforms().get_adc(isample,ichannel);
                                        waveform_adc.push_back((Double_t)adc);
                                    }
                                    Double_t my_baseline = get_baseline( waveform_adc );

                                    if (do_template)
                                    {
                                        if ( charge >= -25000 && ch_charge < -20000 )
                                        {
                                            std::vector<Double_t> temp_vector;
                                            for (uint16_t isample = 0; isample < waveform_number_of_samples; isample++)
                                            {
                                                temp_vector.push_back( waveform_adc[isample] - baseline );
                                            }
                                            update_temp_vector( template_vectors, temp_vector, OM_ID );
                                        }

                                    } else {
                                        Int_t n_try = 60;
                                        std::vector<Double_t> mf_output;
                                        Double_t norm_temp = sqrt( get_inner_product( template_vectors[OM_ID], template_vectors[OM_ID] ) );

                                        for (int i = 0; i < n_try; ++i)
                                        {
                                            std::vector<Double_t> temp_vector;

                                            for (uint16_t isample = ch_peak_cell - 30 - n_try/2 + i; isample < ch_peak_cell + 100 - n_try/2 + i; isample++)
                                            {
                                                Double_t volts = (Double_t)waveform_adc[isample] - my_baseline;
                                                temp_vector.push_back( volts );
                                            }

                                            Double_t norm_test = sqrt( get_inner_product( temp_vector, temp_vector ) );
                                            Double_t mf = get_inner_product( temp_vector, template_vectors[OM_ID] )/( norm_temp * norm_test );

                                            mf_output.push_back(mf);

                                            if ( chosen_OM == 1000 ){} else if (OM_ID == chosen_OM)
                                            {
                                                draw_pulse(template_vectors[OM_ID], temp_vector, i, mf,
                                                           (ch_peak_cell - 30 - n_try/2 + i)/2.56, eventn);
                                            }
                                        }
                                        /*save_hist(mf_output,"Sample window time /ns","shape index","Pulse_time_finder",
                                                    "mf_output_" + std::to_string(OM_ID) + ".png",n_try,
                                                    (ch_peak_cell - 30 - n_try/2)/2.56 ,
                                                    (ch_peak_cell - 30 + n_try/2)/2.56, output_file);*/

                                        if ( chosen_OM == 1000 ){} else if (OM_ID == chosen_OM)
                                        {
                                            draw_waveform(waveform_adc, waveform_number_of_samples, my_baseline, eventn, output_file_name);
                                            return 1;
                                        }

                                        if (mf_output.size() > 0)
                                        {
                                            calo_time = get_pulse_time_mf(mf_output) + (Double_t)ch_peak_cell - 30.0 - (Double_t)n_try / 2.0;
                                            // std::cout << "calo_time: " << calo_time << std::endl;
                                        }
                                        // Int_t other_calo_time = get_max_value(mf_output) + ch_peak_cell - 30 - n_try/2;
                                        // std::cout << "int_calo_time: " << other_calo_time << std::endl;
                                    }
                                    //waveform = temp_vector;

                                } else {}
                                tree.Fill();
                            }
                        }
	                }
	            } //end of channels
            }//end of calohit
            event_num ++;
        }   //end of file
    
        std::clog<<"Events processed : " << rtd_counter<< " entries" << std::endl;
        output_file->cd();
        output_file->Write();
        output_file->Close();

        //write_templates( template_vectors );

    } catch (std::exception & error)
    {
            std::cerr << "[error] " << error.what() << std::endl;
            error_code = EXIT_FAILURE;
    }

    sncabling::terminate();
    return error_code;
}

std::vector<std::vector<Double_t>> get_template_pulses( std::string template_file , Int_t n_temp )
{
    std::vector<std::vector<Double_t>> template_pulses;
    std::cout << std::endl;
    std::cout << "Template file name : " << template_file << std::endl;
    TFile temp_root_file(template_file.c_str(), "READ");
    for (Int_t itemp = 0; itemp < n_temp; itemp++)
    {
        //std::cout << "Template: " << itemp << std::endl;
        if ( itemp == 83 || itemp == 109 || itemp == 201 )
        {
            std::vector<Double_t> temp_vector(130, 0.0);
            template_pulses.push_back(temp_vector);
            continue;
        }
        std::vector<Double_t> temp_vector; // Define a temporary filling vector
        //Get the template histogram from the file
        std::string hist_name = "Template_Ch" + std::to_string(itemp);

        TH1D* template_hist = (TH1D*)temp_root_file.Get(hist_name.c_str());

        for (Int_t ihist = 1; ihist < template_hist->GetEntries() + 1; ihist++)
        {
            temp_vector.push_back(template_hist->GetBinContent(ihist));
            //std::cout << ihist << " : " << temp_vector[ihist-1] << std::endl;
        }
        delete template_hist;
        Double_t norm = sqrt( get_inner_product( temp_vector, temp_vector ) );

        if (norm <= 0)
        {
            std::cout << "Error: Abnormal template pulse" << std::endl;
            exit(1);
        }

        for (int ivec = 0 ; ivec < (Int_t)temp_vector.size() ;  ivec++)
        {
            temp_vector[ivec] = temp_vector[ivec]/norm;
            //std::cout << ivec << " : " << temp_vector[ivec] << std::endl;
        }
        //std::cout << std::endl;
        template_pulses.push_back(temp_vector);
    }
    temp_root_file.Close();

    std::cout << "Success..." << std::endl;

    return template_pulses;
}
Double_t get_inner_product( std::vector<Double_t> &vec1, std::vector<Double_t> &vec2 )
{
    if ( vec1.size() != vec2.size() )
    {
        std::cout << ">>> Length of vectors must be the same for an inner product to be calculated" << std::endl;
        std::cout << ">>> Length of vec1: " << vec1.size() << " != length of vec2: " << vec2.size() << std::endl;
        exit(1);
    }

    Double_t inner_product = 0;
    for ( Int_t i_vec = 0 ; i_vec < (Int_t)vec1.size() ; i_vec++ )
    {
        inner_product += vec1[i_vec]*vec2[i_vec];
    }
    return inner_product;
}
void update_temp_vector( std::vector<std::vector<Double_t>> &template_vectors, std::vector<Double_t> new_vector, Int_t OM_ID )
{
    Int_t temp_length = 130;
    Int_t peak_cell = get_peak_cell( new_vector );

    Int_t lower_edge = 30;
    Int_t higher_edge = 100;

    Double_t my_baseline = get_baseline( new_vector );

    Int_t j = 0;
    for (Int_t i = peak_cell - lower_edge; i < peak_cell + higher_edge; ++i)
    {
        template_vectors[OM_ID][j] += new_vector[i] - my_baseline;
        j++;

        if ( j == temp_length )
        {
            break;
        }
    }
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
void write_templates( std::vector<std::vector<Double_t>> &template_vectors )
{
    TFile* template_root_file = new TFile("templates.root", "RECREATE");
    template_root_file->cd();

    for (Int_t i_temp = 0; i_temp < (Int_t)template_vectors.size(); ++i_temp)
    {
        std::string name = "Template_Ch" + std::to_string(i_temp);
        std::cout << std::endl;
        std::cout << name << std::endl;
        TH1D* hist = new TH1D(name.c_str(), name.c_str(), template_vectors[i_temp].size(), 0, template_vectors[i_temp].size());

        Double_t norm = sqrt( get_inner_product( template_vectors[i_temp], template_vectors[i_temp] ) );
        if ( norm == 0 )
        {
            std::cout << ">>> Normalised template vector : " << i_temp << " is 0" << std::endl;
            continue;
        }

        for (int j_bin = 1; j_bin < (Int_t)template_vectors[i_temp].size() + 1; ++j_bin)
        {
            hist->SetBinContent(j_bin, template_vectors[i_temp][j_bin - 1]/norm);
            std::cout << j_bin - 1 << " " << template_vectors[i_temp][j_bin - 1]/norm << std::endl;
        }
        hist->Write();
        delete hist;
    }
}
Double_t get_baseline( std::vector<Double_t> &vec )
{
    Double_t baseline = 0;
    for ( Int_t i = 0 ; i < 20 ; i++ )
    {
        baseline += vec[i];
    }
    return (Double_t)baseline/20.0;
}
Int_t get_max_value( std::vector<Double_t> &vec )
{
    Double_t temp = 0.0;
    Int_t pos = 0;
    for (int i = 0; i < (Int_t)vec.size(); ++i)
    {
        if ( vec[i] > temp )
        {
            temp = vec[i];
            pos = i;
        }
    }
    return pos;
}
void draw_waveform( std::vector<Double_t> &vec, Int_t n_samples,
        Double_t baseline, EVENTN &eventn, std::string output_directory)
{
    TCanvas* waveform_canvas = new TCanvas();
    waveform_canvas->cd();
    gStyle->SetOptStat(0);
    TH1D* waveform_hist = new TH1D("waveform", "waveform", n_samples, 0, n_samples/2.56);

    for (uint16_t i_sample = 0; i_sample < n_samples; i_sample++)
    {
        Double_t volts = vec[i_sample] - baseline;
        waveform_hist->SetBinContent(i_sample+1, volts/2.048);
    }

    waveform_hist->SetXTitle("Time stamp /ns");
    waveform_hist->SetYTitle("Voltage /mV");
    std::string w_title = std::to_string(eventn.wall) + ":" + std::to_string(eventn.col) +
                          ":" + std::to_string(eventn.row) + " Waveform";
    waveform_hist->SetTitle(w_title.c_str());
    std::string can_name = output_directory + "_" + std::to_string(eventn.ID) + "_" +
                           std::to_string(eventn.wall) + "_" + std::to_string(eventn.col) +
                           "_" + std::to_string(eventn.row) + "_waveform.png";

    waveform_hist->Draw();
    waveform_canvas->SetGrid(true);
    waveform_canvas->Update();
    waveform_canvas->SaveAs(can_name.c_str());

    delete waveform_hist;
    delete waveform_canvas;
}
void draw_pulse( std::vector<Double_t> &temp, std::vector<Double_t> &test, Int_t i,
        Double_t convo, Double_t sample_time, EVENTN &eventn)
{
    TCanvas* my_canvas = new TCanvas();
    my_canvas->cd();

    gStyle->SetOptStat(0);
    std::string test_name = "test_hist_" + std::to_string(i);
    std::string temp_name = "temp_hist_" + std::to_string(eventn.OM_ID);

    TH1D* test_hist = new TH1D(test_name.c_str(), test_name.c_str(), (Int_t)test.size(), 0, (Double_t)test.size()/2.56);
    TH1D* temp_hist = new TH1D(temp_name.c_str(), temp_name.c_str(), (Int_t)test.size(), 0, (Double_t)test.size()/2.56);
    test_hist->SetMaximum(10);
    test_hist->SetMinimum(test_hist->GetMinimum() - 50);

    TLegend* legend = new TLegend(0.7, 0.1, 0.9, 0.2);
    //gStyle->SetLegendBorderSize(0);
    for (int j = 1; j <= (Int_t)test.size(); ++j)
    {
        test_hist->SetBinContent(j, test[j-1]/2.048);
    }
    for (int k = 1; k <= (Int_t)temp.size(); ++k)
    {
        temp_hist->SetBinContent(k, temp[k-1]);
    }

    temp_hist->Scale(test_hist->Integral()/temp_hist->Integral());
    //temp_hist->Sumw2();

    temp_hist->SetLineColor(2);
    test_hist->SetLineColor(1);

    test_hist->SetXTitle("Relative time /ns");
    test_hist->SetYTitle("Voltage /mV");
    std::string title = std::to_string(eventn.wall) + ":" + std::to_string(eventn.col) + ":" + std::to_string(eventn.row)
            + " MF: " + std::to_string(convo) + " FBT:" + std::to_string(sample_time) + " ns";
    test_hist->SetTitle(title.c_str());

    std::string can_name = "mf_output_" + std::to_string(eventn.wall) + "_" + std::to_string(eventn.col) + "_" +
            std::to_string(eventn.row) + "_N" + std::to_string(i) + ".png";
    test_hist->Draw("HIST");
    temp_hist->Draw("HIST SAME C");
    legend->AddEntry(test_hist, "test");
    legend->AddEntry(temp_hist, "template");
    legend->Draw();
    my_canvas->SetGrid(true);
    my_canvas->Update();
    my_canvas->SaveAs(can_name.c_str());

    delete test_hist;
    delete temp_hist;
    delete my_canvas;
    delete legend;
}
void save_hist( std::vector<Double_t> &vec, std::string x_label, std::string y_label, std::string title,
        std::string file_name, Int_t n_bins, Double_t min_bin, Double_t max_bin, TFile* root_file)
{
    TCanvas* new_canvas = new TCanvas();
    new_canvas->cd();
    new_canvas->SetGrid(true);

    TH1D* new_hist = new TH1D(title.c_str(), title.c_str(), n_bins, min_bin, max_bin);

    for (int k = 0; k < vec.size(); ++k) {
        new_hist->SetBinContent(k+1, vec[k]);
    }

    new_hist->SetXTitle(x_label.c_str());
    new_hist->SetYTitle(y_label.c_str());
    new_hist->SetTitle(title.c_str());

    gStyle->SetOptStat(0);
    new_hist->Draw("HIST");
    new_hist->Write();
    new_canvas->SaveAs(file_name.c_str());

    delete new_canvas;
}
Double_t get_pulse_time_mf(std::vector<Double_t> &vec)
{
    Int_t guess_mean = get_max_value(vec);
    Int_t lower_bound = guess_mean-5;
    Int_t upper_bound = guess_mean+5;
    TH1D* hist = new TH1D("h","h",upper_bound-lower_bound, lower_bound, upper_bound);

    for (int i = 0; i < hist->GetNbinsX(); ++i)
    {
        hist->SetBinContent(i+1, vec[lower_bound + i]);
    }
    /* std::cout << "low: " << lower_bound << std::endl;
    std::cout << "hig: " << upper_bound << std::endl;
    std::cout << "NBins: " << hist->GetNbinsX() << std::endl; */

    TF1 fit("fit","[0]*TMath::Gaus(x,[1],[2])",lower_bound,upper_bound);
    fit.SetParNames("A","mu","sigma");

    fit.SetParLimits(0,0,10);
    fit.SetParLimits(1,guess_mean-1,guess_mean+1);
    fit.SetParLimits(2,0,10);
    fit.SetParameters(1,guess_mean,1);

    hist->Fit("fit","Q");

    Double_t chi2       = fit.GetChisquare()/fit.GetNDF();
    Double_t mu         = fit.GetParameter(1);
    Double_t mu_err     = fit.GetParError(1);
    Double_t sigma      = fit.GetParameter(2);
    Double_t sigma_err  = fit.GetParError(2);

    delete hist;

    return mu;
}
std::vector<Double_t> read_energy_coef( std::string filename )
{
    std::ifstream file(filename);
    std::string line;

    std::vector<Double_t> vec;
    for (int i = 0; i < 260 ; ++i)
    {
        vec.push_back(0.0);
    }

    if (!file.good())
    {
        std::cout << "<<<ERROR>>> cannot open configuration file : " << filename << std::endl;
        std::cout << "EXIT" << std::endl;
        exit(1);
    }

    while ( std::getline( file, line ) && !file.eof() )
    {
        if  (line.empty())
        {
            // Empty line so ignore
            continue;
        }else{
            std::vector<std::string> line_vec = split( line, '\t' );

            Int_t OM = std::stoi(line_vec[0]) - 260;
            Double_t coef = std::stod(line_vec[0]);
            std::cout << OM << " " << coef << std::endl;
            vec[OM] = coef;
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

