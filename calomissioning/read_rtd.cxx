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


Double_t get_inner_product( std::vector<Double_t> &vec1, std::vector<Double_t> &vec2 );
std::vector<std::vector<Double_t>> get_template_pulses( std::string template_file , Int_t n_temp );
void update_temp_vector( std::vector<std::vector<Double_t>> &template_vectors, std::vector<Double_t> new_vector, Int_t OM_ID );
Int_t get_peak_cell( std::vector<Double_t> &vec );
void write_templates( std::vector<std::vector<Double_t>> &template_vectors );
Double_t get_baseline( std::vector<Double_t> &vec );
Int_t get_max_value( std::vector<Double_t> &vec );

bool debug = true;


void usage(){

  std::clog<<std::endl;
  std::clog<<"+--------------------------------------------------+"<<std::endl;
  std::clog<<"| SuperNEMO calorimeter commissioning tutorial lv0 |"<<std::endl;
  std::clog<<"+--------------------------------------------------+"<<std::endl;

  std::clog<<"How to : "<<std::endl;
  std::clog<<" "<<std::endl;
  std::clog<<std::endl;


}

// Main program

int main(int argc, char **argv)
{
    sncabling::initialize();
    int error_code = EXIT_SUCCESS;

    std::string input_file_name, output_file_name;

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
	        }
        }
    
        if (input_file_name.length() < 1)
        {
	        std::clog<<"Invalid input file"<<std::endl;
	        return 0;
        }

        std::clog<<"Input file name : "<<input_file_name<<std::endl;

        //std::vector<Double_t> temp(130, 0.0);
        //std::vector<std::vector<Double_t>> template_vectors(260, temp);
        std::vector<std::vector<Double_t>> template_vectors = get_template_pulses( "templates.root", 260 );
        //std::cout << "Initialise template vectors" << std::endl;
        /*for (int j = 0; j < (Int_t)template_vectors.size(); ++j)
        {
            std::cout << std::endl;
            std::cout << "Template " << j << std::endl;
            for (int i = 0; i < (Int_t)template_vectors[j].size(); ++i)
            {
                std::cout << "( " << i << " , " << template_vectors[j][i] << " )" << std::endl;
            }
        }*/

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
        Int_t calo_tdc;
        Int_t run_num;
        Int_t wall_num;
        Int_t trig_id;
        Int_t calo_time;
        std::vector<uint16_t> waveform;

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
	            if (board_num >= 10){ board_num++; }; // convert board_num  from [10-19] to [11-20]
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

		                    /*
		                    if (column == 2)
		                    {
		                        if (row == 6)
		                        {
			                        std::string name = "M:1.2.6 Example Waveform HT == 0 " + std::to_string(event_num);
			                        TH1I* ex_hist = new TH1I(name.c_str(), name.c_str(), 1024, 0, 400);
			                        ex_hist->SetXTitle("Timestamp /ns");
			                        ex_hist->SetYTitle("ADC counts /mV");
			                        uint16_t waveform_number_of_samples = calo_hit.get_waveform_number_of_samples();
			                        for (uint16_t isample = 0; isample < waveform_number_of_samples; isample++)
			                        {
			                            uint16_t adc = calo_hit.get_waveforms().get_adc(isample,ichannel);
			                            ex_hist->SetBinContent(isample + 1, adc);
			                        }
			                        a_file.cd();
			                        TCanvas* c1 = new TCanvas();
			                        std::string canvas_name = "/sps/nemo/scratch/wquinn/PDFs/ev_" + std::to_string(event_num) + "_2_6.pdf";
                                    c1->cd();
                                    c1->SetGrid();
                                    ex_hist->Draw();
                                    gStyle->SetOptStat(0);
                                    c1->SaveAs(canvas_name.c_str());
                                    delete ex_hist;
                                    delete c1;

		                        }
		                        else if ( row == 7)
		                        {
			                        std::string name = "M:1.2.7 Example Waveform HT == 0 " + std::to_string(event_num);
                                    TH1I* ex_hist = new TH1I(name.c_str(), name.c_str(), 1024, 0, 400);
		                            ex_hist->SetXTitle("Timestamp /ns");
			                        ex_hist->SetYTitle("ADC counts /mV");
                                    uint16_t waveform_number_of_samples = calo_hit.get_waveform_number_of_samples();
                                    for (uint16_t isample = 0; isample < waveform_number_of_samples; isample++)
		                            {
                                        uint16_t adc = calo_hit.get_waveforms().get_adc(isample,ichannel);
                                        ex_hist->SetBinContent(isample + 1, adc);
                                    }
                                    a_file.cd();
                                    TCanvas* c1 = new TCanvas();
			                        std::string canvas_name = "/sps/nemo/scratch/wquinn/PDFs/ev_" + std::to_string(event_num) + "_2_7.pdf";
                                    c1->cd();
                                    c1->SetGrid();
                                    ex_hist->Draw();
                                    gStyle->SetOptStat(0);
                                    c1->SaveAs(canvas_name.c_str());
			                        delete ex_hist;
		                            delete c1;
		                        }
		                    }
		                    */
		                    // Select a small charge range to add to template pulses
		                    //std::cout << "OM_ID: " << OM_ID << " charge: " << charge << std::endl;
                            /*if ( charge >= -25000 && ch_charge < -20000 )
                            {
                                std::vector<Double_t> temp_vector;
                                uint16_t waveform_number_of_samples = calo_hit.get_waveform_number_of_samples();
                                for (uint16_t isample = 0; isample < waveform_number_of_samples; isample++)
                                {
                                    uint16_t adc = calo_hit.get_waveforms().get_adc(isample,ichannel);
                                    temp_vector.push_back( (Double_t)adc - baseline);
                                    //std::cout << isample << " : " << adc - baseline << std::endl;
                                }
                                update_temp_vector( template_vectors, temp_vector, OM_ID );
                            }*/

                            if (OM_ID != 124){ continue; }

                            if ( ch_peak_cell > 1024 - 160)
                            {
                                continue;
                            }

                            Int_t n_try = 50;
                            std::vector<Double_t> mf_output;
                            Double_t norm_temp = sqrt( get_inner_product( template_vectors[OM_ID], template_vectors[OM_ID] ) );

                            uint16_t waveform_number_of_samples = calo_hit.get_waveform_number_of_samples();

                            TCanvas* waveform_canvas = new TCanvas();
                            waveform_canvas->cd();
                            gStyle->SetOptStat(0);
                            TH1D* waveform_hist = new TH1D("waveform", "waveform", waveform_number_of_samples, 0, waveform_number_of_samples/2.56);
                            std::vector<Double_t> waveform_adc;

                            for (uint16_t isample = 0; isample < waveform_number_of_samples; isample++)
                            {
                                uint16_t adc = calo_hit.get_waveforms().get_adc(isample,ichannel);
                                waveform_adc.push_back(adc);
                            }

                            Double_t my_baseline = get_baseline( waveform_adc );

                            for (uint16_t isample = 0; isample < waveform_number_of_samples; isample++)
                            {
                                Double_t volts = (Double_t)waveform_adc[isample] - my_baseline;
                                waveform_hist->SetBinContent(isample+1, volts/2.048);
                            }

                            waveform_hist->SetXTitle("Time stamp /ns");
                            waveform_hist->SetYTitle("Voltage /mV");
                            std::string w_title = "M:1:9:7 Waveform";
                            waveform_hist->SetTitle(w_title.c_str());

                            std::string can_name = "waveform.png";
                            waveform_hist->Draw();
                            waveform_canvas->SetGrid(true);
                            waveform_canvas->Update();
                            waveform_canvas->SaveAs(can_name.c_str());

                            delete waveform_hist;
                            delete waveform_canvas;

                            std::string t_name = "template_" + std::to_string(OM_ID);
                            TH1D* temp_hist = new TH1D(t_name.c_str(), t_name.c_str(), 130, 0, 130/2.56);
                            temp_hist->SetLineColor(2);

                            for (int k = 0; k < 130; ++k) {
                                temp_hist->SetBinContent(k+1, template_vectors[OM_ID][k]);
                            }

                            for (int i = 0; i < n_try; ++i)
                            {
                                std::vector<Double_t> temp_vector;

                                for (uint16_t isample = ch_peak_cell - 30 - n_try/2 + i; isample < ch_peak_cell + 100 - n_try/2 + i; isample++)
                                {
                                    //uint16_t adc = calo_hit.get_waveforms().get_adc(isample,ichannel);
                                    Double_t volts = (Double_t)waveform_adc[isample] - my_baseline;
                                    temp_vector.push_back( volts );
                                    //std::cout << isample << " : " << adc - baseline << std::endl;
                                }

                                TCanvas* my_canvas = new TCanvas();
                                my_canvas->cd();

                                gStyle->SetOptStat(0);
                                std::string name = "hist_" + std::to_string(i);

                                TH1D* hist = new TH1D(name.c_str(), name.c_str(), 130, 0, 130/2.56);
                                hist->SetMaximum(10);
                                hist->SetMinimum(-200);
                                TLegend* legend = new TLegend(0.7, 0.1, 0.9, 0.2);
                                //gStyle->SetLegendBorderSize(0);

                                cont = false;
                                for (int j = 1; j <= 130; ++j)
                                {
                                    hist->SetBinContent(j, temp_vector[j-1]/2.048);
                                }
                                temp_hist->Scale(hist->Integral()/temp_hist->Integral());
                                temp_hist->Sumw2();


                                Double_t norm_test = sqrt( get_inner_product( temp_vector, temp_vector ) );
                                Double_t mf = get_inner_product( temp_vector, template_vectors[OM_ID] )/( norm_temp * norm_test );

                                hist->SetLineColor(1);
                                hist->SetXTitle("Relative time /ns");
                                hist->SetYTitle("Voltage /mV");
                                std::string title = "M:1:9:7 MF: " + std::to_string(mf) + " FBT:" + std::to_string((ch_peak_cell - 30 - n_try/2 + i)/2.56) + " ns";
                                hist->SetTitle(title.c_str());

                                std::string can_name = "mf_output_" + std::to_string(i) + ".png";
                                hist->Draw("HIST");
                                legend->AddEntry(hist, "test");
                                legend->AddEntry(temp_hist, "template");
                                my_canvas->SetGrid(true);
                                temp_hist->Draw("HIST SAME C");
                                legend->Draw();
                                my_canvas->Update();
                                my_canvas->SaveAs(can_name.c_str());

                                delete hist;
                                delete my_canvas;
                                delete legend;

                                mf_output.push_back(mf);
                            }

                            TCanvas* new_canvas = new TCanvas();
                            new_canvas->cd();
                            new_canvas->SetGrid(true);

                            TH1D* new_hist = new TH1D("mf_output", "mf_output", n_try,
                                                      (ch_peak_cell - 30 - n_try/2)/2.56 , (ch_peak_cell - 30 + n_try/2)/2.56);

                            for (int k = 0; k < n_try; ++k) {
                                new_hist->SetBinContent(k+1, mf_output[k]);
                            }

                            new_hist->SetXTitle("Waveform time /ns");
                            new_hist->SetYTitle("mf shape  output /mV");
                            new_hist->SetTitle("M:1.9.7 pulse start finder");

                            gStyle->SetOptStat(0);
                            new_hist->Draw("HIST");

                            if (cont) {} else { return 1; }


                            calo_time = get_max_value(mf_output) + ch_peak_cell - 30 - n_try/2;

                            //waveform = temp_vector;

                            tree.Fill();
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
    //std::cout << std::endl;
    //std::cout << "<<< update_temp_vector >>> " << std::endl;
    //std::cout << "OM_ID: " << OM_ID << std::endl;
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

    //std::cout << std::endl;
    //std::cout << "Template " << OM_ID << std::endl;
    /*for (int i = 0; i < (Int_t)template_vectors[OM_ID].size(); ++i)
    {
        std::cout << "( " << i << " , " << template_vectors[OM_ID][i] << " )" << std::endl;
    }*/

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
