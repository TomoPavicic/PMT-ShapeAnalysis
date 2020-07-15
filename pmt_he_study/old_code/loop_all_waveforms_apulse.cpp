
/********************************************************************************
 *                                                                              *
 *                      Afterpulse testing Code                                 *
 *                          March 2019                                          *
 *                                                                              *
 ********************************************************************************/

{
    //get the waveform file:
    TFile* file = new TFile("/Users/adilaislam/helium_pmt/1.4kV_2Ch_March/channel_1/loop_waveforms/all_nonempty_waveforms.root");

    // Make an output file to store histograms with amplitude thresholds 25, 30, 35:
    TFile* histFile = new TFile("afterpulse_numbers_thres_25.root","RECREATE");
    //TFile* histFile = new TFile("afterpulse_numbers_thres_30.root","RECREATE");
    //TFile* histFile = new TFile("afterpulse_numbers_thres_35.root","RECREATE");


    //make the histogram for the afterpulse numbers, timestamp location and amplitude index:
    TH1D* afterpulse_number = new TH1D("num afterpulses","Number of Afterpulses for PMT GA0607 (Amplitude Index Threshold = 25)", 100, 0, 100);
    TH1D* afterpulse_time = new TH1D("time afterpulses","Afterpulse Location Timestamp for PMT GA0607 (Amplitude Index Threshold = 25)", 7000, 0, 7000);
    TH1D* afterpulse_amplitude = new TH1D("amplitude afterpulses","Amplitude Index of the Afterpulses for PMT GA0607 (Amplitude Index Threshold = 25)", 50, 0, 50);

    
    /******************************************************
    
                   GET TEMPLATE WAVEFORM
    
    ********************************************************/
    
    //extract template array file
    ifstream template_array_file;
    template_array_file.open("/Users/adilaislam/helium_pmt/1.4kV_2Ch_March/channel_1/template_files/template_array_March1.4kV_CH1.dat");
    
    //create and fill a vector that contains the elements of the template vector
    std::vector<double> vector_template;
    
    if(template_array_file.is_open())
    {
        string array_data;
        
        // Read in the file storing the template and fill it into the template array
        while (getline(template_array_file, array_data))
        {
            stringstream template_array_file_Readin;
            
            template_array_file_Readin << array_data;
            
            Double_t arrayElement = 0;
            
            template_array_file_Readin >> arrayElement;
            
            vector_template.push_back(arrayElement);
            
        }
        
    }
    else
    {
        cout << "-- Error opening template array file --";
    }
    

    for (int iwaveform = 0; iwaveform < 10000; ++iwaveform)
    {
    
        // generate the name of the i'th waveform:
        char histname_i[20];
        sprintf(histname_i,"waveform_%i",iwaveform);
        cout << "histname = " << histname_i << endl;

        // TH1F* hist = (TH1F*)file->Get("waveform_10");
        TH1F* hist = (TH1F*)file->Get(histname_i);

        //*****************************************************//
        
                    /*GET TEST WAVEFORM*/
        
        //*****************************************************//
        
        //create and a fill a vector with the ADC values
        
        std::vector<double> ADC_values;
        for (int i=0; i< hist->GetEntries(); i++)
        {
            ADC_values.push_back(hist->GetBinContent(i));
            
            //std::cout<< hist->GetBinContent(i) << endl;
        }
        
        // Average the pre-trigger region to get the baseline:
        double baseline_average = 0.;
        double baseline_calculated = 0.;
        double pre_trigger = 590;
        for (int j = 0; j < pre_trigger - 20.; j++)
        
        {
            baseline_average += ADC_values.at(j);
            
        }
        
        baseline_calculated = baseline_average/(pre_trigger - 20.);
        
        //subtract the baseline from the ADC values:
        for (int i = 0; i < hist->GetEntries(); i++)
        {
            ADC_values.at(i) =  ADC_values.at(i) - baseline_calculated;
            
        }
        
        //define the region where the afterpulsing starts and the size of the sweep window:
        int apulse_region_start = 0;
        
        int sweep_window = vector_template.size();
        
        // Calculate the normalisation of the template vector
        double norm1 = 0;
        for (int r = 0; r< vector_template.size(); r++)
        {
            double inter = 0;
            inter = (vector_template.at(r) * vector_template.at(r));
            norm1 = norm1 + inter;
        }
        double template_norm = sqrt(norm1);
        
        int apulse_counter = 0; //afterpulse counter.

        
        //count the number of afterpulses
        //do the calculation for amplitude index again, but define the afterpulse region so that the first main pulse isn't included in the counting.

        int first_i; //first timestamp point you hit the pulse
        bool first_pulse_hit = true; //first time you hit the pulse
        double amp_threshold = 25; //set the threshold for the amplitude index.
        
        
        for (int apulse_i = 0; apulse_i < (int)hist->GetEntries();apulse_i++)
        {
            if (apulse_i >= 1000 + int(sweep_window/2) && apulse_i<= (int)hist->GetEntries() - int(sweep_window/2))
                    
            {
                //create and fill a vector for the test template
                std::vector<double> vector_test;
                for (int j =0; j < sweep_window; j++)
                        
                {
                    int k = apulse_i - int(sweep_window/2);
                    vector_test.push_back(ADC_values.at(k+j));
                }
                    
                        
                //make the match filter, inner product between a template signal vector and a sample test vector from the waveform
                double match_filter = 0;
                for (int r = 0; r < sweep_window ; r++)
                {
                    double inter = 0;
                    //cout << vector_template.at(r) << vector_test.at(r) << endl;
                    inter = (vector_template.at(r) * vector_test.at(r));
                    //cout << "inter "<<inter << endl;
                    match_filter += inter;
                            
                }
                
                //Amplitude Index:
                double Amplitude = match_filter/template_norm;
                
                //start afterpulse counting
                if (first_pulse_hit == true)
                {
                    if (Amplitude > amp_threshold)
                    {
                        apulse_counter++;
                                
                        first_i = apulse_i; //storing the point you first hit the pulse
                                
                        first_pulse_hit = false; //so that it will run this "if statement" once
                                
                        cout << first_i << " " << apulse_counter << endl;
                       
                        //fill the time and amplitude histograms:
                        afterpulse_time->Fill(first_i);
                        afterpulse_amplitude->Fill(Amplitude);

                    }
                }
                        
                else
                {
                    if (apulse_i > (first_i+80)) //you have now skipped 2 pulse widths and are searching for the first point of the next afterpulse.
                    {
                        if (Amplitude > amp_threshold)
                        {
                            apulse_counter++;
                                    
                            first_i = apulse_i; //storing the point you first hit the next pulse
                            cout << first_i << " " << apulse_counter << endl;

                            afterpulse_time->Fill(first_i);
                            afterpulse_amplitude->Fill(Amplitude);

                        }
                                

                    }

                }

            }

        }
        cout << "Total number of afterpulses in " << histname_i << ": " <<apulse_counter<<endl;
        
        // Register the number of afterpulses for this waveform:
        afterpulse_number->Fill(apulse_counter);

    }
    //label histogram axis:
    afterpulse_number->GetXaxis()->SetTitle("Number of afterpulses");
    afterpulse_number->GetYaxis()->SetTitle("Frequency");

    afterpulse_time->GetXaxis()->SetTitle("Timestamp of afterpulse (ns)");
    afterpulse_time->GetYaxis()->SetTitle("Frequency");
        
    afterpulse_amplitude->GetXaxis()->SetTitle("Amplitude Index of afterpulse");
    afterpulse_amplitude->GetYaxis()->SetTitle("Frequency");

    //write the histogram files:
    afterpulse_number->Write();
    afterpulse_time->Write();
    afterpulse_amplitude->Write();

    histFile->cd();
    histFile->Close();

    cout << " -- " << endl;
    cout << " -- " << endl;
        
    cout << "-- Code terminated --" << endl;
        
    cout << " -- " << endl;
    cout << " -- " << endl;
        
        
}
