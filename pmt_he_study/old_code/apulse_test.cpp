
/********************************************************************************
 *                                                                              *
 *                      Afterpulse testing Code                                 *
 *                          March 2019                                          *
 *                                                                              *
 ********************************************************************************/

//get the waveform file
{
    TFile* file = new TFile("/Users/willquinn/Documents/PhD/PMT_Permeation_Project/ROOT_files/waveforms/22.11_waveforms.root");
    
    //get the waveform from the file
    TH1F* hist = (TH1F*)file->Get("waveform_75941");

    TGraph* shape_graph = new TGraph();
    TGraph* amplitude_graph = new TGraph();
    TGraph* template_graph = new TGraph();
    TGraph* shape_v_amp_graph = new TGraph();

    //*****************************************************//
        
                /*GET TEST WAVEFORM*/
        
    //*****************************************************//
        
    //create and a fill a vector with the ADC values
        
    std::vector<double> ADC_values;
    
        // Average the pre-trigger region to get the baseline:

    double baseline_average = 0.;
    double baseline_calculated = 0.;
    double pre_trigger = 150;
    for (int j = 0; j < pre_trigger - 20.; j++)
    {
        //cout << hist->GetBinContent(j) << endl;
        baseline_average += hist->GetBinContent(j);

    }

    baseline_calculated = baseline_average/(pre_trigger - 20.);

    for (int i = 0; i < hist->GetEntries(); i++)
    {
        ADC_values.push_back(hist->GetBinContent(i) - baseline_calculated);

    }

        //*****************************************************//

                    /*GET TEMPLATE WAVEFORM*/

        //*****************************************************//

    //extract template array file
    TFile* template_file = new TFile("/Users/willquinn/Documents/PhD/PMT_Permeation_Project/ROOT_files/A1400_B1400_t1445_templates.root","READ");

    //create and fill a vector that contains the elements of the template vector
    std::vector<double> vector_template;


    TH1D* template_hist = (TH1D*)template_file->Get("Template_Waveform_Channel1_A1400_B1400_t1445");

    for (int i = 2000; i < 2000 + template_hist->GetEntries(); i++)
    {
        ADC_values.at(i) += template_hist->GetBinContent(i-2000) * 0.0005;
    }

    for(int m = 0; m < template_hist->GetEntries();m++)
    {
        double arrayElement = template_hist->GetBinContent(m);

        vector_template.push_back(arrayElement);

        template_graph->SetPoint(m,m,arrayElement);

            //cout << arrayElement << endl;
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
         norm1 += inter;
    }
    double template_norm = sqrt(norm1);

    int apulse_counter = 0;

    //counter to fill in the graphs
    int graph_i = 0;

    //do the sweep
    for (int apulse_i = 0; apulse_i < (int)hist->GetEntries();apulse_i++)
    {
        if (apulse_i >= apulse_region_start + int(sweep_window/2) && apulse_i<= (int)hist->GetEntries() - int(sweep_window/2))
        {
            //create and fill a vector for the test template
            std::vector<double> vector_test;
            for (int j =0; j < sweep_window; j++)
            {
                int k = apulse_i - int(sweep_window/2);
                vector_test.push_back(ADC_values.at(k+j));
            }

            // Calculate the normalisation of the test vector
            double norm2 = 0;
            for (int r = 0; r < sweep_window ; r++)
            {
                double inter = 0;
                inter = (vector_test.at(r) * vector_test.at(r));
                norm2 += inter;
            }
            double test_norm = sqrt(norm2);


            //make the match filter, inner product between a template signal vector and a sample test vector from the waveform
            double match_filter = 0;
            for (int r = 0; r < sweep_window ; r++)
            {
                double inter = 0;
                //cout << vector_template.at(r) << " " << vector_test.at(r) << endl;
                inter = (vector_template.at(r) * vector_test.at(r));
                //cout << "inter "<<inter << endl;
                match_filter += inter;
            }

            //cout<< "match filter"<< match_filter << endl;

            // work out the shape and amplitude index
            double Shape = match_filter/(template_norm*test_norm);
            //cout<< "Shape " << Shape << endl;

            double Amplitude = match_filter/template_norm;

            //fill the graphs

            shape_graph->SetPoint(graph_i, graph_i, Shape);

            amplitude_graph->SetPoint(graph_i, graph_i, Amplitude);

            shape_v_amp_graph->SetPoint(graph_i, Amplitude, Shape);

            //cout << "Time " << graph_i << " Amplitude " << Amplitude << " Shape " << Shape << endl;

            graph_i++;

        }

    }

    //count the number of afterpulses
    //do the calculation for amplitude index again, but define the afterpulse region so that the first main pulse isn't included in the counting.

    int first_i; //first timestamp point you hit the pulse
    bool first_pulse_hit = true; //first time you hit the pulse
    double amp_threshold = 5;

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
                }
            }

            else
            {
                if (apulse_i > (first_i+80)) //you have now skipped 2 pulse widths and are searching for the first point of the next afterpulse.
                    {
                        if (Amplitude > amp_threshold)
                        {
                            apulse_counter++;

                            first_i = apulse_i; //storing the point you first the pulse
                            cout << first_i << " " << apulse_counter << endl;

                        }

                    }

            }

        }

    }

    TGraph* waveform_graph = new TGraph();
    for(int i = 0;i < ADC_values.size();i++)
    {
        waveform_graph->SetPoint(i,i,ADC_values.at(i));
    }

            //**************************************************************//

                                //Draw the graphs//

            //**************************************************************//


            //////////////////////////////////////////////////////////////////

            TCanvas* c1 = new TCanvas();
            amplitude_graph->SetTitle("Amplitude Index vs. Timestamp for Waveform 3 PMT GA0607");
            amplitude_graph->SetMarkerStyle(31);
            amplitude_graph->SetMarkerColor(4);
            amplitude_graph->GetXaxis()->SetTitle("Timestamp /ns");
            amplitude_graph->GetYaxis()->SetTitle("Amplitude Index");
            amplitude_graph->GetYaxis()->SetRange(0,1000);
            //amplitude_graph->Write();
            amplitude_graph->Draw("AP");
            c1->SetLogy();
            c1->SetGrid();

//            ///////////////////////////////////////////////////////////////////////

            TCanvas* c2 = new TCanvas();
            shape_graph->SetTitle("Shape Index vs. Timestamp for Waveform 3 PMT GA0607");
            shape_graph->SetMarkerStyle(31);
            shape_graph->SetMarkerColor(2);
            shape_graph->GetXaxis()->SetTitle("Timestamp/ns");
            shape_graph->GetYaxis()->SetTitle("Shape Index");
            //shape_graph->Write();
            shape_graph->Draw("AP");
            c2->SetGrid();
//
//
//            ///////////////////////////////////////////////////////////////////////
//
////
            TCanvas* c3 = new TCanvas();
            template_graph->SetTitle("Template Shape Waveform for PMT GA0607");
            template_graph->SetMarkerStyle(31);
            template_graph->SetMarkerColor(6);
            template_graph->GetXaxis()->SetTitle("Timestamp /ns");
            template_graph->GetYaxis()->SetTitle("ADC-count /mV");
            //template_graph->Write();
            template_graph->Draw("AP");
            c3->SetGrid();
//

            ///////////////////////////////////////////////////////////////////////


            TCanvas* c4 = new TCanvas();
            //shape_v_amp_graph->Write();
            shape_v_amp_graph->GetXaxis()->SetTitle("Amplitude Index");
            shape_v_amp_graph->SetTitle("Amplitide vs. Indices Shape");
            shape_v_amp_graph->GetYaxis()->SetTitle("Shape Index");
            shape_v_amp_graph->Draw("AP*");
            c4->SetGrid();
    //

            ///////////////////////////////////////////////////////////////////////

             TCanvas *c5 = new TCanvas();
             waveform_graph->GetXaxis()->SetTitle("Timestamp /ns");
            waveform_graph->SetTitle("Waveform");
            waveform_graph->GetYaxis()->SetTitle("ADC counts");
            waveform_graph->Draw("AP*");
            c5->SetGrid();

/*
            TCanvas *c5 = new TCanvas("c5","3 Graphs",700,900);
//
            auto *p2 = new TPad("p2","p3",0.,0.,1.,0.3);
            p2->Draw();
            p2->SetTopMargin(0.001);
            p2->SetBottomMargin(0.3);
            p2->SetGrid();

            auto *p1 = new TPad("p1","p1",0.,0.3,1.,1.);
            p1->Draw();
            p1->SetBottomMargin(0.001);
            p1->cd();
            p1->SetGrid();
//
//
            amplitude_graph->SetLineColor(4);
            amplitude_graph->SetMarkerColor(1);
            amplitude_graph->SetMarkerStyle(8);
            amplitude_graph->SetMarkerSize(0.8);
            amplitude_graph->SetTitle("");

            amplitude_graph->SetTitle("Amplitude vs Timestamp");
//
            amplitude_graph->GetXaxis()->SetTitle("Timestamp /ns");
            amplitude_graph->GetXaxis()->CenterTitle();
//
            shape_graph->GetYaxis()->SetTitle("Shape Index vs Timestamp");
            amplitude_graph->GetYaxis()->CenterTitle();
            amplitude_graph->GetYaxis()->SetRangeUser(-1,350);
            amplitude_graph->GetYaxis()->SetTitleOffset(1.5);
            amplitude_graph->GetXaxis()->SetTickSize(0.);
//
            amplitude_graph->Draw("P");
//
            hist->SetMarkerColor(kBlue);
            hist->Draw();
//
            TLegend *leg = new TLegend(0.15,0.75,0.5,0.85);
//
            leg->AddEntry(shape_graph,"shape index","lp");
            leg->AddEntry(hist,"waveform","lp");
            leg->Draw();
//
            p2->cd();
            amplitude_graph->Draw("AL");
//
//
*/
            ///////////////////////////////////////////////////////////////////////////////

//
            //overlayed graphs
            gROOT->ForceStyle(0);
            c6 = new TCanvas("c6","Shape and Amplitide Index vs Timestamp Waveform 2 PMT GA0607",200,10,700,2000);
            TPad *pad = new TPad("pad","",0,0,1,1);
            pad->SetGrid();
            pad->Draw();
            pad->cd();


            // draw a frame to define the range
            TH1F *hr = c6->DrawFrame(0,-1000,ADC_values.size(),0);
            hr->SetTitle("Amplitude Index vs PMT Waveform");
            hr->SetXTitle("Timestamp /ns");
            hr->SetYTitle("ADC Counts");
            pad->GetFrame()->SetBorderSize(12);

            // create first graph

            //waveform_graph->SetMarkerColor(kBlue);
            //waveform_graph->SetMarkerStyle(21);
            waveform_graph->Draw("Po");


            //create a transparent pad drawn on top of the main pad
            c6->cd();
            TPad *overlay = new TPad("overlay","",0,0,1,1);
            overlay->SetFillStyle(0);
            overlay->SetFillColor(0);
            overlay->SetFrameFillStyle(0);
            overlay->Draw("FA");
            overlay->cd();
            amplitude_graph->SetMarkerColor(kRed);
            amplitude_graph->SetMarkerStyle(20);
            amplitude_graph->SetName("shape_graph");
            Double_t xmin = 0;
            Double_t ymin = 0;
            Double_t xmax = ADC_values.size();
            Double_t ymax = 200;
            TH1F *hframe = overlay->DrawFrame(xmin,ymin,xmax,ymax);
            hframe->GetXaxis()->SetLabelOffset(99);
            hframe->GetYaxis()->SetTickLength(0);
            hframe->GetYaxis()->SetLabelOffset(99);
            amplitude_graph->Draw("P*");


            //Draw an axis on the right side
            TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,510,"+L");
            axis->SetLineColor(kRed);
            axis->SetLabelColor(kRed);
            axis->SetTitleColor(kRed);
            axis->SetTitle("Amplitude Index");
            axis->Draw();



            /////////////////////////////////////////////////////////////////////////


            cout << " -- " << endl;
            cout << " -- " << endl;
        
            cout << "-- Code terminated --" << endl;
        
            cout << " -- " << endl;
            cout << " -- " << endl;
        
    }
