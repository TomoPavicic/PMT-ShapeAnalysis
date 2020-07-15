{
    for(int file_num = 180;file_num<181;file_num++)
    {
        char filename[100];
        sprintf(filename,"/Users/willquinn/Documents/PhD/SuperNEMO/SNEMO_ComData_Analysis/MAIN_WALL/ROOT_files/Average_Waveforms_run%i.root",file_num);
        TFile* rootfile = new TFile(filename,"READ");
        if (rootfile == nullptr)
        {
            //cout << "Error: can't open file"<<endl;
        }else{
            for(int slot_num = 0; slot_num<20;slot_num++)
            {
                char Slot_Num[10];
                sprintf(Slot_Num,"Slot%i",slot_num);
                TDirectory* Slot = (TDirectory*)rootfile->GetDirectory(Slot_Num);
                if (Slot != nullptr)
                {
                    cout << "Slot: "<<slot_num<<endl;
                    for(int channel_num = 0; channel_num<20;channel_num++)
                    {
                        char Channel_Num[15];
                        sprintf(Channel_Num,"Channel%i",channel_num);
                        TDirectory* Channel = (TDirectory*)rootfile->GetDirectory(Slot_Num)->GetDirectory(Channel_Num);

                        if (Channel != nullptr)
                        {
                            cout << "Channel: "<<channel_num<<endl;
                            char Shape_Number[30];
                            sprintf(Shape_Number,"amp_shape%i_%i",slot_num,channel_num);
                            TH2F* hist = (TH2F*)rootfile->GetDirectory(Slot_Num)->GetDirectory(Channel_Num)->Get(Shape_Number);

                            if(hist != nullptr)
                            {
                                cout<<"hist: "<<Shape_Number<<endl;
                                char CanvasPDFName[100];
                                sprintf(CanvasPDFName,"/Users/willquinn/Documents/PhD/SuperNEMO/SNEMO_ComData_Analysis/MAIN_WALL/PDFs/run%i/shape_vs_amplitude/shape_vs_amplitude_%i_%i_%i.pdf",file_num,file_num,slot_num,channel_num);
                                TCanvas* c1 = new TCanvas(CanvasPDFName);
                                hist->Draw("colz");
                                gStyle->SetOptStat(0);
                                c1->Update();
                                c1->SaveAs(CanvasPDFName,"pdf");
                                delete c1;
                            }
                            delete hist;
                        }
                        delete Channel;
                    }
                    delete Slot;
                }
            }
        }
        rootfile->Close();
        delete rootfile;
    }
}