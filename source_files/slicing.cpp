#ifndef SLICING_CPP
#define SLICING_CPP

#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/slicing.h"

using namespace xsec;

void Slices ( TH2D *h_unfolded, TH2D *h_true, TH2D *h_reco, const char file[1024] ){

    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    // Take the bin edges to be the title
    int x_bins = h_unfolded->GetNbinsX(); // Cos theta
    int y_bins = h_unfolded->GetNbinsY(); // Tmu

    TCanvas *c_Tmu   = new TCanvas ( "c_Tmu", "", 800, 600 );
   
    TLegend *leg_T   = new TLegend( 0.12, 0.68, 0.32, 0.88 );

    TH1D *h_Tmu      = new TH1D ( "h_Tmu", "", x_bins, -1, 1 );
    TH1D *h_Tmu_true = new TH1D ( "h_Tmu_true", "", x_bins, -1, 1 );
    TH1D *h_Tmu_reco = new TH1D ( "h_Tmu_reco", "", x_bins, -1, 1 );
    
    leg_T->AddEntry( h_Tmu, " Unfolded ", "l" );
    leg_T->AddEntry( h_Tmu_true, " True ", "l" );
    leg_T->AddEntry( h_Tmu_reco, " Reco ", "l" );
    
    for ( int i = 1; i <= y_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_T;
        double up_edge_T;

        low_edge_T = h_unfolded->GetYaxis()->GetBinLowEdge(i);
        up_edge_T = h_unfolded->GetYaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        stringstream conv;
        conv.clear();

        string title;
        title.clear();

        char file_name[1024];

        conv << setprecision(4) << file << "comp_Tmu_slice_" << low_edge_T << ";" << up_edge_T << ".png";
        title = conv.str();
        
        strcpy( file_name, title.c_str() );

        // For the title of the histogram
        stringstream hist;
        hist.clear();

        char hist_name[1024];
        string temp1;

        hist << setprecision(4) << "T_{#mu} slice: " << low_edge_T << "_" << up_edge_T;
        temp1 = hist.str();

        strcpy( hist_name, temp1.c_str() );

        // Fill the histogram
        for ( int j = 1; j <= x_bins; ++j ){
            h_Tmu->SetBinContent( j, h_unfolded->GetBinContent(j, i) );
            h_Tmu_true->SetBinContent( j, h_true->GetBinContent(j, i) );
            h_Tmu_reco->SetBinContent( j, h_reco->GetBinContent(j, i) );
        }

        double unfo_max = h_Tmu->GetBinContent(h_Tmu->GetMaximumBin());
        double true_max = h_Tmu_true->GetBinContent(h_Tmu_true->GetMaximumBin());
        double reco_max = h_Tmu_reco->GetBinContent(h_Tmu_reco->GetMaximumBin());

        double max = 1.1 * TMath::Max(unfo_max, TMath::Max(true_max, reco_max)); 

        h_Tmu->Draw();
        h_Tmu_true->Draw("same");
        h_Tmu_reco->Draw("same");
        h_Tmu->SetTitle(hist_name);
        h_Tmu->GetYaxis()->SetRangeUser(0,max);
        h_Tmu->GetXaxis()->SetTitle("cos#theta_{#mu}");   
        h_Tmu->GetYaxis()->SetTitle("Event rate");   
        h_Tmu->SetLineColor( kRed + 2 );
        h_Tmu_true->SetLineColor( kGreen + 2 );
        h_Tmu_reco->SetLineColor( kBlue + 2 );
        h_Tmu_reco->SetLineStyle( 7 );

        h_Tmu->SetTitleOffset(1.45, "Y");
        h_Tmu->SetStats(kFALSE);

        leg_T->Draw();
        c_Tmu->SaveAs(file_name);

    } 
   
    delete h_Tmu;
    delete h_Tmu_true;
    delete h_Tmu_reco;

    delete c_Tmu;
    
    delete leg_T;

    TCanvas *c_cosmu   = new TCanvas ( "c_cosmu", "", 800, 600 );
    
    TLegend *leg_c     = new TLegend( 0.68, 0.68, 0.88, 0.88 );

    TH1D *h_cosmu      = new TH1D ( "h_cosmu", "", y_bins, 0.2, 2 );
    TH1D *h_cosmu_true = new TH1D ( "h_cosmu_true", "", y_bins, 0.2, 2 );
    TH1D *h_cosmu_reco = new TH1D ( "h_cosmu_reco", "", y_bins, 0.2, 2 );
     
    leg_c->AddEntry( h_cosmu, " Unfolded ", "l" );
    leg_c->AddEntry( h_cosmu_true, " True ", "l" );
    leg_c->AddEntry( h_cosmu_reco, " Reco ", "l" );
    
    // Cos theta mu slices
    for ( int i = 1; i <= x_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_cos;
        double up_edge_cos;

        low_edge_cos = h_unfolded->GetXaxis()->GetBinLowEdge(i);
        up_edge_cos = h_unfolded->GetXaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        stringstream conv1;
        conv1.clear();

        string title1;
        title1.clear();

        char file_name1[1024];

        conv1 << setprecision(4) << file << "comp_cos_thetamu_slice_" << low_edge_cos << ";" << up_edge_cos << ".png";
        title1 = conv1.str();
        
        strcpy( file_name1, title1.c_str() );

        // For the title of the histogram
        stringstream hist1;
        hist1.clear();

        char hist_name1[1024];
        string temp2;

        hist1 << setprecision(4) << "cos#theta_{#mu} slice: " << low_edge_cos << "," << up_edge_cos;
        temp2 = hist1.str();

        strcpy( hist_name1, temp2.c_str() );

        // Fill the histogram
        for ( int j = 1; j <= y_bins; ++j ){
            h_cosmu->SetBinContent( j, h_unfolded->GetBinContent(i, j) );
            h_cosmu_true->SetBinContent( j, h_true->GetBinContent(i, j) );
            h_cosmu_reco->SetBinContent( j, h_reco->GetBinContent(i, j) );
        }

        double unfo_max_c = h_cosmu->GetBinContent(h_cosmu->GetMaximumBin());
        double true_max_c = h_cosmu_true->GetBinContent(h_cosmu_true->GetMaximumBin());
        double reco_max_c = h_cosmu_reco->GetBinContent(h_cosmu_reco->GetMaximumBin());

        double max_c = 1.1 * TMath::Max(unfo_max_c, TMath::Max(true_max_c, reco_max_c)); 
        
        h_cosmu->Draw();
        h_cosmu_true->Draw("same");
        h_cosmu_reco->Draw("same");
        h_cosmu->SetTitle(hist_name1);
        h_cosmu->GetYaxis()->SetRangeUser(0,max_c);
        h_cosmu->GetXaxis()->SetTitle("T_{#mu}");   
        h_cosmu->GetYaxis()->SetTitle("Event rate");   
        h_cosmu->SetLineColor( kRed + 2 );
        h_cosmu_true->SetLineColor( kGreen + 2 );
        h_cosmu_reco->SetLineColor( kBlue + 2 );
        h_cosmu_reco->SetLineStyle( 7 );
        h_cosmu->SetTitleOffset(1.45, "Y");
        h_cosmu->SetStats(kFALSE);

        leg_c->Draw();
        c_cosmu->SaveAs(file_name1);

    } 
    
    delete h_cosmu;
    delete h_cosmu_true;
    delete h_cosmu_reco;

    delete c_cosmu;

    delete leg_c;

}

void Slices ( TH2D *h_ddxsec, const char file[1024] ){

    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    // Take the bin edges to be the title
    int x_bins = h_ddxsec->GetNbinsX(); // Cos theta
    int y_bins = h_ddxsec->GetNbinsY(); // Tmu

    TCanvas *c_Tmu   = new TCanvas ( "c_Tmu", "", 800, 600 );
   
    TH1D *h_Tmu      = new TH1D ( "h_Tmu", "", x_bins, -1, 1 );
    
    for ( int i = 1; i <= y_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_T;
        double up_edge_T;

        low_edge_T = h_ddxsec->GetYaxis()->GetBinLowEdge(i);
        up_edge_T = h_ddxsec->GetYaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        stringstream conv;
        conv.clear();

        string title;
        title.clear();

        char file_name[1024];

        conv << setprecision(4) << file << "xsec_Tmu_slice_" << low_edge_T << ";" << up_edge_T << ".png";
        title = conv.str();
        
        strcpy( file_name, title.c_str() );

        // For the title of the histogram
        stringstream hist;
        hist.clear();

        char hist_name[1024];
        string temp1;

        hist << setprecision(4) << "T_{#mu} slice: " << low_edge_T << "_" << up_edge_T;
        temp1 = hist.str();

        strcpy( hist_name, temp1.c_str() );

        // Fill the histogram
        for ( int j = 1; j <= x_bins; ++j ){
            h_Tmu->SetBinContent( j, h_ddxsec->GetBinContent(j, i) );
            h_Tmu->SetBinError( j, h_ddxsec->GetBinError(j, i) ); 
        }

        double max = 1.2 * h_Tmu->GetBinContent(h_Tmu->GetMaximumBin());

        h_Tmu->Draw("e1x0");
        h_Tmu->SetTitle(hist_name);
        h_Tmu->GetYaxis()->SetRangeUser(0,max);
        h_Tmu->GetXaxis()->SetTitle("cos#theta_{#mu}");   
        h_Tmu->GetYaxis()->SetTitle("CC0#pi, d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]");

        h_Tmu->SetTitleOffset(1.2, "Y");
        h_Tmu->SetStats(kFALSE);

        c_Tmu->SaveAs(file_name);

    } 
   
    delete h_Tmu;

    delete c_Tmu;

    TCanvas *c_cosmu   = new TCanvas ( "c_cosmu", "", 800, 600 );
    
    TH1D *h_cosmu      = new TH1D ( "h_cosmu", "", y_bins, 0.2, 2 );
     
    // Cos theta mu slices
    for ( int i = 1; i <= x_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_cos;
        double up_edge_cos;

        low_edge_cos = h_ddxsec->GetXaxis()->GetBinLowEdge(i);
        up_edge_cos = h_ddxsec->GetXaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        stringstream conv1;
        conv1.clear();

        string title1;
        title1.clear();

        char file_name1[1024];

        conv1 << setprecision(4) << file << "xsec_cos_thetamu_slice_" << low_edge_cos << ";" << up_edge_cos << ".png";
        title1 = conv1.str();
        
        strcpy( file_name1, title1.c_str() );

        // For the title of the histogram
        stringstream hist1;
        hist1.clear();

        char hist_name1[1024];
        string temp2;

        hist1 << setprecision(4) << "cos#theta_{#mu} slice: " << low_edge_cos << "," << up_edge_cos;
        temp2 = hist1.str();

        strcpy( hist_name1, temp2.c_str() );

        // Fill the histogram
        for ( int j = 1; j <= y_bins; ++j ){
            h_cosmu->SetBinContent( j, h_ddxsec->GetBinContent(i, j) );
            h_cosmu->SetBinError( j, h_ddxsec->GetBinError(i, j) ); 
        }

        double max_c = 1.2 * h_cosmu->GetBinContent(h_cosmu->GetMaximumBin());
        
        h_cosmu->Draw("e1x0");
        h_cosmu->SetTitle(hist_name1);
        h_cosmu->GetYaxis()->SetRangeUser(0,max_c);
        h_cosmu->GetXaxis()->SetTitle("T_{#mu}");   
        h_cosmu->GetYaxis()->SetTitle("CC0#pi, d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]");
        h_cosmu->SetTitleOffset(1.2, "Y");
        h_cosmu->SetStats(kFALSE);

        c_cosmu->SaveAs(file_name1);

    } 
    
    delete h_cosmu;

    delete c_cosmu;

}

void Characterisation ( TH2D *h_smeared, std::vector< Particle > primary, std::vector< Event > reco_events,  const char file_path[1024] ){

    int n_entries  = reco_events.size();

    // Count the number of events with each type of physical process
    int qel   = 0;
    int mec   = 0;
    int res   = 0;
    int dis   = 0;
    int coh   = 0;
    int singlepi   = 0;
    int other = 0;
    
    // Find the physical process of the event and enumerate according to 
    //  QE           : 0
    //  MEC          : 1
    //  RES          : 2
    //  DIS          : 3
    //  COH          : 4
    //  Non RES, 1pi : 5
    //  Other        : 6
    for ( int i = 0; i < n_entries; ++i ){
    
        if( reco_events[i].GetPhysicalProcess() == 0 ){
            ++qel; 
        }
        else if( reco_events[i].GetPhysicalProcess() == 1 ){
            ++mec; 
        }
        else if( reco_events[i].GetPhysicalProcess() == 2 ){
            ++res; 
        }
        else if( reco_events[i].GetPhysicalProcess() == 3 ){
            ++dis; 
        }
        else if( reco_events[i].GetPhysicalProcess() == 4 ){
            ++coh; 
        }
        else if( reco_events[i].GetPhysicalProcess() == 5 ){
            ++singlepi; 
        }
        else if( reco_events[i].GetPhysicalProcess() == 6 ){
            ++other; 
        }
    }

    std::cout << " QEL : " << qel << std::endl;
    std::cout << " MEC : " << mec << std::endl;
    std::cout << " RES : " << res << std::endl;
    std::cout << " DIS : " << dis << std::endl;
    std::cout << " COH : " << coh << std::endl;
    std::cout << " 1pi : " << singlepi << std::endl;
    std::cout << " Oth : " << other << std::endl;
    
    vector< TH1D* > signal_h;
    
    TCanvas *c_Tmu         = new TCanvas ( "c_Tmu", "", 800, 600 );
    TLegend *leg_T         = new TLegend( 0.16, 0.66, 0.30, 0.88 );

    // Take the bin edges to be the title
    int x_bins = h_smeared->GetNbinsX(); // Cos theta
    int y_bins = h_smeared->GetNbinsY(); // Tmu
   
    TH1D *h_Tmu_muccqe     = new TH1D ( "h_Tmu_muccqe",    "", x_bins, -1, 1 );
    TH1D *h_Tmu_muccmec    = new TH1D ( "h_Tmu_muccmec",   "", x_bins, -1, 1 );
    TH1D *h_Tmu_muccres    = new TH1D ( "h_Tmu_muccres",   "", x_bins, -1, 1 );
    TH1D *h_Tmu_muccdis    = new TH1D ( "h_Tmu_muccdis",   "", x_bins, -1, 1 );
    TH1D *h_Tmu_mucccoh    = new TH1D ( "h_Tmu_mucccoh",   "", x_bins, -1, 1 );
    TH1D *h_Tmu_mucc1pi    = new TH1D ( "h_Tmu_mucc1pi",   "", x_bins, -1, 1 );
    TH1D *h_Tmu_muccother  = new TH1D ( "h_Tmu_muccother", "", x_bins, -1, 1 );
    
    // For any counters that don't = 0, initialise the histogram and legend entry
    if( qel > 0 ){
        
        leg_T->AddEntry( h_Tmu_muccqe,     " #nu_{#mu} CC QE ", "f" );
        signal_h.push_back(h_Tmu_muccqe); 
    }
    if( mec > 0 ){

        leg_T->AddEntry( h_Tmu_muccmec,    " #nu_{#mu} CC MEC ", "f" );
        signal_h.push_back(h_Tmu_muccmec);
    }
    if( res > 0 ){

        leg_T->AddEntry( h_Tmu_muccres, " #nu_{#mu} CC RES ", "f" );
        signal_h.push_back(h_Tmu_muccres);
    }
    if( dis > 0 ){

        leg_T->AddEntry( h_Tmu_muccdis, " #nu_{#mu} CC DIS ", "f" );
        signal_h.push_back(h_Tmu_muccdis);
    }
    if( coh > 0 ){

        leg_T->AddEntry( h_Tmu_mucccoh, " #nu_{#mu} CC COH ", "f" );
        signal_h.push_back(h_Tmu_mucccoh);
    }
    if( singlepi > 0 ){    
        
        leg_T->AddEntry( h_Tmu_mucc1pi,    " #nu_{#mu} CC non-RES 1#pi ", "f" );
        signal_h.push_back(h_Tmu_mucc1pi);
    }
    if( other > 0 ){
        
        leg_T->AddEntry( h_Tmu_muccother,  " #nu_{#mu} CC Other ", "f" );
        signal_h.push_back(h_Tmu_muccother);
    }

    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    for ( int i = 1; i <= y_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_T;
        double up_edge_T;

        low_edge_T = h_smeared->GetYaxis()->GetBinLowEdge(i);
        up_edge_T  = h_smeared->GetYaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        // Filenames for signal histograms
        stringstream conv;
        conv.clear();

        string title;
        title.clear();

        char file_name[1024];

        conv << setprecision(4) << file_path << "_pre_FSI_Tmu_slice_" << low_edge_T << "_" << up_edge_T << ".png";
        title = conv.str();
        
        strcpy( file_name, title.c_str() );

        // For the title of the histogram
        stringstream hist;
        hist.clear();

        char hist_name[1024];
        string temp;

        hist << setprecision(4) << "T_{#mu} in [" << low_edge_T << "," << up_edge_T << "], GeV ";
        temp = hist.str();

        strcpy( hist_name, temp.c_str() );

        for( unsigned int i = 0; i < signal_h.size(); ++i ){
            signal_h[i]->Reset("M");
        }  

        for ( int k = 0; k < n_entries; ++k ){
      
            Event ev    = reco_events[k];            
            int process = ev.GetPhysicalProcess();
            
            Particle p  = primary[k];
            double cth  = p.GetCosSmeared();
            double T    = p.GetTSmeared();


            // Muon neutrino in this bin
            if ( T >= low_edge_T 
              && T < up_edge_T
              && ev.GetNeutrinoFlavour() == 14
              && ev.GetIsCC() ){
                
                // CCQE
                if ( process == 0  ) {
                    h_Tmu_muccqe->Fill( cth,1 );
                }
                
                // CCMEC
                else if ( process == 1 ){
                    h_Tmu_muccmec->Fill( cth, 1 );
                }
                
                // CCRES
                else if ( process == 2 ){
                    h_Tmu_muccres->Fill( cth, 1 );
                }
                
                // CCDIS
                else if ( process == 3 ){
                    h_Tmu_muccdis->Fill( cth, 1 );
                }
                
                // CCCOH
                else if ( process == 4 ){
                    h_Tmu_mucccoh->Fill( cth, 1 );
                }
                
                // CC non-RES, 1Pi
                else if ( process == 5 ){
                    h_Tmu_mucc1pi->Fill( cth, 1 );
                }
                
                // CCOther
                else if ( process == 6 ){
                    h_Tmu_muccother->Fill( cth, 1 );
                }
            }    
        }

        int pal[12];
        pal[0]  = kRed;
        pal[1]  = kGreen + 2;
        pal[2]  = kOrange + 7;
        pal[3]  = kBlue;
        pal[4]  = kMagenta + 1;
        pal[5]  = kCyan + 2;
        pal[6]  = kYellow + 1;
        pal[7]  = kAzure + 2;
        pal[8]  = kViolet + 2;
        pal[9]  = kTeal - 1;
        pal[10] = kPink + 2;
        pal[11] = kSpring + 2;
        
        gStyle->SetPalette( 12, pal );
        gStyle->SetHatchesLineWidth( 1 );
        gStyle->SetHatchesSpacing( 0.5 );
        gStyle->SetTitleOffset(1.5, "Y");
        gStyle->SetOptStat( 0 );

        double norm_Tmu = 0;
        
        THStack *hsT         = new THStack("hsT","pre-FSI ");
        
        for ( unsigned int i = 0; i < signal_h.size(); ++i ){
            norm_Tmu += signal_h[i]->Integral();
        }

        for ( unsigned int i = 0; i < signal_h.size(); ++i ){

            signal_h[i]->SetFillColor(pal[i]);
            signal_h[i]->SetLineColor(pal[i]);

            signal_h[i]->SetLineWidth(1.5);
            signal_h[i]->Scale(1/norm_Tmu);

            hsT->Add(signal_h[i]);
            
        }
        
        // Fill the stacks 
        // Signal
        hsT->Draw();
        
        hsT->SetTitle(hist_name);
        hsT->GetXaxis()->SetTitle("cos#theta_{#mu}");   
        hsT->GetYaxis()->SetTitle("Normalised event rate");   
        
        leg_T->Draw();
        c_Tmu->SaveAs(file_name);

        delete hsT;
    } 

    for ( unsigned int i = 0; i < signal_h.size(); ++i ){
        delete signal_h[i];
    }
    delete c_Tmu;
    
    delete leg_T;
   /*
    // ------------------------------------------------------------------------------
    //                              Cos slices
    // ------------------------------------------------------------------------------
    
    vector< TH1D* > signal_c_h;
    
    TCanvas *c_cosmu = new TCanvas ( "c_cosmu", "", 800, 600 );
    
    TLegend *leg_cos = new TLegend( 0.72, 0.66, 0.86, 0.88 );
 
    TH1D *h_cosmu_muccqe     = new TH1D ( "h_cosmu_muccqe",     "", y_bins, 0.2, 2 );
    TH1D *h_cosmu_muccmec    = new TH1D ( "h_cosmu_muccmec",    "", y_bins, 0.2, 2 );
    TH1D *h_cosmu_muccres1pi = new TH1D ( "h_cosmu_muccres1pi", "", y_bins, 0.2, 2 );
    TH1D *h_cosmu_mucc1pi    = new TH1D ( "h_cosmu_mucc1pi",    "", y_bins, 0.2, 2 );
    TH1D *h_cosmu_muccother  = new TH1D ( "h_cosmu_muccother",  "", y_bins, 0.2, 2 );
    TH1D *h_cosmu_mubarcc    = new TH1D ( "h_cosmu_mubarcc",    "", y_bins, 0.2, 2 );
    TH1D *h_cosmu_mubarnc    = new TH1D ( "h_cosmu_mubarnc",    "", y_bins, 0.2, 2 );
    
    signal_c_h.push_back(h_cosmu_muccqe);
    signal_c_h.push_back(h_cosmu_muccmec);
    signal_c_h.push_back(h_cosmu_muccres1pi);
    signal_c_h.push_back(h_cosmu_mucc1pi);
    signal_c_h.push_back(h_cosmu_muccother);
    signal_c_h.push_back(h_cosmu_mubarcc);
    signal_c_h.push_back(h_cosmu_mubarnc);

    leg_cos->AddEntry( h_cosmu_muccother,  " #nu_{#mu} CCOther ", "f" );
    //leg_cos->AddEntry( h_cosmu_mucc1pi,    " #nu_{#mu} CC non-RES 1#pi ", "f" );
    //leg_cos->AddEntry( h_cosmu_muccres1pi, " #nu_{#mu} CCRES ", "f" );
    leg_cos->AddEntry( h_cosmu_muccmec,    " #nu_{#mu} CCMEC ", "f" );
    leg_cos->AddEntry( h_cosmu_muccqe,     " #nu_{#mu} CCQE ", "f" );
    
    //leg_cos->AddEntry( h_cosmu_mubarcc,    " #bar{ #nu_{#mu} } CC ", "f" );
    //leg_cos->AddEntry( h_cosmu_mubarnc,    " #bar{ #nu_{#mu} } NC ", "f" );
   
    for ( int i = 1; i <= x_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_cos;
        double up_edge_cos;

        low_edge_cos = h_smeared->GetXaxis()->GetBinLowEdge(i);
        up_edge_cos  = h_smeared->GetXaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        // Filenames for signal histograms
        stringstream conv_c;
        conv_c.clear();

        string title_c;
        title_c.clear();

        char file_name_c[1024];

        conv_c << setprecision(4) << "working_dir/categorisation/" << pr_path << "/" << pr_path << "_pre_FSI_cosmu_slice_" << low_edge_cos << "_" << up_edge_cos << ".png";
        title_c = conv_c.str();
        
        strcpy( file_name_c, title_c.c_str() );

        // For the title of the histogram
        stringstream hist_c;
        hist_c.clear();

        char hist_name_c[1024];
        string temp_c;

        hist_c << setprecision(4) << "cos#theta_{#mu} in [" << low_edge_cos << "," << up_edge_cos << "]";
        temp_c = hist_c.str();

        strcpy( hist_name_c, temp_c.c_str() );


        for( unsigned int i = 0; i < signal_c_h.size(); ++i ){
            signal_c_h[i]->Reset("M");
        }  

        THStack *hscos = new THStack("hscos","pre_FSI");
            
        for ( int k = 0; k < n_entries; ++k ){
       
            info->GetEntry(k);

            int ev     = l_ev->GetValue();
            double cth = l_cth->GetValue();
            double T   = l_T->GetValue();
            int neu    = l_neu->GetValue();
            int cc     = l_cc->GetValue();
            int nc     = l_nc->GetValue();
            int qel    = l_qel->GetValue();
            int res    = l_res->GetValue();
            int mec    = l_mec->GetValue();
            int npip   = l_npip->GetValue();
            int npim   = l_npim->GetValue();
            int npi0   = l_npi0->GetValue();

            // Muon neutrino in this bin
            if ( cth >= low_edge_cos 
              && cth < up_edge_cos
              && neu == 14
              && cc ){
                
                // CCQE
                if ( qel ) {
                    h_cosmu_muccqe->Fill( T, 1 );
                }
                
                // CCMEC
                else if ( mec ){
                    h_cosmu_muccmec->Fill( T, 1 );
                }
                
                // CCRES, 1Pi
                else if ( res 
                     && ( npip + npim == 1 ) ){
                    h_cosmu_muccres1pi->Fill( T, 1 );
                }
                // CC non-RES, 1Pi
                else if ( !res
                     && ( npip + npim == 1 ) ){
                    h_cosmu_mucc1pi->Fill( T, 1 );
                }
                
                // CCOther
                else{
                    h_cosmu_muccother->Fill( T, 1 );
                }
            }    
            else if ( cth >= low_edge_cos
                   && cth < up_edge_cos
                   && neu == -14 ){
                
                // Numubar CC
                if ( cc ){
                    h_cosmu_mubarcc->Fill( T, 1 );
                }

                // Numubar NC
                else if ( nc ){
                    h_cosmu_mubarnc->Fill( T, 1 );
                }
            } 
        }
        
        int pal[12];
        pal[0]  = kRed;
        pal[1]  = kGreen + 2;
        pal[2]  = kOrange + 7;
        pal[3]  = kBlue;
        pal[4]  = kMagenta + 1;
        pal[5]  = kCyan + 2;
        pal[6]  = kYellow + 1;
        pal[7]  = kAzure + 2;
        pal[8]  = kViolet + 2;
        pal[9]  = kTeal - 1;
        pal[10] = kPink + 2;
        pal[11] = kSpring + 2;
        
        gStyle->SetPalette( 12, pal );
        gStyle->SetHatchesLineWidth( 1 );
        gStyle->SetHatchesSpacing( 0.5 );
        gStyle->SetTitleOffset(1.5, "Y");
        gStyle->SetOptStat(0);

        double norm_cosmu = 0;
       
        for ( unsigned int i = 0; i < signal_c_h.size() - 2; ++i ){
            norm_cosmu += signal_c_h[i]->Integral();
        }
        
        for ( unsigned int i = 0; i < signal_c_h.size() - 2; ++i ){

            signal_c_h[i]->SetFillColor(pal[i]);
            signal_c_h[i]->SetLineColor(pal[i]);

            signal_c_h[i]->SetLineWidth(1.5);
            signal_c_h[i]->Scale(1/norm_cosmu);

            hscos->Add(signal_c_h[i]);
        }
        
        // Fill the stacks 
        // Signal
        hscos->Draw();
        
        hscos->SetTitle(hist_name_c);
        hscos->GetXaxis()->SetTitle("T_{#mu} [GeV]");   
        hscos->GetYaxis()->SetTitle("Normalised event rate");   
    
        leg_cos->Draw();
        c_cosmu->SaveAs(file_name_c);

        delete hscos;
    } 

    for ( unsigned int i = 0; i < signal_c_h.size(); ++i ){
        delete signal_c_h[i];
    }
    
    delete c_cosmu;
    
    delete leg_cos;
*/
}

#endif
