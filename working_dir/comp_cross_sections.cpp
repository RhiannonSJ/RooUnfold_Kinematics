#include <iostream>
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/working_dir/roo_unfold.h"

void Slices( TH2D *h_all, TH2D *h_1p, TH2D *h_2p, TH2D *h_3p, TH2D *h_a, TH2D *h_1, TH2D *h_2, TH2D *h_3, double norm );

void comp_cross_secs(){

    // ===========================================================================
    // READ IN FILE
    // ===========================================================================

    TFile f_hists("/hepstore/rjones/Exercises/RooUnfold_Kinematics/working_dir/selection_ddxsec.root");
    if(f_hists.IsZombie()){
        std::cerr << " Error opening file " << std::endl;
        exit(1);
    }
    else{
        std::cout << "======================== Histogram file open ========================" << std::endl;
    }
 
    TFile f_hists2("/hepstore/rjones/Exercises/RooUnfold_Kinematics/working_dir/flux_selection_ddxsec.root");
    if(f_hists2.IsZombie()){
        std::cerr << " Error opening file " << std::endl;
        exit(1);
    }
    else{
        std::cout << "======================== Histogram file open ========================" << std::endl;
    }
 
    //==============================================================================
    // Reading in the flux root file 
    //==============================================================================
    TFile fflux("/hepstore/rjones/Exercises/Fluxes/sbn_FHC_flux_hist.root");
    if(fflux.IsZombie()){
        std::cerr << " Error opening file " << std::endl;
        exit(1);
    }
    else{
        std::cout << "======================== SBND flux file open ========================" << std::endl;
    }

    // Get the histograms 
    
    TH1D *h_flux = (TH1D*) fflux.Get("h_numu_110m");
   
    TH2D *h_all = (TH2D*) f_hists.Get("all_signal");
    TH2D *h_1p  = (TH2D*) f_hists.Get("1_p");
    TH2D *h_2p  = (TH2D*) f_hists.Get("2_p");
    TH2D *h_3p  = (TH2D*) f_hists.Get("3_p");

    TH2D *h_all2 = (TH2D*) f_hists2.Get("all_signal");
    TH2D *h_1p2  = (TH2D*) f_hists2.Get("1_p");
    TH2D *h_2p2  = (TH2D*) f_hists2.Get("2_p");
    TH2D *h_3p2  = (TH2D*) f_hists2.Get("3_p");

    // Normalization variables
    double cos_bins = h_all->GetXaxis()->GetBinWidth(1);
    double Tmu_bins = h_all->GetYaxis()->GetBinWidth(1); // GeV
    // Units of flux: #/m^2/50MeV/10^6 POT
    double flux_int = h_flux->Integral();
    double Na       = 6.022e23;  //
    double M_fid    = 112000; // (kg / m^3) * m^3
    double A_ar     = 0.039948;  // kg / mol
    double tot_tgt  = ( Na * M_fid ) / ( A_ar ); // number of target particles  
    double POT      = 6.6e20;    // POT
    double scaling  = 5e8;      // Sort units from factor of 1 / ( integrated flux * Nt * bin widths * POT )
    double cm_conv  = 1e-38;     // take out factor of 10e-38 cm^2 

    // Playing with flux systematic
    double err_flux = 0.05; // 5%

    double scalar_norm = scaling / ( cm_conv * cos_bins * Tmu_bins * flux_int * tot_tgt * POT );

    TH2D *h_ddxsec_all = new TH2D( *h_all );
    
    TH2D *h_ddxsec_1p = new TH2D( *h_1p );
    
    TH2D *h_ddxsec_2p = new TH2D( *h_2p );
    
    TH2D *h_ddxsec_3p = new TH2D( *h_3p );
    
    // Slice the histograms
    Slices( h_all2, h_1p2, h_2p2, h_3p2, h_ddxsec_all, h_ddxsec_1p, h_ddxsec_2p, h_ddxsec_3p, scalar_norm );



}

void Slices( TH2D *h_all, TH2D *h_1p, TH2D *h_2p, TH2D *h_3p, TH2D *h_a, TH2D *h_1, TH2D *h_2, TH2D *h_3, double norm ){

    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    // Take the bin edges to be the title
    int x_bins = h_1->GetNbinsX(); // Cos theta
    int y_bins = h_1->GetNbinsY(); // Tmu

    TCanvas *c_Tmu   = new TCanvas ( "c_Tmu", "", 800, 600 );
   
    TLegend *leg_T   = new TLegend( 0.12, 0.68, 0.38, 0.88 );

    TH1D *h_Tmu_a  = new TH1D ( "h_Tmu_a", "", x_bins, -1, 1 );
    TH1D *h_Tmu_1  = new TH1D ( "h_Tmu_1", "", x_bins, -1, 1 );
    TH1D *h_Tmu_2  = new TH1D ( "h_Tmu_2", "", x_bins, -1, 1 );
    TH1D *h_Tmu_3  = new TH1D ( "h_Tmu_3", "", x_bins, -1, 1 );
    TH1D *h_Tmu_a2 = new TH1D ( "h_Tmu_a2", "", x_bins, -1, 1 );
    TH1D *h_Tmu_12 = new TH1D ( "h_Tmu_12", "", x_bins, -1, 1 );
    TH1D *h_Tmu_22 = new TH1D ( "h_Tmu_22", "", x_bins, -1, 1 );
    TH1D *h_Tmu_32 = new TH1D ( "h_Tmu_32", "", x_bins, -1, 1 );
    
    leg_T->AddEntry( h_Tmu_a, " CC 0#pi ", "l" );
    leg_T->AddEntry( h_Tmu_1, " CC 0#pi, 1p ", "l" );
    leg_T->AddEntry( h_Tmu_2, " CC 0#pi, 2p ", "l" );
    leg_T->AddEntry( h_Tmu_3, " CC 0#pi, 3p ", "l" );
    
    for ( int i = 1; i <= y_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_T;
        double up_edge_T;

        low_edge_T = h_a->GetYaxis()->GetBinLowEdge(i);
        up_edge_T = h_a->GetYaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        std::stringstream conv;
        conv.clear();

        std::string title;
        title.clear();

        char file_name[1024];

        conv << std::setprecision(4) << "/hepstore/rjones/Exercises/RooUnfold_Kinematics/working_dir/test_xsec/Tmu_slice_" << low_edge_T << ";" << up_edge_T << ".png";
        title = conv.str();
        
        strcpy( file_name, title.c_str() );

        // For the title of the histogram
        std::stringstream hist;
        hist.clear();

        char hist_name[1024];
        std::string temp1;

        hist << std::setprecision(4) << "T_{#mu} in [" << low_edge_T << "," << up_edge_T << "], GeV";
        temp1 = hist.str();

        strcpy( hist_name, temp1.c_str() );

        // Fill the histogram
        for ( int j = 1; j <= x_bins; ++j ){
            h_Tmu_a->SetBinContent( j, h_a->GetBinContent(j, i) );
            h_Tmu_1->SetBinContent( j, h_1->GetBinContent(j, i) );
            h_Tmu_2->SetBinContent( j, h_2->GetBinContent(j, i) );
            h_Tmu_3->SetBinContent( j, h_3->GetBinContent(j, i) );
            h_Tmu_a2->SetBinContent( j, h_all->GetBinContent(j, i) );
            h_Tmu_12->SetBinContent( j, h_1p->GetBinContent(j, i) );
            h_Tmu_22->SetBinContent( j, h_2p->GetBinContent(j, i) );
            h_Tmu_32->SetBinContent( j, h_3p->GetBinContent(j, i) );
            
            h_Tmu_a->SetBinError( j, h_a->GetBinError(j, i) );
            h_Tmu_1->SetBinError( j, h_1->GetBinError(j, i) );
            h_Tmu_2->SetBinError( j, h_2->GetBinError(j, i) );
            h_Tmu_3->SetBinError( j, h_3->GetBinError(j, i) );
            h_Tmu_a2->SetBinError( j, h_all->GetBinError(j, i) );
            h_Tmu_12->SetBinError( j, h_1p->GetBinError(j, i) );
            h_Tmu_22->SetBinError( j, h_2p->GetBinError(j, i) );
            h_Tmu_32->SetBinError( j, h_3p->GetBinError(j, i) );
        }

        
        h_Tmu_a->Scale( norm );
        h_Tmu_1->Scale( norm );
        h_Tmu_2->Scale( norm );
        h_Tmu_3->Scale( norm );
        

        double max_a = h_Tmu_a->GetBinContent(h_Tmu_a->GetMaximumBin());
        double max_1 = h_Tmu_1->GetBinContent(h_Tmu_1->GetMaximumBin());
        double max_2 = h_Tmu_2->GetBinContent(h_Tmu_2->GetMaximumBin());
        double max_3 = h_Tmu_3->GetBinContent(h_Tmu_3->GetMaximumBin());

        double max = 1.15 * TMath::Max(max_1, TMath::Max(max_2, TMath::Max( max_3, max_a ))); 

        h_Tmu_a->Draw("e1x0");
        h_Tmu_1->Draw("e1x0same");
        h_Tmu_2->Draw("e1x0same");
        h_Tmu_3->Draw("e1x0same");
        h_Tmu_a2->Draw("e1x0same");
        h_Tmu_12->Draw("e1x0same");
        h_Tmu_22->Draw("e1x0same");
        h_Tmu_32->Draw("e1x0same");
        h_Tmu_a->Draw("hist same");
        h_Tmu_1->Draw("hist same");
        h_Tmu_2->Draw("hist same");
        h_Tmu_3->Draw("hist same");
        h_Tmu_a->SetTitle(hist_name);
        h_Tmu_a->GetYaxis()->SetRangeUser(0,max);
        h_Tmu_a->GetXaxis()->SetTitle("cos#theta_{#mu}");   
        h_Tmu_a->GetYaxis()->SetTitle("CC0#pi, d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]");   
        h_Tmu_a->SetLineColor( kGreen + 2 );
        h_Tmu_1->SetLineColor( kBlack );
        h_Tmu_2->SetLineColor( kMagenta + 2 );
        h_Tmu_3->SetLineColor( kBlue + 2 );
        h_Tmu_a2->SetLineColor( kGreen );
        h_Tmu_12->SetLineColor( kBlack );
        h_Tmu_22->SetLineColor( kMagenta );
        h_Tmu_32->SetLineColor( kBlue );

        h_Tmu_a->SetTitleOffset(1.4, "Y");
        h_Tmu_a->SetStats(kFALSE);

        leg_T->Draw();
        c_Tmu->SaveAs(file_name);

    } 
   
    delete h_Tmu_a;
    delete h_Tmu_1;
    delete h_Tmu_2;
    delete h_Tmu_3;

    delete h_Tmu_a2;
    delete h_Tmu_12;
    delete h_Tmu_22;
    delete h_Tmu_32;

    delete c_Tmu;
    
    delete leg_T;

    TCanvas *c_cosmu   = new TCanvas ( "c_cosmu", "", 800, 600 );
    
    TLegend *leg_c     = new TLegend( 0.62, 0.68, 0.88, 0.88 );

    TH1D *h_cosmu_a = new TH1D ( "h_cosmu_a", "", y_bins, 0, 2 );
    TH1D *h_cosmu_1 = new TH1D ( "h_cosmu_1", "", y_bins, 0, 2 );
    TH1D *h_cosmu_2 = new TH1D ( "h_cosmu_2", "", y_bins, 0, 2 );
    TH1D *h_cosmu_3 = new TH1D ( "h_cosmu_3", "", y_bins, 0, 2 );
    TH1D *h_cosmu_a2 = new TH1D ( "h_cosmu_a2", "", y_bins, 0, 2 );
    TH1D *h_cosmu_12 = new TH1D ( "h_cosmu_12", "", y_bins, 0, 2 );
    TH1D *h_cosmu_22 = new TH1D ( "h_cosmu_22", "", y_bins, 0, 2 );
    TH1D *h_cosmu_32 = new TH1D ( "h_cosmu_32", "", y_bins, 0, 2 );
     
    leg_c->AddEntry( h_cosmu_a, " CC 0#pi ", "l" );
    leg_c->AddEntry( h_cosmu_1, " CC 0#pi, 1p ", "l" );
    leg_c->AddEntry( h_cosmu_2, " CC 0#pi, 2p ", "l" );
    leg_c->AddEntry( h_cosmu_3, " CC 0#pi, 3p ", "l" );
    
    // Cos theta mu slices
    for ( int i = 1; i <= x_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_cos;
        double up_edge_cos;

        low_edge_cos = h_1->GetXaxis()->GetBinLowEdge(i);
        up_edge_cos = h_1->GetXaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        std::stringstream conv1;
        conv1.clear();

        std::string title1;
        title1.clear();

        char file_name1[1024];

        conv1 << std::setprecision(4) << "/hepstore/rjones/Exercises/RooUnfold_Kinematics/working_dir/test_xsec/cos_thetamu_slice_" << low_edge_cos << ";" << up_edge_cos << ".png";
        title1 = conv1.str();
        
        strcpy( file_name1, title1.c_str() );

        // For the title of the histogram
        std::stringstream hist1;
        hist1.clear();

        char hist_name1[1024];
        std::string temp2;

        hist1 << std::setprecision(4) << "cos#theta_{#mu} in [ " << low_edge_cos << "," << up_edge_cos << "]";
        temp2 = hist1.str();

        strcpy( hist_name1, temp2.c_str() );

        // Fill the histogram
        for ( int j = 1; j <= y_bins; ++j ){
            h_cosmu_a->SetBinContent( j, h_a->GetBinContent(i, j) );
            h_cosmu_1->SetBinContent( j, h_1->GetBinContent(i, j) );
            h_cosmu_2->SetBinContent( j, h_2->GetBinContent(i, j) );
            h_cosmu_3->SetBinContent( j, h_3->GetBinContent(i, j) );
            h_cosmu_a2->SetBinContent( j, h_all->GetBinContent(i, j) );
            h_cosmu_12->SetBinContent( j, h_1p->GetBinContent(i, j) );
            h_cosmu_22->SetBinContent( j, h_2p->GetBinContent(i, j) );
            h_cosmu_32->SetBinContent( j, h_3p->GetBinContent(i, j) );
            
            h_cosmu_a->SetBinError( j, h_a->GetBinError(i, j) );
            h_cosmu_1->SetBinError( j, h_1->GetBinError(i, j) );
            h_cosmu_2->SetBinError( j, h_2->GetBinError(i, j) );
            h_cosmu_3->SetBinError( j, h_3->GetBinError(i, j) );
            h_cosmu_a2->SetBinError( j, h_all->GetBinError(i, j) );
            h_cosmu_12->SetBinError( j, h_1p->GetBinError(i, j) );
            h_cosmu_22->SetBinError( j, h_2p->GetBinError(i, j) );
            h_cosmu_32->SetBinError( j, h_3p->GetBinError(i, j) );
        }

        h_cosmu_a->Scale( norm ); 
        h_cosmu_1->Scale( norm ); 
        h_cosmu_2->Scale( norm ); 
        h_cosmu_3->Scale( norm );
        
        double max_a_c = h_cosmu_a->GetBinContent(h_cosmu_a->GetMaximumBin());
        double max_1_c = h_cosmu_1->GetBinContent(h_cosmu_1->GetMaximumBin());
        double max_2_c = h_cosmu_2->GetBinContent(h_cosmu_2->GetMaximumBin());
        double max_3_c = h_cosmu_3->GetBinContent(h_cosmu_3->GetMaximumBin());

        double max_c = 1.15 * TMath::Max(max_1_c, TMath::Max(max_2_c, TMath::Max(max_3_c, max_a_c))); 
        

        h_cosmu_a->Draw(" hist ");
        h_cosmu_1->Draw("hist same");
        h_cosmu_2->Draw("hist same");
        h_cosmu_3->Draw("hist same");
        h_cosmu_a->Draw("e1x0same");
        h_cosmu_1->Draw("e1x0same");
        h_cosmu_2->Draw("e1x0same");
        h_cosmu_3->Draw("e1x0same");
        h_cosmu_a2->Draw("e1x0same");
        h_cosmu_12->Draw("e1x0same");
        h_cosmu_22->Draw("e1x0same");
        h_cosmu_32->Draw("e1x0same");
        h_cosmu_a->SetTitle(hist_name1);
        h_cosmu_a->GetYaxis()->SetRangeUser(0,max_c);
        h_cosmu_a->GetXaxis()->SetTitle("T_{#mu} [GeV]");   
        h_cosmu_a->GetYaxis()->SetTitle("CC0#pi, d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]");   
        h_cosmu_a->SetLineColor( kGreen + 2 );
        h_cosmu_1->SetLineColor( kBlack );
        h_cosmu_2->SetLineColor( kMagenta + 2 );
        h_cosmu_3->SetLineColor( kBlue + 2 );
        h_cosmu_a2->SetLineColor( kGreen );
        h_cosmu_12->SetLineColor( kBlack  );
        h_cosmu_22->SetLineColor( kMagenta  );
        h_cosmu_32->SetLineColor( kBlue  );
        
        h_cosmu_a->SetTitleOffset(1.4, "Y");
        h_cosmu_a->SetStats(kFALSE);

        leg_c->Draw();
        c_cosmu->SaveAs(file_name1);

    } 
    
    delete h_cosmu_a;
    delete h_cosmu_1;
    delete h_cosmu_2;
    delete h_cosmu_3;

    delete c_cosmu;

    delete leg_c;
}
