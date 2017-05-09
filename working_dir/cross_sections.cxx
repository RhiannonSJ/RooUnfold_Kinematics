#include <iostream>
using std::cout;
using std::endl;

#include "roo_unfold.h"

#include "../src/RooUnfoldResponse.h"
#include "../src/RooUnfoldBayes.h"


void GetTruth( TTree *tree, int n_protons, std::vector<double> & truth_T, std::vector<double> & truth_cos, std::vector<bool> & truth_detectable, std::vector<double> & impur_T, std::vector<double> & impur_cos );

void Smear( const std::vector<double> & truth_T, const std::vector<double> & truth_cos, std::vector<double> & smear_T, std::vector<double> & smear_cos ); 

void GetResponse( const std::vector<double> & truth_T, const std::vector<double> & truth_cos, const std::vector<bool> & truth_detectable, const std::vector<double> & smear_T, const std::vector<double> & smear_cos, const std::vector<double> & smear_impur_T, const std::vector<double> & smear_impur_cos, RooUnfoldResponse & response); 

void Slices( TH2D *h_unfolded, TH2D *h_true, TH2D *h_reco, const char n_pr[1024] );

void cross_sections() {
    //==============================================================================
    // Reading in the event root files 
    //==============================================================================
    TFile f_test("/hepstore/rjones/Exercises/Flavours/Default+MEC/sbnd/1M/gntp.10000.gst.root");
    if(f_test.IsZombie()){
        std::cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "======================== Default + MEC event file open ========================" << endl;
    }
 
    TFile f_train("/hepstore/rjones/Exercises/Flavours/Default/sbnd/1M/gntp.10000.gst.root");
    if(f_train.IsZombie()){
        std::cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "=========================== Default event file open ===========================" << endl;
    }

    //==============================================================================
    // Reading in the flux root file 
    //==============================================================================
    TFile fflux("/hepstore/rjones/Exercises/Fluxes/sbn_FHC_flux_hist.root");
    if(fflux.IsZombie()){
        std::cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "============================ SBND flux file open ==============================" << endl;
    }

    /*
    //==============================================================================
    // Reading in the efficiency root file 
    //==============================================================================
    TFile feff("/hepstore/rjones/Exercises/RooUnfold_Kinematics/working_dir/efficiency.root");
    if(feff.IsZombie()){
        std::cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "========================== Efficiency flux file open ==========================" << endl;
    }
    
    //==============================================================================
    // Get the efficiency flux histogram
    //==============================================================================
    
    TH2D *h_eff = (TH2D*) feff.Get("h_eff");

    */

    //==============================================================================
    // Get the SBND flux histogram
    //==============================================================================
    
    TH1D *h_flux = (TH1D*) fflux.Get("h_numu_110m");

    TTree *gst_train = (TTree*) f_train.Get("gst");
    TTree *gst_test  = (TTree*) f_test.Get("gst");

    
    //==============================================================================
    // Get the truth vectors 
    //==============================================================================
    std::vector<double> truth_T_train; 
    std::vector<double> truth_cos_train; 
    std::vector<bool>   truth_detectable_train;
    std::vector<double> impur_T_train; 
    std::vector<double> impur_cos_train; 

    GetTruth( gst_train, 2, truth_T_train, truth_cos_train, truth_detectable_train, impur_T_train, impur_cos_train );

    std::vector<double> truth_T_test; 
    std::vector<double> truth_cos_test; 
    std::vector<bool>   truth_detectable_test;
    std::vector<double> impur_T_test; 
    std::vector<double> impur_cos_test; 

    GetTruth( gst_test, 2,  truth_T_test, truth_cos_test, truth_detectable_test, impur_T_test, impur_cos_test );
    
    //==============================================================================
    // Smearing 
    //==============================================================================
    
    std::vector<double> smear_T_train; 
    std::vector<double> smear_cos_train; 
    std::vector<double> impur_smear_T_train; 
    std::vector<double> impur_smear_cos_train; 
    
    Smear( truth_T_train, truth_cos_train, smear_T_train, smear_cos_train );
    Smear( impur_T_train, impur_cos_train, impur_smear_T_train, impur_smear_cos_train );
    
    std::vector<double> smear_T_test; 
    std::vector<double> smear_cos_test; 
    std::vector<double> impur_smear_T_test; 
    std::vector<double> impur_smear_cos_test; 
    
    Smear( truth_T_test, truth_cos_test, smear_T_test, smear_cos_test );
    Smear( impur_T_test, impur_cos_test, impur_smear_T_test, impur_smear_cos_test );
    
    //==============================================================================
    // Train response matrix
    //==============================================================================

    TH2D *h_true_train  = new TH2D( "h_true_train", " true ", 20, -1, 1, 18, 0, 2 );     
    TH2D *h_reco_train  = new TH2D( "h_reco_train", " reco ", 20, -1, 1, 18, 0, 2 );     

    RooUnfoldResponse response( h_reco_train, h_true_train );

    GetResponse( truth_T_train, truth_cos_train, truth_detectable_train, smear_T_train, smear_cos_train, impur_smear_T_train, impur_smear_cos_train, response); 
            
    //==============================================================================
    // Test unfolding
    //==============================================================================
   
    TCanvas * c = new TCanvas();

    TH2D *h_true_test  = new TH2D( "h_true_test", " true ", 20, -1, 1, 18, 0, 2 );     
    TH2D *h_cut_test   = new TH2D( "h_cut_test",  " cut  ", 20, -1, 1, 18, 0, 2 );     
    TH2D *h_reco_test  = new TH2D( "h_reco_test", " reco ", 20, -1, 1, 18, 0, 2 );     
    
    for ( unsigned int i = 0; i < truth_T_test.size(); ++i ) {
        h_true_test->Fill( truth_cos_test[i], truth_T_test[i]);

        if ( truth_detectable_test[i] ) {
            h_reco_test->Fill( smear_cos_test[i], smear_T_test[i]);
            h_cut_test->Fill( smear_cos_test[i], smear_T_test[i]);
        }
    }

    for ( unsigned int i = 0; i < impur_smear_T_test.size(); ++i ) {
        h_reco_test->Fill( impur_smear_cos_test[i], impur_smear_T_test[i]);
    }

    RooUnfoldBayes unfold( &response, h_reco_test, 1 );
    TH2D *h_unfold_test =  (TH2D*) unfold.Hreco();

    gStyle->SetPalette(55);

    // TRUE
    h_true_test->SetStats(kFALSE);
    h_true_test->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_true_test->GetYaxis()->SetTitle("T_{#mu}");
    h_true_test->SetTitle("True #mu kinematics");
    h_true_test->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/2_p/2D_true.png" );
    
    // CUT
    h_cut_test->SetStats(kFALSE);
    h_cut_test->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_cut_test->GetYaxis()->SetTitle("T_{#mu}");
    h_cut_test->SetTitle("#mu kinematics with cuts and smearing");
    h_cut_test->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/2_p/2D_cuts.png" );
    
    // RECO
    h_reco_test->SetStats(kFALSE);
    h_reco_test->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_reco_test->GetYaxis()->SetTitle("T_{#mu}");
    h_reco_test->SetTitle("Reco kinematics");
    h_reco_test->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/2_p/2D_reco.png" );
    
    // UNFOLDED
    h_unfold_test->SetStats(kFALSE);
    h_unfold_test->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_unfold_test->GetYaxis()->SetTitle("T_{#mu}");
    h_unfold_test->SetTitle("Unfolded #mu kinematics");
    h_unfold_test->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/2_p/2D_unfolded.png" );
   
    // COMPARISON
    TH2D *h_comp = new TH2D( *h_unfold_test );
    h_comp->Add( h_true_test, -1 );
    h_comp->Divide( h_true_test );
    h_comp->Scale( 100 );

    h_comp->SetStats(kFALSE);
    h_comp->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_comp->GetYaxis()->SetTitle("T_{#mu}");
    h_comp->SetTitle("CC0#pi, bin content difference between true and unfolded");
    h_comp->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/2_p/true_unfolding_comp.png" );

    // EFFICIENCY
    TH2D *h_eff = new TH2D( *h_cut_test );
    h_eff->Divide( h_true_test );

    h_eff->SetStats(kFALSE);
    h_eff->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_eff->GetYaxis()->SetTitle("T_{#mu}");
    h_eff->SetTitle("Efficiency");
    h_eff->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/2_p/efficiency.png" );


    // CROSS-SECTIONS
  
    // Normalization variables
    double cos_bins = h_unfold_test->GetXaxis()->GetBinWidth(1);
    double Tmu_bins = h_unfold_test->GetYaxis()->GetBinWidth(1);
    double flux_int = h_flux->Integral();
    double Na       = 6.63e34;
    double M_fid    = 1.016e8; // grams
    double A_ar     = 39.948;
    double tot_tgt  = ( Na * M_fid ) / ( A_ar ); // given in flux unit definition  
    double barns    = 1e38;

    double scalar_norm = ( barns ) / ( cos_bins * Tmu_bins * flux_int * tot_tgt );
    
    TH2D *h_ddxsec = new TH2D( *h_unfold_test );
    h_ddxsec->Scale( scalar_norm );
    h_ddxsec->Divide( h_eff );

    h_ddxsec->SetStats(kFALSE);
    h_ddxsec->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_ddxsec->GetYaxis()->SetTitle("T_{#mu}");
    h_ddxsec->SetTitle("CC0#pi, d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]");
    h_ddxsec->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/2_p/2D_unfolding_ddxsec.png" );

    // SLICING
    TH2D *h_true_ddxsec = new TH2D( *h_true_test );
    h_true_ddxsec->Scale( scalar_norm );
    h_true_ddxsec->Divide( h_eff );

    TH2D *h_reco_ddxsec = new TH2D( *h_reco_test );
    h_reco_ddxsec->Scale( scalar_norm );
    h_reco_ddxsec->Divide( h_eff );

    Slices( h_ddxsec, h_true_ddxsec, h_reco_ddxsec, "2_p" );

    delete h_true_train;
    delete h_reco_train;
    delete h_true_test;
    delete h_reco_test;
    delete h_cut_test;
    delete h_comp;
    delete h_eff;
    delete h_true_ddxsec;
    delete h_reco_ddxsec; 
    delete h_ddxsec;
}


void GetTruth( TTree *tree, int n_protons, std::vector<double> & truth_T, std::vector<double> & truth_cos, std::vector<bool> & truth_detectable, std::vector<double> & impur_T, std::vector<double> & impur_cos ) {

    TBranch *b_Ef    = tree->GetBranch("Ef");
    TBranch *b_El    = tree->GetBranch("El");
    TBranch *b_pl    = tree->GetBranch("pl");
    TBranch *b_cthl  = tree->GetBranch("cthl");
    TBranch *b_cthf  = tree->GetBranch("cthf");
    TBranch *b_nfpi0 = tree->GetBranch("nfpi0");
    TBranch *b_nfpip = tree->GetBranch("nfpip");
    TBranch *b_nfpim = tree->GetBranch("nfpim");
    TBranch *b_nfp   = tree->GetBranch("nfp");
    TBranch *b_cc    = tree->GetBranch("cc");
    TBranch *b_nc    = tree->GetBranch("nc");
    TBranch *b_nf    = tree->GetBranch("nf");
    TBranch *b_fspl  = tree->GetBranch("fspl");
    TBranch *b_pdgf  = tree->GetBranch("pdgf");

    // Variables 
    double m_mu = 0.10566; // Muon mass, GeV
    double m_pi = 0.13957; // Charged pion mass, GeV
    double m_pr = 0.93827; // Proton mass, GeV
   
    double T_thresh_mu = 0.05;
    double T_thresh_pi = 0.05;
    double T_thresh_pr = 0.05;

    int impurity = 5;

    int n_entries = tree->GetEntries();

    // Fill the vectors
    for ( int i = 0; i < n_entries; ++i ){
    
        tree->GetEntry(i);
    
        int nf    = b_nf->GetLeaf("nf")->GetValue();
        int fspl  = b_fspl->GetLeaf("fspl")->GetValue();
        int pdgf  = b_pdgf->GetLeaf("pdgf")->GetValue();
        int nfpi0 = b_nfpi0->GetLeaf("nfpi0")->GetValue();
        int nfpip = b_nfpip->GetLeaf("nfpip")->GetValue();
        int nfpim = b_nfpim->GetLeaf("nfpim")->GetValue();
        int nfp   = b_nfp->GetLeaf("nfp")->GetValue();
        int cc    = b_cc->GetLeaf("cc")->GetValue();
        int nc    = b_nc->GetLeaf("nc")->GetValue();

        double e_mu   = b_El->GetLeaf("El")->GetValue();
        double cos_mu = b_cthl->GetLeaf("cthl")->GetValue();
        
        int n_detectable_pr = 0;

        // Count how many protons are above 50MeV
        for ( int j = 0; j < nf; ++j ) {
 
            b_pdgf->GetEntry(i);
            b_Ef->GetEntry(i);
 
            int pdgf      = b_pdgf->GetLeaf("pdgf")->GetValue(j);
            double e_pr   = b_Ef->GetLeaf("Ef")->GetValue(j);
 
 
            // Calculate the kinetic energy of the pions
            if (pdgf == 2212 ){
 
                double T_pr = e_pr - m_pr;
                if( T_pr > T_thresh_pr ){
                    n_detectable_pr++;
                }
            }
        }

        if ( ( n_protons != -1 && nfp == n_protons ) || ( n_protons == -1 ) ){
            // If there are no pions fill the kinetic energy and cos theta with the muon energy   
            if ( cc == 1 && nfpip + nfpim + nfpi0 == 0 ){
                // Calculate the kinetic energy for muons
                if ( fspl == 13 ){
 
                    // Energy of the final state primary lepton
                    double T_mu = e_mu - m_mu;
 
                    bool isDetectable = false;
                    if ( n_protons == -1 ) {
                        isDetectable = (T_mu > T_thresh_mu);
                    }
                    else{
                        isDetectable = (T_mu > T_thresh_mu && n_protons == n_detectable_pr);
                    }

                    truth_T.push_back(T_mu);
                    truth_cos.push_back(cos_mu);
                    truth_detectable.push_back( isDetectable );
                }
            }

            if ( nc == 1 && nfpip + nfpim == 1 && nfpi0 == 0 ){
 
                // For all the final state hadronic particles, get their pdg code
                for ( int j = 0; j < nf; ++j ) {
 
                    b_pdgf->GetEntry(i);
                    b_cthf->GetEntry(i);
                    b_Ef->GetEntry(i);
 
                    int pdgf      = b_pdgf->GetLeaf("pdgf")->GetValue(j);
                    double e_pi   = b_Ef->GetLeaf("Ef")->GetValue(j);
                    double cos_pi = b_cthf->GetLeaf("cthf")->GetValue(j);
 
 
                    // Calculate the kinetic energy of the pions
                    if (pdgf == 211 || pdgf == -211 ){
 
                        double T_pi = e_pi - m_pi;
 
                        bool isDetectable = false;
                        if ( n_protons == -1 ) {
                            isDetectable = (T_pi > T_thresh_pi);
                        }
                        else{
                            isDetectable = (T_pi > T_thresh_pi && n_protons == n_detectable_pr);
                        }

                        if ( isDetectable ) {
                            int random;
                            random = rand() % impurity + 1;

                            if( random == 1 ){
                                impur_T.push_back(T_pi);
                                impur_cos.push_back(cos_pi);
                            }
                        }
                    }
                }
            }
        }
    }
}

void Smear( const std::vector<double> & truth_T, const std::vector<double> & truth_cos, std::vector<double> & smear_T, std::vector<double> & smear_cos ) {

    // Initiate the random number generation
    ROOT::Math::GSLRngMT *_random_gen = new ROOT::Math::GSLRngMT;
    _random_gen->Initialize();
    _random_gen->SetSeed( time( NULL ) );
 
    int n_values = truth_T.size();

    if ( truth_T.size() != truth_cos.size() ){
        std::cerr << " Vectors must be the same length " << endl;
        exit(1);
    }
    // Event by event, generate Tmu_prime and Tpi_prime: lognormal
    // Then find thetamu_prime and thetapi_prime: gaussian
    for ( int i = 0; i < n_values; ++i ){
 
        // -------------------------------------------------------
        //                   Kinetic energy
        // -------------------------------------------------------
        // Calculate the mean and sigma for the LogNormal function
        //     zeta  = TMath::Log( m * ( 1 / sqrt( 1 + ( var / pow( m, 2 ) ) ) ) );
        //     sigma = sqrt( log( 1 + ( var / pow( m, 2 ) ) ) );
        //     m     = expectation value = Tl
        //     var   = variance = s.d.^2 = ( Tl * 0.1 ) ^ 2
 
        double var_mu     = TMath::Power( truth_T[i] * 0.1, 2 );
        double sigma_mu   = TMath::Sqrt( TMath::Log( 1 + ( var_mu / TMath::Power(truth_T[i], 2 ) ) ) );
        double zeta_mu    = TMath::Log( truth_T[i] * ( 1 / TMath::Sqrt( 1 + ( var_mu / TMath::Power( truth_T[i], 2 ) ) ) ) );
        double lognorm_mu = _random_gen->LogNormal( zeta_mu, sigma_mu );
     
        // -------------------------------------------------------
        //                  Cos theta
        // -------------------------------------------------------
     
        // Calculate the mean and sigma for the LogNormal function
        //      theta = acos(costtheta)
        //      var   = 5 degrees
     
        double sd_thetamu    = TMath::Pi() / 36; // 5 degrees
        double gaus_theta    = TMath::ACos( truth_cos[i] ) + _random_gen->Gaussian( sd_thetamu );
        double gaus_costheta = TMath::Cos( gaus_theta );
     
        smear_T.push_back(lognorm_mu);
        smear_cos.push_back(gaus_costheta);

    }
}

void GetResponse( const std::vector<double> & truth_T, const std::vector<double> & truth_cos, const std::vector<bool> & truth_detectable, const std::vector<double> & smear_T, const std::vector<double> & smear_cos, const std::vector<double> & smear_impur_T, const std::vector<double> & smear_impur_cos, RooUnfoldResponse & response) {

    for ( unsigned int i = 0; i < truth_T.size(); ++i ) {
        if ( truth_detectable[i] ) {
            response.Fill( smear_cos[i], smear_T[i], truth_cos[i], truth_T[i] ); 
        }
        else{
            response.Miss( truth_cos[i], truth_T[i] );
        }
    }

    for ( unsigned int i = 0; i < smear_impur_T.size(); ++i ) {
        response.Fake( smear_impur_cos[i], smear_impur_T[i] );
    }
} 

void Slices ( TH2D *h_unfolded, TH2D *h_true, TH2D *h_reco, const char n_pr[1024] ){

    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    // Take the bin edges to be the title
    int x_bins = h_unfolded->GetNbinsX(); // Cos theta
    int y_bins = h_unfolded->GetNbinsY(); // Tmu

    TCanvas *c_Tmu   = new TCanvas ( "c_Tmu", "", 800, 600 );
   
    TLegend *leg_T   = new TLegend( 0.12, 0.78, 0.28, 0.88 );

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

        conv << setprecision(4) << "working_dir/unfolded_distributions/" << n_pr << "/Tmu_slice_" << low_edge_T << ";" << up_edge_T << ".png";
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
        h_Tmu->GetYaxis()->SetTitle("d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]");   
        h_Tmu->SetLineColor( kRed + 2 );
        h_Tmu_true->SetLineColor( kGreen + 2 );
        h_Tmu_reco->SetLineColor( kBlue + 2 );
        h_Tmu_reco->SetLineStyle( 7 );

        h_Tmu->SetTitleOffset(1.1, "Y");
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
    
    TLegend *leg_c     = new TLegend( 0.72, 0.78, 0.88, 0.88 );

    TH1D *h_cosmu      = new TH1D ( "h_cosmu", "", y_bins, 0, 2 );
    TH1D *h_cosmu_true = new TH1D ( "h_cosmu_true", "", y_bins, 0, 2 );
    TH1D *h_cosmu_reco = new TH1D ( "h_cosmu_reco", "", y_bins, 0, 2 );
     
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

        conv1 << setprecision(4) << "working_dir/unfolded_distributions/" << n_pr << "/cos_thetamu_slice_" << low_edge_cos << ";" << up_edge_cos << ".png";
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
        h_cosmu->GetYaxis()->SetTitle("d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]");   
        h_cosmu->SetLineColor( kRed + 2 );
        h_cosmu_true->SetLineColor( kGreen + 2 );
        h_cosmu_reco->SetLineColor( kBlue + 2 );
        h_cosmu_reco->SetLineStyle( 7 );
        h_cosmu->SetTitleOffset(1.1, "Y");
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

