//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#include <iostream>
using std::cout;
using std::endl;

#include "roo_unfold.h"

#include "../src/RooUnfoldResponse.h"
#include "../src/RooUnfoldBayes.h"

//==============================================================================
// Slice the unfolded histogram
//==============================================================================

void FillVect ( TTree *event, 
                int n_protons,
                vector< int >    &impur,
                vector< double > &Tmu_vect, 
                vector< double > &cosmu_vect,
                vector< double > &Tpi_vect, 
                vector< double > &cospi_vect ){

    TBranch *b_Ef    = event->GetBranch("Ef");
    TBranch *b_El    = event->GetBranch("El");
    TBranch *b_pl    = event->GetBranch("pl");
    TBranch *b_cthl  = event->GetBranch("cthl");
    TBranch *b_cthf  = event->GetBranch("cthf");
    TBranch *b_nfpi0 = event->GetBranch("nfpi0");
    TBranch *b_nfpip = event->GetBranch("nfpip");
    TBranch *b_nfpim = event->GetBranch("nfpim");
    TBranch *b_nfp   = event->GetBranch("nfp");
    TBranch *b_cc    = event->GetBranch("cc");
    TBranch *b_nc    = event->GetBranch("nc");
    TBranch *b_nf    = event->GetBranch("nf");
    TBranch *b_fspl  = event->GetBranch("fspl");
    TBranch *b_pdgf  = event->GetBranch("pdgf");

    // Variables 
    double m_mu = 0.10566; // Muon mass, GeV
    double m_pi = 0.13957; // Charged pion mass, GeV
    double m_pr = 0.93827; // Proton mass, GeV
    
    int n_entries = event->GetEntries();

    // Fill the vectors
    for ( int i = 0; i < n_entries; ++i ){
    
        event->GetEntry(i);
    
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
        
        double T_mu;

        int p_above_thresh = 0;

        // Count how many protons are above 50MeV
        for ( int j = 0; j < nf; ++j ) {
 
            b_pdgf->GetEntry(i);
            b_cthf->GetEntry(i);
            b_Ef->GetEntry(i);
 
            int pdgf      = b_pdgf->GetLeaf("pdgf")->GetValue(j);
            double e_pr   = b_Ef->GetLeaf("Ef")->GetValue(j);
 
 
            // Calculate the kinetic energy of the pions
            if (pdgf == 2212 ){
 
                double T_pr = e_pr - m_pr;
                if( T_pr > 0.05 ){
                    p_above_thresh++;
                }
            }
        }


        if ( n_protons != -1 ){
            // If there are no pions fill the kinetic energy and cos theta with the muon energy   
            if ( p_above_thresh == n_protons && cc == 1 && nfpip + nfpim + nfpi0 == 0 ){
                // Calculate the kinetic energy for muons
                if ( fspl == 13 ){
 
                    // Energy of the final state primary lepton
                    T_mu = e_mu - m_mu;
 
                    Tmu_vect.push_back(T_mu);
                    cosmu_vect.push_back(cos_mu);
 
                }
         
                // If the final state primary is a lepton, push back a number that will
                // be removed in the cuts later
                else if ( fspl == 11 ){
                    Tmu_vect.push_back(-99999);
                    cosmu_vect.push_back(-99999);
                }
 
                else{
                    Tmu_vect.push_back(-99999);
                    cosmu_vect.push_back(-99999);
                }
            }
            else{
                Tmu_vect.push_back(-99999);
                cosmu_vect.push_back(-99999);
            }

            if ( p_above_thresh == n_protons && nc == 1 && nfpip + nfpim == 1 && nfpi0 == 0 ){
 
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
 
                    Tpi_vect.push_back(T_pi);
                    cospi_vect.push_back(cos_pi);
            
                    }
                }
            }
            else{ 
                Tpi_vect.push_back(-99999);
                cospi_vect.push_back(-99999);
            }

            // Apply the impurity cut
            if ( p_above_thresh == n_protons && nc == 1 && nfpip + nfpim == 1 && nfpi0 == 0 ){
                int random;
                random = rand() % 5 + 1;

                // If the random number = 5 then push back a 1 onto the vector
                // If the random number < 5 push back a 0
                if( random == 5 ){
                    impur.push_back(1);
                }
                else{
                    impur.push_back(0);
                }
            }
            else{
                impur.push_back(0);
            } 
        }
        else{
            // If there are no pions fill the kinetic energy and cos theta with the muon energy   
            if ( cc == 1 && nfpip + nfpim + nfpi0 == 0 ){
                // Calculate the kinetic energy for muons
                if ( fspl == 13 ){
 
                    // Energy of the final state primary lepton
                    T_mu = e_mu - m_mu;
 
                    Tmu_vect.push_back(T_mu);
                    cosmu_vect.push_back(cos_mu);
 
                }
         
                // If the final state primary is a lepton, push back a number that will
                // be removed in the cuts later
                else if ( fspl == 11 ){
                    Tmu_vect.push_back(-99999);
                    cosmu_vect.push_back(-99999);
                }
 
                else{
                    Tmu_vect.push_back(-99999);
                    cosmu_vect.push_back(-99999);
                }
            }
            else{
                Tmu_vect.push_back(-99999);
                cosmu_vect.push_back(-99999);
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
 
                    Tpi_vect.push_back(T_pi);
                    cospi_vect.push_back(cos_pi);
            
                    }  
                }
            }
            else{ 
                Tpi_vect.push_back(-99999);
                cospi_vect.push_back(-99999);
            }

            // Apply the impurity cut
            if ( nc == 1 && nfpip + nfpim == 1 && nfpi0 == 0 ){
                int random;
                random = rand() % 5 + 1;

                // If the random number = 5 then push back a 1 onto the vector
                // If the random number < 5 push back a 0
                if( random == 5 ){
                    impur.push_back(1);
                }
                else{
                    impur.push_back(0);
                }
            }
            else{
                impur.push_back(0);
            } 
        }
    }
}

//==============================================================================
// Smear the kinetic energy and cos thetas
//==============================================================================

void Smearing( std::vector< double > T_mu, 
               std::vector< double > cos_mu, 
               std::vector< double > &T_mu_sm, 
               std::vector< double > &cos_mu_sm ){

    // Initiate the random number generation
    ROOT::Math::GSLRngMT *_random_gen = new ROOT::Math::GSLRngMT;
    _random_gen->Initialize();
    _random_gen->SetSeed( time( NULL ) );
 
    int n_values = T_mu.size();

    if ( T_mu.size() != cos_mu.size() ){
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
 
        if ( T_mu[i] > 0.05 ){
            double var_mu     = TMath::Power( T_mu[i] * 0.1, 2 );
            double sigma_mu   = TMath::Sqrt( TMath::Log( 1 + ( var_mu / TMath::Power(T_mu[i], 2 ) ) ) );
            double zeta_mu    = TMath::Log( T_mu[i] * ( 1 / TMath::Sqrt( 1 + ( var_mu / TMath::Power( T_mu[i], 2 ) ) ) ) );
            double lognorm_mu = _random_gen->LogNormal( zeta_mu, sigma_mu );
     
            // -------------------------------------------------------
            //                  Cos theta
            // -------------------------------------------------------
     
            // Calculate the mean and sigma for the LogNormal function
            //      theta = acos(costtheta)
            //      var   = 5 degrees
     
            double sd_thetamu    = TMath::Pi() / 36; // 5 degrees
            double gaus_theta    = TMath::ACos( cos_mu[i] ) + _random_gen->Gaussian( sd_thetamu );
            double gaus_costheta = TMath::Cos( gaus_theta );
     
            T_mu_sm.push_back(lognorm_mu);
            cos_mu_sm.push_back(gaus_costheta);
        }
        else{
            T_mu_sm.push_back(-99999);
            cos_mu_sm.push_back(-99999);
        }
    }
}

//==============================================================================
// Slice the unfolded histogram
//==============================================================================

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


//==============================================================================
// The main function
//==============================================================================
void unfolding_with_categorization() { 
    
    //==============================================================================
    // Reading in the event root files 
    //==============================================================================
    TFile f("/hepstore/rjones/Exercises/Flavours/Default+MEC/sbnd/1M/gntp.10000.gst.root");
    if(f.IsZombie()){
        std::cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "======================== Default + MEC event file open ========================" << endl;
    }

    TFile f1("/hepstore/rjones/Exercises/Flavours/Default/sbnd/1M/gntp.10000.gst.root");
    if(f.IsZombie()){
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
    if(f.IsZombie()){
        std::cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "============================ SBND flux file open ==============================" << endl;
    }

    //==============================================================================
    // Reading in the efficiency root file 
    //==============================================================================
    TFile feff("/hepstore/rjones/Exercises/Kinematics_Unfolding/working_dir/efficiency.root");
    if(f.IsZombie()){
        std::cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "========================== Efficiency flux file open ==========================" << endl;
    }
    
    //==============================================================================
    // Get the SBND flux histogram
    //==============================================================================
    
    TH1D *h_flux = (TH1D*) fflux.Get("h_numu_110m");

    //==============================================================================
    // Get the efficiency flux histogram
    //==============================================================================
    
    TH2D *h_eff = (TH2D*) feff.Get("h_eff");

    //==============================================================================
    // Get everything from the event tree
    //==============================================================================

    TTree *gst = (TTree*) f.Get("gst");
    TTree *gst_train = (TTree*) f1.Get("gst");

    TBranch *b_Ef    = gst->GetBranch("Ef");
    TBranch *b_El    = gst->GetBranch("El");
    TBranch *b_pl    = gst->GetBranch("pl");
    TBranch *b_cthl  = gst->GetBranch("cthl");
    TBranch *b_cthf  = gst->GetBranch("cthf");
    TBranch *b_nfpi0 = gst->GetBranch("nfpi0");
    TBranch *b_nfpip = gst->GetBranch("nfpip");
    TBranch *b_nfpim = gst->GetBranch("nfpim");
    TBranch *b_nfp   = gst->GetBranch("nfp");
    TBranch *b_cc    = gst->GetBranch("cc");
    TBranch *b_nc    = gst->GetBranch("nc");
    TBranch *b_nf    = gst->GetBranch("nf");
    TBranch *b_fspl  = gst->GetBranch("fspl");
    TBranch *b_pdgf  = gst->GetBranch("pdgf");

    TBranch *b_tr_Ef    = gst_train->GetBranch("Ef");
    TBranch *b_tr_El    = gst_train->GetBranch("El");
    TBranch *b_tr_pl    = gst_train->GetBranch("pl");
    TBranch *b_tr_cthl  = gst_train->GetBranch("cthl");
    TBranch *b_tr_cthf  = gst_train->GetBranch("cthf");
    TBranch *b_tr_nfpi0 = gst_train->GetBranch("nfpi0");
    TBranch *b_tr_nfpip = gst_train->GetBranch("nfpip");
    TBranch *b_tr_nfpim = gst_train->GetBranch("nfpim");
    TBranch *b_tr_nfp   = gst_train->GetBranch("nfp");
    TBranch *b_tr_cc    = gst_train->GetBranch("cc");
    TBranch *b_tr_nc    = gst_train->GetBranch("nc");
    TBranch *b_tr_nf    = gst_train->GetBranch("nf");
    TBranch *b_tr_fspl  = gst_train->GetBranch("fspl");
    TBranch *b_tr_pdgf  = gst_train->GetBranch("pdgf");

    // Define and fill vectors for the muon energy, momentum and opening angle
    // Define an empty vector to hold the reconstructed energy

    std::vector< double > T_mu_tr_vect;
    std::vector< double > T_pi_tr_vect;
    std::vector< double > cos_mu_tr_vect;
    std::vector< double > cos_pi_tr_vect;

    std::vector< double > T_mu_vect;
    std::vector< double > T_pi_vect;
    std::vector< double > cos_mu_vect;
    std::vector< double > cos_pi_vect;

    std::vector< int > Impurity_tr;
    std::vector< int > Impurity;
   
    std::vector< double > T_mu_tr_vect_1;
    std::vector< double > T_pi_tr_vect_1;
    std::vector< double > cos_mu_tr_vect_1;
    std::vector< double > cos_pi_tr_vect_1;

    std::vector< double > T_mu_vect_1;
    std::vector< double > T_pi_vect_1;
    std::vector< double > cos_mu_vect_1;
    std::vector< double > cos_pi_vect_1;
     
    std::vector< int > Impurity_tr_1;
    std::vector< int > Impurity_1;
    
    std::vector< double > T_mu_tr_vect_2;
    std::vector< double > T_pi_tr_vect_2;
    std::vector< double > cos_mu_tr_vect_2;
    std::vector< double > cos_pi_tr_vect_2;

    std::vector< double > T_mu_vect_2;
    std::vector< double > T_pi_vect_2;
    std::vector< double > cos_mu_vect_2;
    std::vector< double > cos_pi_vect_2;
     
    std::vector< int > Impurity_tr_2;
    std::vector< int > Impurity_2;

    std::vector< double > T_mu_tr_vect_3;
    std::vector< double > T_pi_tr_vect_3;
    std::vector< double > cos_mu_tr_vect_3;
    std::vector< double > cos_pi_tr_vect_3;

    std::vector< double > T_mu_vect_3;
    std::vector< double > T_pi_vect_3;
    std::vector< double > cos_mu_vect_3;
    std::vector< double > cos_pi_vect_3;
     
    std::vector< int > Impurity_tr_3;
    std::vector< int > Impurity_3;

    int n_entries = gst->GetEntries();
    
    //==============================================================================
    // Fill the vectors with only signal for chosen selecetions
    //==============================================================================

    // Use n_protons == -1 to account for them all 
    // Use the same sample to train all histograms
    FillVect( gst_train, -1, Impurity_tr, T_mu_tr_vect, cos_mu_tr_vect, T_pi_tr_vect, cos_pi_tr_vect );
   
    // Use the same sample to train 1p histograms
    FillVect( gst_train, 1, Impurity_tr_1, T_mu_tr_vect_1, cos_mu_tr_vect_1, T_pi_tr_vect_1, cos_pi_tr_vect_1 );

    // Use the same sample to train 2p histograms
    FillVect( gst_train, 2, Impurity_tr_2, T_mu_tr_vect_2, cos_mu_tr_vect_2, T_pi_tr_vect_2, cos_pi_tr_vect_2 );

    // Use the same sample to train 3p histograms
    FillVect( gst_train, 3, Impurity_tr_3, T_mu_tr_vect_3, cos_mu_tr_vect_3, T_pi_tr_vect_3, cos_pi_tr_vect_3 );

    // Fill individual selection criteria
    // All number of protons as signal
    FillVect( gst, -1, Impurity, T_mu_vect, cos_mu_vect, T_pi_vect, cos_pi_vect );
    
    // 1 proton signal
    FillVect( gst, 1, Impurity_1, T_mu_vect_1, cos_mu_vect_1, T_pi_vect_1, cos_pi_vect_1 );
    
    // 2 protons signal
    FillVect( gst, 2, Impurity_2, T_mu_vect_2, cos_mu_vect_2, T_pi_vect_2, cos_pi_vect_2 );
    
    // 3 protons signal
    FillVect( gst, 3, Impurity_3, T_mu_vect_3, cos_mu_vect_3, T_pi_vect_3, cos_pi_vect_3 );

    //==============================================================================
    // Smearing
    //==============================================================================

    std::vector< double > T_mu_prime;
    std::vector< double > cos_mu_prime;
    
    std::vector< double > T_mu_prime_1;
    std::vector< double > cos_mu_prime_1;
    
    std::vector< double > T_mu_prime_2;
    std::vector< double > cos_mu_prime_2;
    
    std::vector< double > T_mu_prime_3;
    std::vector< double > cos_mu_prime_3;
    
    std::vector< double > T_pi_prime;
    std::vector< double > cos_pi_prime;
    
    std::vector< double > T_pi_prime_1;
    std::vector< double > cos_pi_prime_1;
    
    std::vector< double > T_pi_prime_2;
    std::vector< double > cos_pi_prime_2;
    
    std::vector< double > T_pi_prime_3;
    std::vector< double > cos_pi_prime_3;
    
    // Fill the smeared T and cos vectors 
    Smearing( T_mu_vect, cos_mu_vect, T_mu_prime, cos_mu_prime );  
    Smearing( T_mu_vect_1, cos_mu_vect_1, T_mu_prime_1, cos_mu_prime_1 );
    Smearing( T_mu_vect_2, cos_mu_vect_2, T_mu_prime_2, cos_mu_prime_2 );
    Smearing( T_mu_vect_3, cos_mu_vect_3, T_mu_prime_3, cos_mu_prime_3 );
    
    Smearing( T_pi_vect, cos_pi_vect, T_pi_prime, cos_pi_prime );  
    Smearing( T_pi_vect_1, cos_pi_vect_1, T_pi_prime_1, cos_pi_prime_1 );
    Smearing( T_pi_vect_2, cos_pi_vect_2, T_pi_prime_2, cos_pi_prime_2 );
    Smearing( T_pi_vect_3, cos_pi_vect_3, T_pi_prime_3, cos_pi_prime_3 );
    
    //==============================================================================
    // Define the histograms
    //==============================================================================
    
    TCanvas *c = new TCanvas( "c", "Unfolding", 800, 600 );
   
    gStyle->SetPalette(55);

    // Fill the smeared, training and to-be-unfolded truth histograms
    TH2D *h_sm_all = new TH2D( "h_sm_all", " CC0#pi, smeared ", 20, -1, 1, 18, 0, 2 );    
    TH2D *h_sm_1p  = new TH2D( "h_sm_1p", " CC0#pi 1p, smeared ", 20, -1, 1, 18, 0, 2 );    
    TH2D *h_sm_2p  = new TH2D( "h_sm_2p", " CC0#pi 2p, smeared ", 20, -1, 1, 18, 0, 2 );    
    TH2D *h_sm_3p  = new TH2D( "h_sm_3p", " CC0#pi 3p, smeared ", 20, -1, 1, 18, 0, 2 );    
  
    // Training
    TH2D *h_un_tr  = new TH2D( "h_un_tr", " CC0#pi, truth ", 20, -1, 1, 18, 0, 2 );     
    gst_train->Draw("( El - 0.10566 ):cthl>>h_un_tr","fspl == 13 && cc && (nfpip + nfpim + nfpi0  == 0)","colz");
 
    TH2D *h_un_tr_1 = new TH2D( "h_un_tr_1", " CC0#pi 1p, truth ", 20, -1, 1, 18, 0, 2 );     
    gst_train->Draw("( El - 0.10566 ):cthl>>h_un_tr_1"," nfp == 1 && fspl == 13 && cc && (nfpip + nfpim + nfpi0  == 0)","colz");
 
    TH2D *h_un_tr_2 = new TH2D( "h_un_tr_2", " CC0#pi 2p, truth ", 20, -1, 1, 18, 0, 2 );     
    gst_train->Draw("( El - 0.10566 ):cthl>>h_un_tr_2"," nfp == 2 && fspl == 13 && cc && (nfpip + nfpim + nfpi0  == 0)","colz");
 
    TH2D *h_un_tr_3 = new TH2D( "h_un_tr_3", " CC0#pi 3p, truth ", 20, -1, 1, 18, 0, 2 );     
    gst_train->Draw("( El - 0.10566 ):cthl>>h_un_tr_1"," nfp == 3 && fspl == 13 && cc && (nfpip + nfpim + nfpi0  == 0)","colz");
 
    // True histograms
    TH2D *h_un_all = new TH2D( "h_un_all", " CC0#pi, truth ", 20, -1, 1, 18, 0, 2 );     
    gst->Draw("( El - 0.10566 ):cthl>>h_un_all","fspl == 13 && cc && (nfpip + nfpim + nfpi0  == 0)","colz");
 
    h_un_all->SetStats(kFALSE);
    h_un_all->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_un_all->GetYaxis()->SetTitle("T_{#mu}");
    c->SetRightMargin(0.13);
    c->SaveAs("working_dir/unfolded_distributions/all_signal/true_kinematic_distribution.png");

    TH2D *h_un_1p  = new TH2D( "h_un_1p", " CC0#pi 1p, truth ", 20, -1, 1, 18, 0, 2 );     
    gst->Draw("( El - 0.10566 ):cthl>>h_un_1p","nfp == 1 && fspl == 13 && cc && (nfpip + nfpim + nfpi0  == 0)","colz");
 
    h_un_1p->SetStats(kFALSE);
    h_un_1p->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_un_1p->GetYaxis()->SetTitle("T_{#mu}");
    c->SetRightMargin(0.13);
    c->SaveAs("working_dir/unfolded_distributions/1_p/true_kinematic_distribution.png");

    TH2D *h_un_2p  = new TH2D( "h_un_2p", " CC0#pi 2p, truth ", 20, -1, 1, 18, 0, 2 );     
    gst->Draw("( El - 0.10566 ):cthl>>h_un_2p","nfp == 2 && fspl == 13 && cc && (nfpip + nfpim + nfpi0  == 0)","colz");
 
    h_un_2p->SetStats(kFALSE);
    h_un_2p->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_un_2p->GetYaxis()->SetTitle("T_{#mu}");
    c->SetRightMargin(0.13);
    c->SaveAs("working_dir/unfolded_distributions/2_p/true_kinematic_distribution.png");

    TH2D *h_un_3p  = new TH2D( "h_un_3p", " CC0#pi 3p, truth ", 20, -1, 1, 18, 0, 2 );     
    gst->Draw("( El - 0.10566 ):cthl>>h_un_3p","nfp == 3 && fspl == 13 && cc && (nfpip + nfpim + nfpi0  == 0)","colz");
 
    h_un_3p->SetStats(kFALSE);
    h_un_3p->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_un_3p->GetYaxis()->SetTitle("T_{#mu}");
    c->SetRightMargin(0.13);
    c->SaveAs("working_dir/unfolded_distributions/3_p/true_kinematic_distribution.png");

    //==============================================================================
    // Train the unfolding algorithm and get the response matrix
    //==============================================================================

    cout << "==================================== TRAIN ====================================" << endl;
    
    // Define the binning for the response matrix
    RooUnfoldResponse response ( h_un_tr, h_un_tr );
    RooUnfoldResponse response_1 ( h_un_tr_1, h_un_tr_1 );
    RooUnfoldResponse response_2 ( h_un_tr_2, h_un_tr_2 );
    RooUnfoldResponse response_3 ( h_un_tr_3, h_un_tr_3 );
    
    for ( int i = 0; i < n_entries; ++i ){
 
        gst->GetEntry(i);
 
        int cc    = b_cc->GetLeaf( "cc" )->GetValue();
        int nc    = b_nc->GetLeaf( "nc" )->GetValue();
        int fspl  = b_fspl->GetLeaf( "fspl" )->GetValue();
        int nfpip = b_nfpip->GetLeaf( "nfpip" )->GetValue();
        int nfpim = b_nfpim->GetLeaf( "nfpim" )->GetValue();
        int nfpi0 = b_nfpi0->GetLeaf( "nfpi0" )->GetValue();
        int nfp   = b_nfp->GetLeaf( "nfp" )->GetValue();
 
        // Fill the response matrix for all signal
        if ( fspl == 13
          && cc
          && Impurity[i] == 0
          && ( nfpip + nfpim + nfpi0 == 0 )
          && T_mu_prime[i] > 0.05 ){
 
            response.Fill( cos_mu_prime[i], T_mu_prime[i], cos_mu_vect[i], T_mu_vect[i]  );       
        }
        else if ( Impurity[i] == 1
               && T_pi_prime[i] > 0.05 ){
 
            response.Fill( cos_pi_prime[i], T_pi_prime[i], cos_pi_vect[i], T_pi_vect[i] );
        }

        // Fill the response matrix for 1p
        if ( fspl == 13
          && cc
          && Impurity_1[i] == 0
          && ( nfpip + nfpim + nfpi0 == 0 )
          && T_mu_prime_1[i] > 0.05 ){
 
            response_1.Fill( cos_mu_prime_1[i], T_mu_prime_1[i], cos_mu_vect_1[i], T_mu_vect_1[i]  );       
        }
        else if ( Impurity_1[i] == 1
                  && T_pi_prime_1[i] > 0.05 ){
 
            response_1.Fill( cos_pi_prime_1[i], T_pi_prime_1[i], cos_pi_vect_1[i], T_pi_vect_1[i] );
        
        }
        
        // Fill the response matrix for 2p
        if ( fspl == 13
          && cc
          && Impurity_2[i] == 0
          && ( nfpip + nfpim + nfpi0 == 0 )
          && T_mu_prime_2[i] > 0.05 ){
 
            response_2.Fill( cos_mu_prime_2[i], T_mu_prime_2[i], cos_mu_vect_2[i], T_mu_vect_2[i]  );       
        }
        else if ( Impurity_2[i] == 1
                  && T_pi_prime_2[i] > 0.05 ){
 
            response_2.Fill( cos_pi_prime_2[i], T_pi_prime_2[i], cos_pi_vect_2[i], T_pi_vect_2[i] );
        
        }
        
        // Fill the response matrix for 3p
        if ( fspl == 13
          && cc
          && Impurity_3[i] == 0
          && ( nfpip + nfpim + nfpi0 == 0 )
          && T_mu_prime_3[i] > 0.05 ){
 
            response_3.Fill( cos_mu_prime_3[i], T_mu_prime_3[i], cos_mu_vect_3[i], T_mu_vect_3[i]  );       
        }
        else if ( Impurity_3[i] == 1
                  && T_pi_prime_3[i] > 0.05 ){
 
            response_3.Fill( cos_pi_prime_3[i], T_pi_prime_3[i], cos_pi_vect_3[i], T_pi_vect_3[i] );
        
        }
    }
    
    //==============================================================================
    // Fill the true histogram for comparison with the unfolded histogram 
    //==============================================================================

    for ( int i = 0; i < n_entries; ++i ){
 
        gst->GetEntry(i);
 
        int cc    = b_cc->GetLeaf( "cc" )->GetValue();
        int nc    = b_nc->GetLeaf( "nc" )->GetValue();
        int fspl  = b_fspl->GetLeaf( "fspl" )->GetValue();
        int nfpip = b_nfpip->GetLeaf( "nfpip" )->GetValue();
        int nfpim = b_nfpim->GetLeaf( "nfpim" )->GetValue();
        int nfpi0 = b_nfpi0->GetLeaf( "nfpi0" )->GetValue();
        int nfp   = b_nfp->GetLeaf( "nfp" )->GetValue();
 
        if ( fspl == 13
          && cc
          && Impurity[i] == 0
          && ( nfpip + nfpim + nfpi0 == 0 )
          && T_mu_prime[i] > 0.05 ){
 
            h_sm_all->Fill(cos_mu_prime[i], T_mu_prime[i]);
        }
        else if ( Impurity[i] == 1
               && T_pi_prime[i] > 0.05 ){
 
            h_sm_all->Fill(cos_pi_prime[i], T_pi_prime[i]);
        }
        
        if ( fspl == 13
          && cc
          && Impurity_1[i] == 0
          && ( nfpip + nfpim + nfpi0 == 0 )
          && T_mu_prime_1[i] > 0.05 ){
 
            h_sm_1p->Fill(cos_mu_prime_1[i], T_mu_prime_1[i]);
        }
        else if ( Impurity_1[i] == 1
                  && T_pi_prime_1[i] > 0.05 ){
 
            h_sm_1p->Fill(cos_pi_prime_1[i], T_pi_prime_1[i]);
        }
        
        if ( fspl == 13
          && cc
          && Impurity_2[i] == 0
          && ( nfpip + nfpim + nfpi0 == 0 )
          && T_mu_prime_2[i] > 0.05 ){
 
            h_sm_2p->Fill( cos_mu_prime_2[i], T_mu_prime_2[i] );       
        }
        else if ( Impurity_2[i] == 1
                  && T_pi_prime_2[i] > 0.05 ){
 
            h_sm_2p->Fill( cos_pi_prime_2[i], T_pi_prime_2[i] );
        
        }
        
        if ( fspl == 13
          && cc
          && Impurity_3[i] == 0
          && ( nfpip + nfpim + nfpi0 == 0 )
          && T_mu_prime_3[i] > 0.05 ){
 
            h_sm_3p->Fill( cos_mu_prime_3[i], T_mu_prime_3[i] );       
        }
        else if ( Impurity_3[i] == 1
                  && T_pi_prime_3[i] > 0.05 ){
 
            h_sm_3p->Fill( cos_pi_prime_3[i], T_pi_prime_3[i] );
        
        }

    }
    
    //==============================================================================
    // Unfold the first half of the smeared statistics
    //==============================================================================

    cout << "==================================== UNFOLD ===================================" << endl;

    RooUnfoldBayes    unfold_all  ( &response, h_sm_all, 1 );
    RooUnfoldBayes    unfold_1p   ( &response_1, h_sm_1p, 1 );
    RooUnfoldBayes    unfold_2p   ( &response_2, h_sm_2p, 1 );
    RooUnfoldBayes    unfold_3p   ( &response_3, h_sm_3p, 1 );

    // Unfold
    // Histogram output
    TH2D *hunfold_all =  (TH2D*) unfold_all.Hreco();
    TH2D *hunfold_1p  =  (TH2D*) unfold_1p.Hreco();
    TH2D *hunfold_2p  =  (TH2D*) unfold_2p.Hreco();
    TH2D *hunfold_3p  =  (TH2D*) unfold_3p.Hreco();
    
    // Vector and covariance matrix output
    //TVectorD unfolded_dist = unfold.Vreco();
    //TMatrixD unfolded_errs = unfold.Ereco;

    unfold_all.PrintTable (cout, hunfold_all);
   
    hunfold_all->SetStats(kFALSE);
    hunfold_all->GetXaxis()->SetTitle("cos#theta_{#mu}");
    hunfold_all->GetYaxis()->SetTitle("T_{#mu}");
    hunfold_all->SetTitle("Unfolded #mu kinematics");
    hunfold_all->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/all_signal/2D_unfolding.png" );

    hunfold_1p->SetStats(kFALSE);
    hunfold_1p->GetXaxis()->SetTitle("cos#theta_{#mu}");
    hunfold_1p->GetYaxis()->SetTitle("T_{#mu}");
    hunfold_1p->SetTitle("Unfolded #mu kinematics");
    hunfold_1p->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/1_p/2D_unfolding.png" );
    
    hunfold_2p->SetStats(kFALSE);
    hunfold_2p->GetXaxis()->SetTitle("cos#theta_{#mu}");
    hunfold_2p->GetYaxis()->SetTitle("T_{#mu}");
    hunfold_2p->SetTitle("Unfolded #mu kinematics");
    hunfold_2p->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/2_p/2D_unfolding.png" );
    
    hunfold_3p->SetStats(kFALSE);
    hunfold_3p->GetXaxis()->SetTitle("cos#theta_{#mu}");
    hunfold_3p->GetYaxis()->SetTitle("T_{#mu}");
    hunfold_3p->SetTitle("Unfolded #mu kinematics");
    hunfold_3p->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/3_p/2D_unfolding.png" );
    
    //==============================================================================
    // Fill the comparison histogram and draw
    //==============================================================================
   
    TH2D *h_comp = new TH2D( *h_un_all );
    h_comp->Add( hunfold_all, -1 );
    h_comp->Divide( h_un_all );
    h_comp->Scale( 100 );

    h_comp->SetStats(kFALSE);
    h_comp->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_comp->GetYaxis()->SetTitle("T_{#mu}");
    h_comp->SetTitle("CC0#pi, bin content difference between true and unfolded");
    h_comp->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/all_signal/true_unfolding_comp.png" );

    TH2D *h_comp_1 = new TH2D( *h_un_1p );
    h_comp_1->Add( hunfold_1p, -1 );
    h_comp_1->Divide( h_un_1p );
    h_comp_1->Scale( 100 );

    h_comp_1->SetStats(kFALSE);
    h_comp_1->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_comp_1->GetYaxis()->SetTitle("T_{#mu}");
    h_comp_1->SetTitle("CC0#pi 1p, bin content difference between true and unfolded");
    h_comp_1->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/1_p/true_unfolding_comp.png" );

    TH2D *h_comp_2 = new TH2D( *h_un_2p );
    h_comp_2->Add( hunfold_2p, -1 );
    h_comp_2->Divide( h_un_2p );
    h_comp_2->Scale( 100 );

    h_comp_2->SetStats(kFALSE);
    h_comp_2->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_comp_2->GetYaxis()->SetTitle("T_{#mu}");
    h_comp_2->SetTitle("CC0#pi 2p, bin content difference between true and unfolded");
    h_comp_2->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/2_p/true_unfolding_comp.png" );

    TH2D *h_comp_3 = new TH2D( *h_un_3p );
    h_comp_3->Add( hunfold_3p, -1 );
    h_comp_3->Divide( h_un_3p );
    h_comp_3->Scale( 100 );

    h_comp_3->SetStats(kFALSE);
    h_comp_3->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_comp_3->GetYaxis()->SetTitle("T_{#mu}");
    h_comp_3->SetTitle("CC0#pi 3p, bin content difference between true and unfolded");
    h_comp_3->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/3_p/true_unfolding_comp.png" );

    //==============================================================================
    // Fill the ddxsec histogram and draw
    //==============================================================================
  
    // Normalization variables
    double cos_bins = hunfold_all->GetXaxis()->GetBinWidth(1);
    double Tmu_bins = hunfold_all->GetYaxis()->GetBinWidth(1);
    double flux_int = h_flux->Integral();
    double Na       = 6.63e34;
    double M_fid    = 1.016e8; // grams
    double A_ar     = 39.948;
    double tot_tgt  = ( Na * M_fid ) / ( A_ar ); // given in flux unit definition  
    double barns    = 1e38;

    double scalar_norm = ( barns ) / ( cos_bins * Tmu_bins * flux_int * tot_tgt );
    
    TH2D *h_ddxsec = new TH2D( *hunfold_all );
    h_ddxsec->Scale( scalar_norm );
    h_ddxsec->Divide( h_eff );

    h_ddxsec->SetStats(kFALSE);
    h_ddxsec->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_ddxsec->GetYaxis()->SetTitle("T_{#mu}");
    h_ddxsec->SetTitle("CC0#pi, d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]");
    h_ddxsec->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/all_signal/2D_unfolding_ddxsec.png" );

    TH2D *h_ddxsec_1 = new TH2D( *hunfold_1p );
    h_ddxsec_1->Scale( scalar_norm );
    h_ddxsec_1->Divide( h_eff );

    h_ddxsec_1->SetStats(kFALSE);
    h_ddxsec_1->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_ddxsec_1->GetYaxis()->SetTitle("T_{#mu}");
    h_ddxsec_1->SetTitle("CC0#pi 1p, d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]");
    h_ddxsec->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/1_p/2D_unfolding_ddxsec.png" );

    TH2D *h_ddxsec_2 = new TH2D( *hunfold_2p );
    h_ddxsec_2->Scale( scalar_norm );
    h_ddxsec_2->Divide( h_eff );

    h_ddxsec_2->SetStats(kFALSE);
    h_ddxsec_2->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_ddxsec_2->GetYaxis()->SetTitle("T_{#mu}");
    h_ddxsec_2->SetTitle("CC0#pi 2p, d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]");
    h_ddxsec->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/2_p/2D_unfolding_ddxsec.png" );

    TH2D *h_ddxsec_3 = new TH2D( *hunfold_3p );
    h_ddxsec_3->Scale( scalar_norm );
    h_ddxsec_3->Divide( h_eff );

    h_ddxsec_3->SetStats(kFALSE);
    h_ddxsec_3->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_ddxsec_3->GetYaxis()->SetTitle("T_{#mu}");
    h_ddxsec_3->SetTitle("CC0#pi 3p, d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]");
    h_ddxsec->Draw("colz");
    
    c->SetRightMargin(0.13);
    c->SaveAs( "working_dir/unfolded_distributions/3_p/2D_unfolding_ddxsec.png" );

    //==============================================================================
    // Slices and comparisons
    //==============================================================================

    // Normalise the true and reco histogram
    h_un_all->Scale( scalar_norm );
    h_un_all->Divide( h_eff );

    h_sm_all->Scale( scalar_norm );
    h_sm_all->Divide( h_eff );

    h_un_1p->Scale( scalar_norm );
    h_un_1p->Divide( h_eff );

    h_sm_1p->Scale( scalar_norm );
    h_sm_1p->Divide( h_eff );

    h_un_2p->Scale( scalar_norm );
    h_un_2p->Divide( h_eff );

    h_sm_2p->Scale( scalar_norm );
    h_sm_2p->Divide( h_eff );

    h_un_3p->Scale( scalar_norm );
    h_un_3p->Divide( h_eff );

    h_sm_3p->Scale( scalar_norm );
    h_sm_3p->Divide( h_eff );

    Slices( h_ddxsec, h_un_all, h_sm_all, "all_signal" );
    Slices( h_ddxsec_1, h_un_1p, h_sm_1p, "1_p" );
    Slices( h_ddxsec_2, h_un_2p, h_sm_2p, "2_p" );
    Slices( h_ddxsec_3, h_un_3p, h_sm_3p, "3_p" );
    
    //==============================================================================
    // Delete pointers
    //==============================================================================
    
    delete hunfold_all;
    delete hunfold_1p;
    delete hunfold_2p;
    delete hunfold_3p;

    delete h_comp;
    delete h_comp_1;
    delete h_comp_2;
    delete h_comp_3;

    delete h_ddxsec;
    delete h_ddxsec_1;
    delete h_ddxsec_2;
    delete h_ddxsec_3;

    delete h_sm_all;
    delete h_sm_1p;
    delete h_sm_2p;
    delete h_sm_3p;

    delete h_un_tr;
    delete h_un_tr_1;
    delete h_un_tr_2;
    delete h_un_tr_3;

    delete h_un_all;
    delete h_un_1p;
    delete h_un_2p;
    delete h_un_3p;
    
    delete c;
    
} 
