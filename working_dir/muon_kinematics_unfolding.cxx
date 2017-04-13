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
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"

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
 
        double var_mu     = TMath::Power( T_mu[i] * 0.1, 2 );
        double sigma_mu   = TMath::Sqrt( TMath::Log( 1 + ( var_mu / TMath::Power( T_mu[i], 2 ) ) ) );
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
 
        if ( T_mu[i] > 0.05 ){
            T_mu_sm.push_back(lognorm_mu);
            cos_mu_sm.push_back(gaus_costheta);
 
        }

    }
}


//==============================================================================
// The main function
//==============================================================================
void muon_kinematics_unfolding() { 
    //==============================================================================
    // Reading in the root file 
    //==============================================================================
    TFile f("/hepstore/rjones/Exercises/Flavours/Default+MEC/sbnd/1M/gntp.10000.gst.root");
    if(f.IsZombie()){
        std::cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "============================== Default + MEC open =============================" << endl;
    }

    //==============================================================================
    // Get everything from the tree
    //==============================================================================

    TTree *gst = (TTree*) f.Get("gst");

    TBranch *b_Ef    = gst->GetBranch("Ef");
    TBranch *b_Ev    = gst->GetBranch("Ev");
    TBranch *b_El    = gst->GetBranch("El");
    TBranch *b_pl    = gst->GetBranch("pl");
    TBranch *b_cthl  = gst->GetBranch("cthl");
    TBranch *b_cthf  = gst->GetBranch("cthf");
    TBranch *b_nfpi0 = gst->GetBranch("nfpi0");
    TBranch *b_nfpip = gst->GetBranch("nfpip");
    TBranch *b_nfpim = gst->GetBranch("nfpim");
    TBranch *b_cc    = gst->GetBranch("cc");
    TBranch *b_nc    = gst->GetBranch("nc");
    TBranch *b_nf    = gst->GetBranch("nf");
    TBranch *b_pdgf  = gst->GetBranch("pdgf");
    TBranch *b_fspl  = gst->GetBranch("fspl");

    // Variables 
    double m_mu = 0.10566; // Muon mass, GeV
    double m_pi = 0.13957; // Charged pion mass, GeV

    // Define and fill vectors for the muon energy, momentum and opening angle
    // Define an empty vector to hold the reconstructed energy
    std::vector< double > T_mu_vect;
    std::vector< double > T_pi_vect;
    std::vector< double > cos_mu_vect;
    std::vector< double > cos_pi_vect;

    std::vector< int > Impurity;
    
    int n_entries = gst->GetEntries();

    // Fill the vectors
    for ( int i = 0; i < n_entries; ++i ){
    
        gst->GetEntry(i);
    
        int nf    = b_nf->GetLeaf("nf")->GetValue();
        int fspl  = b_fspl->GetLeaf("fspl")->GetValue();
        int nfpi0 = b_nfpi0->GetLeaf("nfpi0")->GetValue();
        int nfpip = b_nfpip->GetLeaf("nfpip")->GetValue();
        int nfpim = b_nfpim->GetLeaf("nfpim")->GetValue();
 
        double e_mu   = b_El->GetLeaf("El")->GetValue();
        double cos_mu = b_cthl->GetLeaf("cthl")->GetValue();
 
        double T_mu;
 
        // Calculate the kinetic energy for muons
        if ( fspl == 13 ){
 
            // Energy of the final state primary lepton
            T_mu = e_mu - m_mu;
 
            T_mu_vect.push_back(T_mu);
            cos_mu_vect.push_back(cos_mu);
 
        }
        // If the final state primary is a lepton, push back a number that will
        // be removed in the cuts later
        else if ( fspl == 11 ){
            T_mu_vect.push_back(-99999);
            cos_mu_vect.push_back(-99999);
        }
        else{
            T_mu_vect.push_back(-99999);
            cos_mu_vect.push_back(-99999);
        }

        if ( nfpip + nfpim == 1 && nfpi0 == 0 ){
 
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
 
                T_pi_vect.push_back(T_pi);
                cos_pi_vect.push_back(cos_pi);
        
                }
            }
        }
        else{
            T_pi_vect.push_back(-99999);
            cos_pi_vect.push_back(-99999);
        }

        // Apply the impurity cut
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 && nfpip + nfpim == 1 && nfpi0 == 0 ){
            int random;
            random = rand() % 5 + 1;

            // If the random number = 5 then push back a 1 onto the vector
            // If the random number < 5 push back a 0
            if( random == 5 ){
                Impurity.push_back(1);
            }
            else{
                Impurity.push_back(0);
            }
        }
        else{
            Impurity.push_back(0);
        } 
    }


    //==============================================================================
    // Smearing
    //==============================================================================

    std::vector< double > T_mu_prime;
    std::vector< double > cos_mu_prime;
    
    // Fill the smeared T and cos vectors 
    Smearing( T_mu_vect, cos_mu_vect, T_mu_prime, cos_mu_prime );
    
    
    //==============================================================================
    // Define the histograms
    //==============================================================================
    
    TCanvas *c = new TCanvas( "c", "Unfolding", 800, 600 );
   
    gStyle->SetPalette(55);

    // Fill 2 histograms with the first and last half of the distribution for training 
    // And unfolding separately
    TH2D *h_un = new TH2D( "h_un", " CC0#pi, truth ", 20, -1, 1, 18, 0, 2 );     
    gst->Draw("( El - 0.10566 ):cthl>>h_un","fspl == 13 && cc && (nfpip + nfpim + nfpi0  == 0)","colz");
 
    TH2D *h_un_1 = new TH2D( "h_un_1", " CC0#pi, truth ", 20, -1, 1, 18, 0, 2 );     
    gst->Draw("( El - 0.10566 ):cthl>>h_un_1","iev <= 500000 && fspl == 13 && cc && (nfpip + nfpim + nfpi0  == 0)","colz");
 
    TH2D *h_un_2 = new TH2D( "h_un_2", " CC0#pi, truth ", 20, -1, 1, 18, 0, 2 );     
    gst->Draw("( El - 0.10566 ):cthl>>h_un_2","iev > 500000 && fspl == 13 && cc && (nfpip + nfpim + nfpi0  == 0)","colz");
 
    h_un_2->SetStats(kFALSE);
    h_un_2->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_un_2->GetYaxis()->SetTitle("T_{#mu}");
    c->SetRightMargin(0.13);
    //c->SetLogz();
    c->SaveAs("working_dir/unfolded_distributions/true_kinematic_distribution.png");

    TH2D *h_sm   = new TH2D( "h_sm", " CC0#pi, smeared ", 20, -1, 1, 18, 0, 2 );    
    TH2D *h_sm_1 = new TH2D( "h_sm_1", " CC0#pi, smeared ", 20, -1, 1, 18, 0, 2 );    
    TH2D *h_sm_2 = new TH2D( "h_sm_2", " CC0#pi, smeared ", 20, -1, 1, 18, 0, 2 );    
   
    //==============================================================================
    // Train the unfolding algorithm and get the response matrix
    //==============================================================================

    int n_half = TMath::Floor( double( n_entries ) / 2. );

    cout << "==================================== TRAIN ====================================" << endl;
    
    RooUnfoldResponse response ( h_sm, h_un );
    
    for ( int i = 0; i < n_half; ++i ){
 
        gst->GetEntry(i);
 
        int cc    = b_cc->GetLeaf( "cc" )->GetValue();
        int nc    = b_nc->GetLeaf( "nc" )->GetValue();
        int fspl  = b_fspl->GetLeaf( "fspl" )->GetValue();
        int nfpip = b_nfpip->GetLeaf( "nfpip" )->GetValue();
        int nfpim = b_nfpim->GetLeaf( "nfpim" )->GetValue();
        int nfpi0 = b_nfpi0->GetLeaf( "nfpi0" )->GetValue();
 
        // Set the energy and impurity cuts
        // Filling the signal ntuple
        if ( fspl == 13
          && cc
          && ( nfpip + nfpim + nfpi0 == 0 )
          && T_mu_vect[i] != -99999 ){
 
            // Filling the cuts histogram
            h_sm_1->Fill( cos_mu_prime[i], T_mu_prime[i] );
            response.Fill( cos_mu_prime[i], T_mu_prime[i], cos_mu_vect[i], T_mu_vect[i]  );       
        }
        else if ( Impurity[i] == 1
               && T_pi_vect[i] > 0.05 ){
 
            h_sm_1->Fill( cos_pi_vect[i], T_pi_vect[i] );
            response.Fill( cos_pi_vect[i], T_pi_vect[i], cos_mu_vect[i], T_mu_vect[i] );
        
        }

    }
    
    //==============================================================================
    // Fill the true histogram for comparison with the unfolded histogram 
    //==============================================================================

    for ( int i = n_half; i < n_entries; ++i ){
 
        gst->GetEntry(i);
 
        int cc    = b_cc->GetLeaf( "cc" )->GetValue();
        int nc    = b_nc->GetLeaf( "nc" )->GetValue();
        int fspl  = b_fspl->GetLeaf( "fspl" )->GetValue();
        int nfpip = b_nfpip->GetLeaf( "nfpip" )->GetValue();
        int nfpim = b_nfpim->GetLeaf( "nfpim" )->GetValue();
        int nfpi0 = b_nfpi0->GetLeaf( "nfpi0" )->GetValue();
 
        // Set the energy and impurity cuts
        if ( fspl == 13
          && cc
          && ( nfpip + nfpim + nfpi0 == 0 )
          && T_mu_vect[i] > 0.05 ){
 
            // Filling the cuts histogram
            h_sm_2->Fill(cos_mu_prime[i], T_mu_prime[i]);
        }
        else if ( Impurity[i] == 1
               && T_pi_vect[i] > 0.05 ){
 
            h_sm_2->Fill(cos_pi_vect[i], T_pi_vect[i]);
        }

    }
    
    //==============================================================================
    // Unfold the first half of the smeared statistics
    //==============================================================================

    cout << "==================================== UNFOLD ===================================" << endl;
    RooUnfoldBayes    unfold   ( &response, h_sm_2, 1 ); // Try different numbers of iterations

    // Unfold
    // Histogram output
    TH2D *hUnfold =  (TH2D*) unfold.Hreco();
    
    // Vector and covariance matrix output
    //TVectorD unfolded_dist = unfold.Vreco();
    //TMatrixD unfolded_errs = unfold.Ereco;

    unfold.PrintTable (cout, h_un_2);
   
    // double norm = hUnfold->Integral();

    hUnfold->SetStats(kFALSE);
    hUnfold->GetXaxis()->SetTitle("cos#theta_{#mu}");
    hUnfold->GetYaxis()->SetTitle("T_{#mu}");
    hUnfold->SetTitle("Unfolded #mu kinematics");
    hUnfold->Draw("colz");
    
    c->SetRightMargin(0.13);
    //c->SetLogz();
    c->SaveAs( "working_dir/unfolded_distributions/2D_unfolding.png" );

    //==============================================================================
    // Fill the comparison histogram and draw
    //==============================================================================
   
    TH2D *h_comp = new TH2D( *h_un_2 );
    h_comp->Add( hUnfold, -1 );
    //h_comp->Divide( h_un_2 );
    //h_comp->Scale( 100 );

    h_comp->SetStats(kFALSE);
    h_comp->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_comp->GetYaxis()->SetTitle("T_{#mu}");
    h_comp->SetTitle("CC0#pi, bin content difference between true and unfolded");
    h_comp->Draw("colz");
    
    c->SetRightMargin(0.13);
    //c->SetLogz();
    c->SaveAs( "working_dir/unfolded_distributions/true_unfolding_comp.png" );

    //==============================================================================
    // Delete pointers
    //==============================================================================
    delete hUnfold;

    delete h_comp;

    delete h_sm;
    delete h_sm_1;
    delete h_sm_2;

    delete h_un;
    delete h_un_1;
    delete h_un_2;
    
    delete c;
    
} 
