#include <iostream>
using std::cout;
using std::endl;

#include "roo_unfold.h"

#include "../src/RooUnfoldResponse.h"
#include "../src/RooUnfoldBayes.h"


void GetTruth( TTree *tree, int n_protons, std::vector<double> & truth_T, std::vector<double> & truth_cos, std::vector<bool> & truth_detectable, std::vector<double> & impur_T, std::vector<double> & impur_cos );

void Smear( const std::vector<double> & truth_T, const std::vector<double> & truth_cos, std::vector<double> & smear_T, std::vector<double> & smear_cos ); 



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

    GetTruth( gst_train, -1, truth_T_train, truth_cos_train, truth_detectable_train, impur_T_train, impur_cos_train );

    std::vector<double> truth_T_test; 
    std::vector<double> truth_cos_test; 
    std::vector<bool>   truth_detectable_test;
    std::vector<double> impur_T_test; 
    std::vector<double> impur_cos_test; 

    GetTruth( gst_test, -1,  truth_T_test, truth_cos_test, truth_detectable_test, impur_T_test, impur_cos_test );
    
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
            b_cthf->GetEntry(i);
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
                            impur_T.push_back(T_pi);
                            impur_cos.push_back(cos_pi);
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
