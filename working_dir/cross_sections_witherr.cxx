#include <iostream>
using std::cout;
using std::endl;

#include "roo_unfold.h"

#include "../src/RooUnfoldResponse.h"
#include "../src/RooUnfoldBayes.h"


void GetTruth( TTree *tree, 
               int n_min_bins, 
               int n_max_bins,  
               int n_protons, 
               std::vector<double> & truth_T, 
               std::vector<double> & truth_cos, 
               std::vector<bool>   & truth_detectable, 
               std::vector<double> & impur_T, 
               std::vector<double> & impur_cos );

bool isReconstructed(double tmu, double cos, TRandom *_rand);

bool isContained(double tmu, double cos, TRandom *_rand);

double RangeSmear(double ke, double m, ROOT::Math::GSLRngMT *_random_gen);

double McsSmear(double ke, ROOT::Math::GSLRngMT *_random_gen);

void Smear( const double &mass,
            const std::vector<double> & truth_T, 
            const std::vector<double> & truth_cos, 
            std::vector<double>       & smear_T, 
            std::vector<double>       & smear_cos );

void GetResponse( const std::vector<double> & truth_T, 
                  const std::vector<double> & truth_cos,
                  const std::vector<bool>   & truth_detectable, 
                  const std::vector<double> & smear_T, 
                  const std::vector<double> & smear_cos, 
                  const std::vector<double> & smear_impur_T, 
                  const std::vector<double> & smear_impur_cos, 
                  RooUnfoldResponse         & response); 

void Slices( TH2D *h_unfolded, TH2D *h_true, TH2D *h_reco, const char n_pr[1024] );

void cross_sections_witherr( const int &n_protons, const char pr_path[1024] ) {
    //==============================================================================
    // Reading in the event root files 
    //==============================================================================
    TFile f_test("/hepstore/rjones/Exercises/Flavours/G16_01b/sbnd/1M/gntp.10000.gst.root");
    if(f_test.IsZombie()){
        std::cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "=========================== G16_01b event file open ===========================" << endl;
    }
 
    TFile f_train("/hepstore/rjones/Exercises/Flavours/G16_01b_2/sbnd/1M/gntp.10000.gst.root");
    if(f_train.IsZombie()){
        std::cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "========================== G16_01b_2 event file open =========================" << endl;
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

    //==============================================================================
    // Get the SBND flux histogram
    //==============================================================================
    
    TH1D *h_flux = (TH1D*) fflux.Get("h_numu_110m");

    TTree *gst_train = (TTree*) f_train.Get("gst");
    TTree *gst_test  = (TTree*) f_test.Get("gst");

    //==============================================================================
    // Open file to record statistics
    //==============================================================================
    
    char stats_path[1024];
    
    stringstream temp_stats;
    temp_stats.clear();

    string stats;
    stats.clear();

    temp_stats << "working_dir/err_test/" << pr_path << "/stats.txt";
    
    stats = temp_stats.str();    

    strcpy( stats_path, stats.c_str() );

    ofstream s_file;
    s_file.open(stats_path); 

    if ( n_protons != -1 ){
        s_file << " File to record the statistics for signal with " << n_protons << " protons " << endl; 
        s_file << " ======================================================= " << endl; 
    }
    else{
        s_file << " File to record the statistics for signal with all protons " << endl; 
        s_file << " ========================================================= " << endl; 
    }

    //==============================================================================
    // Get the truth vectors 
    //==============================================================================

    int n_half  = gst_test->GetEntries() / 2;
    int n_full  = gst_test->GetEntries();
    int n_train = gst_train->GetEntries();

    std::vector<double> truth_T_train; 
    std::vector<double> truth_cos_train; 
    std::vector<bool>   truth_detectable_train;
    std::vector<double> impur_T_train; 
    std::vector<double> impur_cos_train; 

    GetTruth( gst_train, 0, n_train, n_protons, truth_T_train, truth_cos_train, truth_detectable_train, impur_T_train, impur_cos_train );

    std::vector<double> truth_T_test; 
    std::vector<double> truth_cos_test; 
    std::vector<bool>   truth_detectable_test;
    std::vector<double> impur_T_test; 
    std::vector<double> impur_cos_test; 

    GetTruth( gst_test, 0, n_full, n_protons, truth_T_test, truth_cos_test, truth_detectable_test, impur_T_test, impur_cos_test );
    
    //==============================================================================
    // Smearing 
    //==============================================================================
    
    std::vector<double> smear_T_train; 
    std::vector<double> smear_cos_train; 
    std::vector<double> impur_smear_T_train; 
    std::vector<double> impur_smear_cos_train; 
    double mass_mu = 0.105658;
    double mass_pi = 0.13957;

    Smear( mass_mu, truth_T_train, truth_cos_train, smear_T_train, smear_cos_train );
    Smear( mass_pi, impur_T_train, impur_cos_train, impur_smear_T_train, impur_smear_cos_train );
    
    std::vector<double> smear_T_test; 
    std::vector<double> smear_cos_test; 
    std::vector<double> impur_smear_T_test; 
    std::vector<double> impur_smear_cos_test; 
    
    Smear( mass_mu, truth_T_test, truth_cos_test, smear_T_test, smear_cos_test );
    Smear( mass_pi, impur_T_test, impur_cos_test, impur_smear_T_test, impur_smear_cos_test );
    
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
   
    // Scaling factor for all SBND events
    double scale = 7.56016;
    //double scale = 0.2;


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

    cout << " Reco before : " << h_reco_test->Integral() << endl;
    cout << " True before : " << h_true_test->Integral() << endl << endl;
    
    h_reco_test->Scale(scale);
    h_true_test->Scale(scale);
    h_cut_test->Scale(scale);

    cout << " Reco after  : " << h_reco_test->Integral() << endl;
    cout << " True after  : " << h_true_test->Integral() << endl;

    RooUnfoldBayes unfold( &response, h_reco_test, 1 );
    TH2D *h_unfold_test =  (TH2D*) unfold.Hreco();

    gStyle->SetPalette(55);
    gStyle->SetNumberContours(250);

    // TRUE
    h_true_test->SetStats(kFALSE);
    h_true_test->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_true_test->GetYaxis()->SetTitle("T_{#mu}");
    h_true_test->SetTitle("True #mu kinematics");
    h_true_test->Draw("colz");
   
    // Get the integral to see statistics 
    double true_int = h_true_test->Integral();

    s_file << " True : " << true_int << endl;

    char true_path[1024];
    
    stringstream temp_tr;
    temp_tr.clear();

    string image_tr;
    image_tr.clear();

    temp_tr << "working_dir/err_test/" << pr_path << "/" << pr_path << "_2D_true.png";
    
    image_tr = temp_tr.str();    

    strcpy( true_path, image_tr.c_str() );

    c->SetRightMargin(0.13);
    c->SaveAs( true_path );

    // CUT
    h_cut_test->SetStats(kFALSE);
    h_cut_test->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_cut_test->GetYaxis()->SetTitle("T_{#mu}");
    h_cut_test->SetTitle("#mu kinematics with cuts and smearing");
    h_cut_test->Draw("colz");
    
    // Get the integral to see statistics 
    double cut_int = h_cut_test->Integral();

    s_file << " Cut  : " << cut_int << endl;

    char cut_path[1024];
    
    stringstream temp_ct;
    temp_ct.clear();

    string image_ct;
    image_ct.clear();

    temp_ct << "working_dir/err_test/" << pr_path << "/" << pr_path << "_2D_cut.png";
    
    image_ct = temp_ct.str();    

    strcpy( cut_path, image_ct.c_str() );

    c->SetRightMargin(0.13);
    c->SaveAs( cut_path );
    
    // RECO
    h_reco_test->SetStats(kFALSE);
    h_reco_test->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_reco_test->GetYaxis()->SetTitle("T_{#mu}");
    h_reco_test->SetTitle("Data-like kinematics");
    h_reco_test->Draw("colz");
    
    // Get the integral to see statistics 
    double reco_int = h_reco_test->Integral();

    s_file << " Reco : " << reco_int << endl;

    char reco_path[1024];
    
    stringstream temp_re;
    temp_re.clear();

    string image_re;
    image_re.clear();

    temp_re << "working_dir/err_test/" << pr_path << "/" << pr_path << "_2D_reco.png";
    
    image_re = temp_re.str();    

    strcpy( reco_path, image_re.c_str() );

    c->SetRightMargin(0.13);
    c->SaveAs( reco_path );
    
    // UNFOLDED
    h_unfold_test->SetStats(kFALSE);
    h_unfold_test->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_unfold_test->GetYaxis()->SetTitle("T_{#mu}");
    h_unfold_test->SetTitle("Unfolded #mu kinematics");
    h_unfold_test->Draw("colz");
    
    // WRITE TO FILE
    TFile file("working_dir/selection_ddxsec.root", "UPDATE"); 
    h_unfold_test->Write(pr_path);
    file.Close();

    // Get the integral to see statistics 
    double unfold_int = h_unfold_test->Integral();

    s_file << " Unf  : " << unfold_int << endl;

    char unfold_path[1024];
    
    stringstream temp_unf;
    temp_unf.clear();

    string image_unf;
    image_unf.clear();

    temp_unf << "working_dir/err_test/" << pr_path << "/" << pr_path << "_2D_unfolded.png";
    
    image_unf = temp_unf.str();    

    strcpy( unfold_path, image_unf.c_str() );

    c->SetRightMargin(0.13);
    c->SaveAs( unfold_path );
    
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
    
    char comp_path[1024];
    
    stringstream temp_cmp;
    temp_cmp.clear();

    string image_cmp;
    image_cmp.clear();

    temp_cmp << "working_dir/err_test/" << pr_path << "/" << pr_path << "_true_unfold_comp.png";
    
    image_cmp = temp_cmp.str();    

    strcpy( comp_path, image_cmp.c_str() );

    c->SetRightMargin(0.13);
    c->SaveAs( comp_path );

    // EFFICIENCY
    TH2D *h_eff = new TH2D( *h_cut_test );
    h_eff->Divide( h_true_test );

    h_eff->SetStats(kFALSE);
    h_eff->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_eff->GetYaxis()->SetTitle("T_{#mu}");
    h_eff->SetTitle("Efficiency");
    h_eff->Draw("colz");
    
    char eff_path[1024];
    
    stringstream temp_eff;
    temp_eff.clear();

    string image_eff;
    image_eff.clear();

    temp_eff << "working_dir/err_test/" << pr_path <<  "/" << pr_path <<"_efficiency.png";
    
    image_eff = temp_eff.str();    

    strcpy( eff_path, image_eff.c_str() );

    c->SetRightMargin(0.13);
    c->SaveAs( eff_path );
    
    // CROSS-SECTIONS
  
    // Normalization variables
    double cos_bins = h_unfold_test->GetXaxis()->GetBinWidth(1);
    double Tmu_bins = h_unfold_test->GetYaxis()->GetBinWidth(1); // GeV
    // Units of flux: #/m^2/50MeV/10^6 POT
    double flux_int = h_flux->Integral();
    double Na       = 6.022e23;  //
    double M_fid    = 112000; // (kg / m^3) * m^3
    double A_ar     = 0.039948;  // kg / mol
    double tot_tgt  = ( Na * M_fid ) / ( A_ar ); // number of target particles  
    double POT      = 6.6e20;    // POT
    double scaling  = 5e8;      // Sort units from factor of 1 / ( integrated flux * Nt * bin widths * POT )
    double cm_conv  = 1e-38;     // take out factor of 10e-38 cm^2 
    double units_scale = 10e48; // Scale for units 

    // Playing with flux systematic
    double err_flux = 0.15; // 15%

    double scalar_norm = scaling / ( cm_conv * cos_bins * Tmu_bins * flux_int * tot_tgt * POT );

    TH2D *h_ddxsec = new TH2D( *h_unfold_test );
    h_ddxsec->Scale( scalar_norm );

    int n_x = h_ddxsec->GetNbinsX();
    int n_y = h_ddxsec->GetNbinsY();

    ofstream cov_file;
    cov_file.open("working_dir/cov_matrix.txt"); 

    cov_file << " Covariance matrix of flux systematics " << endl;
    cov_file << "=======================================" << endl << endl;

    // Linearise the 2D histogram

    TCanvas *c_lin = new TCanvas();

    TH1D *h_ddxsec_linear = new TH1D("h_ddxsec_linear", "", 360, 0, 360);

    for ( int i = 0; i < n_y; ++i ){
        
        for ( int j = 0; j < n_x; ++j ){
            
            int k = Map(j,i,n_x,n_y) + 1;
            
            h_ddxsec_linear->Fill( k , h_ddxsec->GetBinContent(j,i));        

        }
    }

    int n_lin_bins = h_ddxsec_linear->GetNbinsX();
        
    std::vector< std::vector< double > > flux_cov;

    std::vector< double > diagonals;

    // Fill covariance matrix

    for ( int i = 0; i < n_lin_bins; ++i ){

        vector< double > temp;

        for ( int j = 0; j < n_lin_bins; ++j ){
            
            double cov = h_ddxsec->GetBinContent(i+1) * h_ddxsec->GetBinContent(j+1) * err_flux * err_flux * (1/flux_int) * (1/flux_int);

            temp.push_back(cov);
        
            cov_file << setprecision(4) << setw(4) << cov << ", ";

            if ( i == j ){
            
                diagonals.push_back(cov);

            }

        }

        flux_cov.push_back(temp);

        cov_file << endl;
    }

    // Set the error on the 1D bins using the diagonal elements of the covariance matrix
    
    for ( int i = 0; i < n_lin_bins; ++i ){
        
        h_ddxsec_linear->SetBinError(i, diagonals[i]);
    
    }

    h_ddxsec_linear->Draw("hist");
    h_ddxsec_linear->Draw("e1x0same");

    c_lin->SaveAs("working_dir/err_test/linearised_ddxsec.png");

    delete h_ddxsec_linear; 

    cov_file.close();
    
    h_ddxsec->SetStats(kFALSE);
    h_ddxsec->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_ddxsec->GetYaxis()->SetTitle("T_{#mu}");
    h_ddxsec->SetTitle("CC0#pi, d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]");
    h_ddxsec->Draw("colz");
      
    cout << " Scale : " << scaling << endl;
    cout << " Cos   : " << cos_bins << endl;
    cout << " T     : " << Tmu_bins << endl;
    cout << " Flux  : " << flux_int << endl;
    cout << " Target: " << tot_tgt << endl;
    cout << " POT   : " << POT << endl << endl;
    cout << " Scaling factor : " << scalar_norm << endl;
    cout << " Integrated tot : " << h_unfold_test->Integral() << endl << endl;

    char ddxsec_path[1024];
    
    stringstream temp_ddx;
    temp_ddx.clear();

    string image_ddx;
    image_ddx.clear();

    temp_ddx << "working_dir/err_test/" << pr_path << "/" << pr_path << "_2D_unfolding_ddxsec.png";
    
    image_ddx = temp_ddx.str();    

    strcpy( ddxsec_path, image_ddx.c_str() );

    c->SetRightMargin(0.13);
    c->SaveAs( ddxsec_path );
   
    // SLICING
    TH2D *h_true_ddxsec = new TH2D( *h_true_test );
    h_true_ddxsec->Scale( scalar_norm );

    TH2D *h_reco_ddxsec = new TH2D( *h_reco_test );
    h_reco_ddxsec->Scale( scalar_norm );

    Slices( h_unfold_test, h_true_test, h_reco_test, pr_path );

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


void GetTruth( TTree *tree, 
               int n_min_bins, 
               int n_max_bins, 
               int n_protons, 
               std::vector<double> & truth_T, 
               std::vector<double> & truth_cos, 
               std::vector<bool>   & truth_detectable, 
               std::vector<double> & impur_T, 
               std::vector<double> & impur_cos ) {

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
    for ( int i = n_min_bins; i < n_max_bins; ++i ){
    
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

        TRandom *random = new TRandom();

        if ( ( n_protons != -1 && nfp == n_protons ) || ( n_protons == -1 ) ){
            // If there are no pions fill the kinetic energy and cos theta with the muon energy   
            if ( cc == 1 && nfpip + nfpim + nfpi0 == 0 ){
                // Calculate the kinetic energy for muons
                if ( fspl == 13 ){
 
                    // Energy of the final state primary lepton
                    double T_mu = e_mu - m_mu;
 
                    bool isDetectable = false;
                    if ( n_protons == -1 ) {
                        isDetectable = ( isReconstructed( T_mu, cos_mu, random ) );
                    }
                    else{
                        isDetectable = ( isReconstructed( T_mu, cos_mu, random ) && n_protons == n_detectable_pr);
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
    /*
        else if ( ( n_protons != -1 && nfp != n_protons ) ){
            // If there are no pions fill the kinetic energy and cos theta with the muon energy   
            if ( cc == 1 && nfpip + nfpim + nfpi0 == 0 ){
                // Calculate the kinetic energy for muons
                if ( fspl == 13 ){
 
                    // Energy of the final state primary lepton
                    double T_mu = e_mu - m_mu;
 
                    if ( T_mu > T_thresh_mu && nfp == n_detectable_pr ) {
                        impur_T.push_back(T_mu);
                        impur_cos.push_back(cos_mu);
                    }
                }
            }
        }*/
    }
}

bool isReconstructed(double tmu, double cos, TRandom *_rand){

    if(tmu<0.05) return 0;
    
    //CCQE efficiency as function of momentum from Fig10.b of uboone DocDB 7561
    
    double effmom[] = {0.79,0.95,0.975,0.98,0.98,0.99,0.99,0.995,1.0,1.0};
    
    //CCQE efficiency as function of theta from Fig10.c of uboone DocDB 7561
    
    double efftheta[] = {1,1,1,1,1,1,0.995,0.985,1,0.985,0.995,1,1,1,1,0.995,0.995,0.995,0.99,0.98,0.97,0.97,0.975,0.98,0.975,0.97,0.95,0.945,0.94,0.92,0.87,0.91};
    
    //Find position of muon momentum in momentum efficiency array
    
    int posmom = 9; 
    
    for(int i=0; i<10; i++){
    
      if(tmu<0.1*(i+1)) {posmom = i; break;}
    
    }
    
    //Find position of muon angle in angle efficiency array
    
    int postheta = 0; 
    
    for(int j=0; j<32; j++){
    
      if(TMath::ACos(cos)<0.09817*(j+1)) {postheta = j; break;}
    
    }
    
    //Don't have correlations between angle and momentum efficiencies, take the smallest value to be conservative
    
    double efficiency = 1; 
    
    if(effmom[posmom]<efftheta[postheta]) efficiency = effmom[posmom];
    
    else efficiency = efftheta[postheta];
    
    //Generate random number between 0 and 1
    
    double prob = _rand->Uniform(0,1);
    
    if(prob<=efficiency) return 1;
    
    else return 0;

}

bool isContained(double tmu, double cos, TRandom *_rand){

      // Values from http://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.png

      double momentum[] = {0.047,0.056,0.068,0.085,0.100,0.153,0.176,0.222,0.287,0.392,0.495,0.900,1.101,1.502,2.103}; //GeV
    
      double range[] = {0.007,0.013,0.024,0.047,0.076,0.221,0.304,0.482,0.761,1.236,1.708,3.534,4.415,6.126,8.610}; //m
        
      // Find position of muon momentum in momentum array
    
      int pos = 15;
    
      for(int i=0; i<15; i++){
    
        if(tmu<momentum[i]) {pos = i; break;}
    
      }
    
      //Use linear interpolation to calculate expected range in LAr
    
      double rmu = 0;
    
      if(pos==15) rmu = 9;
    
      else if(pos==0) rmu = 0.007;
    
      else{
    
        rmu = range[pos-1]+(tmu-momentum[pos-1])*(range[pos]-range[pos-1])/(momentum[pos]-momentum[pos-1]);
    
      }
    
      //Generate a random position inside detector
    
      double zstart = _rand->Uniform(0,5.312);
    
      double ystart = _rand->Uniform(0,4.294);
    
      double xstart = _rand->Uniform(0,4.312);
    
      double phi = _rand->Uniform(0,6.283);
    
      double theta = TMath::ACos(cos);
    
      //Calculate end position
    
      double zend = zstart + rmu*cos;
    
      double yend = ystart + rmu*TMath::Sin(theta)*TMath::Sin(phi);
    
      double xend = xstart + rmu*TMath::Sin(theta)*TMath::Cos(phi);
    
      //If range * cos thetamu <= z position
    
      if(xend>=0 && xend<=4.312 && yend>=0 && yend<=4.294 && zend>=0 && zend<=5.312) return 1;
    
      else return 0;

}

double RangeSmear(double ke, double m, ROOT::Math::GSLRngMT *_random_gen){

      //For contained muons use range based bias and resolution
    
      //Values from Fig 12 of https://arxiv.org.png/1703.06187.png
    
      double bias[] = {-0.0035,-0.0059,-0.0047,-0.0059,-0.0035,-0.0029,-0.0076,-0.0059,0.0006};
    
      double resolution[] = {0.017,0.021,0.023,0.026,0.025,0.030,0.030,0.040,0.032};


      double e = ke + m;

      int pos = 8;

      for(int i=0; i<9; i++){

        if(e<(0.33+0.186*(i+1))) {pos = i; break;}

    }

    /* The lognormal version
    double var_mu     = TMath::Power( ke * resolution[pos], 2 );
    double sigma_mu   = TMath::Sqrt( TMath::Log( 1 + ( var_mu / TMath::Power(ke, 2 ) ) ) );
    double zeta_mu    = TMath::Log( ke * ( 1 / TMath::Sqrt( 1 + ( var_mu / TMath::Power( ke, 2 ) ) ) ) );
    double ke_smear   = ke + _random_gen->LogNormal( zeta_mu, sigma_mu )+bias[pos]*ke;
      */       

    double e_smear = e + _random_gen->Gaussian(resolution[pos]*e)+bias[pos]*e;
    
    double ke_smear = e_smear - m;

    if(ke_smear<0) ke_smear = 0;  return ke_smear;

}

double McsSmear(double ke, ROOT::Math::GSLRngMT *_random_gen){

    //For exiting muons use multiple coulomb scattering bias and resolution

    //Values from Fig 5 of https://arxiv.org.png/1703.06187.png

    double bias[] = {0.0273,0.0409,0.0352,0.0250,0.0227,0.0068,0.0364,0.0273,0.0227};//0.0409

    double resolution[] = {0.127,0.145,0.143,0.141,0.164,0.177,0.250,0.266,0.341};//0.145

    int pos = 8;

    for(int i=0; i<9; i++){

        if(ke<(0.34+0.41*(i+1))) {pos = i; break;}

    }

    /* The lognormal version
    double var_mu     = TMath::Power( ke * resolution[pos], 2 );
    double sigma_mu   = TMath::Sqrt( TMath::Log( 1 + ( var_mu / TMath::Power(ke, 2 ) ) ) );
    double zeta_mu    = TMath::Log( ke * ( 1 / TMath::Sqrt( 1 + ( var_mu / TMath::Power( ke, 2 ) ) ) ) );
    double lognorm_mu = _random_gen->LogNormal( zeta_mu, sigma_mu );
    double ke_smear   = ke + lognorm_mu + bias[pos]*ke;
      */       

    double var = _random_gen->Gaussian(resolution[pos]*ke);

    double ke_smear = ke + var + bias[pos]*ke;

    if(ke_smear<0) ke_smear = 0;  return ke_smear;

}

void Smear( const double              & mass, 
            const std::vector<double> & truth_T, 
            const std::vector<double> & truth_cos, 
            std::vector<double>       & smear_T, 
            std::vector<double>       & smear_cos ){

    // Initiate the random number generation
    ROOT::Math::GSLRngMT *_random_gen = new ROOT::Math::GSLRngMT;
    _random_gen->Initialize();
    _random_gen->SetSeed( time( NULL ) );

    // Random number for position in the detector
    TRandom *random = new TRandom();

    int n_values = truth_T.size();

    if ( truth_T.size() != truth_cos.size() ){
        std::cerr << " Vectors must be the same length " << endl;
        exit(1);
    }

    // Event by event, generate Tmu_prime and Tpi_prime: lognormal
    // Then find thetamu_prime and thetapi_prime: gaussian
    for ( int i = 0; i < n_values; ++i ){
 
        // Find out if the event was reconstructed
        //if ( isReconstructed( truth_T[i], truth_cos[i], random ) ){

            // Find out if the event was contained
            // If contained: use Range smear
            //               smear cos as normal
            // If exiting  : use MCS smear
            //               smear cos as normal
            if ( isContained(truth_T[i], truth_cos[i], random ) ){
                
                // -------------------------------------------------------
                //                   Kinetic energy
                // -------------------------------------------------------
                // Calculate the mean and sigma for the LogNormal function
                //     zeta  = TMath::Log( m * ( 1 / sqrt( 1 + ( var / pow( m, 2 ) ) ) ) );
                //     sigma = sqrt( log( 1 + ( var / pow( m, 2 ) ) ) );
                //     m     = expectation value = Tl
                //     var   = variance = s.d.^2 = ( Tl * 0.1 ) ^ 2
                double range_smeared_mu = RangeSmear( truth_T[i], mass, _random_gen );
    
                smear_T.push_back(range_smeared_mu);
                    
    
                // -------------------------------------------------------
                //                  Cos theta
                // -------------------------------------------------------
             
                // Calculate the mean and sigma for the LogNormal function
                //      theta = acos(costtheta)
                //      var   = 5 degrees
             
                double sd_thetamu    = TMath::Pi() / 36; // 5 degrees
                double gaus_theta    = TMath::ACos( truth_cos[i] ) + _random_gen->Gaussian( sd_thetamu );
                double gaus_costheta = TMath::Cos( gaus_theta );
                 
                smear_cos.push_back(gaus_costheta);
      
            }
            else{
    
                // -------------------------------------------------------
                //                   Kinetic energy
                // -------------------------------------------------------
                // Calculate the mean and sigma for the LogNormal function
                //     zeta  = TMath::Log( m * ( 1 / sqrt( 1 + ( var / pow( m, 2 ) ) ) ) );
                //     sigma = sqrt( log( 1 + ( var / pow( m, 2 ) ) ) );
                //     m     = expectation value = Tl
                //     var   = variance = s.d.^2 = ( Tl * 0.1 ) ^ 2
                double mcs_smeared_mu = McsSmear( truth_T[i], _random_gen );
    
                smear_T.push_back(mcs_smeared_mu);
                 
    
                // -------------------------------------------------------
                //                  Cos theta
                // -------------------------------------------------------
          
                // Calculate the mean and sigma for the LogNormal function
                //      theta = acos(costtheta)
                //      var   = 5 degrees
                 
                double sd_thetamu    = TMath::Pi() / 36; // 5 degrees
                double gaus_theta    = TMath::ACos( truth_cos[i] ) + _random_gen->Gaussian( sd_thetamu );
                double gaus_costheta = TMath::Cos( gaus_theta );
                 
                smear_cos.push_back(gaus_costheta);
            }
        }
   // }
    delete _random_gen;
    delete random;
}

void GetResponse( const std::vector<double> & truth_T, 
                  const std::vector<double> & truth_cos, 
                  const std::vector<bool>   & truth_detectable, 
                  const std::vector<double> & smear_T, 
                  const std::vector<double> & smear_cos, 
                  const std::vector<double> & smear_impur_T, 
                  const std::vector<double> & smear_impur_cos, 
                  RooUnfoldResponse         & response) {

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

        conv << setprecision(4) << "working_dir/err_test/" << n_pr << "/Tmu_slice_" << low_edge_T << ";" << up_edge_T << ".png";
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

        conv1 << setprecision(4) << "working_dir/err_test/" << n_pr << "/cos_thetamu_slice_" << low_edge_cos << ";" << up_edge_cos << ".png";
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

