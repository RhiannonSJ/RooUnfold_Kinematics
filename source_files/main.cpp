/* The Main 
 *
 *  Run all functions to get the cross section distributions
 *  for various final state topologies
 *
 * =============================================
 *
 *  Author: Rhiannon Jones
 *  Date  : July 2017
 *
 * =============================================
 */

#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/all_class_headers.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/unfolding_functions.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/cross_section.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/slicing.h"

using namespace xsec;

int main_func(){
   
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
    //         Get the SBND flux histogram and trees from the root files
    //==============================================================================
    
    TH1D *h_flux = (TH1D*) fflux.Get("h_numu_110m");

    TTree *gst_train = (TTree*) f_train.Get("gst");
    TTree *gst_test  = (TTree*) f_test.Get("gst");
    
    // Typedef for the map
    typedef std::map< std::vector< int >, int > top_map; 

    // Test the filling of the event list and cout some events
    
    // Define the vectors of events
    std::vector< Event > training_events;
    std::vector< Event > testing_events;
    
    // Fill the vectors
    LoadEventList( gst_train, training_events );
    LoadEventList( gst_test,  testing_events );

    // Define the chosen signal and background interactions 
    // For the purpose of testing, use 
    //      signal: CC 0pi, N protons  
    //      backgrounds: NC 1pi+
    //                 : NC 1pi-
    // Maps
    top_map signal_map;
   
    std::vector< int > cc_0pi_mu;
    std::vector< int > cc_0pi_pi;
    
    cc_0pi_mu.push_back( 13 );
    cc_0pi_pi.push_back( 211 );
    cc_0pi_pi.push_back(-211 );
    cc_0pi_pi.push_back( 111 );
    
    signal_map.insert( std::make_pair( cc_0pi_mu, 1 ) );
    signal_map.insert( std::make_pair( cc_0pi_pi, 0 ) );
    
    // Signal 
    Interaction signal( signal_map, 13, true );

    // Vector of interactions for the background
    std::vector< Interaction > backgrounds;

    top_map background_map_nc1pip;
    top_map background_map_nc1pim;
    
    std::vector< int > nc_1pi_pip;
    std::vector< int > nc_1pi_pim;
    std::vector< int > nc_1pi_pi0;
    
    nc_1pi_pip.push_back( 211 );
    nc_1pi_pim.push_back(-211 );
    nc_1pi_pi0.push_back( 111 );
    
    background_map_nc1pip.insert( std::make_pair( nc_1pi_pip, 1 ) );
    background_map_nc1pip.insert( std::make_pair( nc_1pi_pim, 0 ) );
    background_map_nc1pip.insert( std::make_pair( nc_1pi_pi0, 0 ) );

    background_map_nc1pim.insert( std::make_pair( nc_1pi_pip, 0 ) );
    background_map_nc1pim.insert( std::make_pair( nc_1pi_pim, 1 ) );
    background_map_nc1pim.insert( std::make_pair( nc_1pi_pi0, 0 ) );

    // Background interactions
    Interaction background_1( background_map_nc1pip, 211, false, 0.2 );
    Interaction background_2( background_map_nc1pim, -211, false, 0.2 );

    backgrounds.push_back( background_1 );
    backgrounds.push_back( background_2 );
    
    // Initialise histograms
    // Training
    TH2D *h_true_train = new TH2D( "h_true_train", " true ", 20, -1, 1, 18, 0.2, 2 );
    TH2D *h_reco_train = new TH2D( "h_reco_train", " reco ", 20, -1, 1, 18, 0.2, 2 );

    // Testing
    TH2D *h_true_test = new TH2D( "h_true_test", " true ", 20, -1, 1, 18, 0.2, 2 );
    TH2D *h_reco_test = new TH2D( "h_reco_test", " reco ", 20, -1, 1, 18, 0.2, 2 );

    // Get the response matrix
    RooUnfoldResponse response( h_reco_train, h_true_train );

    GetResponse( training_events, signal, backgrounds, response );
    
    // Fill the testing histograms
    GetTrueRecoHists( training_events, signal, backgrounds, h_true_test, h_reco_test );

    // Scaling factor for full POT and G16_01b
    double scale_POT = 7.56016;
    
    h_true_test->Scale(scale_POT);
    h_reco_test->Scale(scale_POT);
    
    // Unfold 
    RooUnfoldBayes unfold( &response, h_reco_test, 1 );

    // Get unfolded histogram with errors
    TH2D *h_unfold_test = (TH2D*) unfold.Hreco((RooUnfold::ErrorTreatment)2);
   
    // Cross-section
    double scalar_norm = GetCrossSectionScale( h_unfold_test, h_flux );
    TH2D *h_ddxsec = new TH2D ( *h_unfold_test );
    h_ddxsec->Scale( scalar_norm );

    // Get the ddxsec errors
    // Linearised histogram
    TH1D *h_ddxsec_linear = new TH1D("h_ddxsec_linear", "", 360, -0.5, 359.5);
    
    // File path
    const char file_path[1024] = "/hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_0pi/all_signal/";
    // Set 10% flux error on cross-section
    Set2DErrors( h_ddxsec, h_unfold_test, h_ddxsec_linear, 0.1, file_path, scalar_norm );

    // Get the slices
    Slices( h_unfold_test, h_true_test, h_reco_test, file_path );
    Slices( h_ddxsec, file_path );

    // Get the reconstructed event list
    std::vector< Event >    reco_events;
    std::vector< Particle > primaries;
    
    GetRecoEventList( testing_events, signal, backgrounds, primaries, reco_events );
    Characterisation( h_reco_test, primaries, reco_events, file_path );
    
    // Drawing
    // Canvas
    TCanvas *c_test = new TCanvas();

    // File path for saving CC_0Pi, n proton distributions 
    //   /hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_0pi/all_signal

    // Set the information for each histogram and then save them
    Set2DHistInfo( h_true_test, "cos#theta_{#mu}", "T_{#mu}", "True #mu kinematics", "colz" );
   
    c_test->SaveAs( "/hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_0pi/all_signal/cc0pi_true_2D.png" );

    Set2DHistInfo( h_reco_test, "cos#theta_{#mu}", "T_{#mu}", "Data-like kinematics", "colz" );
    
    c_test->SaveAs( "/hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_0pi/all_signal/cc0pi_reco_2D.png" );

    Set2DHistInfo( h_unfold_test, "cos#theta_{#mu}", "T_{#mu}", "Unfolded kinematics", "colz" );
    
    c_test->SaveAs( "/hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_0pi/all_signal/cc0pi_unfolded_2D.png" );

    Set2DHistInfo( h_ddxsec, "cos#theta_{#mu}", "T_{#mu}", "CC0#pi, d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]", "colz" );
    
    c_test->SaveAs( "/hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_0pi/all_signal/cc0pi_ddxsec_2D.png" );

    delete c_test;
    delete h_true_train;
    delete h_reco_train;
    delete h_true_test;
    delete h_reco_test;
    delete h_unfold_test;
    delete h_ddxsec;
    
    //==============================================================================
    //         Things previously used to test stages of the process
    //==============================================================================
    
    /*
    std::cout << " Size of training list : " << training_events.size() << std::endl;
    std::cout << " Size of testing list  : " << testing_events.size() << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << " Random events          : " << std::endl;
    std::cout << testing_events[100] << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << testing_events[10110] << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << testing_events[12] << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << testing_events[33152] << std::endl;
    std::cout << "-------------------------------------" << std::endl;
*/

    /*
     *          QUICK TEST
    // CC 0pi, 1pi event
    Event ev; 

    ev.Add( Particle( 2212, 0.05, 0.8 ) );
    ev.Add( Particle( 13,   0.05, 0.8 ) );

    // Map to hold chosen signal
    top_map signal_cc_0pi;

    std::vector< int > cc_0pi_mu;
    std::vector< int > cc_0pi_p;
    std::vector< int > cc_0pi_pi;
    
    cc_0pi_mu.push_back( 13 );
    cc_0pi_p.push_back( 2212 );
    cc_0pi_pi.push_back( 211 );
    cc_0pi_pi.push_back(-211 );
    cc_0pi_pi.push_back( 111 );
    
    signal_cc_0pi.insert( std::make_pair( cc_0pi_mu, 1 ) );
    signal_cc_0pi.insert( std::make_pair( cc_0pi_p,  1 ) );
    signal_cc_0pi.insert( std::make_pair( cc_0pi_pi, 0 ) );

    std::cout << "--------------------------------" << std::endl;
    std::cout << "--------- True CC 0pi ----------" << std::endl;
    std::cout << ev << std::endl;
    std::cout << ( ev.CheckIfTrue( true, signal_cc_0pi ) ? " True signal " : " Not true signal " ) << std::endl;
    std::cout << ( ev.CheckIfReconstructed( signal_cc_0pi ) ? " Reco signal " : " Not reco signal " ) << std::endl;
    std::cout << "--------------------------------" << std::endl;
    */

    return 0;

}
