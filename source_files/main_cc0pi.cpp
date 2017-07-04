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
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/run_all_per_signal.h"

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
    top_map signal_map_all;
    top_map signal_map_1p;
    top_map signal_map_2p;
    top_map signal_map_3p;
   
    std::vector< int > cc_0pi_mu;
    std::vector< int > cc_0pi_pi;
    std::vector< int > cc_0pi_p;
    
    cc_0pi_mu.push_back( 13 );
    cc_0pi_pi.push_back( 211 );
    cc_0pi_pi.push_back(-211 );
    cc_0pi_pi.push_back( 111 );
    cc_0pi_p.push_back( 2212 );
    
    signal_map_all.insert( std::make_pair( cc_0pi_mu, 1 ) );
    signal_map_all.insert( std::make_pair( cc_0pi_pi, 0 ) );
    
    signal_map_1p.insert( std::make_pair( cc_0pi_mu, 1 ) );
    signal_map_1p.insert( std::make_pair( cc_0pi_pi, 0 ) );
    signal_map_1p.insert( std::make_pair( cc_0pi_p,  1 ) );
    
    signal_map_2p.insert( std::make_pair( cc_0pi_mu, 1 ) );
    signal_map_2p.insert( std::make_pair( cc_0pi_pi, 0 ) );
    signal_map_2p.insert( std::make_pair( cc_0pi_p,  2 ) );
    
    signal_map_3p.insert( std::make_pair( cc_0pi_mu, 1 ) );
    signal_map_3p.insert( std::make_pair( cc_0pi_pi, 0 ) );
    signal_map_3p.insert( std::make_pair( cc_0pi_p,  3 ) );
    
    // Signal all
    Interaction signal_all( signal_map_all, 13, true );
    Interaction signal_1p(  signal_map_1p,  13, true );
    Interaction signal_2p(  signal_map_2p,  13, true );
    Interaction signal_3p(  signal_map_3p,  13, true );

    // Vector of interactions for the background
    std::vector< Interaction > backgrounds_all;
    std::vector< Interaction > backgrounds_1p;
    std::vector< Interaction > backgrounds_2p;
    std::vector< Interaction > backgrounds_3p;

    top_map background_map_nc1pip_all;
    top_map background_map_nc1pim_all;
    
    top_map background_map_nc1pip_1p;
    top_map background_map_nc1pim_1p;
    
    top_map background_map_nc1pip_2p;
    top_map background_map_nc1pim_2p;
    
    top_map background_map_nc1pip_3p;
    top_map background_map_nc1pim_3p;
    
    std::vector< int > nc_1pi_pip;
    std::vector< int > nc_1pi_pim;
    std::vector< int > nc_1pi_pi0;
    std::vector< int > nc_1pi_p;
    
    nc_1pi_pip.push_back( 211 );
    nc_1pi_pim.push_back(-211 );
    nc_1pi_pi0.push_back( 111 );
    nc_1pi_p.push_back( 2212 );
    
    background_map_nc1pip_all.insert( std::make_pair( nc_1pi_pip, 1 ) );
    background_map_nc1pip_all.insert( std::make_pair( nc_1pi_pim, 0 ) );
    background_map_nc1pip_all.insert( std::make_pair( nc_1pi_pi0, 0 ) );

    background_map_nc1pim_all.insert( std::make_pair( nc_1pi_pip, 0 ) );
    background_map_nc1pim_all.insert( std::make_pair( nc_1pi_pim, 1 ) );
    background_map_nc1pim_all.insert( std::make_pair( nc_1pi_pi0, 0 ) );

    background_map_nc1pip_1p.insert( std::make_pair( nc_1pi_pip, 1 ) );
    background_map_nc1pip_1p.insert( std::make_pair( nc_1pi_pim, 0 ) );
    background_map_nc1pip_1p.insert( std::make_pair( nc_1pi_pi0, 0 ) );
    background_map_nc1pip_1p.insert( std::make_pair( nc_1pi_p,   1 ) );

    background_map_nc1pim_1p.insert( std::make_pair( nc_1pi_pip, 0 ) );
    background_map_nc1pim_1p.insert( std::make_pair( nc_1pi_pim, 1 ) );
    background_map_nc1pim_1p.insert( std::make_pair( nc_1pi_pi0, 0 ) );
    background_map_nc1pim_1p.insert( std::make_pair( nc_1pi_p,   1 ) );
    
    background_map_nc1pip_2p.insert( std::make_pair( nc_1pi_pip, 1 ) );
    background_map_nc1pip_2p.insert( std::make_pair( nc_1pi_pim, 0 ) );
    background_map_nc1pip_2p.insert( std::make_pair( nc_1pi_pi0, 0 ) );
    background_map_nc1pip_2p.insert( std::make_pair( nc_1pi_p,   2 ) );

    background_map_nc1pim_2p.insert( std::make_pair( nc_1pi_pip, 0 ) );
    background_map_nc1pim_2p.insert( std::make_pair( nc_1pi_pim, 1 ) );
    background_map_nc1pim_2p.insert( std::make_pair( nc_1pi_pi0, 0 ) );
    background_map_nc1pim_2p.insert( std::make_pair( nc_1pi_p,   2 ) );
    
    background_map_nc1pip_3p.insert( std::make_pair( nc_1pi_pip, 1 ) );
    background_map_nc1pip_3p.insert( std::make_pair( nc_1pi_pim, 0 ) );
    background_map_nc1pip_3p.insert( std::make_pair( nc_1pi_pi0, 0 ) );
    background_map_nc1pip_3p.insert( std::make_pair( nc_1pi_p,   3 ) );

    background_map_nc1pim_3p.insert( std::make_pair( nc_1pi_pip, 0 ) );
    background_map_nc1pim_3p.insert( std::make_pair( nc_1pi_pim, 1 ) );
    background_map_nc1pim_3p.insert( std::make_pair( nc_1pi_pi0, 0 ) );
    background_map_nc1pim_3p.insert( std::make_pair( nc_1pi_p,   3 ) );
    
    // Background interactions
    Interaction background_1( background_map_nc1pip_all, 211, false, 0.2 );
    Interaction background_2( background_map_nc1pim_all, -211, false, 0.2 );

    Interaction background_1_1p( background_map_nc1pip_1p, 211, false, 0.2 );
    Interaction background_2_1p( background_map_nc1pim_1p, -211, false, 0.2 );

    Interaction background_1_2p( background_map_nc1pip_2p, 211, false, 0.2 );
    Interaction background_2_2p( background_map_nc1pim_2p, -211, false, 0.2 );

    Interaction background_1_3p( background_map_nc1pip_3p, 211, false, 0.2 );
    Interaction background_2_3p( background_map_nc1pim_3p, -211, false, 0.2 );

    backgrounds_all.push_back( background_1 );
    backgrounds_all.push_back( background_2 );
   
    backgrounds_1p.push_back( background_1_1p );
    backgrounds_1p.push_back( background_2_1p );
   
    backgrounds_2p.push_back( background_1_2p );
    backgrounds_2p.push_back( background_2_2p );
   
    backgrounds_3p.push_back( background_1_3p );
    backgrounds_3p.push_back( background_2_3p );
   
    // File path
    const char file_path[1024]  = "/hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_0pi/all_signal/";
    const char file_path1[1024] = "/hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_0pi/1_p/";
    const char file_path2[1024] = "/hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_0pi/2_p/";
    const char file_path3[1024] = "/hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_0pi/3_p/";
   
    TH2D *h_ddxsec_all = RunAll( h_flux, training_events, testing_events, signal_all, backgrounds_all, file_path );
    TH2D *h_ddxsec_1p  = RunAll( h_flux, training_events, testing_events, signal_1p,  backgrounds_1p,  file_path1 );
    TH2D *h_ddxsec_2p  = RunAll( h_flux, training_events, testing_events, signal_2p,  backgrounds_2p,  file_path2 );
    TH2D *h_ddxsec_3p  = RunAll( h_flux, training_events, testing_events, signal_3p,  backgrounds_3p,  file_path3 );

    std::vector< std::string > names;
    names.push_back( "CC 0#pi, Np" );
    names.push_back( "CC 0#pi, 1p" );
    names.push_back( "CC 0#pi, 2p" );
    names.push_back( "CC 0#pi, 3p" );

    std::vector< TH2D* > ddxsec;
    ddxsec.push_back( h_ddxsec_all );
    ddxsec.push_back( h_ddxsec_1p );
    ddxsec.push_back( h_ddxsec_2p );
    ddxsec.push_back( h_ddxsec_3p );
    
    SignalComparison( ddxsec, "/hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_0pi/comparison/", names ); 

    delete h_ddxsec_all;
    delete h_ddxsec_1p;
    delete h_ddxsec_2p;
    delete h_ddxsec_3p;
    
    return 0;

}
