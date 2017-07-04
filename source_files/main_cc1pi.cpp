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
    //      signal: CC 1pi, N protons  
    //            : CC 1pip, N protons
    //            : CC 1pim, N protons
    //
    //      backgrounds: NC 2pi+, N protons
    //                 : NC 2pi-, N protons
    //                 : NC 2pi,  N protons                 
    // Maps
    top_map signal_map_all;
    top_map signal_map_pip;
    top_map signal_map_pim;
   
    std::vector< int > cc_0pi_mu;
    std::vector< int > cc_0pi_pi;
    std::vector< int > cc_0pi_pip;
    std::vector< int > cc_0pi_pim;
    
    cc_0pi_mu.push_back(   13 );
    cc_0pi_pi.push_back(  211 );
    cc_0pi_pi.push_back( -211 );
    cc_0pi_pip.push_back( 211 );
    cc_0pi_pim.push_back(-211 );
    
    signal_map_all.insert( std::make_pair( cc_0pi_mu,  1 ) );
    signal_map_all.insert( std::make_pair( cc_0pi_pi,  1 ) );
    
    signal_map_pip.insert( std::make_pair( cc_0pi_mu,  1 ) );
    signal_map_pip.insert( std::make_pair( cc_0pi_pip, 1 ) );
    signal_map_pip.insert( std::make_pair( cc_0pi_pim, 0 ) );
    
    signal_map_pim.insert( std::make_pair( cc_0pi_mu,  1 ) );
    signal_map_pim.insert( std::make_pair( cc_0pi_pip, 0 ) );
    signal_map_pim.insert( std::make_pair( cc_0pi_pim, 1 ) );
    
    // Signal all
    Interaction signal_all(  signal_map_all,  13, true );
    Interaction signal_pip(  signal_map_pip,  13, true );
    Interaction signal_pim(  signal_map_pim,  13, true );

    // Vector of interactions for the background
    std::vector< Interaction > backgrounds_all;
    std::vector< Interaction > backgrounds_pip;
    std::vector< Interaction > backgrounds_pim;

    top_map background_map_nc2pip;
    top_map background_map_nc2pim;
    top_map background_map_nc1pipm;
    
    std::vector< int > nc_1pi_pip;
    std::vector< int > nc_1pi_pim;
    std::vector< int > nc_1pi_pi0;
    
    nc_1pi_pip.push_back( 211 );
    nc_1pi_pim.push_back(-211 );
    nc_1pi_pi0.push_back( 111 );
    
    background_map_nc2pip.insert( std::make_pair( nc_1pi_pip, 2 ) );
    background_map_nc2pip.insert( std::make_pair( nc_1pi_pim, 0 ) );
    background_map_nc2pip.insert( std::make_pair( nc_1pi_pi0, 0 ) );

    background_map_nc2pim.insert( std::make_pair( nc_1pi_pip, 0 ) );
    background_map_nc2pim.insert( std::make_pair( nc_1pi_pim, 2 ) );
    background_map_nc2pim.insert( std::make_pair( nc_1pi_pi0, 0 ) );

    background_map_nc1pipm.insert( std::make_pair( nc_1pi_pip, 1 ) );
    background_map_nc1pipm.insert( std::make_pair( nc_1pi_pim, 1 ) );
    background_map_nc1pipm.insert( std::make_pair( nc_1pi_pi0, 0 ) );

    // Background interactions
    Interaction background_1( background_map_nc2pip,   211, false, 0.2 );
    Interaction background_2( background_map_nc2pim,  -211, false, 0.2 );
    Interaction background_3( background_map_nc1pipm,  211, false, 0.2 );
    Interaction background_4( background_map_nc1pipm, -211, false, 0.2 );

    backgrounds_all.push_back( background_1 );
    backgrounds_all.push_back( background_2 );
    backgrounds_all.push_back( background_3 );
    backgrounds_all.push_back( background_4 );
   
    backgrounds_pip.push_back( background_1 );
    backgrounds_pip.push_back( background_3 );
    backgrounds_pip.push_back( background_4 );
   
    backgrounds_pim.push_back( background_2 );
    backgrounds_pim.push_back( background_3 );
    backgrounds_pim.push_back( background_4 );
   
    // File path
    const char file_path[1024]  = "/hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_1pi/all/";
    const char file_path1[1024] = "/hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_1pi/pip/";
    const char file_path2[1024] = "/hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_1pi/pim/";
   
    TH2D *h_ddxsec_all = RunAll( h_flux, training_events, testing_events, signal_all, backgrounds_all,  file_path );
    TH2D *h_ddxsec_pip = RunAll( h_flux, training_events, testing_events, signal_pip, backgrounds_pip, file_path1 );
    TH2D *h_ddxsec_pim = RunAll( h_flux, training_events, testing_events, signal_pim, backgrounds_pim, file_path2 );

    std::vector< std::string > names;
    names.push_back( "CC 1#pi^{#pm}" );
    names.push_back( "CC 1#pi^{+}" );
    names.push_back( "CC 1#pi^{-}" );

    std::vector< TH2D* > ddxsec;
    ddxsec.push_back( h_ddxsec_all );
    ddxsec.push_back( h_ddxsec_pip );
    ddxsec.push_back( h_ddxsec_pim );
    
    SignalComparison( ddxsec, "/hepstore/rjones/Exercises/RooUnfold_Kinematics/unfolding_general/distributions/cc_1pi/comparison/", names ); 

    delete h_ddxsec_all;
    delete h_ddxsec_pip;
    delete h_ddxsec_pim;
    
    return 0;

}
