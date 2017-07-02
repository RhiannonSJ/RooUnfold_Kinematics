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
    // Get the SBND flux histogram and trees from the root files
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
