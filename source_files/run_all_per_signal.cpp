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

#ifndef RUN_ALL_PER_SIGNAL_CPP
#define RUN_ALL_PER_SIGNAL_CPP

#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/run_all_per_signal.h"

using namespace xsec;

TH2D *RunAll( TH1D *h_flux, 
              const std::vector< Event >       &training_events,
              const std::vector< Event >       &testing_events,
              const Interaction                &signal, 
              const std::vector< Interaction > &backgrounds, 
              const char file_path[1024] ){
    
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

    // Get the filename for each 2D histogram
    std::stringstream conv;
    conv.clear();
    
    std::string title;
    title.clear();
    
    char file_name1[1024];
    
    conv << std::setprecision(4) << file_path << "cc0pi_true_2D.png";
    title = conv.str();
            
    strcpy( file_name1, title.c_str() );
    
    Set2DHistInfo( h_true_test, "cos#theta_{#mu}", "T_{#mu}", "True #mu kinematics", "colz" );
   
    c_test->SaveAs( file_name1 );

    conv.clear();
    title.clear();
    
    conv.str(std::string());
    
    conv << std::setprecision(4) << file_path << "cc0pi_reco_2D.png";
    title = conv.str();
    
    char file_name2[1024];
    strcpy( file_name2, title.c_str() );

    Set2DHistInfo( h_reco_test, "cos#theta_{#mu}", "T_{#mu}", "Data-like kinematics", "colz" );
    
    c_test->SaveAs( file_name2 );

    conv.clear();
    title.clear();
    
    conv.str(std::string());

    conv << std::setprecision(4) << file_path << "cc0pi_unfolded_2D.png";
    title = conv.str();
    
    char file_name3[1024];
    strcpy( file_name3, title.c_str() );

    Set2DHistInfo( h_unfold_test, "cos#theta_{#mu}", "T_{#mu}", "Unfolded kinematics", "colz" );
    
    c_test->SaveAs( file_name3 );

    conv.clear();
    title.clear();
    
    conv.str(std::string());

    conv << std::setprecision(4) << file_path << "cc0pi_ddxsec_2D.png";
    title = conv.str();
    
    char file_name4[1024];
    strcpy( file_name4, title.c_str() );

    Set2DHistInfo( h_ddxsec, "cos#theta_{#mu}", "T_{#mu}", "CC0#pi, d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]", "colz" );
    
    c_test->SaveAs( file_name4 );

    delete c_test;
    delete h_true_train;
    delete h_reco_train;
    delete h_true_test;
    delete h_reco_test;
    delete h_unfold_test;
    delete h_ddxsec_linear;
    
    return h_ddxsec;

}

#endif
