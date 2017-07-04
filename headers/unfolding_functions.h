#ifndef UNFOLDING_FUNCTIONS_H
#define UNFOLDING_FUNCTIONS_H

//using std::cout;
//using std::endl;


#include "../src/RooUnfoldResponse.h"
#include "../src/RooUnfoldBayes.h"

#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/roo_unfold.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/all_class_headers.h"
using namespace xsec;


typedef std::map< std::vector< int >, int > top_map; 

// Load the list of events from the gst tree and push onto a vector of events
void LoadEventList( TTree *gst, std::vector< Event > &event_list );


// Get the response matrix for the chosen topology
void GetResponse( const std::vector< Event >       & event_list,
                  const Interaction                & interaction,
                  const std::vector< Interaction > & background,
                  RooUnfoldResponse                & response );
                             
// Get the true and reconstructed histograms
void GetTrueRecoHists( const std::vector< Event >       & event_list,
                       const Interaction                & interaction,
                       const std::vector< Interaction > & background,
                       TH2D *true_hist,
                       TH2D *reco_hist );

// Get event list of reconstructed events
void GetRecoEventList( const std::vector< Event >       & event_list,
                       const Interaction                & interaction,
                       const std::vector< Interaction > & background,
                       std::vector< Particle >          & primary_list,
                       std::vector< Event >             & reco_event_list );

// Set info for 2D histograms
void Set2DHistInfo( TH2D *hist, const char x_axis[1024], const char y_axis[1024], const char title[1024], const char draw_opt[1024] );

#endif


