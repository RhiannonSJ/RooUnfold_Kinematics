#ifndef SLICING_H
#define SLICING_H

#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/roo_unfold.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/cov_func.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/all_class_headers.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/unfolding_functions.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/cross_section.h"

using namespace xsec;

typedef std::map< std::vector< int >, int > top_map; 

// Unfolding, reco and truth comparison slices
void Slices ( TH2D *h_unfolded, TH2D *h_true, TH2D *h_reco, const char file[1024] );

// Cross-section slices with errors
void Slices ( TH2D *h_ddxsec, const char file[1024] );

// Get the characterised slices
void Characterisation ( TH2D *h_smeared, std::vector< Particle > primary, std::vector< Event > reco_events,  const char file_path[1024] );

#endif


