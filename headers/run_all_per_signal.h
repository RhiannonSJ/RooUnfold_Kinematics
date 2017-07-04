#ifndef RUN_ALL_PER_SIGNAL_H
#define RUN_ALL_PER_SIGNAL_H


#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/roo_unfold.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/cov_func.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/all_class_headers.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/unfolding_functions.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/cross_section.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/slicing.h"

using namespace xsec;

typedef std::map< std::vector< int >, int > top_map; 

TH2D *RunAll( TH1D *h_flux,
              const std::vector< Event >       &training_events,
              const std::vector< Event >       &testing_events,
              const Interaction                &signal, 
              const std::vector< Interaction > &backgrounds, 
              const char file_path[1024] );

#endif


