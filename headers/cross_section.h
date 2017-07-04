#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include "../src/RooUnfoldResponse.h"
#include "../src/RooUnfoldBayes.h"

#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/roo_unfold.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/cov_func.h"
#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/all_class_headers.h"

using namespace xsec;

typedef std::map< std::vector< int >, int > top_map; 

// Cross section scaling calculation
double GetCrossSectionScale( TH2D *unfold, TH1D *flux );

// Linearisation function
void GetLinearisatedHistogram( TH2D *xsec, TH1D *linear);

// Get covariance matrix
TMatrixDSym GetCovariance( TH1D *linear, double flux_err, const char file[1024] );

// Set errors on the linear histogram
void SetLinearErrors( TMatrixDSym covariance, TH1D *linear, const char file[1024] );

// Set errors on the 2D histogram from the covariance matrix
void SetDDXSecErrors( TH1D *linear, TH2D *xsec, TMatrixDSym covariance );

// Set total errors due to flux and statistics
void SetFluxStatsErrors( TH2D *unfold, TH2D *xsec, double norm );

// Get errors using linearisation
void Set2DErrors( TH2D *xsec, TH2D *unfold, TH1D *linear, double flux_err, const char file[1024], double norm );

#endif


