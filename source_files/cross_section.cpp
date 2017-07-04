#ifndef CROSS_SECTION_CPP
#define CROSS_SECTION_CPP

#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/cross_section.h"
using namespace xsec;

double GetCrossSectionScale( TH2D *unfold, TH1D *flux ){

    // Normalization variables
    double cos_bins = unfold->GetXaxis()->GetBinWidth(1);
    double Tmu_bins = unfold->GetYaxis()->GetBinWidth(1); // GeV
    // Units of flux: #/m^2/50MeV/10^6 POT
    double flux_int = flux->Integral("width");
    double Na       = 6.022e23;  //
    double M_fid    = 112000; // (kg / m^3) * m^3
    double A_ar     = 0.039948;  // kg / mol
    double tot_tgt  = ( Na * M_fid ) / ( A_ar ); // number of target particles  
    double POT      = 6.6e20;    // POT
    double scaling  = 5e8;      // Sort units from factor of 1 / ( integrated flux * Nt * bin widths * POT )
    double cm_conv  = 1e-38;     // take out factor of 10e-38 cm^2 
    double units_scale = 10e48; // Scale for units 

    double scalar_norm = scaling / ( cm_conv * cos_bins * Tmu_bins * flux_int * tot_tgt * POT );
    
    return scalar_norm;

}

// Linearisation function
void GetLinearisatedHistogram( TH2D *xsec, TH1D *linear){

    int n_x = xsec->GetNbinsX();
    int n_y = xsec->GetNbinsY();
    
    for ( int i = 0; i < n_y; ++i ){
        
        for ( int j = 0; j < n_x; ++j ){
            
            int k = Map(j,i,n_x,n_y) + 1;
            
            linear->SetBinContent( k , xsec->GetBinContent(j+1,i+1));        

        }
    }
}

// Get covariance matrix
TMatrixDSym GetCovariance( TH1D *linear, double flux_err, const char file[1024] ){

    int n_lin_bins = linear->GetNbinsX();

    TMatrixDSym covariance(n_lin_bins);

    // Fill covariance matrix
    for ( int i = 0; i < n_lin_bins; ++i ){

        for ( int j = 0; j < n_lin_bins; ++j ){
            
            covariance[i][j] = linear->GetBinContent(i+1) * linear->GetBinContent(j+1) * flux_err * flux_err;

        }
    }

    // For the file name
    stringstream name;
    name.clear();

    char file_path[1024];
    string temp;

    name << file << "covariance_matrix.root";
    temp = name.str();

    strcpy( file_path, temp.c_str() );
    
    TFile file_cov( file_path, "RECREATE");

    covariance.Write("covariance");

    TCanvas * c_cov = new TCanvas();

	TH2D h_cov( "h_data_cov", "Covariance;bin;bin",
  		         n_lin_bins , -0.5, n_lin_bins - 0.5,
				 n_lin_bins , -0.5, n_lin_bins - 0.5 )  ;


   	    for ( int i = 0; i < n_lin_bins; ++i ) {
   		    for ( int j = 0 ; j < n_lin_bins; ++j ) {
   		        h_cov.Fill( i, j, covariance[i][j] ) ;
   		    }
   		}

    // For the file name
    stringstream name1;
    name1.clear();

    char file_path1[1024];
    string temp1;

    name1 << file << "covariance_plot.png";
    temp1 = name1.str();

    strcpy( file_path1, temp1.c_str() );
    
    h_cov.Write("h_cov");
    h_cov.SetStats(kFALSE);
    h_cov.Draw("colz");
    c_cov->SaveAs( file_path1 );

    return covariance;

}

// Set errors due to flux on the linearised histogram
void SetLinearErrors( TMatrixDSym covariance, TH1D *linear, const char file[1024] ){

    // Set the error on the 1D bins using the diagonal elements of the covariance matrix
    
    int n_lin_bins = linear->GetNbinsX();

    for ( int i = 1; i <= n_lin_bins; ++i ){
        
        linear->SetBinError(i, TMath::Sqrt(covariance[i-1][i-1]));

    }

    // For the file name
    stringstream name;
    name.clear();

    char file_path[1024];
    string temp;

    name << file << "linearised_ddxsec.png";
    temp = name.str();

    strcpy( file_path, temp.c_str() );
    
    // Save the linearised histogram with errors
    TCanvas *c_lin = new TCanvas();
    linear->Draw("e1x0");
    c_lin->SaveAs( file_path );
    delete c_lin;

}

// Set errors on the 2D histogram using the linearised histogram
void SetDDXSecErrors( TH1D *linear, TH2D *xsec, TMatrixDSym covariance ){

    int n_x        = xsec->GetNbinsX();
    int n_y        = xsec->GetNbinsY();
    int n_lin_bins = linear->GetNbinsX();
    
    for ( int i = 0; i < n_lin_bins; ++i ){
        
        int x = GetXBin(i, n_lin_bins, n_x);
        int y = GetYBin(i, n_lin_bins, n_x);

        xsec->SetBinError(x+1,y+1, TMath::Sqrt(covariance[i][i]));

    } 
}

void SetFluxStatsErrors( TH2D *unfold, TH2D *xsec, double norm ){

    int n_x = xsec->GetNbinsX();
    int n_y = xsec->GetNbinsY();
    
    for ( int i = 0; i < n_y; ++i ){
        
        for ( int j = 0; j < n_x; ++j ){

            double flux_error_squared = xsec->GetBinError( j, i );
            double statistical_error  = unfold->GetBinError( j, i );

            double total_error = flux_error_squared + ( statistical_error * statistical_error * norm * norm );

            xsec->SetBinError( j, i, TMath::Sqrt( total_error ) );

        }
    }
}

// Get errors using linearisation
void Set2DErrors( TH2D *xsec, TH2D *unfold, TH1D *linear, double flux_err, const char file[1024], double norm ){

    // Linearisation function
    GetLinearisatedHistogram( xsec, linear);

    // Get covariance matrix
    TMatrixDSym covariance = GetCovariance( linear, flux_err, file );

    // Set errors on the linear histogram
    SetLinearErrors( covariance, linear, file );

    // Set errors on the 2D histogram from the covariance matrix
    SetDDXSecErrors( linear, xsec, covariance );

    // Set the total error contributions from flux and statistics
    SetFluxStatsErrors( unfold, xsec, norm );
}

#endif
