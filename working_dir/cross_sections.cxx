#include <iostream>
using std::cout;
using std::endl;

#include "roo_unfold.h"

#include "../src/RooUnfoldResponse.h"
#include "../src/RooUnfoldBayes.h"


void GetTruth( TTree *tree, int n_protons, std::vector<double> & truth_T, std::vector<double> & truth_cos, std::vector<bool> & truth_detectable, std::vector<double> & impur_T, std::vector<double> & imput_cos );

// void Smear( const std::vector<double> & truth_T, const std::vector<double> & truth_cos, std::vector<double> & smear_T, std::vector<double> & smear_cos ); 

void cross_sections() {
    //==============================================================================
    // Reading in the event root files 
    //==============================================================================
    TFile f_test("/hepstore/rjones/Exercises/Flavours/Default+MEC/sbnd/1M/gntp.10000.gst.root");
    if(f_test.IsZombie()){
        std::cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "======================== Default + MEC event file open ========================" << endl;
    }
 
    TFile f_train("/hepstore/rjones/Exercises/Flavours/Default/sbnd/1M/gntp.10000.gst.root");
    if(f_train.IsZombie()){
        std::cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "=========================== Default event file open ===========================" << endl;
    }

    
    //==============================================================================
    // Get the truth vectors 
    //==============================================================================
    std::vector<double> truth_T_train; 
    std::vector<double> truth_cos_train; 
    std::vector<bool>   truth_detectable_train;
    std::vector<double> impur_T_train; 
    std::vector<double> impur_cos_train; 

    GetTruth( truth_T_train, truth_cos_train, truth_detectable_train, impur_T_train, impur_cos_train );

    std::vector<double> truth_T_test; 
    std::vector<double> truth_cos_test; 
    std::vector<bool>   truth_detectable_test;
    std::vector<double> impur_T_test; 
    std::vector<double> impur_cos_test; 

    GetTruth( truth_T_test, truth_cos_test, truth_detectable_test, impur_T_test, impur_cos_test );

}
