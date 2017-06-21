/*
 *--------------------------------------------------------------
 *
 * Author: Rhiannon Jones
 * Date  : June 2017
 *
 *
*/
#ifndef COV_FUNC
#define COV_FUNC

#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TObjArray.h"

using namespace std;

// -------------------------------------------------------------------------


int Map( unsigned int i, unsigned int j, unsigned int x_size, unsigned int y_size );


int GetXBin( unsigned int i, unsigned int size, unsigned int x_size );


int GetYBin( unsigned int i, unsigned int size, unsigned int x_size );


// -------------------------------------------------------------------------

#endif
