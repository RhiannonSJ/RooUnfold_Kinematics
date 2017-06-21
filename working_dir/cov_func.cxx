/*
 * The functions used throughout the plotting of MC info
 *
 *--------------------------------------------------------------
 *
 * Author: Rhiannon Jones
 * Date  : June 2017
 *
*/

#include "cov_func.h"

using namespace std; 


// -------------------------------------------------------------------------


int Map( unsigned int i, unsigned int j, unsigned int x_size, unsigned int y_size ){

    if ( i >= x_size ) return -1;  
    if ( j >= y_size ) return -1;  

    return i + j * x_size;

}

int GetXBin( unsigned int i, unsigned int size, unsigned int x_size ){

    if ( i >= size ) return -1;

    return i % x_size; 

}


int GetYBin( unsigned int i, unsigned int size, unsigned int x_size ){

    if ( i >= size ) return -1;

    return i / x_size; 

}


// -------------------------------------------------------------------------

