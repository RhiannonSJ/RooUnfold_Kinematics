/*
 *--------------------------------------------------------------
 *
 * Author: Rhiannon Jones
 * Date  : June 2017
 *
 *
*/
#ifndef COV_FUNC_H
#define COV_FUNC_H

namespace xsec{

// -------------------------------------------------------------------------


int Map( unsigned int i, unsigned int j, unsigned int x_size, unsigned int y_size );


int GetXBin( unsigned int i, unsigned int size, unsigned int x_size );


int GetYBin( unsigned int i, unsigned int size, unsigned int x_size );


// -------------------------------------------------------------------------
} // namespace
#endif
