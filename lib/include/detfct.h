/*
 * Library to compute determinant for GF2 matrices
 *
 * After importing the library please redefine N with the used size
 */

#ifndef DETFCT_H
#define DETFCT_H

#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <bitset>

#define DETFCT_MAX_SIZE 512

namespace GF2_Utils {
	/*
	   Triangularise M sur place donc M originale detruite.
	   Determinant est egale aux produits des elements de la diagonale de la triangularisation.
	   Des qu'un element nul est rencontre sur la diagonale alors retourne zero.
	   */
	bool det_b(std::vector<std::bitset<DETFCT_MAX_SIZE>> M, long msize);
}

#endif
