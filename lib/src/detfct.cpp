/*
 * Library to compute determinant for GF2 matrices
 *
 * After importing the library please redefine N with the used size
 */

#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <bitset>

#include "detfct.h"

namespace GF2_Utils {
	long posit_table[DETFCT_MAX_SIZE];//On conserve les position des indexes de rangee sur lesquelles operees.

	long count_effective_length_posit_table;//Pour calculer le nombre d'entrees non nulles pour une colonne.

	/*
	   Triangularise M sur place donc M originale detruite.
	   Determinant est egale aux produits des elements de la diagonale de la triangularisation.
	   Des qu'un element nul est rencontre sur la diagonale alors retourne zero.
	   */

	bool det_b(std::vector<std::bitset<DETFCT_MAX_SIZE>> M, long msize)
	{
		if (msize > DETFCT_MAX_SIZE) 
			throw std::runtime_error("Matrix is too big to compute determinant");

		for(long c = 0; c < msize; c++)
		{
			count_effective_length_posit_table=0;
			for(long r = c; r < msize; r++)
			{
				if(M[r].test(c))//fouille les elements non nuls en-dessous de M[c][c]
				{
					posit_table[count_effective_length_posit_table]=r;//enregistre leurs positions
					count_effective_length_posit_table++;//conteur nombre elts non nulle
				}
			}
			if(count_effective_length_posit_table==0)//rien a faire, assurement M[c][c]=0
			{
				return 0;
			}
			else
			{
				for(long row_index = 1; row_index < count_effective_length_posit_table ; row_index++)//boucle n'est pas effective si pos.size() == 1 et seulement le swap des rangee doit etre fait
				{
					M[posit_table[row_index]] ^= M[posit_table[0]];

				}
				std::swap(M[c],M[posit_table[0]]);//Apres le swap, maintenant element courant de la diagonale de la triangularisation est non nulle
			}
		}

		//Jusqu'a ici, tous les elements de la diagonale de la matrice triangularisee sont non nuls. Il ne reste qu'a verifier l'entree dans le coin inferieur droit c'est-a-dire simplement retourne sa valeur.
		return M[msize-1][msize-1];

	}
}
