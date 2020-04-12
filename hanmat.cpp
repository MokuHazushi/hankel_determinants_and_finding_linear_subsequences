/***** BEGINNING OF DISCLAIMER *****
       
       This code is provided without warranty of any kind, either express or implied. Use at your own risk.
       
       The use of this is done at your own discretion and risk and with agreement that you will be solely responsible for any damage to your computer system or loss of data that results from the use of this code. You are solely responsible for adequate protection and backup of the data and equipment used in connection with this code, and we will not be liable for any damages that you may suffer in connection with using, modifying or distributing any of this code. No advice or information, whether oral or written, obtained by you from us or from this website shall create any warranty for this code.
       
       We make makes no warranty that:
       
       1) the code will meet your requirements
       2) the code will be uninterrupted, timely, secure or error-free
       3) the results that may be obtained from the use of the software will be effective, accurate or reliable
       4) the quality of the code will meet your expectations
       5) any errors in the code obtained from us will be corrected. 
       
       The code and its documentation made available on this website:
       
       1) could include technical or other mistakes, inaccuracies or typographical errors. We may make changes to the code or documentation made available on its web site at any time without prior-notice.
       2) may be out of date, and we make no commitment to update such materials. 
       
       We assume no responsibility for errors or omissions in the code or documentation available from its web site.
       
       In no event shall we be liable to you or any third parties for any special, punitive, incidental, indirect or consequential damages of any kind, or any damages whatsoever, including, without limitation, those resulting from loss of use, data or profits, and on any theory of liability, arising out of or in connection with the use of this code. 
       
       ***** END OF DISCLAIMER *****/

/*
  
  AUTHOR: Claude Gravel
  
  VERSION: 1 (March 2020)
  
  TECHNICAL REFERENCE: (arxiv) "Finding linearly generated subsequences" by Claude Gravel and Daniel Panario

  NOTES:
  
  This code relies on NTL library which depends on GMPL.
  For information about NTL (Number Theory Library) see: https://www.shoup.net/ntl/
  For information about GMPL (GNU Multiple Precision Arithmetic Library) see: https://gmplib.org/
  This code is single thread. The pthread option is necessary however for NTL.

  COMPILE AND LINK COMMAND: "g++ -Wall -Wpedantic -g -O5 -I<NTL directory>/include -I<GMPL directory>/include hanmat.cpp -o hanmat -std=c++11 -L<NTL library>/lib -lntl -L<GMPL libray>/lib -lgmp -lm -lpthread"

  RUN COMMAND EXAMPLE: "./hanmat" or "./hanmat > results.txt" or "./hanmat > results.txt &2>&1 &"
  
*/


#include <iomanip>
#include <iostream>
#include <fstream>
#include <random>
#include <time.h>
#include <strings.h>
#include <sys/resource.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/GF2.h>
#include <NTL/vec_vec_GF2.h>
#include <NTL/mat_GF2.h>

/**** ***** ***** ***** *****/

const int n=31;//length data vector
const int max_dim=(n/2)*( n%2 == 0) + ((n+1)/2)*( n%2 == 1);
const double ratio_size=1.0/2.0;
const int len_C_gen = (int)floor((double)n*ratio_size);//length of generating vector

/*
  Possible locations of the hidden linear subsequence.
  Case [1] is representative of a hard case.
  Case [2] is representative of an easy case.
*/

//const int left_index = (int)(floor(7*len_C_gen));//where the linear subsequence begins, case [1]
//const int right_index = (int)(floor(9*len_C_gen));//where the linear subsequence ends, case [1]
const int left_index = len_C_gen;//where the linear subsequence begins, case [2]
const int right_index = n;//where the linear subsequence ends, case [2]

/**** ***** ***** ***** *****/

const int max_len_side_grid=7;//Research required to increase value here...

/**** ***** ***** ***** *****/

const bool RUN_NAIVE_FLAG = true;
//If RUN_NAIVE_FLAG is false, then skip naive (thus resulting in nan, warning of mismatch between naive and new methods unless the last part is modified appropriately.)
//RUN_NAIVE_FLAG can be false when the new method is wanted solely and no need to compare methods.

const bool PRINT_COUNTS_FOR_DIRECT_COMPUTATION = false;
// PRINT_COUNTS_FOR_DIRECT_COMPUTATION = true prints counts of the determinants that we were not able to catch with our theory and that required explicit computations.

const bool PRINT_TRIA_TABLE=true;
// To print triangular table(s) out.

/**** ***** ***** ***** *****/

NTL::Vec<NTL::GF2> C_gen;//generation vector
NTL::Vec<NTL::GF2> V0;//data vector to be filled with a linear subsequence
NTL::Mat<NTL::GF2> M;
NTL::Mat<bool> flags_M;
NTL::Mat<NTL::GF2> new_M;

/**** ***** ***** ***** *****/

int ct_naive;
int ct_sqfill;
int ct_nsew;
int *ct_det_tab;
int *ct_grid;
int *ct_size_explicit;

/**** ***** ***** ***** *****/

//NTL::Mat<NTL::GF2> *AspC = new NTL::Mat<NTL::GF2> [max_dim+1];
NTL::Mat<NTL::GF2> *AspTL = new NTL::Mat<NTL::GF2> [max_dim+1];
//NTL::Mat<NTL::GF2> *AspTR = new NTL::Mat<NTL::GF2> [max_dim+1];
//NTL::Mat<NTL::GF2> *AspBL = new NTL::Mat<NTL::GF2> [max_dim+1];
//No need of the list of Bottom Right matrices, i.e. AspBR, since it contains the unknown d_{i,j}.

/**** ***** ***** ***** *****/

double get_ms_res_time(void)//This is the most accurate measure of time I could think. Anything better?
{
  
  double current_time;
  
  struct rusage ruse;
  getrusage(RUSAGE_THREAD, &ruse);
  current_time = ((double) ruse.ru_utime.tv_sec * 1000.0) + (((double) ruse.ru_utime.tv_usec) / 1000.0); //resultat en millisecondes
  
  return current_time;
}

/**** ***** ***** ***** *****/

/*
  DESCRIPTION
  print_ntl_gf2_mat prints on the standard output M with max_dim rows, and n columns.
  out_width1 is to control the number of caracters for the index of the row.
*/
void print_ntl_gf2_mat(NTL::Mat<NTL::GF2> M,int max_dim,int n,int out_width1)
{
  
  std::cout << "\n";
  
  std::cout << std::setw(out_width1) << 0 << " : ";
  for(int c0=0;c0 < n;c0++)
    {
      std::cout << M[0][c0] << " ";
    }
  std::cout << "\n";
  for(int r0=1;r0<max_dim+1;r0++)
    {
      std::cout << std::setw(out_width1) << r0 << " : ";
      
      for(int c0=0;c0 < r0-1;c0++)
	{
	  std::cout << "  ";
	}
      
      for(int c0=r0-1;c0 <= n-r0;c0++)
	{
	  std::cout << M[r0][c0] << " ";
	}
      std::cout << "\n";
    }
}

/**** ***** ***** ***** *****/

/*
  DESCRIPTION:
  Verifies that M1 and M2 (with the same size) are *not* the same.
  So it returns false is the tables are the same and true otherwise.
*/
bool chk_triangular_tables_not_the_same(NTL::Mat<NTL::GF2> M1, NTL::Mat<NTL::GF2> M2,int max_dim,int n)
{
  /*
    return false if tables are same
    return true if tables not the same
  */
  
  bool chk=false;
  for(int c0=0;c0 < n;c0++)
    {
      chk = (M1[0][c0]!=M2[0][c0]);
      if(chk){return chk;}
    }
  
  for(int r0=1;r0<max_dim+1;r0++)
    {
      
      for(int c0=r0-1;c0 <= n-r0;c0++)
	{
	  chk = (M1[r0][c0]!=M2[r0][c0]);
	  if(chk){return chk;}
	}
      
    }
  
  return chk;
}

/**** ***** ***** ***** *****/

/*
  DESCRIPTION:
  Fill the triangular table with zeros when the top part of a square is detected.
  By top part, for position indices i_0, i_1 such that j-1 < i_0<i_1 <n-j, and j>=2,
  If a run of non-zero elements is detected from i_0-1 to to i_1+1 on level j-2,
  and a run of zero elements is detected from i_0 to i_1 on level j-1, then fill the table
  below accordingly.
*/

void s_f(int j)
{
  /*
    Use knowledge of table with rows 0-th to (j-1)-th included to build square of zeros
  
    StartAtLeft & EndInMiddle
    StartAtLeft & EndAtRight
    StartInMiddle & EndInMiddle
    StartInMiddle & EndAtRight
  */
  
  /*
    (...)                                   ...                 ...                                         (...)
    
          (j-2,j-3) |*| (j-2,j-2) (j-2,j-1) ....................... (j-2, n-j) (j-2,n-j+1) |*| (j,n-j+2) < n-(j-2)+1=n-j+3
    
                        (j-1,j-2) (j-1,j-1) ....................... (j-1, n-j) (j-1,n-j+1)               < n-(j-1)+1=n-j+2
    
    --------------------------->  (j  ,j-1) ....................... (j  , n-j)                           < n-(j  )+1=n-j+1

  */
  
  int i = j-1;
 
  while(i<=n-j)//using knowledge of (j-2)th, (j-1)th rows, stop when right position = n-j+1
    {
      int lp = 0;//index of the left non-zero element
      int rp = 0;//index of the right non-zero element
     
      bool chkL1 = (i==j-1) && (new_M[j-2][i-1]!=0) && (new_M[j-1][i-1]==0) && (new_M[j-2][i]!=0) && (new_M[j-1][i]==0);

      bool chkL2 = (new_M[j-2][i-1]!=0) && (new_M[j-2][i]!=0) && (new_M[j-1][i-1]!=0) && (new_M[j-1][i]==0);
    
      if(chkL1||chkL2)
	{
	  lp = i;
	  int t1 = lp+1;
	  
	  while(t1 <= n-j)
	    {
	      bool chkR1 = (new_M[j-2][t1-1]!=0) && (new_M[j-2][t1]!=0) && (new_M[j-1][t1-1]==0) && (new_M[j-1][t1]!=0);
	      bool chkR2 = (t1==n-j) && (new_M[j-2][t1]!=0) && (new_M[j-1][t1]==0) && (new_M[j-2][t1+1]!=0) && (new_M[j-1][t1+1]==0);
	      
	      if(chkR1||chkR2)
		{
		  rp = t1-1;
		  break;
		}
	      else
		{
		  t1 = t1 + 1;
		}
	    }
	  if(rp>lp)
	    {
	      
	      int stop_dim;
	      if(max_dim<=(j+rp-lp-1))
		{
		  stop_dim = max_dim;
		}
	      else
		{
		  stop_dim = (j+rp-lp-1);
		}
	     	      
	      for(int x_j=j;x_j<stop_dim;x_j++)//inside square excluding boundaries, square of zeros top row = j-1 which has been completed, thus begin at j
		{
		  int x_i_l,x_i_r;//extreme possible values of i as the filling process goes from top part to bottom part

		  if(lp>=x_j-1)
		    {
		      x_i_l = lp;
		    }
		  else
		    {
		      x_i_l = x_j-1;
		    }
		  if(rp<=n-x_j)
		    {
		      x_i_r = rp;
		    }
		  else
		    {
		      x_i_r = n-x_j;
		    }
		  for(int x_i=x_i_l;x_i<=x_i_r;x_i++)
		    {
		      new_M[x_j][x_i]=0;
		      flags_M[x_j][x_i]=1;
		      
		      ct_sqfill+=1;
		      
		    }
		}
	      i = rp;
	    }
	  else
	    {
	      i = i + 1;
	    }
	}
      else
	{
	  i = i+1;
	}
    }      
}

/**** ***** ***** ***** *****/

NTL::GF2 solve_eq_for_lower_corner(int j, int i, int q)
{
 
  /*
    q is the minimum length of the grid required to solve for d_{i,j}.
    The Main matrix has size therefore (q+1) X (q+1).
    The Top/Bottom & Left/Right matrices used to solve for d_{i,j} has size q X q.
  */
  
  NTL::GF2 ret = (NTL::GF2)0;
  
  NTL::Mat<NTL::GF2> T;
  T.SetDims(q,q);//temporary matrix of size q X q for a minor obtained from the (q+1) X (q+1) Main matrix
  
  for(int g=0;g<q;g++)//index mineure
    {
      for(int h=0;h<q+1;h++)//index rangee de la mineure
	{
	  if(g<h)
	    {
	      for(int k=0;k<q;k++)//index colonne de la mineure
		{
		  T[h-1][k]=new_M[j-2*q+h+k][i-h+k];
		}
	    }
	  else if (g>h)
	    {
	      for(int k=0;k<q;k++)
		{
		  T[h][k]  =new_M[j-2*q+h+k][i-h+k];
		} 
	    }
	  //else it is not part of the expansion of the determinant
	}
      ret = ret+new_M[j-q+g][i+q-g]*NTL::determinant(T);
    } 
  return ret;
}

/**** ***** ***** ***** *****/

bool chk_conds_for_solvability(int j, int i, int *effective_length)
{
  bool ret=false;
  
  for(int w1=2;w1<=max_len_side_grid;w1++)//w1 is the length of the side which is growing so that the main matrix has size (w1+1)x(w1+1)
    {
      for(int rx=0;rx<w1;rx++)
	{
	  for(int cx=0;cx<w1;cx++)
	    {
	      //The size for any of the Top/Bottom & Left/Right matrices is w1 x w1 since the main matrix has size (w1+1)x(w1+1)

	      AspTL[w1][rx][cx]=new_M[j-2*w1+rx+cx][i-rx+cx];
	      //AspTR[w1][rx][cx]=new_M[j-2*w1+(rx+1)+cx][i-rx+(cx+1)];
	      //AspBL[w1][rx][cx]=new_M[j-2*w1+(rx+1)+cx][i-rx+(cx-1)];
	      
	    }
	  /*
	  for(int rx=0;rx<w1-1;rx++)
	    {
	      for(int cx=0;cx<w1-1;cx++)
		{
		  // The size of the sub centered matrix is (w1-1) X (w1-1)
		  AspC[w1-1][rx][cx]=new_M[j-2*(w1-1)+rx+cx][i-rx+cx];
		}
	    }
	  */
	}

      //if ( (NTL::determinant(AspC[w1-1])==0) && (NTL::determinant(AspTL[w1])==0) && ( (NTL::determinant(AspTR[w1])==0) || (NTL::determinant(AspBL[w1])==0) ) )
    
      //if ( (NTL::determinant(AspC[w1-1])==0) && (NTL::determinant(AspTL[w1])!=0) && (new_M[j-w1][i-1]*new_M[j-w1][i+1]+new_M[j-w1-1][i]*new_M[j-w1+1][i]==0))//this is right for side length of the grid <=7;
      
      //if ( (NTL::determinant(AspTL[w1])!=0) && (new_M[j-w1][i-1]*new_M[j-w1][i+1]+new_M[j-w1-1][i]*new_M[j-w1+1][i]==0) )/******/
      if ( (NTL::determinant(AspTL[w1])!=0) && (new_M[j-w1][i]==0) )/******/
	{
	  ret=true;
	  *effective_length=w1;//length of a side of the grid so the main square matrix of size (w1+1) X (w1+1)	  
	  break;
	}
    }
  return ret;
}

/**** ***** ***** ***** *****/

void d_c_s(int j)
{
  int i = j-1;

  while(i<n-j+1)
    {
      int effective_len;
      if(!flags_M[j][i])
	{
	  /************/
	  
	  if( (j>=2) && (new_M[j-2][i]==(NTL::GF2)1) )//North-South-East-West
	    {
	      
	      new_M[j][i]=new_M[j-1][i-1]*new_M[j-1][i+1]+new_M[j-1][i];
	      flags_M[j][i]=true;
	      ct_nsew+=1;
	      
	    }
	  else if ( (j>=2*max_len_side_grid) && chk_conds_for_solvability(j,i,&effective_len) )
	    {
	      new_M[j][i] = solve_eq_for_lower_corner(j,i,effective_len);
	      flags_M[j][i] = true;
	      ct_grid[effective_len]+=1;
	    }
	  else
	    {
	      int t_x=1;
	      
	      while(t_x<=j)
		{
		  if(new_M[j-t_x][i]==0)//start check from j-1 to 0 (new_M[0][i] == 1 by definition) ( level of sequence )
		    {
		      t_x+=1;
		    }
		  else
		    {
		      break;//when here, new_M[j-t_x][i]!=0 AND new_M[j-1, j-2, ..., j-(t_x-1)][i]==0
		    }
		}
	      //Note that t_x can never be 2 actually for it is handled by the North-South-East-West 
	      if(t_x==1)//explicit computation
		{
		  
		  NTL::Mat<NTL::GF2> tmp;
		  tmp.SetDims(j,j);
		  
		  for(int r = 0;r<j ; r++)
		    {
		      for(int c = 0; c<j ; c++)
			{
			  tmp[c][r] = new_M[t_x][i-r+c];//==V0[i-j+r+c+1];
			}
		    }
		  ct_det_tab[0]+=1;
		  new_M[j][i] = NTL::determinant(tmp);
		  flags_M[j][i]=true;
		  ct_size_explicit[j]+=1;
		}
	      else//computation by cross identities
		{
		  NTL::Mat<NTL::GF2> tmp;
		  tmp.SetDims(t_x,t_x);
		 
		  for(int r=0;r<t_x;r++)
		    {
		      for(int c=0;c<t_x;c++)
			{
			  tmp[r][c]=new_M[j-(t_x-1)][i+c-r];
			}
		    }
		  new_M[j][i] = NTL::determinant(tmp);
		  flags_M[j][i]=true;
		  ct_det_tab[t_x-1]+=1;//we use info at level t_x-1 above j, new_M[j-t_x][i]!=0 AND new_M[j-1, j-2, ..., j-(t_x-1)][i]==0, we therefore increase counter indexed by t_x-1
		}
	    }
	}
      i=i+1;
    }
}

/**** ***** ***** ***** *****/

int main(void)
{
  std::random_device r_dev;
  std::mt19937_64 mt(r_dev());
  std::bernoulli_distribution dist(0.5);
  //std::uniform_real_distribution<double> dist(0, 1);
  
  std::cout.precision(6);
  std::cout.setf( std::ios::fixed, std::ios::floatfield );

  if(right_index>n){std::cerr << "\n\nOut of bound\n\n"; exit(-1);}

  C_gen.SetLength(len_C_gen);
  
  for(unsigned int i1=0;i1<C_gen.length();i1++)
    {
      C_gen[i1]=(NTL::GF2)(dist(mt));
    }
  
  C_gen[C_gen.length()-1]=(NTL::GF2)1;//make sure the generating vector is nonzero
  
  V0.SetLength(n);
  for(int x1=0;x1<n;x1++)
    {
      V0[x1]=(NTL::GF2)(dist(mt));
    }
  
  for(int x1=left_index ;  x1 < right_index; x1++)
    {
      NTL::GF2 tmpsum=(NTL::GF2)0;
      for(int y1 = 0; y1 < C_gen.length()-1; y1++)
	{
	  tmpsum = tmpsum + (C_gen[ C_gen.length()-2-y1 ]*V0[ x1 -1 - y1 ]);//scalar product computation --- linear subsequence
	}
      V0[x1]=(NTL::GF2)(tmpsum);//scalar product value
    }
  
  /**************************************************************/

  std::cout << "\n\n";
  std::cout << "Original sequence (is identical to row 1 of triangular table), length = " << V0.length() << ", sequence:\n     ";
  for(unsigned int l1 = 0; l1<V0.length(); l1++)
    {
      std::cout << V0[l1] << " ";
    }
  
  std::cout << "\n\n";
  std::cout << "Generating vector, length = " << C_gen.length() << ", vector:\n     ";
  for(unsigned int l1 = 0; l1<C_gen.length(); l1++)
    {
      std::cout << C_gen[l1] << " ";
    }
  
  std::cout << "\n";
  
  /**************************************************************/
 
  M.SetDims(max_dim+1,n);//number of rows = max_dim + 1 (an extra one required for initializing the top part)

  for(int r1=0;r1<max_dim+1;r1++)
    {
      for(int c1=0;c1<n;c1++)
	{
	  M[r1][c1]=(NTL::GF2)0;
	}
    }

  double tnaive0 = get_ms_res_time();
  for(int c1=0;c1<n;c1++)
    {
      M[0][c1]=(NTL::GF2)1;
    }
  for(int c1=0;c1<n;c1++)
    {
      M[1][c1]=V0[c1];
    }
  
  ct_naive=0;
 
  if(RUN_NAIVE_FLAG)
    {
      for(int j = 2; j<max_dim+1; j++)
	{
	  NTL::Mat<NTL::GF2> tmp;
	  tmp.SetDims(j,j);
	  
	  for(int r1=0;r1<j;r1++)
	    {
	      for(int c1=0;c1<j;c1++)
		{
		  tmp[r1][c1]=0;
		}
	    }
	  
	  for(int i = j-1; i<n-j+1; i++)
	    {
	      for( int r = 0; r<j ; r++)
		{
		  for( int c = 0; c<j ; c++)
		    {
		      tmp[r][c] = V0[(i+1-j+r+c)];
		    }
		}
	      ct_naive+=1;
	      M[j][i] = NTL::determinant(tmp);
	    }
	}
    }
  
  double tnaive1 = get_ms_res_time();
  double time_diff_naive = tnaive1-tnaive0;
  
  /**************************************************************/
  
  flags_M.SetDims(max_dim+1,n);
  
  new_M.SetDims(max_dim+1,n);

  AspTL[0].SetDims(1,1);
  //AspC[0].SetDims(1,1);
  //AspTR[0].SetDims(1,1);
  //AspBL[0].SetDims(1,1);

  AspTL[0][0][0]=0;//init but never needed
  //AspC[0][0][0]=0;//init but never needed
  //AspTR[0][0][0]=0;//init but never needed
  //AspBL[0][0][0]=0;//init but never needed
  
  for(int t=1;t<max_dim+1;t++)
    {
      AspTL[t].SetDims(t,t);
      //AspC[t].SetDims(t,t);
      //AspTR[t].SetDims(t,t);
      //AspBL[t].SetDims(t,t);
    }
  
  ct_sqfill=0;
  ct_nsew=0;
  ct_det_tab = new int [n];
  bzero(ct_det_tab,n*sizeof(int));
  ct_grid = new int [n];
  bzero(ct_grid,n*sizeof(int));
  ct_size_explicit = new int [max_dim];
  bzero(ct_size_explicit,max_dim*sizeof(int));
  
  double tnew0 = get_ms_res_time();
  
  for(int i1=0;i1<max_dim+1;i1++)
    {
      for(int i2=0;i2<n;i2++)
	{
	  new_M[i1][i2]=(NTL::GF2)0;
	  flags_M[i1][i2]=false;
	}
    }
 
  for(int i2=0;i2<n;i2++)//j=0,1
    {
      new_M[0][i2]=(NTL::GF2)1;
      flags_M[0][i2]=true;//flag as completed
      new_M[1][i2]=(NTL::GF2)V0[i2];
      flags_M[1][i2]=true;
    }

  for(int j = 2; j<max_dim+1; j++)//levels < j are completed
    {
      s_f(j);//make squares using full knowledge for previous rows 0, 1,...,j-1 inclusively
      d_c_s(j);//fill j-th row for determinants not found yet
    }
  
  double tnew1 = get_ms_res_time();
  double time_diff_new = tnew1-tnew0;
    
  if(chk_triangular_tables_not_the_same(M,new_M,max_dim,n))
    {
      std::cout<<"\n\n*****There are mismatches between naive and new dynamic methods. Exit now.*****\n";
      exit(-1);
    }

  int number_of_entries = ((n%2)==0)*((n/2)-1)*(n/2)+((n%2)==1)*(n/2)*(n-1)/2;
  int out_width0 = (int)floor(1.0+(log(n)/log(10)));
  int out_width1 = (int)floor(1.0+(log(number_of_entries)/log(10)));
  
  std::cout << "\n[A." << std::setw(out_width0) <<0<< "]  length of a sequence              = " << n;
  std::cout << "\n[A." << std::setw(out_width0) <<1<< "]  length of generating vector       = " << len_C_gen;
  std::cout << "\n[A." << std::setw(out_width0) <<2<< "]  left index linear subsequence     = " << left_index;
  std::cout << "\n[A." << std::setw(out_width0) <<3<< "]  right index linear subsequene     = " << right_index;
  std::cout << "\n[A." << std::setw(out_width0) <<4<< "]  number of entries of a table      = " << number_of_entries << " (analytic count)\n";

  std::cout << "\n[B." << std::setw(out_width0) <<0<< "]  time new                          = " << time_diff_new << " ms";
  std::cout << "\n[B." << std::setw(out_width0) <<1<< "]  time naive                        = " << time_diff_naive << " ms";
  std::cout << "\n[B." << std::setw(out_width0) <<2<< "]  ratio B.0/B.1                     = " << time_diff_new/time_diff_naive << "\n";

  long sum_ct_det_tab=ct_det_tab[0];//direct eval
  long sum_ct_grid=0;
 
  for(int q1=1;q1<n;q1++)
    {
      if(ct_det_tab[q1]!=0)
	{
	  sum_ct_det_tab+=ct_det_tab[q1];
	  std::cout << "\n[C." << std::setw(out_width0) << q1 << "]  counts cross level " << std::setw(out_width1) << q1 << "          = " << ct_det_tab[q1];
	}
      
    }

  std::cout << "\n";
  for(int q1=0;q1<n;q1++)
    {
      if(ct_grid[q1]!=0)
	{
	  sum_ct_grid+=ct_grid[q1];
	  std::cout << "\n[D." << std::setw(out_width0) << q1 << "]  counts inter grid "<<q1<<"by"<<q1 << std::setw(out_width1) << "            = " << ct_grid[q1];
	}
      
    }

  std::cout << "\n\n[E." << std::setw(out_width0) <<0<< "]  counts North-South-East-West      = " << ct_nsew ;
  std::cout << "\n\n[F." << std::setw(out_width0) <<0<< "]  counts square filling             = " << ct_sqfill;
  std::cout << "\n\n[G." << std::setw(out_width0) <<0<< "]  counts direct evaluations         = " << ct_det_tab[0];

  std::cout << "\n\ntotal sum of counts                        = " << ct_nsew+sum_ct_det_tab+ct_sqfill+sum_ct_grid << "\n\n";

 
  if(PRINT_TRIA_TABLE)
    {
      
      int out_width2 = (int)floor(1.0+(log(max_dim)/log(10)));//output number of max character for row indices --- alignment outputs
      print_ntl_gf2_mat(new_M,max_dim,n,out_width2);
      //print_ntl_gf2_mat(M,max_dim,n,out_width2);
      
    }
 
  if(PRINT_COUNTS_FOR_DIRECT_COMPUTATION)
    {
      
      std::cout << "max possible size = " << max_dim << "\n";
      long sum_ct_size=0;

      for(int z1=0;z1<max_dim;z1++)
	{
	  if(ct_size_explicit[z1]!=0)
	    {
	      sum_ct_size+=ct_size_explicit[z1];
	      std::cout << z1 << " --> " << ct_size_explicit[z1] << "\n";
	    }
	}
      if(sum_ct_size!=ct_det_tab[0])
	{
	  std::cout << "This can't happen. Something wrong.\n";
	}
    }

  delete [] ct_grid;
  delete [] ct_det_tab;
  delete [] ct_size_explicit;
  
  return EXIT_SUCCESS;
}
