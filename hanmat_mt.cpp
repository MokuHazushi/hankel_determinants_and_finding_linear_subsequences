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
  
  AUTHORS: Claude Gravel and Bastien Rigault
  
  VERSION: 2 (April 2020)
  
  TECHNICAL REFERENCE: (arxiv) "Finding linearly generated subsequences" by Claude Gravel and Daniel Panario

  NOTES:
  
  This code relies on NTL library which depends on GMPL.
  For information about NTL (Number Theory Library) see: https://www.shoup.net/ntl/
  For information about GMPL (GNU Multiple Precision Arithmetic Library) see: https://gmplib.org/

 
*/

#include <iomanip>
#include <iostream>
#include <fstream>
#include <random>
#include <ctime>
#include <chrono>
#include <sys/resource.h>
#include <thread>
#include <mutex>

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/GF2.h>
#include <NTL/vec_vec_GF2.h>
#include <NTL/mat_GF2.h>

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

const int n=512;//length data vector
const int max_dim=(n/2)*( n%2 == 0) + ((n+1)/2)*( n%2 == 1);
const double ratio_size=1.0/16.0;
const int len_C_gen = (int)round((double)n*ratio_size);//length of generating vector
const int max_number_of_thread = 8; // Maximum number of thread that can be created for each j-th row

/*
  Possible locations of the hidden linear subsequence.
  Case [1] is representative of a hard case.
  Case [2] is representative of an easy case.
*/

const int left_index = (int)(floor(7*len_C_gen));//where the linear subsequence begins, case [1]
const int right_index = (int)(floor(9*len_C_gen));//where the linear subsequence ends, case [1]
//const int left_index = len_C_gen;//where the linear subsequence begins, case [2]
//const int right_index = n;//where the linear subsequence ends, case [2]

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

const int max_len_side_grid=7;//Research required to increase value here...

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

const bool RUN_NAIVE_FLAG = true;
//If RUN_NAIVE_FLAG is false, then skip naive (thus resulting in nan, warning of mismatch between naive and new methods unless the last part is modified appropriately.)
//RUN_NAIVE_FLAG can be false when the new method is wanted solely and no need to compare methods.

const bool PRINT_COUNTS_FOR_DIRECT_COMPUTATION = false;
// PRINT_COUNTS_FOR_DIRECT_COMPUTATION = true prints counts of the determinants that we were not able to catch with our theory and that required explicit computations.

const bool PRINT_TRIA_TABLE=true;
// To print triangular table(s) out.

bool MULTI_THREAD_ON=false;

enum Thread_Allocation_Strategy {
  UNIFORM_ALLOCATION,
  OPTIMIZED_ALLOCATION
};

Thread_Allocation_Strategy THREAD_ALLOCATION_STRATEGY = UNIFORM_ALLOCATION;

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

NTL::Vec<NTL::GF2> C_gen;//generation vector
NTL::Vec<NTL::GF2> V0;//data vector to be filled with a linear subsequence

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

struct experiment {
  NTL::Mat<NTL::GF2> M;
  NTL::Mat<bool> flags_M; // Not used by trivial approach
  NTL::Mat<NTL::GF2> AspTL[max_dim+1]; // Not used by trivial approach
};

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

double time_triv_st,time_triv_mt,time_fast_st,time_fast_mt,time_fast_mt_opt;// for timing

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

/* PRINT FUNCTIONS */

void print_stats()
{
  int number_of_entries = ((n%2)==0)*((n/2)-1)*(n/2)+((n%2)==1)*(n/2)*(n-1)/2;
  int out_width0 = (int)floor(1.0+(log(n)/log(10)));
  
  std::cout << "\n[A." << std::setw(out_width0) <<0<< "]  length of a sequence                  = " << n;
  std::cout << "\n[A." << std::setw(out_width0) <<1<< "]  length of generating vector           = " << len_C_gen;
  std::cout << "\n[A." << std::setw(out_width0) <<2<< "]  left index linear subsequence         = " << left_index;
  std::cout << "\n[A." << std::setw(out_width0) <<3<< "]  right index linear subsequene         = " << right_index;
  std::cout << "\n[A." << std::setw(out_width0) <<4<< "]  number of entries of a table          = " << number_of_entries << " (analytic count)\n";

  
  std::cout << "\n[B." << std::setw(out_width0) <<0<< "]  Wall time trivial algo single thread                       = " << time_triv_st << " ms";
  std::cout << "\n[B." << std::setw(out_width0) <<1<< "]  Wall time trivial algo multi thread                        = " << time_triv_mt << " ms";
  std::cout << "\n[B." << std::setw(out_width0) <<2<< "]  Wall time faster algo single thread                        = " << time_fast_st << " ms";
  std::cout << "\n[B." << std::setw(out_width0) <<3<< "]  Wall time faster algo multi thread uniform allocation      = " << time_fast_mt << " ms";
  std::cout << "\n[B." << std::setw(out_width0) <<4<< "]  Wall time faster algo multi thread optimized allocation    = " << time_fast_mt_opt << " ms \n\n";
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

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
    std::cout << M[0][c0] << " ";

  std::cout << "\n";
  for(int r0=1;r0<max_dim+1;r0++)
  {
    std::cout << std::setw(out_width1) << r0 << " : ";
    for(int c0=0;c0 < r0-1;c0++)
      std::cout << "  ";
      
    for(int c0=r0-1;c0 <= n-r0;c0++)
      std::cout << M[r0][c0] << " ";
    std::cout << "\n";
  }
}

void how_many_differences(NTL::Mat<NTL::GF2> M,int max_dim,int n,int & nb_differences)
{
  
  
  nb_differences=0;
  for(int c0=0;c0 < n;c0++)
    {
      if(M[0][c0]==NTL::GF2(1))
	{
	  nb_differences++;
	}
    }

  for(int r0=1;r0<max_dim+1;r0++)
  {
      
    for(int c0=r0-1;c0 <= n-r0;c0++)
      {
	if(M[r0][c0]==NTL::GF2(1))
	  {
	    nb_differences++;
	  }
      }
  }
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

/* UTILS FUNCTIONS */
void init_generating_vector()
{
  std::random_device r_dev;
  std::mt19937_64 mt(r_dev());
  std::bernoulli_distribution dist(0.5);

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
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void init_experiment(struct experiment & experiment)
{
  experiment.M.SetDims(max_dim+1, n); //number of rows = max_dim + 1 (an extra one required for initializing the top part)
  experiment.flags_M.SetDims(max_dim+1, n);
  experiment.AspTL[0].SetDims(1, 1);
  
  experiment.AspTL[0][0][0]=0;//init but never needed
    
  for(int r1=0;r1<max_dim+1;r1++)
    {
      for(int c1=0;c1<n;c1++)
	{
	  experiment.M[r1][c1]=NTL::GF2(0);
	  experiment.flags_M[r1][c1]=false;
	}
    }
  
  for(int c1=0;c1<n;c1++)
    {
      experiment.M[0][c1]=NTL::GF2(1);
      experiment.M[1][c1]=V0[c1];
      experiment.flags_M[0][c1]=true;
      experiment.flags_M[1][c1]=true;
    }
  
  for(int t=1;t<max_dim+1;t++)
    {
      experiment.AspTL[t].SetDims(t,t);
    }
  
  for(int i2=0;i2<n;i2++)//j=0,1
    
    {
      experiment.M[0][i2]=NTL::GF2(1);
      experiment.flags_M[0][i2]=true;//flag as completed
      experiment.M[1][i2]=V0[i2];
      experiment.flags_M[1][i2]=true;
    }
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

/*
  DESCRIPTION:
  Verifies that M1 and M2 (with the same size) are *not* the same.
  So it returns false is the tables are the same and true otherwise.
*/
bool chk_triangular_tables_not_the_same(NTL::Mat<NTL::GF2> M1, NTL::Mat<NTL::GF2> M2, int max_dim, int n)
{
  /*
    return false if tables are same
    return true if tables not the same
  */
  
  bool chk=false;
  for(int c0=0;c0 < n;c0++)
    {
      chk = (M1[0][c0]!=M2[0][c0]);
      if(chk)
	{
	  return chk;
	}
    }
  
  for(int r0=1;r0<max_dim+1;r0++)
    {
      for(int c0=r0-1;c0 <= n-r0;c0++)
	{
	  chk = (M1[r0][c0]!=M2[r0][c0]);
	  if(chk)
	    {
	      return chk;
	    }
	}
    }
  return chk;
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

/* ALGORITHM FUNCTIONS */

/*
  DESCRIPTION
  partition_thread_load for a given j-th row, determines each thread load (start index and range)
  Returned value is a list of pair, pair.first corresponds to start index and pair.second the range.
  The size of the list might be smaller than the max_number_of_thread value in cases when the number of
  determinant to compute is too small.
  Returns an empty list if the computation of the j-th row should be handle by one (current) thread
*/
std::vector<std::pair<int,int>> partition_thread_load(int j)
{
  std::vector<std::pair<int, int>> result;
  int start_index = j-1;
  int end_index = n-j+1;
  //Est-ce que end_index est inclusif ? n = 5, donc positions sont 0,1,2,3,4
  //Si j = 2, alors positions sont 1,2,3
  //et start_index=1, end_index = 5-2+1=4, donc exclusif
  // nombre de determinant = 4 - 1 = 3, ok
  
  int nb_determinant = end_index - start_index;
  
  // Uniform share strategy
  int remaining_determinant = nb_determinant;
  int determinant_per_thread = nb_determinant / max_number_of_thread;
  int start = start_index;
  
  if (determinant_per_thread == 0)
    return result;
  
  while (remaining_determinant > 0)
    {
      if ((remaining_determinant - determinant_per_thread) > 0)
	{
	  result.push_back(std::pair<int, int>(start, determinant_per_thread));
	  remaining_determinant -= (determinant_per_thread+1);
	  start += determinant_per_thread + 1;
	}
      else
	{
	  result.push_back(std::pair<int, int>(start, remaining_determinant));
	  remaining_determinant = 0;
	}
    }
  
  return result;
}

/*
  DESCRIPTION
  partition_thread_load_optimized for a given j-th row, determines each thread load (start index and range)
  Returned value is a list of pair, pair.first corresponds to start index and pair.second the range.
  The size of the list is independant from the max_number_of_thread parameter. The max_number_of_thread parameter
  only give a maximum boundary for the number of determinant that can be allocated to a thread.
  Returns an empty list if the computation of the j-th row should be handle by one (current) thread
*/
std::vector<std::pair<int,int>> partition_thread_load_optimized(int j, NTL::Mat<bool> & flags_M)
{
  std::vector<std::pair<int, int>> result;
  int start_index = j-1;
  int end_index = n-j+1;
  int nb_determinant = end_index - start_index;
  int determinant_per_thread = nb_determinant / max_number_of_thread;
  
  if (determinant_per_thread == 0) // Allocate only one thread (current one)
    return result;
  
  // Optimized allocation, skip flagged determinant
  int start = start_index;
  int nb_determinant_allocated = 0;
  for (int i=start_index; i<=end_index; i++)// Verif <= versus < ? 
    {
      if (flags_M[j][i])
	{
	  if (nb_determinant_allocated > 0)
	    {
	      // Allocate a thread
	      result.push_back(std::pair<int, int>(start, nb_determinant_allocated));
	      start = i;
	      nb_determinant_allocated = 0;
	    }
	  continue;
	}
      
      nb_determinant_allocated++;
      if (nb_determinant_allocated >= determinant_per_thread)
	{
	  // Allocate a thread
	  result.push_back(std::pair<int, int>(start, nb_determinant_allocated));
	  start = i;
	  nb_determinant_allocated = 0;
	}
    }
  if (nb_determinant_allocated != 0)
    {
      // Allocate a thread
      result.push_back(std::pair<int, int>(start, nb_determinant_allocated));
    }
  
  return result;
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

/*
  DESCRIPTION:
  Fill the triangular table with zeros when the top part of a square is detected.
  By top part, for position indices i_0, i_1 such that j-1 < i_0<i_1 <n-j, and j>=2,
  If a run of non-zero elements is detected from i_0-1 to to i_1+1 on level j-2,
  and a run of zero elements is detected from i_0 to i_1 on level j-1, then fill the table
  below accordingly.
*/

void s_f(struct experiment & experiment, int j)
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
  //struct experiment & experiment = const_cast<struct experiment &>(exp);
  
  while(i<n-j+1)//using knowledge of (j-2)th, (j-1)th rows, stop when right position = n-j+1 excluded
    {
      int lp = 0;//index of the left non-zero element
      int rp = 0;//index of the right non-zero element
      
      bool chkL1 = (i==j-1) && (experiment.M[j-2][i-1]!=0) && (experiment.M[j-1][i-1]==0) && (experiment.M[j-2][i]!=0) && (experiment.M[j-1][i]==0);
      
      bool chkL2 = (experiment.M[j-2][i-1]!=0) && (experiment.M[j-2][i]!=0) && (experiment.M[j-1][i-1]!=0) && (experiment.M[j-1][i]==0);
      
      if(chkL1||chkL2)
	{
	  lp = i;
	  int t1 = lp+1;
	  
	  while(t1 <= n-j)
	    {
	      bool chkR1 = (experiment.M[j-2][t1-1]!=0) && (experiment.M[j-2][t1]!=0) && (experiment.M[j-1][t1-1]==0) && (experiment.M[j-1][t1]!=0);
	      bool chkR2 = (t1==n-j) && (experiment.M[j-2][t1]!=0) && (experiment.M[j-1][t1]==0) && (experiment.M[j-2][t1+1]!=0) && (experiment.M[j-1][t1+1]==0);
	      
	      if(chkR1||chkR2)
		{
		  rp = t1-1;
		  break;
		}
	      else
		t1 = t1 + 1;
	    }
	  if(rp>lp)
	    {
	      int stop_dim;
	      if(max_dim<=(j+rp-lp-1))
		stop_dim = max_dim;
	      else
		stop_dim = (j+rp-lp-1);
	      
	      for(int x_j=j;x_j<stop_dim;x_j++)//inside square excluding boundaries, square of zeros top row = j-1 which has been completed, thus begin at j
		{
		  int x_i_l,x_i_r;//extreme possible values of i as the filling process goes from top part to bottom part
		  
		  if(lp>=x_j-1)
		    x_i_l = lp;
		  else
		    x_i_l = x_j-1;
		  
		  if(rp<=n-x_j)
		    x_i_r = rp;
		  else
		    x_i_r = n-x_j;
		  
		  for(int x_i=x_i_l;x_i<=x_i_r;x_i++)
		    {
		      //assign_GF2(experiment.M, NTL::GF2(0), x_j, x_i);
		      experiment.M[x_j][x_i]=NTL::GF2(0);
		      experiment.flags_M[x_j][x_i]=1;
		    }
		}
	      i = rp;
	    }
	  else
	    i = i + 1;
	}
      else
	i = i+1;
    }
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

NTL::GF2 solve_eq_for_lower_corner(struct experiment & experiment, int j, int i, int q)
{
  
  /*
    q is the minimum length of the grid required to solve for d_{i,j}.
    The Main matrix has size therefore (q+1) X (q+1).
    The Top/Bottom & Left/Right matrices used to solve for d_{i,j} has size q X q.
  */
  
 
  NTL::GF2 ret = NTL::GF2(0);
  
  NTL::Mat<NTL::GF2> T;
  T.SetDims(q,q);//temporary matrix of size q X q for a minor obtained from the (q+1) X (q+1) Main matrix
  
  for(int g=0;g<q;g++)//index mineure
    {
      for(int h=0;h<q+1;h++)//index rangee de la mineure
	{
	  if(g<h)
	    {
	      for(int k=0;k<q;k++)//index colonne de la mineure
		T[h-1][k]=experiment.M[j-2*q+h+k][i-h+k];
	    }
	  else if (g>h)
	    {
	      for(int k=0;k<q;k++)
		T[h][k]=experiment.M[j-2*q+h+k][i-h+k];
	    }
	  //else it is not part of the expansion of the determinant
	}
      ret = ret+experiment.M[j-q+g][i+q-g]*NTL::determinant(T);
    }
  return ret;
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

bool chk_conds_for_solvability(struct experiment & experiment, int j, int i, int & effective_length)
{
 
  bool ret=false;
  
  for(int w1=2;w1<=max_len_side_grid;w1++)//w1 is the length of the side which is growing so that the main matrix has size (w1+1)x(w1+1)
    {
      for(int rx=0;rx<w1;rx++)
	{
	  for(int cx=0;cx<w1;cx++)
	    {
	      experiment.AspTL[w1][rx][cx] = experiment.M[j-2*w1+rx+cx][i-rx+cx];
	    }
	}
      
      if ( (NTL::determinant(experiment.AspTL[w1])!=0) && (experiment.M[j-w1][i]==0) )/******/
	{
	  ret=true;
	  effective_length=w1;//length of a side of the grid so the main square matrix of size (w1+1) X (w1+1)	  
	  break;
	}
    }
  return ret;
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void d_c_s(struct experiment & experiment, int j, int start, int range)
{
 
  int i = start;
  
  while(i <= start+range)
    {
      int effective_len;
      
      if(!experiment.flags_M[j][i])
	{
  	  /************/
	  
	  if( (j>=2) && (experiment.M[j-2][i]==NTL::GF2(1)) )//North-South-East-West
	    {
	      experiment.M[j][i] = experiment.M[j-1][i-1]*experiment.M[j-1][i+1]+experiment.M[j-1][i];
	      experiment.flags_M[j][i]=true;
	    }
	  else if ( (j>=2*max_len_side_grid) && chk_conds_for_solvability(experiment, j, i, effective_len) )
	    {
	      experiment.M[j][i] = solve_eq_for_lower_corner(experiment, j, i, effective_len);
	      experiment.flags_M[j][i] = true;
	    }
	  else
	    {
	      int t_x=1;
	      
	      while(t_x<=j)
		{
		  if(experiment.M[j-t_x][i]==0)//start check from j-1 to 0 (experiment.M[0][i] == 1 by definition) ( level of sequence )
		    t_x+=1;
		  else
		    break;//when here, experiment.M[j-t_x][i]!=0 AND experiment.M[j-1, j-2, ..., j-(t_x-1)][i]==0
		}
	      //Note that t_x can never be 2 actually for it is handled by the North-South-East-West 
	      if(t_x==1)//explicit computation
		{
		  NTL::Mat<NTL::GF2> tmp;
		  tmp.SetDims(j,j);
		  
		  for(int r = 0;r<j ; r++)
		    {
		      for(int c = 0; c<j ; c++)
			tmp[c][r] = experiment.M[t_x][i-r+c];//==V0[i-j+r+c+1];
		    }
		  
		  experiment.M[j][i] = NTL::determinant(tmp);
		  experiment.flags_M[j][i]=true;
		}
	      else//computation by cross identities
		{
		  NTL::Mat<NTL::GF2> tmp;
		  tmp.SetDims(t_x,t_x);
		  
		  for(int r=0;r<t_x;r++)
		    {
		      for(int c=0;c<t_x;c++)
			tmp[r][c]=experiment.M[j-(t_x-1)][i+c-r];
		    }
		  
		  experiment.M[j][i] = NTL::determinant(tmp);
		  experiment.flags_M[j][i]=true;
		}
	    }
	}
      i=i+1;
    }
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void solve_j_trivial(struct experiment & experiment, int j, int start, int range)
{
  //struct experiment & experiment = const_cast<struct experiment&>(exp);
  
  NTL::Mat<NTL::GF2> tmp;
  tmp.SetDims(j,j);
  
  for(int r1=0;r1<j;r1++)
    {
      for(int c1=0;c1<j;c1++)
	tmp[r1][c1]=0;
    }
    
  for(int i = start; i<=start+range; i++)
    {
      for( int r = 0; r<j ; r++)
	{
	  for( int c = 0; c<j ; c++)
	    {
	      tmp[r][c] = V0[(i+1-j+r+c)];
	    }
	}
      experiment.M[j][i] = NTL::determinant(tmp);
    }
   
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void solve_trivial(struct experiment & experiment)
{
  if (!MULTI_THREAD_ON)
    {
      // Single thread solving
      for (int j=2; j<max_dim+1; j++)
	solve_j_trivial(experiment, j, j-1, n-2*j+2);

      return;
    }
  
  // Multi thread solving
  for (int j=2; j<max_dim+1; j++)
    {
      std::vector<std::pair<int, int>> thread_load = partition_thread_load(j);
    
      if (thread_load.empty())
      {
	return;
	//solve_j_trivial(experiment, j, j-1, n-2*j+2);
	//continue;
      }
      
      std::thread* thread_pool = new std::thread[thread_load.size()];
      for (size_t k=0; k<thread_load.size(); k++)
	thread_pool[k] = std::thread( &solve_j_trivial , std::ref(experiment) , j, thread_load[k].first, thread_load[k].second);
      for (size_t k=0; k<thread_load.size(); k++)
	thread_pool[k].join();
      
      
      delete[] thread_pool;
    }
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void solve_j_fast(struct experiment & exp, int j, int start, int range)
{
  s_f(exp, j);
  d_c_s(exp, j, start, range);
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void solve_fast(struct experiment & experiment)
{
  if (!MULTI_THREAD_ON)
    {
      // Single thread solving
      for (int j=2; j<max_dim+1; j++)
	solve_j_fast(experiment, j, j-1, n-2*j+2);
      return;
    }
  
  // Multi thread solving
  for (int j=2; j<max_dim+1; j++)
    {
      std::vector<std::pair<int, int>> thread_load;
      switch (THREAD_ALLOCATION_STRATEGY)
	{
	case UNIFORM_ALLOCATION:
	  thread_load = partition_thread_load(j);

	  /*
	  for(long hh = 0 ; hh<thread_load.size();hh++)
	    {
	      std::cout << "[" << thread_load[hh].first << " " << thread_load[hh].second << ")" << ", ";
	    }
	  std::cout << "\n";
	  */
	  break;
	case OPTIMIZED_ALLOCATION:
	  thread_load = partition_thread_load_optimized(j, experiment.flags_M);
	  break;
	}
      
      if (thread_load.empty())
	{
	  solve_j_fast(experiment, j, j-1, n-2*j+2);
	  continue;
	}
      
      std::thread* thread_pool = new std::thread[thread_load.size()];
      for (size_t k=0; k<thread_load.size(); k++)
	thread_pool[k] = std::thread(&solve_j_fast, std::ref(experiment), j, thread_load[k].first, thread_load[k].second);
      for (size_t k=0; k<thread_load.size(); k++)
	thread_pool[k].join();
      
      delete[] thread_pool;
    }
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

int main(void)
{
  
  std::cout.precision(6);
  std::cout.setf( std::ios::fixed, std::ios::floatfield );
  
  if(right_index>n){std::cerr << "\n\nOut of bound\n\n"; exit(-1);}
  
  // Initialize generating vector
  init_generating_vector();
  std::cout << "\n\n";
  std::cout << "Original sequence (is identical to row 1 of triangular table), length = " << V0.length() << ", sequence:\n     ";
  for(unsigned int l1 = 0; l1<V0.length(); l1++)
    std::cout << V0[l1] << " ";
  
  std::cout << "\n\n";
  std::cout << "Generating vector, length = " << C_gen.length() << ", vector:\n     ";
  for(unsigned int l1 = 0; l1<C_gen.length(); l1++)
    std::cout << C_gen[l1] << " ";
  
  std::cout << "\n";
  
  // Initialize experiments
  struct experiment 
    trivial_experiment, 
    trivial_experiment_mt, 
    fast_experiment, 
    fast_experiment_mt,
    fast_experiment_mt_opt;
  
  init_experiment(trivial_experiment);
  init_experiment(trivial_experiment_mt);
  init_experiment(fast_experiment);
  init_experiment(fast_experiment_mt);
  init_experiment(fast_experiment_mt_opt);
  
  // Run experiments
  MULTI_THREAD_ON = false;
 
  auto t_start = std::chrono::high_resolution_clock::now();
  solve_trivial(trivial_experiment);
  auto t_end = std::chrono::high_resolution_clock::now();
  time_triv_st = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  
  t_start = std::chrono::high_resolution_clock::now();
  solve_fast(fast_experiment);
  t_end = std::chrono::high_resolution_clock::now();
  time_fast_st = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  
  
  MULTI_THREAD_ON = true;
  THREAD_ALLOCATION_STRATEGY = UNIFORM_ALLOCATION;  
  
  t_start = std::chrono::high_resolution_clock::now();
  solve_trivial(trivial_experiment_mt);
  t_end = std::chrono::high_resolution_clock::now();
  time_triv_mt = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  
  t_start = std::chrono::high_resolution_clock::now();
  solve_fast(fast_experiment_mt);
  t_end = std::chrono::high_resolution_clock::now();
  time_fast_mt =  std::chrono::duration<double, std::milli>(t_end-t_start).count();

  THREAD_ALLOCATION_STRATEGY = OPTIMIZED_ALLOCATION;
  t_start = std::chrono::high_resolution_clock::now();
  solve_fast(fast_experiment_mt_opt);
  t_end = std::chrono::high_resolution_clock::now();
  time_fast_mt_opt =  std::chrono::duration<double, std::milli>(t_end-t_start).count();


  // Prints
  print_stats();

  int nb_differences;
  // Check results
  if (chk_triangular_tables_not_the_same(trivial_experiment.M, fast_experiment.M, max_dim, n))
    {
      how_many_differences(trivial_experiment.M+fast_experiment.M,max_dim,n,nb_differences);
      std::cout << "\n\n*****Mismatches between trivial single threaded and fast single threaded.*****\n";
      std::cout << "#####Number of mismatches = " << nb_differences << "\n";
      std::cout << "EXIT - big problem\n";
      exit(-1);
    }
  
  if (chk_triangular_tables_not_the_same(trivial_experiment.M, trivial_experiment_mt.M, max_dim, n))
    {
      how_many_differences(trivial_experiment.M+trivial_experiment_mt.M,max_dim,n,nb_differences);
     
      std::cout << "\n\n*****Mismatches trivial single threaded and trivial multi threaded.*****\n";
      std::cout << "#####Number of mismatches = " << nb_differences << "\n";
      //exit(-1);
    }
  
  if (chk_triangular_tables_not_the_same(trivial_experiment.M, fast_experiment_mt.M, max_dim, n))
    {
       how_many_differences(trivial_experiment.M+fast_experiment_mt.M,max_dim,n,nb_differences);
      
       std::cout << "\n\n*****Mismatche trivial single threaded and fast multi threaded.*****\n";
       std::cout << "#####Number of mismatches = " << nb_differences << "\n";
       
       //exit(-1);
    }
  
  if (chk_triangular_tables_not_the_same(trivial_experiment.M, fast_experiment_mt_opt.M, max_dim, n))
    {
      how_many_differences(trivial_experiment.M+fast_experiment_mt_opt.M,max_dim,n,nb_differences);
     
      std::cout << "\n\n*****Mismatches trivial and fast optimized multi threaded.*****\n";
      std::cout << "#####Number of mismatches = " << nb_differences << "\n";

      //exit(-1);
    }

  
 
  return 0;
}