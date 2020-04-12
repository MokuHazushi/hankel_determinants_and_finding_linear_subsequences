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
  
  COMPILE AND LINK COMMAND: "g++ -Wall -Wpedantic -g -O5 -I<NTL directory>/include -I<GMPL directory>/include hanmat_sampling.cpp -o hanmat_sampling -std=c++11 -L<NTL library>/lib -lntl -L<GMPL libray>/lib -lgmp -lm -lpthread"
  
  EXAMPLE OF RUN COMMAND: "./hanmat_sampling > res.txt &2>&1 &"
  
  DESCRIPTION:
  It averages the times and counts for hanmat code over a sample split across NUM_THREADS threads.
  Each thread has SAMPLE_SIZE_PER_THREAD sequences to process for time and count evaluations.
  The total sample size is thus NUM_THREADS*SAMPLE_SIZE_PER_THREAD.

  PARAMETERS:
  Number of threads, and sample size per thread can be changed just a few line below.
  Length of sequences, left and right index of linear subsequences can be changed below.
  Length of the generating vector can be changed as well.
  See not far below for a comment that mentions that the parameters can be changed.

  OUTPUT
  Print statistics on the standard output.

  NOTE
  This is almost a cut and paste of hanmat.cpp code adapted in a multithread setting in order to sample.
  Report to hanmat.cpp for more details about modules in this code.
  
*/

#include <iomanip>
#include <iostream>
#include <fstream>
#include <random>
#include <thread>
#include <sys/resource.h>

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/GF2.h>
#include <NTL/vec_vec_GF2.h>
#include <NTL/mat_GF2.h>

/***** ***** ***** ***** *****/
//SOME PARAMETERS CAN BE CHANGED HERE:

const int NUM_THREADS = 507;
const int SAMPLE_SIZE_PER_THREAD = 20;

const int n=4096;
const double ratio_size=1.0/16.0;
const int len_C_gen = (int)floor((double)n*ratio_size);//length of generating vector
//const int left_index = len_C_gen;//where the linear subsequence begins
//const int right_index = n;//where the linear subsequence ends
const int left_index = (int)(floor(7*len_C_gen));//where the linear subsequence begins
const int right_index = (int)(floor(9*len_C_gen));//where the linear subsequence ends

/***** ***** ***** ***** *****/

NTL::Vec<NTL::GF2> C_gen[NUM_THREADS];//generation vector
NTL::Vec<NTL::GF2> V0[NUM_THREADS];
NTL::Mat<NTL::GF2> M[NUM_THREADS];
NTL::Mat<bool> flags_M[NUM_THREADS];
NTL::Mat<NTL::GF2> new_M[NUM_THREADS];

//NTL::Mat<NTL::GF2> *AspC[NUM_THREADS];
NTL::Mat<NTL::GF2> *AspTL[NUM_THREADS];

const int max_dim=(n/2)*( n%2 == 0) + ((n+1)/2)*( n%2 == 1);
const int max_len_side_grid=7;
double timing[NUM_THREADS][2];
double sum_timing[2];

int *ct_det_tab[NUM_THREADS];
int *ct_grid[NUM_THREADS];
int ct_sqfill[NUM_THREADS];
int ct_nsew[NUM_THREADS];

double *sum_ct_det_tab;
double *sum_ct_grid;
double sum_ct_sqfill;
double sum_ct_nsew;

double get_ms_res_time(void)
{ 
 
  double current_time;

  struct rusage ruse;
  getrusage(RUSAGE_THREAD, &ruse);
  current_time = ((double) ruse.ru_utime.tv_sec * 1000.0) + (((double) ruse.ru_utime.tv_usec) / 1000.0); //resultat en millisecondes   

  return current_time;
}

bool chk_triangular_tables_not_the_same(NTL::Mat<NTL::GF2> M1, NTL::Mat<NTL::GF2> M2,int max_dim,int n)
{
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


void s_f(int index, int j)
{
  int i = j-1;
 
  while(i<=n-j)
    {
      int lp = 0;
      int rp = 0;
     
      bool chkL1 = (i==j-1) && (new_M[index][j-2][i-1]!=0) && (new_M[index][j-1][i-1]==0) && (new_M[index][j-2][i]!=0) && (new_M[index][j-1][i]==0);

      bool chkL2 = (new_M[index][j-2][i-1]!=0) && (new_M[index][j-2][i]!=0) && (new_M[index][j-1][i-1]!=0) && (new_M[index][j-1][i]==0);
    
      if(chkL1||chkL2)
	{
	  lp = i;
	  int t1 = lp+1;
	  
	  while(t1 <= n-j)
	    {
	      bool chkR1 = (new_M[index][j-2][t1-1]!=0) && (new_M[index][j-2][t1]!=0) && (new_M[index][j-1][t1-1]==0) && (new_M[index][j-1][t1]!=0);
	      bool chkR2 = (t1==n-j) && (new_M[index][j-2][t1]!=0) && (new_M[index][j-1][t1]==0) && (new_M[index][j-2][t1+1]!=0) && (new_M[index][j-1][t1+1]==0);
	      
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
		  stop_dim=max_dim;
		}
	      else
		{
		  stop_dim=j+rp-lp-1;
		}
	      
	      for(int x_j=j;x_j<stop_dim;x_j++)
		{
		  int x_i_l,x_i_r;

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
		      new_M[index][x_j][x_i]=0;
		      flags_M[index][x_j][x_i]=1;
		      
		      ct_sqfill[index]+=1;
		      
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

NTL::GF2 solve_eq_for_lower_corner(int index, int j, int i, int q)
{
  NTL::GF2 ret = (NTL::GF2)0;
  
  NTL::Mat<NTL::GF2> T;
  T.SetDims(q,q);
  
  for(int g=0;g<q;g++)
    {
      for(int h=0;h<q+1;h++)
	{
	  if(g<h)
	    {
	      for(int k=0;k<q;k++)
		{
		  T[h-1][k]=new_M[index][j-2*q+h+k][i-h+k];
		}
	    }
	  else if (g>h)
	    {
	      for(int k=0;k<q;k++)
		{
		  T[h][k]  =new_M[index][j-2*q+h+k][i-h+k];
		} 
	    }
	}
      ret = ret+new_M[index][j-q+g][i+q-g]*NTL::determinant(T);
    } 
  return ret;
}

int effective_length[NUM_THREADS];

bool chk_conds_for_solvability(int index, int j, int i)
{
  bool ret=false;
  
  for(int w1=2;w1<=max_len_side_grid;w1++)
    {
      for(int rx=0;rx<w1;rx++)
	{
	  for(int cx=0;cx<w1;cx++)
	    {
	      AspTL[index][w1][rx][cx]=new_M[index][j-2*w1+rx+cx][i-rx+cx];
	    }
	  /*
	  for(int rx=0;rx<w1-1;rx++)
	    {
	      for(int cx=0;cx<w1-1;cx++)
		{
		  AspC[index][w1-1][rx][cx]=new_M[index][j-2*(w1-1)+rx+cx][i-rx+cx];
		}
	    }
	  */
	}

      if ( (NTL::determinant(AspTL[index][w1])!=0) && (new_M[index][j-w1][i]==0))
	{
	  ret=true;
	  effective_length[index]=w1;	  
	  break;
	}
    }
  return ret;
}


void d_c_s(int index, int j)
{
  int i = j-1;

  while(i<n-j+1)
    {
      if(!flags_M[index][j][i])
	{
	  if( (j>=2) && (new_M[index][j-2][i]==(NTL::GF2)1) )
	    {
	      
	      new_M[index][j][i]=new_M[index][j-1][i-1]*new_M[index][j-1][i+1]+new_M[index][j-1][i];
	      flags_M[index][j][i]=true;
	      ct_nsew[index]+=1;
	      
	    }
	  else if ( (j>=2*max_len_side_grid) && chk_conds_for_solvability(index,j,i) )
	    {
	      new_M[index][j][i] = solve_eq_for_lower_corner(index,j,i,effective_length[index]);
	      flags_M[index][j][i] = true;
	      ct_grid[index][effective_length[index]]+=1;
	    }
	  else
	    {
	      int t_x=1;
	      
	      while(t_x<=j)
		{
		  if(new_M[index][j-t_x][i]==0)
		    {
		      t_x+=1;
		    }
		  else
		    {
		      break;
		    }
		}
	    
	      if(t_x==1)
		{
		  
		  NTL::Mat<NTL::GF2> tmp;
		  tmp.SetDims(j,j);
		  
		  for(int r = 0;r<j ; r++)
		    {
		      for(int c = 0; c<j ; c++)
			{
			  tmp[c][r] = new_M[index][t_x][i-r+c];
			}
		    }
		  ct_det_tab[index][0]+=1;
		  new_M[index][j][i] = NTL::determinant(tmp);
		  flags_M[index][j][i]=true;
		}
	      else
		{
		  NTL::Mat<NTL::GF2> tmp;
		  tmp.SetDims(t_x,t_x);
		  
		  for(int r=0;r<t_x;r++)
		    {
		      for(int c=0;c<t_x;c++)
			{
			  tmp[r][c]=new_M[index][j-(t_x-1)][i+c-r];
			}
		    }
		  new_M[index][j][i] = NTL::determinant(tmp);
		  flags_M[index][j][i]=true;
		  ct_det_tab[index][t_x-1]+=1;
		}
	    }
	}
      i=i+1;
    }
}



void timing_estimate(int index)
{
  std::random_device r_dev;
  std::mt19937_64 mt(r_dev());
  std::bernoulli_distribution dist(0.5);
  
  if(right_index>n){std::cerr << "\n\nOut of bound\n\n"; exit(-1);}

  /*
    Begin alloc of datum
  */
  C_gen[index].SetLength(len_C_gen);
  V0[index].SetLength(n);
  M[index].SetDims(max_dim+1,n);
  
  //AspC[index] = new NTL::Mat<NTL::GF2> [max_dim+1];
  AspTL[index] = new NTL::Mat<NTL::GF2> [max_dim+1];
  
  flags_M[index].SetDims(max_dim+1,n);
  new_M[index].SetDims(max_dim+1,n);
  
  //AspC[index][0].SetDims(1,1);
  AspTL[index][0].SetDims(1,1);
  for(int t=1;t<max_dim+1;t++)
    {
      //AspC[index][t].SetDims(t,t);
      AspTL[index][t].SetDims(t,t);
    }
  
  ct_grid[index]=new int [n];
  ct_det_tab[index]=new int [n];
  /*
    End of alloc datum
  */

  /*
    Begin init datum
  */
  timing[index][0]=0.0;
  timing[index][1]=0.0;
  bzero(ct_grid[index],n*sizeof(int));
  ct_nsew[index]=0;
  ct_sqfill[index]=0;
  bzero(ct_det_tab[index],n*sizeof(int));
  /*
    End init datum
  */
 
  /*
    Now each thread does its part of sampling, it is a stratified sampling.
  */
  for(int rep=0;rep < SAMPLE_SIZE_PER_THREAD; rep++)
    {
      
      for(unsigned int i1=0;i1<C_gen[index].length();i1++)
	{
	  C_gen[index][i1]=(NTL::GF2)(dist(mt));
	}
  
      C_gen[index][C_gen[index].length()-1]=(NTL::GF2)1;
      
      for(int x1=0;x1<n;x1++)
	{
	  V0[index][x1]=(NTL::GF2)(dist(mt));
	}
  
      for(int x1=left_index ;  x1 < right_index; x1++)
	{
	  NTL::GF2 tmpsum=(NTL::GF2)0;
	  for(int y1 = 0; y1 < C_gen[index].length()-1; y1++)
	    {
	      tmpsum = tmpsum + (C_gen[index][ C_gen[index].length()-2-y1 ]*V0[index][ x1 -1 - y1 ]);
	    }
	  V0[index][x1]=(NTL::GF2)(tmpsum);
	}
      
      /**************************************************************/
      
      for(int r1=0;r1<max_dim+1;r1++)
	{
	  for(int c1=0;c1<n;c1++)
	    {
	      M[index][r1][c1]=(NTL::GF2)0;
	    }
	}
      
      double tnaive0 = get_ms_res_time();
      for(int c1=0;c1<n;c1++)
	{
	  M[index][0][c1]=(NTL::GF2)1;
	}
      for(int c1=0;c1<n;c1++)
	{
	  M[index][1][c1]=V0[index][c1];
	}
      
      bool RUN_NAIVE = true;
      if(RUN_NAIVE){
	for(int j = 2; j<max_dim+1; j++)
	  {
	    //std::cout << j << "\n";
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
			tmp[r][c] = V0[index][(i+1-j+r+c)];
		      }
		  }
		M[index][j][i] = NTL::determinant(tmp);
	      }
	  }
      }
    
      double tnaive1 = get_ms_res_time();
      double time_diff = tnaive1-tnaive0;
      
      timing[index][0] += time_diff;
      /**************************************************************/
  
      double tnew0 = get_ms_res_time();
      
      for(int i1=0;i1<max_dim+1;i1++)
	{
	  for(int i2=0;i2<n;i2++)
	    {
	      new_M[index][i1][i2]=(NTL::GF2)0;
	      flags_M[index][i1][i2]=false;
	    }
	}
      
      for(int i2=0;i2<n;i2++)//j=0,1
	{
	  new_M[index][0][i2]=(NTL::GF2)1;
	  flags_M[index][0][i2]=true;
	  new_M[index][1][i2]=(NTL::GF2)V0[index][i2];
	  flags_M[index][1][i2]=true;
	}
     
      for(int j = 2; j<max_dim+1; j++)
	{
	  s_f(index,j);
	  d_c_s(index,j);
	}
      
      double tnew1 = get_ms_res_time();
      double time_diff_new = tnew1-tnew0;
      
      timing[index][1]+=time_diff_new;
      
      if(chk_triangular_tables_not_the_same(new_M[index], M[index],max_dim, n))
	{
	  std::cerr << "problems with the mathematical thoery - exit now\n"; exit(-1);
	}
    }
 
}

int main(void)
{
  std::cout.precision(6);
  std::cout.setf( std::ios::fixed, std::ios::floatfield );
  
  std::thread th[NUM_THREADS];
  
  for(int i1 = 0; i1 < NUM_THREADS;i1 = i1+1)
    {
      th[i1] = std::thread(timing_estimate,i1);
    }
  
  for(int i1=0;i1< NUM_THREADS;i1=i1+1)
    {
      th[i1].join();
    }
  
  sum_timing[0]=0.0;
  sum_timing[1]=0.0;
  
  sum_ct_sqfill=0.0;
  sum_ct_det_tab=new double[n];
  bzero(sum_ct_det_tab,n*sizeof(double));
  sum_ct_grid=new double[n];
  bzero(sum_ct_grid,n*sizeof(double));
  sum_ct_nsew=0.0;
   
  for(int i1=0;i1<NUM_THREADS;i1=i1+1)
    {
      
      sum_timing[0]+=timing[i1][0];
      sum_timing[1]+=timing[i1][1];
      sum_ct_sqfill+=ct_sqfill[i1];
      sum_ct_nsew+=ct_nsew[i1];
    }
  for(int e2=0;e2<n;e2++)
    {
      for(int i1=0;i1<NUM_THREADS;i1=i1+1)
	{
	  sum_ct_grid[e2]+=ct_grid[i1][e2];
	  sum_ct_det_tab[e2]+=ct_det_tab[i1][e2];
	}
    }

  int number_of_entries = ((n%2)==0)*((n/2)-1)*(n/2)+((n%2)==1)*(n/2)*(n-1)/2;
  int out_width0 = (int)floor(1.0+(log(n)/log(10)));
  int out_width1 = (int)floor(1.0+(log(number_of_entries)/log(10)));
  
  std::cout << "\n[A." << std::setw(out_width0) <<0<< "]  sample size                       = " << SAMPLE_SIZE_PER_THREAD*NUM_THREADS;
  std::cout << "\n[A." << std::setw(out_width0) <<1<< "]  length of a sequence              = " << n;
  std::cout << "\n[A." << std::setw(out_width0) <<2<< "]  length of generating vector       = " << len_C_gen;
  std::cout << "\n[A." << std::setw(out_width0) <<3<< "]  left index linear subsequence     = " << left_index;
  std::cout << "\n[A." << std::setw(out_width0) <<4<< "]  right index linear subsequene     = " << right_index;
  std::cout << "\n[A." << std::setw(out_width0) <<5<< "]  number of entries of a table      = " << number_of_entries << " (analytic count)\n";

  std::cout << "\n[B." << std::setw(out_width0) <<0<< "]  time new                          = " << sum_timing[1]/(SAMPLE_SIZE_PER_THREAD*NUM_THREADS) << " ms";
  std::cout << "\n[B." << std::setw(out_width0) <<1<< "]  time naive                        = " << sum_timing[0]/(SAMPLE_SIZE_PER_THREAD*NUM_THREADS) << " ms";
  std::cout << "\n[B." << std::setw(out_width0) <<2<< "]  ratio B.0/B.1                     = " << sum_timing[1]/sum_timing[0] << "\n";


  for(int q1=1;q1<n;q1++)
    {
      if(sum_ct_det_tab[q1]!=0)
	{
	  std::cout << "\n[C." << std::setw(out_width0) << q1 << "]  counts cross level " << std::setw(out_width1) << q1 << "          = " << sum_ct_det_tab[q1]/(SAMPLE_SIZE_PER_THREAD*NUM_THREADS);
	}
    }
  std::cout << "\n";
  for(int q1=0;q1<n;q1++)
    {
      if(sum_ct_grid[q1]!=0)
	{
	  std::cout << "\n[D." << std::setw(out_width0) << q1 << "]  counts inter grid "<<q1<<"by"<<q1 << std::setw(out_width1) << "            = " << sum_ct_grid[q1]/(SAMPLE_SIZE_PER_THREAD*NUM_THREADS);
	}
      
    }
  
  double tmp_zozodet = 0.0;
  double tmp_zozogrid = 0.0;
  for(int q1=0;q1<n;q1++)
    {
      tmp_zozogrid+=sum_ct_grid[q1];
      tmp_zozodet+=sum_ct_det_tab[q1];
    }
  std::cout << "\n\n[E.0]  counts North-South-East-West        = " << sum_ct_nsew/(SAMPLE_SIZE_PER_THREAD*NUM_THREADS);
  std::cout << "\n\n[F.0]  counts square filling               = " << sum_ct_sqfill/(SAMPLE_SIZE_PER_THREAD*NUM_THREADS);
  std::cout << "\n\n[G.0]  counts direct evaluations           = " << sum_ct_det_tab[0]/(SAMPLE_SIZE_PER_THREAD*NUM_THREADS);
  std::cout << "\n\ntotal sum of counts                        = " << ((double)(tmp_zozodet+tmp_zozogrid+sum_ct_sqfill+sum_ct_nsew)/(double)(SAMPLE_SIZE_PER_THREAD*NUM_THREADS));

  
  delete [] sum_ct_grid;
  delete [] sum_ct_det_tab;
  for(int i1=0;i1<NUM_THREADS;i1++)
    {
      delete [] ct_grid[i1];
      delete [] ct_det_tab[i1];
      //delete [] AspC[i1];
      delete [] AspTL[i1];
    }
 
  std::cout << "\n\n";
  return EXIT_SUCCESS;
}
