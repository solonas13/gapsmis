/**
    GapsMis: flexible sequence alignment with a bounded number of gaps.
    Copyright (C) 2013 Solon P. Pissis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"
#include "output.h"

int main ( int argc, char ** argv)
 {
   struct TSwitch  sw;
   char   * out_file;

   unsigned int MAXnumgaps;	//input argument		
   double gap_open_pen;		//input argument
   double gap_extend_pen;	//input argument
   unsigned int scoring_matrix; //input argument
   unsigned int MAXgap;	

   double MAXscore;		//to be computed
   unsigned int MINnumgaps;	//to be computed		
   
   unsigned int start;		//where to start backtracing
   unsigned int * gaps_pos;	//position of the gap(s)
   unsigned int * gaps_len;	//len of the gap(s)
   unsigned int * where;	//where is the gap(s): text [1] or pattern [2]
   
   struct TSeq * t;		//text t
   unsigned int n; 		//length of text
   struct TSeq * p;		//pattern p
   unsigned int m; 		//length of pattern
   
   double ***       G; 		//dynamic programming matrix
   int *** 	 H; 		//backtracing matrix
   	
   unsigned int swap;           //swap the text and the pattern in case m < n        
   unsigned int i, j;

   /* checks the arguments */
   i = decode_switches ( argc, argv, &sw );

   if ( i < 5 || ! sw . seq_a || ! sw . seq_b ) 
    {
      usage ();
      return ( 1 );
    }
   else 
    {
      gap_open_pen   = - sw . gap_open_pen;	//the penalties should have a negative value
      gap_extend_pen = - sw . gap_extend_pen;
      out_file       =   sw . out_file;

      if ( ! strcmp ( "EDNAFULL", sw . matrix ) )       scoring_matrix = 0;
      else if ( ! strcmp ( "EBLOSUM62", sw . matrix ) ) scoring_matrix = 1;
      else
       {
         fprintf ( stderr, "Error: scoring matrix argument should be `EDNAFULL' for nucleotide sequences or `EBLOSUM62' for protein sequences!!!\n" );
         return ( 1 );
       }
    }
   
   /* reads the text */
   t = read_fasta_file ( sw . seq_a );
   if ( ! t )
    {
      fprintf (stderr, "Error: cannot read file %s!!!\n", sw . seq_a );
      return ( 1 );
    }
   if ( ! t -> header ) t -> header = strdup ( "Seq A" );
   
   /* reads the pattern */
   p = read_fasta_file ( sw . seq_b );
   if ( ! p )
    {
      fprintf( stderr, "Error: cannot read file %s!!!\n", sw . seq_b );
      return ( 1 );
    }
   if ( ! p -> header ) p -> header = strdup ( "Seq B" );
   
   /* calculate text's and pattern's length */
   n = strlen ( t -> data );
   m = strlen ( p -> data );
   
   /* checks the lengths of text and pattern and swaps if needed */
   if ( m > n )
    {
      swap_txt_pat ( &t, &n, &p, &m );
      swap = 1;
    }
   else
      swap = 0;
   
   /* checks the max num of gaps: MAXnumgaps < n */
   MAXnumgaps =  ( sw . max_num_gaps <= -1 ) ?  2 : sw . max_num_gaps;
   if( MAXnumgaps >= n )
    {
      fprintf ( stderr, "Error: the max gap length should be less than the length of the text!!!\n" );
      return ( 1 );
    }

   MAXgap =  ( sw . max_gap <= -1 ) ?  n - 1 : sw . max_gap;
   if( MAXgap >= n )
    {
      fprintf ( stderr, "Error: the max gap length should be less than the length of the text!!!\n" );
      return ( 1 );
    }
   
   /* 3d dynamic memory allocation for matrices G and H*/
   if( ( G = ( double *** ) malloc ( ( MAXnumgaps ) * sizeof( double ** ) ) ) == NULL )
    {
      fprintf( stderr, "G could not be allocated\n" );
      return 0;
    } 
	
   for ( i = 0; i < MAXnumgaps; i ++ )
    {
      if( ( G[i] = ( double ** ) malloc ( ( n + 1 ) * sizeof( double * ) ) ) == NULL )
       {
         fprintf( stderr, "G could not be allocated\n" );
         return 0;
       }
 
      if( ( G[i][0] = ( double * ) calloc ( ( n + 1 ) * ( m + 1 ), sizeof( double ) ) ) == NULL )
       {
	  fprintf( stderr, "G could not be allocated\n" );
	  return 0;
       } 

      for ( j = 1; j < n + 1; j++ )
        G[i][j] = ( void * ) G[i][0] + j * ( m + 1 ) * sizeof( double );
    }

   if( ( H = ( int *** ) malloc ( ( MAXnumgaps ) * sizeof( int ** ) ) ) == NULL )
    {
      fprintf( stderr, "H could not be allocated\n" );
      return 0;
    } 
	
   for ( i = 0; i < MAXnumgaps; i ++ )
    {
      if( ( H[i] = ( int ** ) malloc ( ( n + 1 ) * sizeof( int * ) ) ) == NULL )
       {
         fprintf( stderr, "H could not be allocated\n" );
         return 0;
       }
 
      if( ( H[i][0] = ( int * ) calloc ( ( n + 1 ) * ( m + 1 ), sizeof( int ) ) ) == NULL )
       {
	  fprintf( stderr, "H could not be allocated\n" );
	  return 0;
       } 

      for ( j = 1; j < n + 1; j++ )
        H[i][j] = ( void* ) H[i][0] + j * ( m + 1 ) * sizeof( int );
    }

   /* dynamic programming algorithm */
   if ( MAXgap == n - 1 )
    {
      if ( ! ( dp_algorithm( G, MAXnumgaps, H, t -> data, n, p -> data, m, scoring_matrix, gap_open_pen, gap_extend_pen ) ) )
       {
         fprintf ( stderr, "Error: dp_algorithm() failed!!!\n" );
         return ( 1 );	
       }
    }
   else
    {
      if ( ! ( dp_algorithm_pruned( G, MAXnumgaps, H, t -> data, n, p -> data, m, scoring_matrix, gap_open_pen, gap_extend_pen, MAXgap ) ) )
       {
         fprintf ( stderr, "Error: dp_algorithm_pruned() failed!!!\n" );
         return ( 1 );	
       }
    }

   #if 1
   int s;
   for ( s = 0;  s < MAXnumgaps; s++ )
        {

                for(i = 0; i < n+1; i++)                //Matrix G output
                {
                        for(j = 0; j < m+1; j++)
                        {
                                fprintf( stderr,"%d\t", (int) G[s][i][j] );
                        }
                        fprintf(stderr,"\n");
                }

                fprintf(stderr,"\n");

                #if 1
                for(i = 0; i < n+1; i++)                //Matrix H output
                {
                        for(j = 0; j < m+1; j++)
                        {
                                fprintf( stderr,"%d\t", H[s][i][j]);
                        }
                        fprintf(stderr,"\n");
                }

                fprintf(stderr,"\n");
                #endif
        }
   #endif

   MINnumgaps = 0;	//to be computed		
   /* finds the optimal alignment based on the matrices scores */
   opt_solution ( G, MAXnumgaps, n, m, &MAXscore, &start, &MINnumgaps );

   if( ( gaps_pos = ( unsigned int * ) calloc ( MINnumgaps, sizeof( unsigned int ) ) ) == NULL )
    {
      fprintf( stderr, "gaps_pos could not be allocated\n" );
      return 0;
    }

   if( ( gaps_len = ( unsigned int * ) calloc ( MINnumgaps, sizeof( unsigned int ) ) ) == NULL )
    {
      fprintf( stderr, "gaps_pos could not be allocated\n" );
      return 0;
    }

   if( ( where = ( unsigned int * ) calloc ( MINnumgaps, sizeof( unsigned int ) ) ) == NULL )
    {
      fprintf( stderr, "where could not be allocated\n" );
      return 0;
    }
 
   /* computes the position of the gap */
   backtracing ( H[MINnumgaps - 1], m, n, start, gaps_pos, MINnumgaps, gaps_len, where );

   /* outputs the results */
   if ( ! ( results( out_file, t, n, p, m, MAXscore, gaps_pos, MINnumgaps, gaps_len, where, swap, scoring_matrix, gap_open_pen, gap_extend_pen ) ) )
    {
      fprintf(stderr, "Error: results() failed!!!\n");
      return ( 1 );	
    }

   for ( i = 0;  i < MAXnumgaps; i++ )
    {
      free ( G[i][0] );
      free ( G[i] );
      free ( H[i][0] );
      free ( H[i] );
    }
   free ( G );
   free ( H );
   free ( gaps_pos );
   free ( gaps_len );
   free ( where );
   free ( t -> header );
   free ( p -> header );
   free ( t );
   free ( p );
   free ( sw . out_file );
   free ( sw . matrix );
   return ( 0 );
 }



