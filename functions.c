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
#include <float.h>
#include <getopt.h>
#include "functions.h"
#include "EDNAFULL.h"
#include "EBLOSUM62.h"

/* long options for command-line switches */
static struct option long_options[] =
 {
   { "sequence-a",         required_argument, NULL, 'a' },
   { "sequence-b",         required_argument, NULL, 'b' },
   { "gap-open-penalty",   required_argument, NULL, 'g' },
   { "gap-extend-penalty", required_argument, NULL, 'e' },
   { "output-file",        required_argument, NULL, 'o' },
   { "data-file",          required_argument, NULL, 'd' },
   { "help",               no_argument,       NULL, 'h' },
   { "max-num-gaps",       required_argument, NULL, 'l' },
   { "max-gap",            required_argument, NULL, 'm' },
   { "local",              required_argument, NULL, 'L' },
   { NULL,                 0,                 NULL, 0   }
 };

/* Decode the input switches */
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
 {
   int          oi;
   int          opt;
   double       val;
   char       * ep;

   /* initialisation */
   sw -> seq_a          = NULL;
   sw -> seq_b          = NULL;
   sw -> gap_open_pen   = 10.0;
   sw -> gap_extend_pen = 0.5;
   sw -> max_gap        = -1;
   sw -> max_num_gaps   = -1;
   sw -> max_gap        = -1;
   sw -> out_file       = ( char * ) malloc ( 15 * sizeof ( char ) );
   sw -> matrix         = ( char * ) malloc ( 15 * sizeof ( char ) );
   sw -> L              = 0;
   strcpy ( sw -> out_file, "gapsmis.out" );
   strcpy ( sw -> matrix, "EDNAFULL" );

   while ( ( opt = getopt_long ( argc, argv, "a:b:g:e:o:d:l:m:L:h", long_options, &oi ) ) != - 1 )
    {
      switch ( opt )
       {
         case 'a':
           sw -> seq_a = optarg;
           break;
         
         case 'b':
           sw -> seq_b = optarg;
           break;

         case 'o':
           free ( sw -> out_file );
           sw -> out_file = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> out_file, optarg );
           break;

         case 'd':
           free ( sw -> matrix );
           sw -> matrix = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> matrix, optarg );
           break;

         case 'g':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> gap_open_pen = val;
           break;

         case 'h':
           return ( 0 );

         case 'e':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> gap_extend_pen = val;
           break;

	 case 'l':
           val = strtol ( optarg, &ep, 10 );
           if ( ! ep || ep == optarg )
            {
              return ( 0 );
            }
           sw -> max_num_gaps = val;
           break;
	 case 'm':
           val = strtol ( optarg, &ep, 10 );
           if ( ! ep || ep == optarg )
            {
              return ( 0 );
            }
           sw -> max_gap = val;
           break;
	 case 'L':
           val = strtol ( optarg, &ep, 10 );
           if ( ! ep || ep == optarg )
            {
              return ( 0 );
            }
           sw -> L = val;
           break;
       }
    }

   return ( optind );
 }

/* Read a fasta file (with our without label) and return a structure containing
   the label and the data
   Parameters:
     Input
*/
struct TSeq * read_fasta_file ( const char * szReadsFile )
 {
   struct TSeq        * seq;
   char               * buf = NULL;
   char                 tmp[BUFFER_SIZE];
   FILE               * fd;

   int                  iMaxFileSize = 0;
   int                  iCurFileSize = 0;
   int                  i, j, n;

   /* opens FASTA file */
   if ( ! ( fd = fopen ( szReadsFile, "r" ) ) )
    {
      fprintf ( stderr, "Error: cannot open %s FASTA file!!!\n", szReadsFile );
      return ( NULL );
    }

   /* reads the contents of the FASTA file into a buffer */
   while ( fgets ( tmp, BUFFER_SIZE, fd ) )
    {
      if ( iMaxFileSize - iCurFileSize <= strlen ( tmp ) )
       {
         buf = ( char * ) realloc ( buf, ( iMaxFileSize + MAX_SIZE ) * sizeof ( char ) );
         if ( ! buf )
          {
            fprintf ( stderr, "Error: not enough memory to read file %s!!!\n", szReadsFile );
            fclose ( fd );
            free ( buf );
            return ( NULL );
          }
         iMaxFileSize = ( iMaxFileSize + MAX_SIZE ) * sizeof ( char );
       }
      buf[iCurFileSize] = '\0';
      buf = strcat ( buf, tmp );
      iCurFileSize += strlen ( tmp );
    }

   if ( ! buf || ! strlen ( buf ) )
    {
      fprintf ( stderr, "Error: file %s is empty!!!\n", szReadsFile );
      fclose ( fd );
      free ( buf );
      return ( NULL );
    }

   fclose ( fd );

   n = strlen ( buf );

   i = 0; j = 0;

   /* allocate memory for the placeholder structure */
   seq  = ( struct TSeq * ) calloc ( 1, sizeof ( struct TSeq ) );

   /* in case it is a FASTA file, locate the description header */
   if ( buf[0] == '>' ) 
    {
      while ( buf[i] != '\n' && buf[i] != '\0' ) ++ i;

      if ( buf[i] == '\0' || i < 2 )
       {
         fprintf ( stderr, "Error: FASTA file %s does not contain any data!!!\n", szReadsFile );
         free ( buf );
         free ( seq );
         return ( NULL );
       }

      /* copy the header in the placeholder structure */
      seq -> header = ( char * ) malloc ( ( i ) * sizeof ( char ) );
      strncpy ( seq -> header, buf + sizeof ( char ), i - 1);
      seq -> header[i - 1] = '\0';
    }

   /* reads the data */
   for ( ; i < n; ++ i )
    {
      if ( buf[i] == '\n' || buf[i] == '\t' || buf[i] == ' ' ) continue;
      else buf[j++] = buf[i];
    }
   buf[j] = '\0';

   if ( ! strlen ( buf ) )
    {
      fprintf ( stderr, "Error: FASTA file %s does not contain any data!!!\n", szReadsFile );
      free ( buf );
      free ( seq -> header );
      free ( seq );
      return ( NULL );
    }
   
   /* reduce the allocated memory to the exact amount */
   buf = realloc ( buf, ( j + 1 ) * sizeof ( char ) );     
   seq -> data = buf;

   return ( seq );
 }


/*
The dynamic programming algorithm for calculating matrices G and H
*/
unsigned int dp_algorithm ( double *** G, unsigned int l, int *** H, char * t, unsigned int n, char * p, unsigned int m, unsigned int matrix, double gap_open_penalty, double gap_extend_penalty )
{
	int i;
	int j;
	int s;
	double matching_score = 0;

        for( s = 0; s < l; s++ )		//Initialisations
		for( i = 0; i < n + 1 ; i++ )
		{
			G[s][i][0] = total_scoring( i, 0, gap_open_penalty, gap_extend_penalty );
			H[s][i][0] = i;
		}

        for ( s = 0; s < l; s++ )
		for( j = 0; j < m + 1 ; j++ )
		{
			G[s][0][j] = total_scoring( j, 0, gap_open_penalty, gap_extend_penalty );
			H[s][0][j] = -j;
		}

	for( i = 1; i < n + 1 ; i++ )		//Compute matrix G[0] and H[0] 
		for( j = 1; j < m + 1; j++ )
		{
			matching_score = ( matrix ? pro_delta( t[i - 1], p[j - 1] ) : nuc_delta( t[i - 1], p[j - 1] ) ) ;
			if ( matching_score == ERR )
				return 0;
			if( i < j )
			{
				double u = G[0][i-1][j-1] + matching_score;
				double v = total_scoring( j - i, G[0][i][i], gap_open_penalty, gap_extend_penalty );
				G[0][i][j] = max ( u, v );

				if ( v > u )
					H[0][i][j] = - ( j - i );	
				else
					H[0][i][j] = 0;
			}
			else if ( i > j )
			{
				double u = G[0][i-1][j-1] + matching_score;
				double v = total_scoring( i - j, G[0][j][j], gap_open_penalty, gap_extend_penalty );
				G[0][i][j] = max ( u, v );

				if ( v > u )
					H[0][i][j] = ( i - j );		
				else
					H[0][i][j] = 0;

			}
			else if ( i == j )
			{
				G[0][i][j] = G[0][i-1][j-1] + matching_score;
				H[0][i][j] = 0;
			}
		}

        for ( s = 1;  s < l; s++ )      //Compute matrix H[s] and G[s], for all s = 1,...,l-1
        {
                int * minhval;   //store minimum from edge of matrix
                if( ( minhval = ( int * ) calloc( ( n + 1 ), sizeof( int ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: memory could not be allocated!!!\n");
                        return ( 0 );
                }
                for( j = 1; j < m + 1; j++ )
                {
                        int r = 0;
                        for( i = 1; i < n + 1; i++ )
                        {
				matching_score = ( matrix ? pro_delta( t[i - 1], p[j - 1] ) : nuc_delta( t[i - 1], p[j - 1] ) ) ;
				if ( matching_score == ERR )
					return 0;

                                unsigned int update_ver = 0;
                                if( G[s - 1][r][j] < G[s - 1][i][j] )
                                {
                                        r = i;
                                        update_ver = 1;
                                }

                                unsigned int update_hor = 0;
                                if( G[s - 1][i][minhval[i]] < G[s - 1][i][j] )
                                {
                                        minhval[i] = j;
                                        update_hor = 1;
                                }

				double u = total_scoring( j - minhval[i], G[s - 1][i][minhval[i]], gap_open_penalty, gap_extend_penalty );
				double v = total_scoring( i - r, G[s - 1][r][j], gap_open_penalty, gap_extend_penalty );
				double w = G[s][i - 1][j - 1] + matching_score;

                                if( u > w )
                                {
                                        if( v > u )
                                        {
                                                G[s][i][j] = v;
                                                if ( ! ( update_ver ) )
                                                        H[s][i][j] = i - r;
                                                else
                                                        H[s][i][j] = H[s - 1][i][j];
                                        }
                                        else
                                        {
                                                G[s][i][j] = u;
                                                if ( ! ( update_hor ) )
                                                        H[s][i][j] = - ( j - minhval[i] );
                                                else
                                                        H[s][i][j] = H[s - 1][i][j];
                                        }
                                }
                                else
                                {
                                        if( v > w )
                                        {
                                                G[s][i][j] = v;
                                                if ( ! ( update_ver ) )
                                                        H[s][i][j] = i - r;
                                                else
                                                        H[s][i][j] = H[s - 1][i][j];
                                        }
                                        else
                                        {
                                                G[s][i][j] = w;
                                                H[s][i][j] = 0; //there is no gap
                                        }
                                }
                        }
                }
                free ( minhval );
        }

	return 1;
}

/*
The dynamic programming algorithm for calculating pruned matrices G and H
*/
unsigned int dp_algorithm_pruned ( double *** G, unsigned int l, int *** H, char * t, unsigned int n, char * p, unsigned int m, unsigned int matrix, double gap_open_penalty, double gap_extend_penalty, unsigned int MAXgap )
{
	int i;
	int j;
	int s;
	double matching_score = 0;
	unsigned int i_max = min ( n, m + MAXgap );

        for( s = 0; s < l; s++ )		//Initialisations
		for( i = 0; i < n + 1 ; i++ )
		{
			G[s][i][0] = total_scoring( i, 0, gap_open_penalty, gap_extend_penalty );
			H[s][i][0] = i;
		}

        for ( s = 0; s < l; s++ )
		for( j = 0; j < m + 1 ; j++ )
		{
			G[s][0][j] = total_scoring( j, 0, gap_open_penalty, gap_extend_penalty );
			H[s][0][j] = -j;
		}

	for( i = 1; i < i_max + 1; i++)
	{
		unsigned int j_min = max ( 1, ( int ) ( i - MAXgap ));
                unsigned int j_max = min ( m, ( int ) ( i + MAXgap ));
		for( j = j_min; j <= j_max; j++ )
		{
			matching_score = ( matrix ? pro_delta( t[i - 1], p[j - 1] ) : nuc_delta( t[i - 1], p[j - 1] ) ) ;
			if ( matching_score == ERR )
				return 0;
			if( i < j )
			{
				double u = G[0][i-1][j-1] + matching_score;
				double v = total_scoring( j - i, G[0][i][i], gap_open_penalty, gap_extend_penalty );
				G[0][i][j] = max ( u, v );

				if ( v > u )
					H[0][i][j] = - ( j - i );	
				else
					H[0][i][j] = 0;
			}
			else if ( i > j )
			{
				double u = G[0][i-1][j-1] + matching_score;
				double v = total_scoring( i - j, G[0][j][j], gap_open_penalty, gap_extend_penalty );
				G[0][i][j] = max ( u, v );

				if ( v > u )
					H[0][i][j] = ( i - j );		
				else
					H[0][i][j] = 0;

			}
			else if ( i == j )
			{
				G[0][i][j] = G[0][i-1][j-1] + matching_score;
				H[0][i][j] = 0;
			}
		}
	}

        for ( s = 1;  s < l; s++ )      //Compute matrix H[s] and G[s], for all s = 1,...,l-1
        {
                int * minvval;   //store minimum from edge of matrix
                if( ( minvval = ( int * ) calloc( ( m + 1 ), sizeof( int ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: memory could not be allocated!!!\n");
                        return ( 0 );
                }
                for( i = 1; i < i_max + 1; i++ )
                {
                        int r = 0;
			unsigned int j_min = max ( 1, ( int ) ( i - MAXgap ));
                	unsigned int j_max = min ( m, ( int ) ( i + MAXgap ));
			for( j = j_min; j <= j_max; j++ )
                        {
				matching_score = ( matrix ? pro_delta( t[i - 1], p[j - 1] ) : nuc_delta( t[i - 1], p[j - 1] ) ) ;
				if ( matching_score == ERR )
					return 0;

                                unsigned int update_hor = 0;
                                if( G[s - 1][i][r] < G[s - 1][i][j] )
                                {
                                        r = j;
                                        update_hor = 1;
                                }

                                unsigned int update_ver = 0;
                                if( G[s - 1][minvval[j]][j] < G[s - 1][i][j] )
                                {
                                        minvval[j] = i;
                                        update_ver = 1;
                                }

				double u = total_scoring( i - minvval[j], G[s - 1][minvval[j]][j], gap_open_penalty, gap_extend_penalty );
				double v = total_scoring( j - r, G[s - 1][i][r], gap_open_penalty, gap_extend_penalty );
				double w = G[s][i - 1][j - 1] + matching_score;

                                if( u > w )
                                {
                                        if( v > u )
                                        {
                                                G[s][i][j] = v;
                                                if ( ! ( update_hor ) )
                                                        H[s][i][j] = - ( j - r );
                                                else
                                                        H[s][i][j] = H[s - 1][i][j];
                                        }
                                        else
                                        {
                                                G[s][i][j] = u;
                                                if ( ! ( update_ver ) )
                                                        H[s][i][j] = + ( i - minvval[j] );
                                                else
                                                        H[s][i][j] = H[s - 1][i][j];
                                        }
                                }
                                else
                                {
                                        if( v > w )
                                        {
                                                G[s][i][j] = v;
                                                if ( ! ( update_hor ) )
                                                        H[s][i][j] = - ( j - r );
                                                else
                                                        H[s][i][j] = H[s - 1][i][j];
                                        }
                                        else
                                        {
                                                G[s][i][j] = w;
                                                H[s][i][j] = 0; //there is no gap
                                        }
                                }
                        }
                }
                free ( minvval );
        }

	return 1;
}

/*
The dynamic programming algorithm for calculating matrices G and H
*/
unsigned int dp_algorithm_lcl ( double *** G, unsigned int l, int *** H, char * t, unsigned int n, char * p, unsigned int m, unsigned int matrix, double gap_open_penalty, double gap_extend_penalty )
{
	int i;
	int j;
	int s;
	double matching_score = 0;

        for( s = 0; s < l; s++ )		//Initialisations
		for( i = 0; i < n + 1 ; i++ )
		{
			G[s][i][0] = 0;
			H[s][i][0] = i;
		}

        for ( s = 0; s < l; s++ )
		for( j = 0; j < m + 1 ; j++ )
		{
			G[s][0][j] = 0;
			H[s][0][j] = -j;
		}

	for( i = 1; i < n + 1 ; i++ )		//Compute matrix G[0] and H[0] 
		for( j = 1; j < m + 1; j++ )
		{
			matching_score = ( matrix ? pro_delta( t[i - 1], p[j - 1] ) : nuc_delta( t[i - 1], p[j - 1] ) ) ;
			if ( matching_score == ERR )
				return 0;
			if( i < j )
			{
				double u = G[0][i-1][j-1] + matching_score;
				double v = total_scoring( j - i, G[0][i][i], gap_open_penalty, gap_extend_penalty );
				G[0][i][j] = max ( max ( u, v ), 0 );

				if ( v > u )
					H[0][i][j] = - ( j - i );	
				else
					H[0][i][j] = 0;
			}
			else if ( i > j )
			{
				double u = G[0][i-1][j-1] + matching_score;
				double v = total_scoring( i - j, G[0][j][j], gap_open_penalty, gap_extend_penalty );
				G[0][i][j] = max ( max ( u, v ), 0 );

				if ( v > u )
					H[0][i][j] = ( i - j );		
				else
					H[0][i][j] = 0;

			}
			else if ( i == j )
			{
				G[0][i][j] = max ( G[0][i-1][j-1] + matching_score, 0 );
				H[0][i][j] = 0;
			}
		}

        for ( s = 1;  s < l; s++ )      //Compute matrix H[s] and G[s], for all s = 1,...,l-1
        {
                int * minhval;   //store minimum from edge of matrix
                if( ( minhval = ( int * ) calloc( ( n + 1 ), sizeof( int ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: memory could not be allocated!!!\n");
                        return ( 0 );
                }
                for( j = 1; j < m + 1; j++ )
                {
                        int r = 0;
                        for( i = 1; i < n + 1; i++ )
                        {
				matching_score = ( matrix ? pro_delta( t[i - 1], p[j - 1] ) : nuc_delta( t[i - 1], p[j - 1] ) ) ;
				if ( matching_score == ERR )
					return 0;

                                unsigned int update_ver = 0;
                                if( G[s - 1][r][j] < G[s - 1][i][j] )
                                {
                                        r = i;
                                        update_ver = 1;
                                }

                                unsigned int update_hor = 0;
                                if( G[s - 1][i][minhval[i]] < G[s - 1][i][j] )
                                {
                                        minhval[i] = j;
                                        update_hor = 1;
                                }

				double u = total_scoring( j - minhval[i], G[s - 1][i][minhval[i]], gap_open_penalty, gap_extend_penalty );
				double v = total_scoring( i - r, G[s - 1][r][j], gap_open_penalty, gap_extend_penalty );
				double w = G[s][i - 1][j - 1] + matching_score;

                                if( u > w )
                                {
                                        if( v > u )
                                        {
                                                G[s][i][j] = max ( v, 0 );
                                                if ( ! ( update_ver ) )
                                                        H[s][i][j] = i - r;
                                                else
                                                        H[s][i][j] = H[s - 1][i][j];
                                        }
                                        else
                                        {
                                                G[s][i][j] = max ( u, 0 );
                                                if ( ! ( update_hor ) )
                                                        H[s][i][j] = - ( j - minhval[i] );
                                                else
                                                        H[s][i][j] = H[s - 1][i][j];
                                        }
                                }
                                else
                                {
                                        if( v > w )
                                        {
                                                G[s][i][j] = max ( v, 0 );
                                                if ( ! ( update_ver ) )
                                                        H[s][i][j] = i - r;
                                                else
                                                        H[s][i][j] = H[s - 1][i][j];
                                        }
                                        else
                                        {
                                                G[s][i][j] = max ( w, 0 );
                                                H[s][i][j] = 0; //there is no gap
                                        }
                                }
                        }
                }
                free ( minhval );
        }

	return 1;
}

/*
The dynamic programming algorithm for calculating pruned matrices G and H
*/
unsigned int dp_algorithm_pruned_lcl ( double *** G, unsigned int l, int *** H, char * t, unsigned int n, char * p, unsigned int m, unsigned int matrix, double gap_open_penalty, double gap_extend_penalty, unsigned int MAXgap )
{
	int i;
	int j;
	int s;
	double matching_score = 0;
	unsigned int i_max = min ( n, m + MAXgap );

        for( s = 0; s < l; s++ )		//Initialisations
		for( i = 0; i < n + 1 ; i++ )
		{
			G[s][i][0] = 0;
			H[s][i][0] = i;
		}

        for ( s = 0; s < l; s++ )
		for( j = 0; j < m + 1 ; j++ )
		{
			G[s][0][j] = 0;
			H[s][0][j] = -j;
		}

	for( i = 1; i < i_max + 1; i++)
	{
		unsigned int j_min = max ( 1, ( int ) ( i - MAXgap ));
                unsigned int j_max = min ( m, ( int ) ( i + MAXgap ));
		for( j = j_min; j <= j_max; j++ )
		{
			matching_score = ( matrix ? pro_delta( t[i - 1], p[j - 1] ) : nuc_delta( t[i - 1], p[j - 1] ) ) ;
			if ( matching_score == ERR )
				return 0;
			if( i < j )
			{
				double u = G[0][i-1][j-1] + matching_score;
				double v = total_scoring( j - i, G[0][i][i], gap_open_penalty, gap_extend_penalty );
				G[0][i][j] = max ( max ( u, v ), 0 );

				if ( v > u )
					H[0][i][j] = - ( j - i );	
				else
					H[0][i][j] = 0;
			}
			else if ( i > j )
			{
				double u = G[0][i-1][j-1] + matching_score;
				double v = total_scoring( i - j, G[0][j][j], gap_open_penalty, gap_extend_penalty );
				G[0][i][j] = max ( max ( u, v ), 0 );

				if ( v > u )
					H[0][i][j] = ( i - j );		
				else
					H[0][i][j] = 0;

			}
			else if ( i == j )
			{
				G[0][i][j] = max ( G[0][i-1][j-1] + matching_score, 0 );
				H[0][i][j] = 0;
			}
		}
	}

        for ( s = 1;  s < l; s++ )      //Compute matrix H[s] and G[s], for all s = 1,...,l-1
        {
                int * minvval;   //store minimum from edge of matrix
                if( ( minvval = ( int * ) calloc( ( m + 1 ), sizeof( int ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: memory could not be allocated!!!\n");
                        return ( 0 );
                }
                for( i = 1; i < i_max + 1; i++ )
                {
                        int r = 0;
			unsigned int j_min = max ( 1, ( int ) ( i - MAXgap ));
                	unsigned int j_max = min ( m, ( int ) ( i + MAXgap ));
			for( j = j_min; j <= j_max; j++ )
                        {
				matching_score = ( matrix ? pro_delta( t[i - 1], p[j - 1] ) : nuc_delta( t[i - 1], p[j - 1] ) ) ;
				if ( matching_score == ERR )
					return 0;

                                unsigned int update_hor = 0;
                                if( G[s - 1][i][r] < G[s - 1][i][j] )
                                {
                                        r = j;
                                        update_hor = 1;
                                }

                                unsigned int update_ver = 0;
                                if( G[s - 1][minvval[j]][j] < G[s - 1][i][j] )
                                {
                                        minvval[j] = i;
                                        update_ver = 1;
                                }

				double u = total_scoring( i - minvval[j], G[s - 1][minvval[j]][j], gap_open_penalty, gap_extend_penalty );
				double v = total_scoring( j - r, G[s - 1][i][r], gap_open_penalty, gap_extend_penalty );
				double w = G[s][i - 1][j - 1] + matching_score;

                                if( u > w )
                                {
                                        if( v > u )
                                        {
                                                G[s][i][j] = max ( v, 0 );
                                                if ( ! ( update_hor ) )
                                                        H[s][i][j] = - ( j - r );
                                                else
                                                        H[s][i][j] = H[s - 1][i][j];
                                        }
                                        else
                                        {
                                                G[s][i][j] = max ( u, 0 );
                                                if ( ! ( update_ver ) )
                                                        H[s][i][j] = + ( i - minvval[j] );
                                                else
                                                        H[s][i][j] = H[s - 1][i][j];
                                        }
                                }
                                else
                                {
                                        if( v > w )
                                        {
                                                G[s][i][j] = max ( v, 0 );
                                                if ( ! ( update_hor ) )
                                                        H[s][i][j] = - ( j - r );
                                                else
                                                        H[s][i][j] = H[s - 1][i][j];
                                        }
                                        else
                                        {
                                                G[s][i][j] = max ( w, 0 );
                                                H[s][i][j] = 0; //there is no gap
                                        }
                                }
                        }
                }
                free ( minvval );
        }

	return 1;
}

/* Returns the score for matching character a and b based on EDNAFULL matrix */
int nuc_delta ( char a, char b )
 {
   unsigned int index_a = nuc_char_to_index ( a );
   unsigned int index_b = nuc_char_to_index ( b );

   if ( ( index_a < NUC_SCORING_MATRIX_SIZE ) && ( index_b < NUC_SCORING_MATRIX_SIZE ) )
     return ( EDNAFULL_matrix[ index_a ][ index_b ] );
   else //Error
     return ( ERR );
 }

/* Returns the score for matching character a and b based on EBLOSUM62 matrix */
int pro_delta ( char a, char b )
 {
   unsigned int index_a = pro_char_to_index( a );
   unsigned int index_b = pro_char_to_index( b );

   if ( ( index_a < PRO_SCORING_MATRIX_SIZE ) && ( index_b < PRO_SCORING_MATRIX_SIZE ) )
     return ( EBLOSUM62_matrix[ index_a ][ index_b ] );
   else //Error
     return ( ERR );
 }

/* Returns the index of char a in EDNAFULL matrix */
unsigned int nuc_char_to_index ( char a )
 {
   unsigned int index; 

   switch ( a )
    {
      case 'A':
        index = 0; break;

      case 'T':
        index = 1; break;

      case 'G':
        index = 2; break;

      case 'C':
        index = 3; break;

      case 'S':
        index = 4; break;

      case 'W':
        index = 5; break;

      case 'R':
        index = 6; break;

      case 'Y':
        index = 7; break;

      case 'K':
        index = 8; break;

      case 'M':
        index = 9; break;

      case 'B':
        index = 10; break;

      case 'V':
        index = 11; break;

      case 'H':
        index = 12; break;

      case 'D':
        index = 13; break;

      case 'N':
        index = 14; break;

      default:
        fprintf ( stderr, "Error: unrecognizable character in one of the nucleotide sequences!!!\n" );
        index = ERR; break;
    }
   
   return ( index );
 }

/* Returns the index of char a in EBLOSUM62 matrix */
unsigned int pro_char_to_index ( char a )
 {
   unsigned int index; 

   switch ( a )
    {
      case 'A':
        index = 0; break;

      case 'R':
        index = 1; break;

      case 'N':
        index = 2; break;

      case 'D':
        index = 3; break;

      case 'C':
        index = 4; break;

      case 'Q':
        index = 5; break;

      case 'E':
        index = 6; break;

      case 'G':
        index = 7; break;

      case 'H':
        index = 8; break;

      case 'I':
        index = 9; break;

      case 'L':
        index = 10; break;

      case 'K':
        index = 11; break;

      case 'M':
        index = 12; break;

      case 'F':
        index = 13; break;

      case 'P':
        index = 14; break;

      case 'S':
        index = 15; break;

      case 'T':
        index = 16; break;

      case 'W':
        index = 17; break;

      case 'Y':
        index = 18; break;

      case 'V':
        index = 19; break;

      case 'B':
        index = 20; break;

      case 'Z':
        index = 21; break;

      case 'X':
        index = 22; break;

      case '*':
        index = 23; break;

      default:
        fprintf ( stderr, "Error: unrecognizable character in one of the protein sequences!!!\n" );
        index = ERR; break;
    }
   return ( index );
 }

/*
Computes the optimal alignment using matrix G in O(2*MAXgap+1) time
Note:   double gap_open_penalty, double gap_extend_penalty, double gap_open_offset_penalty are arguments given by the user to represent the gap penalty.
*/
unsigned int opt_solution ( double *** G, unsigned int l, unsigned int n, unsigned int m, double * MAXscore, unsigned int * start, unsigned int * MINnumgaps )
{
        double score = -DBL_MAX;
        unsigned int i, s;

        for ( s = 0; s < l ; s++ )
        {
                for ( i = 0 ; i < n + 1 ; i++ )
                {
                        double temp_score = G[s][i][m];
                        if ( temp_score > score )
                        {
                                score = temp_score;
                                ( * MAXscore ) = score;
                                ( * start ) = i;                //backtrace from cell G[start,m]
                                ( * MINnumgaps ) = s + 1;
                        }
                }
        }
        return 1;
}

/*
Computes the optimal alignment using matrix G in O(2*MAXgap+1) time
Note:   double gap_open_penalty, double gap_extend_penalty, double gap_open_offset_penalty are arguments given by the user to represent the gap penalty.
*/
unsigned int opt_solution_lcl ( double *** G, unsigned int l, unsigned int n, unsigned int m, double * MAXscore, unsigned int * istart, unsigned int * jstart,unsigned int * MINnumgaps )
{
        double score = -DBL_MAX;
        unsigned int i, j ,s;

        for ( s = 0; s < l ; s++ )
        {
                for ( i = 0 ; i < n + 1 ; i++ )
                {
                	for ( j = 0 ; j < m + 1 ; j++ )
                	{
				double temp_score = G[s][i][j];
				if ( temp_score > score )
				{
					score = temp_score;
					( * MAXscore ) = score;
					( * istart ) = i;                //backtrace from cell G[s,i,j]
					( * jstart ) = j;                
					( * MINnumgaps ) = s + 1;
				}
			}
                }
        }
        return 1;
}

#if 0
/* Gives the position of the gap in O(m) time */
unsigned int backtracing ( int ** H, unsigned int m, unsigned int n, unsigned int start, unsigned int * gaps_pos, unsigned int beta, unsigned int * gaps_len, unsigned int * where )
{
	int i, j;

	i = start; j = m; 	//we start backtracing from the last column

	while ( i >= 0 && j >= 0)
	{
		if ( H[i][j] == 0 )
		{
			--i; --j;
		}
		else				
		{
			if ( H[i][j] < 0 )  	//the gap is inserted in the text
			{	
				gaps_pos[(beta)] = i;
				gaps_len[(beta)] = -H[i][j];
				where[(beta)] = 1;
			        j = j + H[i][j] ;
			}
			else			//the gap is inserted in the pattern
			{		
				gaps_pos[(beta)] = j;
				gaps_len[(beta)] = H[i][j];
				where[(beta)] = 2;
			        i = i - H[i][j];
			}
			(beta)++;
		}
	}
	return 1;
}
#endif

/* Gives the position of the gap in O(m) time */
unsigned int backtracing ( int ** H, unsigned int m, unsigned int n, unsigned int start, unsigned int * gaps_pos, unsigned int l, unsigned int * gaps_len, unsigned int * where )
{
        int i, j, s;

        i = start; j = m; s = 0;        //we start backtracing from the last column

        while ( i >= 0 && j >= 0)
        {
                if ( H[i][j] == 0 )
                {
                        --i; --j;
                }
                else
                {
                        if ( H[i][j] < 0 )      //the gap is inserted in the text
                        {
                                gaps_pos[s] = i;
                                gaps_len[s] = - H[i][j];
                                where[s] = 1;
                                j = j + H[i][j];
                        }
                        else                    //the gap is inserted in the pattern
                        {
                                gaps_pos[s] = j;
                                gaps_len[s] = H[i][j];
                                where[s] = 2;
                                i = i - H[i][j];
                        }
                        s++;
                }
        }
        return 1;
}

/* Gives the position of the gap in O(m) time */
unsigned int backtracing_lcl ( double ** G, unsigned int m, unsigned int n, int ** H, unsigned int istart, unsigned int jstart, unsigned int * gaps_pos, unsigned int l, unsigned int * gaps_len, unsigned int * where, unsigned int * iend, unsigned int * jend )
{
        int i, j, s;

        i = istart; j = jstart; s = 0;        //we start backtracing from the last column
	( * iend ) = i;                
	( * jend ) = j;

        while ( i >= 0 && j >= 0 )
        {
                if ( H[i][j] == 0 )
                {
                        --i; --j;
                }
                else
                {
                        if ( H[i][j] < 0 )      //the gap is inserted in the text
                        {
                                gaps_pos[s] = i;
                                gaps_len[s] = - H[i][j];
                                where[s] = 1;
                                j = j + H[i][j];
                        }
                        else                    //the gap is inserted in the pattern
                        {
                                gaps_pos[s] = j;
                                gaps_len[s] = H[i][j];
                                where[s] = 2;
                                i = i - H[i][j];
                        }
                        s++;
                }
                if ( G[i][j] > 0 )
		{
			( * iend ) = i;                //backtrace up to cell G[s,i,j]
			( * jend ) = j;
		}
		else
			break;
        }
        return 1;
}

/*
Gives the total score of an alignment in constant time
Note: double matrix_score is the value of G[i][m], i.e. the score of an alignment WITHOUT the gap penalties
*/
double total_scoring( unsigned int gap, double matrix_score, double gap_open_penalty, double gap_extend_penalty )
 {
   return ( matrix_score + ( ( gap > 0 ) ? ( gap - 1 ) * gap_extend_penalty + gap_open_penalty : 0 ) );
 }

/* Swaps the text and the pattern in case m > n */
unsigned int swap_txt_pat ( struct TSeq ** seqa, unsigned int * n, struct TSeq ** seqb, unsigned int * m )
 {
   struct TSeq * tmp;

   tmp   = *seqa;
   *seqa = *seqb;
   *seqb = tmp;
   
   SWAP ( *n, *m );
   
   return ( 1 );
 }

