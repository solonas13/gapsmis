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

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdio.h>
#include <getopt.h>
#include "types.h"

#define MAX_SIZE 1024
#define BUFFER_SIZE 128
#define NUC_SCORING_MATRIX_SIZE 15		
#define PRO_SCORING_MATRIX_SIZE 24		
#define ERR 24					//error number returned if char_to_index returns an invalid index

struct TSeq * read_fasta_file ( const char * szReadsFile );

int decode_switches ( int argc, char ** argv, struct TSwitch * );

unsigned int dp_algorithm ( double *** G, unsigned int l, int *** H, char * t, unsigned int n, char * p, unsigned int m, unsigned int matrix, double gap_open_penalty, double gap_extend_penalty );

unsigned int dp_algorithm_pruned ( double *** G, unsigned int l, int *** H, char * t, unsigned int n, char * p, unsigned int m, unsigned int matrix, double gap_open_penalty, double gap_extend_penalty, unsigned int MAXgap );

unsigned int dp_algorithm_lcl ( double *** G, unsigned int l, int *** H, char * t, unsigned int n, char * p, unsigned int m, unsigned int matrix, double gap_open_penalty, double gap_extend_penalty );

int nuc_delta ( char a, char b );

int pro_delta ( char a, char b );

unsigned int nuc_char_to_index ( char a );

unsigned int pro_char_to_index ( char a );

unsigned int opt_solution ( double *** G, unsigned int l, unsigned int n, unsigned int m, double * MAXscore, unsigned int * start, unsigned int * MINnumgaps );

unsigned int opt_solution_lcl ( double *** G, unsigned int l, unsigned int n, unsigned int m, double * MAXscore, unsigned int * istart, unsigned int * jstart,unsigned int * MINnumgaps );

double total_scoring ( unsigned int gap, double current_score, double gap_open_penalty, double gap_extend_penalty );

unsigned int backtracing ( int ** H, unsigned int m, unsigned int n, unsigned int start, unsigned int * gaps_pos, unsigned int l, unsigned int * gaps_len, unsigned int * where );

unsigned int backtracing_lcl ( double ** G, unsigned int m, unsigned int n, int ** H, unsigned int istart, unsigned int jstart, unsigned int * gaps_pos, unsigned int l, unsigned int * gaps_len, unsigned int * where, unsigned int * iend, unsigned int * jend );

unsigned int swap_txt_pat ( struct TSeq ** seqa, unsigned int * n, struct TSeq ** seqb, unsigned int * m );

#endif
