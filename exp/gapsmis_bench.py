#!/usr/bin/python

import os
import random
import re
import sys

def run_gapsmis ( errors, gap_open_pen, gap_ext_pen, num_gaps ):
  os . system ( "../gapsmis -a tmp_a.seq -b tmp_b.seq -o tmp.gapsmis -g " + str ( gap_open_pen ) + " -e " + str ( gap_ext_pen ) + " -l " + str ( num_gaps ) + " -m 30 " )
  return ( process_gapsmis ( errors ) )

def process_gapsmis ( errors ):
  f = open ( 'tmp.gapsmis', 'r' )
  s = f . read ()
  f . close ()
  
  x = re . search ( 'Number of gaps: ([0-9]+)', s )
  y = re . search ( 'Length of gaps: ([0-9]+)', s )
  z = re . search ( 'Number of mismatches: ([0-9]+)', s )

  if int (  x . group ( 1 ) ) > int( errors[0] ): #invalid
    return ( 0 )
						  #else it is valid	
  if int (  y . group ( 1 ) ) <= int( errors[1] ) and int (  z . group ( 1 ) ) <= int( errors[2] ): #and correct
    return ( 1 ) 
  else:
    return ( 2 ) #or else it is valid but incorrect

def debug ( filename, query, target, errors ):
      f = open ( filename, 'a' )
      f . write ( '==BEGIN_OF_INSTANCE==\n' )
      f . write ( 'a: ' + query + '\n' )
      f . write ( 'b: ' + target + '\n' )
      #f . write ( 'Num of gaps: ' +  errors[0] + '\n' )
      #f . write ( 'Len of gaps: ' +  errors[1] + '\n' )
      #f . write ( 'Num of mis: ' +  errors[2] + '\n' )
      f . write ( '==END_OF_INSTANCE==\n\n' )
      f . close ()

def main ( argv = None ):
  if argv is None: argv = sys . argv

  if len ( argv ) != 5:
    print " usage: ./gapsmis_bench.py <gap open penalty [float]> <gap extend penalty [float]> <num of gaps [int]> <num of tests [int]>";
    sys . exit ( 0 )

  gap_open_pen = float ( argv[1] )
  gap_ext_pen  = float ( argv[2] )
  num_gaps = int ( argv[3] )

  match_gapsmis_correct  = 0
  match_gapsmis_valid   = 0

  print "Running tests..."

  fQuery  = open ('queries.fa', 'r')
  fTarget = open ('targets.fa', 'r')
  fErrors = open ('errors.txt', 'r')

  i = 0;
  while 1:
    query = fQuery.readline()
    if not query: break
    i = i + 1
    query = fQuery.readline()

    target = fTarget.readline()
    target = fTarget.readline()

    errors = fErrors.readline()
    errors = errors.split(',')


    f1 = open ( 'tmp_a.seq', 'w' )
    f2 = open ( 'tmp_b.seq', 'w' )
    f1 . write ( query )
    f2 . write ( target )
    f1 . close ()
    f2 . close ()

    gapsmis  = run_gapsmis ( errors, gap_open_pen, gap_ext_pen, num_gaps )
    if gapsmis == 1: match_gapsmis_correct += 1; match_gapsmis_valid += 1;
    elif gapsmis == 2:
      match_gapsmis_valid  += 1;
      debug ( 'gapsmis.inc', query, target, errors )
    else:
      debug ( 'gapsmis.inv', query, target, errors )
    if i == int ( argv[4] ): break

  fQuery.close()
  fTarget.close()
  fErrors.close()
  print "As `Valid', we define the alignment with num of gaps less or equal to the ones simulated."
  print "As `Correct', we define the valid alignment with total gaps length AND total num of mis less or equal to the ones simulated."
  print "Total tests: " + str ( i )
  print "gapsmis"
  print "   Valid: " +  str ( match_gapsmis_valid ) + " (Invalid alignments: gapsmis.inv)"
  print "   Correct:    " +  str ( match_gapsmis_correct ) + " (Valid but not Correct: gapsmis.inc)"

  os . system ( "rm tmp_a.seq" )
  os . system ( "rm tmp_b.seq" )
  os . system ( "rm tmp.gapsmis" )

if __name__ == "__main__":
  sys . exit ( main ( ) )

