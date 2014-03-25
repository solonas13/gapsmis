#!/usr/bin/python

import os
import random
import re
import sys

def run_water ( errors, gap_open_pen, gap_ext_pen ):
  os . system ( "water tmp_a.seq tmp_b.seq -datafile EDNAFULL -gapopen " + str ( gap_open_pen ) + " -gapextend " + str ( gap_ext_pen ) + " -outfile tmp.water 2> /dev/null" );
  return ( process_water ( "tmp.water", errors ) )

def process_water ( filename, errors):
  f = open ( filename, 'r' )
  s = f . read ()
  f . close ()
    
  gap_begin = 0

  #find mismatches
  x = re.findall ("#[=]+.+#[=]+(.+)",s,re.DOTALL)
  x = "".join(x).strip()

  snp_cnt = [i for i in x if i == '.']
  snp_cnt = len (snp_cnt)

  x = re . findall ( "([0-9])+[ ]+([A-Z\-]+)[ ]+[0-9]+", s );
  s1 = '' . join ( [ x[i][1] for i in range ( 0, len ( x ), 2 ) ] )
  s2 = '' . join ( [ x[i][1] for i in range ( 1, len ( x ), 2 ) ] )
  
  if x[0][0] != x[1][0]:
    gap_begin = 1;
  else:
    snp_cnt = snp_cnt + int(x[0][0]) - 1

  dash_list  = re . findall ( '([\-]+)[^\-]', s1 )
  dash_list += re . findall ( '([\-]+)[^\-]', s2 )

  gap_begin_len = 0
  if gap_begin == 1: 
    gap_begin_len = abs ( int ( x[0][0] ) - int ( x[1][0] ) )

  dash_list_len = 0
  if len ( dash_list ) > 0:
    for i in range( len( dash_list ) ):
      dash_list_len += len ( dash_list[i] )

  if ( len ( dash_list ) + gap_begin ) > int( errors[0] ): return ( 0 )      # INVALID : Num of inserted gaps > Num of simulated gaps
  
  total_gap_len = dash_list_len + gap_begin_len #If we are here then it is certainly VALID
  
  if ( total_gap_len <= int( errors[1] ) and snp_cnt <= int ( errors[2] ) ):
    return ( 1 )     # CORRECT
  else:
    return ( 2 )     # VALID

def debug ( filename, query, target):
      f = open ( filename, 'a' )
      f . write ( '==BEGIN_OF_INSTANCE==\n' )
      f . write ( 'a: ' + query + '\n' )
      f . write ( 'b: ' + target + '\n' )
      f . write ( '==END_OF_INSTANCE==\n\n' )
      f . close ()

def main ( argv = None ):
  if argv is None: argv = sys . argv

  if len ( argv ) != 4:
    print " usage: ./water_bench.py <gap open penalty [float]> <gap extend penalty [float]> <num of tests [int]>";
    sys . exit ( 0 )

  gap_open_pen = float ( argv[1] )
  gap_ext_pen  = float ( argv[2] )

  match_water_correct = 0
  match_water_valid  = 0

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

    water = run_water ( errors, gap_open_pen, gap_ext_pen )
    if water == 1: match_water_correct += 1; match_water_valid += 1;
    elif water == 2: 
      match_water_valid  += 1;
      debug ( 'water.inc', query, target)
    else:
      debug ( 'water.inv', query, target)

    if i == float ( argv[3] ): break

  fQuery.close()
  fTarget.close()
  fErrors.close()

  print "As `Valid', we define the alignment with num of gaps less or equal to the ones simulated."
  print "As `Correct', we define the valid alignment with total gaps length AND total num of mis less or equal to the ones simulated."
  print "Total tests: " + str ( i )
  print "water"
  print "   Valid: " +  str ( match_water_valid ) + " (Invalid alignments: water.inv)"
  print "   Correct:    " +  str ( match_water_correct ) + " (Valid but not Correct: water.inc)"

  os . system ( "rm tmp_a.seq" )
  os . system ( "rm tmp_b.seq" )
  os . system ( "rm tmp.water" )

if __name__ == "__main__":
  sys . exit ( main ( ) )

