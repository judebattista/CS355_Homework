#!/usr/bin/env python3
import itertools
import string
import time


# generate a profile matrix from a set of motifs
def gen_profile_matrix( motifs ) :
  bpi = { 'A': 0, 'C': 1, 'G': 2, 'T': 3} # base pair index to profile
  columns = []
  # change motifs into a list of char lists
  for motif in motifs :
    chars = [ ch for ch in motif ]
    columns.append( chars )
  # transpose the list of char lists so we can score columns
  col_mat = [ list(i) for i in zip(*columns) ]
  profile = []
  # count bps in columns and generate profile matrix per column
  l = len(motifs)
  for col in col_mat :
    counts = [0, 0, 0, 0]
    for bp in col :
      counts[bpi[bp]] += 1 # count this bp
    probs = [ counts[i]/float(l) for i in range(0, 4)]
    profile.append(probs)
  return (profile)

# Score a set of motifs by counting the number of mismatches in the motifs
# we want to minimize this score over all the motifs
def score_motifs( motifs ) :
  bpi = { 'A': 0, 'C': 1, 'G': 2, 'T': 3} # base pair index to profile
  columns = []
  # change motifs into a list of char lists
  for motif in motifs :
    chars = [ ch for ch in motif ]
    columns.append( chars )
  # transpose the list of char lists so we can score columns
  col_mat = [ list(i) for i in zip(*columns) ]
  # count bps in each colum and score the set of motifs
  score = 0
  for col in col_mat :
    counts = [0, 0, 0, 0]
    for bp in col :
      counts[bpi[bp]] += 1 # count this bp
    score += sum(counts) - max(counts)
  return score

# score a potential motif pattern against the
# profile array
def score_motif_with_profile(motif, k, profile ) :
  bpi = { 'A': 0, 'C': 1, 'G': 2, 'T': 3} # base pair index to profile
  score = 1.0
  # compute the probability of the motif based on the profile probabilities
  for i, bp in enumerate(motif) :
    score *= profile[i][bpi[bp]]
  return score

# Look for the most probable kmer based on the profile
def find_profile_most_probable_kmer( dna_string, k, profile ) :
  best_score = -1.0
  top_motifs = []
  for i in range(len(dna_string)-k+1) :
    kmer = dna_string[i:i+k]
    score = score_motif_with_profile( kmer, k, profile)
    if ( score > best_score ) :
      top_motifs.append(kmer)
      best_score = score
      best_motif = kmer
  return best_motif

# Look for the most probable k-motifs in the t dna_strings based on observed profile probabilities
def greedy_motif_search( dna_strings, k, t ) :
  best_score = 0.0
  best_motifs = [ dna[0:k] for dna in dna_strings ]  # initial to first k-mer in each string
  print ( "start best_motifs: ", best_motifs )
  for i in range(len(dna_strings[0])-k+1) :
    kmer = dna_strings[0][i:i+k] # get the next potential kmer from the first string
    next_motifs = [ kmer ]  # reset the next_motifs list
    for i in range(1,t) :
      profile = gen_profile_matrix( next_motifs )
      #print "loop kmer: " , kmer, "next_motifs: ", next_motifs, " profile:", profile
      motif_i = find_profile_most_probable_kmer( dna_strings[i], k, profile )
      next_motifs.append(motif_i)
    print('Profile: ', profile)
    print('Motifs: ', next_motifs)
    score_next = score_motifs(next_motifs)
    score_current_best = score_motifs(best_motifs)
    current_motifs = next_motifs[:]
    if score_next < score_current_best :
      best_motifs = next_motifs
    print('Best Motifs: ', best_motifs)
    print('Best score: ', score_current_best)
  print ( " final_best_motifs: ", best_motifs )
  return best_motifs

###################################################################
# Test the funtion
# read input from the file and output results to another file
def run() :
  filein = open("greedyMotifSearch1.txt")
  data = filein.read()
  args = [s.strip() for s in data.splitlines()]
  params = [ int(s) for s in args[0].split() ]
  print ( params )
  k = params[0]  # the length of the motif
  t = params[1]  # the number of strings to find motifs in
  dna_strings = args[1:] # the dna strings to search for a probable motifs

  # call our  function
  outvar = greedy_motif_search( dna_strings, k, t )
  filein.close()

  fileout = open("rosalind_ba2d_out.txt", 'wt')

  for motif in outvar :
    fileout.write( motif )
    fileout.write("\n")
  fileout.close()

run()
