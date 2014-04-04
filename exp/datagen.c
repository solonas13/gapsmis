#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define         GAP_TYPE_INSERTION              0
#define         GAP_TYPE_DELETION               1
#define         GAP_TYPE_SNP                    2
#define         GAP_RATE                        5.7*10E-6
#define         MAX_GAP_SIZE                    27
#define         GAP_INS_RATE                    0.42
#define         GAP_DEL_RATE                    0.58
#define         SNP_RATE                        0.0016
#define         SEQ_LEN_LIMIT_LOWER             500
#define         DNA_ALPH_SIZE                   4
#define         TARGET_READ_OFFSET              60


#define         MIN(x, y) (((x) < (y)) ? (x) : (y))

struct read_info
 {
   char       * header;
   char       * seq;
 };

struct gap_info
 {
   int          pos;
   int          len;
   int          type;
 };

double gap_rate[MAX_GAP_SIZE] = 
 {
   0.319218341, 0.068318203, 0.397876596, 
   0.033174335, 0.012463456, 0.069979997,
   0.005570088, 0.003908294, 0.031020157, 
   0.000830897, 0.001477150, 0.011724881,
   0.003815972, 0.000646253, 0.021541776, 
   0.000092322, 0.000646253, 0.008524388,
   0.002369595, 0.000153870, 0.004462225, 
   0.000923219, 0.000000000, 0.000707801,
   0.000061548, 0.000030774, 0.000461609
 };

char DNA[4] = "ACGT";

struct read_info * read_fasta(const char * filename, int * seq_len)
{
  FILE               * fp;
  int                  size;
  struct read_info   * read;
  char               * x;
  char               * y;
  char               * dst;
  int                  n;
  
  fp = fopen(filename, "r");

  fseek(fp, 0, SEEK_END);
  size = ftell(fp);
  rewind(fp);

  read = (struct read_info *)malloc(sizeof(struct read_info));
  read->header = (char *)malloc((size + 1) * sizeof(char));

  if ( ! ( fread(read->header, size, sizeof(char), fp) > 0 ) )
    fprintf ( stderr, " Error: input file not in FASTA format!" );
  read->header[size] = '\0';

  /* locate the end of read label */
  x    = strchr(read->header, '\n');
  (*x) = '\0';

  /* specify the sequence section */
  read->seq = ++x;

  
  dst = x = strchr(x, '\n');
  if (x) y = strchr(x + 1, '\n');
  memmove(dst, x + 1, y - x - 1);
  dst += y - x - 1;

  while (1)
   {
     x = y;
     y = strchr(x + 1, '\n');
     if (y)
      {
        memmove(dst, x + 1, y - x - 1);
        dst += y - x - 1;
      }
     else
      {
        n = strlen(x + 1);
        memmove(dst, x + 1, n);
        dst += n;
        (*dst) = '\0';
        break;
      }
   }

  *seq_len = strlen(read->seq);

  fclose(fp);

  return (read);
}

int select_indel(int * ins_cnt, int * del_cnt)
{
  int           coin;

  if (!(*ins_cnt))
   {
     --(*del_cnt);
     return (GAP_TYPE_DELETION);
   }
  
  if (!(*del_cnt))
   {
     --(*ins_cnt);
     return (GAP_TYPE_INSERTION);
   }

  coin = rand() % 2;

  if (coin == GAP_TYPE_INSERTION)
   {
     --(*ins_cnt);
   }
  else
   {
     --(*del_cnt);
   }

  return (coin);
}

int select_pos(char * reserved, int gap_type, int seq_size, int gap_size)
{
  int                   pos = 0;
  int                   dirty = 0;
  int                   i;

  assert(gap_type == GAP_TYPE_DELETION || gap_type == GAP_TYPE_INSERTION);

  switch (gap_type)
   {
     case GAP_TYPE_INSERTION:
       do
        {
          pos = rand() % seq_size;
        } while (reserved[pos] == 1);
       reserved[pos] = 1;
       break;
     
     case GAP_TYPE_DELETION:
     do
        {
          dirty = 0;
          pos   = rand() % (seq_size - gap_size);
       
          for (i = 0; i < gap_size; ++ i)
           {
             if (reserved[pos + i])
              {
                dirty = 1;
                break;
              }
           }
        } while (dirty);
       for (i = 0; i < gap_size; ++ i) reserved[pos + i] = 1;
       break;
   }
  return (pos);
}

struct gap_info ** random_gaps(int seq_len, int gaps, int * new_seq_len)
{
  int                   i, j;
  struct gap_info    * gi;
  int                   sum = 0;
  int                   gaps_occ[MAX_GAP_SIZE];
  char                * reserved;
  struct gap_info    ** gap_mapping;
  int                   ins_cnt, del_cnt;
  int                   ins_sum = 0, del_sum = 0;

  printf ("Total gaps: %d\n", gaps);


  /* Calculate the number of gaps of each size for the given sequence length */
  for (i = 0; i < MAX_GAP_SIZE; ++ i)
   {
     gaps_occ[i] = (int)ceil(gap_rate[i] * gaps);
     sum += gaps_occ[i];
   }

  /* marking which positions of the sequence have been already marked with an indel */
  reserved    = (char *)calloc(seq_len, sizeof(char));
  gap_mapping = (struct gap_info **)calloc(seq_len, sizeof(struct gap_info *));

  /* How many of the gaps are insertions and how many deletions */
  ins_cnt = GAP_INS_RATE * gaps;
  del_cnt = gaps - ins_cnt;
  
  /* fill gi with random gap positions */
  for (i = 0; i < MAX_GAP_SIZE; ++ i)
   {
     for (j = 0; j < gaps_occ[i]; ++ j)
      {
        gi       = (struct gap_info *)malloc(sizeof(struct gap_info));
        gi->type = select_indel(&ins_cnt, &del_cnt);
        gi->len  = i + 1;
        gi->pos  = select_pos(reserved, gi->type, seq_len, gi->len);

        if (gi->type == GAP_TYPE_INSERTION)
         {
           ins_sum += i + 1;
         }
        else
         {
           del_sum += i + 1;
         }

        gap_mapping[gi->pos] = gi;
      }
   }

  (*new_seq_len) = seq_len - del_sum + ins_sum;
  free (reserved);

  return (gap_mapping);
}

char * random_snp(const char * seq, int seq_len, int gaps, struct gap_info ** gap_mapping)
{
  char                * reserved;
  char                * snp_mapping;
  char                  c;
  int                   i, j;
  int                   pos;

  reserved    = (char *)calloc(seq_len, sizeof(char));
  snp_mapping = (char *)calloc(seq_len, sizeof(char));

   /* reconstruct reserved */
  for (i = 0; i < seq_len; ++ i)
   {
     if (gap_mapping[i] && gap_mapping[i]->type == GAP_TYPE_DELETION)
      {
        for (j = 0; j < gap_mapping[i]->len; ++ j)
         {
           reserved[i + j] = 1;
         }
      }
   }

  for (i = 0; i < gaps; ++ i)
   {
       do
        {
          pos = rand() % seq_len;
        } while (reserved[pos]);
       reserved[pos] = 1;

       do
        {
          c = DNA[rand() % DNA_ALPH_SIZE];
        } while (c == seq[pos]);
       snp_mapping[pos] = c;
   }
  free (reserved);

  return (snp_mapping);
}

int calculate_length(struct gap_info ** gm, int seq_len, int pos)
{
  int                   i;
  int                   sum = 0;

  i = pos;
  while (i < seq_len)
   {
     if (gm[i])
      {
        if (gm[i]->type == GAP_TYPE_INSERTION) 
         {
           sum += gm[i]->len;
           ++sum;
           ++i;
           continue;
         }
        else  
         {
           sum -= gm[i]->len;
           i += gm[i]->len;
         }
      }
     else 
      {
        sum ++;
        ++ i;
      }
   }
  return (sum);
}

void compute_reads(int cnt_reads, struct gap_info ** gm, char * seq, int seq_len, int read_len, int nr_gaps, const char * sm, int * cnt_gaps)
{
  FILE                * fpr;
  FILE                * fpg;
  FILE                * fpt;
  int                   i, j, k, m;
  char                * reserved;
  int                   pos;
  char                * read;
  char                * target;

  int                   gaps;
  int                   gap_sum;
  int                   snp_sum;

  m = 0;

  read   = (char *)malloc((read_len + 1) * sizeof(char));
  if (!read)
   {
     printf("MEM ERROR!\n");
     exit(1);
   }
  printf("Read of size %d\n", read_len);
  target = (char *)malloc((read_len + TARGET_READ_OFFSET + 1) * sizeof(char));

  fpt = fopen( "targets.fa", "w" );
  fpr = fopen( "queries.fa", "w" );
  fpg = fopen( "errors.txt", "w" );

  reserved = (char *)calloc(seq_len, sizeof(char));

  for (i = 0; i < cnt_reads; ++ i)
   {
     gaps = gap_sum = snp_sum = 0;
     do
      {
        pos = rand() % seq_len;
        
      } while (reserved[pos] || ((pos > seq_len - nr_gaps - 1) && (calculate_length(gm,seq_len,pos) <= read_len)));
     reserved[pos] = 1;
     //printf ("NEW POS at %d read_len: %d seq: %d remain: %d\n", pos, read_len, seq_len, calculate_length(gm,seq_len,pos));

     strncpy(target, seq + pos, read_len + TARGET_READ_OFFSET);
     target[read_len + TARGET_READ_OFFSET] = 0;

     j = 0; m = 0;
     while (j < read_len)
      {
  //      printf ("MISTAKE: seq: %d pos: %d j: %d pos+j: %d m: %d\n", seq_len, pos, j, pos + j, m);
        if (gm[pos + m])
         {
      //     printf ( "HER!\n");
           ++ gaps;
           gap_sum += gm[pos + m]->len;

           if (gm[pos + m]->type == GAP_TYPE_INSERTION)
            {
        //      printf ("ENTER\n");
            //     printf ("gmlen: %d %d %d\n", gm[pos + m]->len, j, rand() % DNA_ALPH_SIZE);
              for (k = 0; k < MIN(gm[pos + m]->len, read_len - j); ++ k)
               {
                 read[j++] = DNA[rand() % DNA_ALPH_SIZE];
               }
              if (j < read_len)
               {
                 if (sm[pos + m])     /* is there an snp */
                  {
                    read[j++] = sm[pos + m];
                    ++snp_sum;
                  }
                 else
                  {
                    read[j++] = seq[pos + m];
                  }
               }
              ++m;
          //    printf ("EXIT\n");
            }
           else
            {
              m += gm[pos + m]->len;
            }
         }
        else
         {
           if (sm[pos + m])   /* is there an snp */
            {
              read[j++] = sm[pos + m++];
              ++snp_sum;
            }
           else
            {
              read[j++] = seq[pos + m++];
            }
         }
      }
     read[read_len] = 0;
     //printf ("%10d read: %s\n", pos, read);
     //m++;
     //if(m==2)  exit(1);

     fprintf(fpt, ">%s_%d\n%s\n", "target", i + 1, target);
     fprintf(fpr, ">%s_%d\n%s\n", "read", i + 1, read);
     fprintf(fpg, "%d,%d,%d\n", gaps, gap_sum, snp_sum);
	
     if ( gaps > 0 )  ( * cnt_gaps ) ++;	
   }

  free(reserved);
  free(read);
  free(target);

  fclose(fpr);
  fclose(fpg);
  fclose(fpt);
  
}

void filter_unknown(char * seq)
 {
   int                  n, i;

   n = strlen(seq);
   
   for (i = 0; i < n; ++ i)
    {
      if (seq[i] == 'N' || seq[i] == 'n')
       {
         seq[i] = DNA[rand() % DNA_ALPH_SIZE];
       }
    }
 }

int main(int argc, char *argv[])
{
  struct read_info    * ri;
  struct gap_info    ** gm;
  char                * sm;
  int                   seq_len;
  int                   new_seq_len;
  int                   gaps;
  int                   cnt_reads;
  int                   read_size;
  int                   cnt_gaps;

  if ( argc != 4)
   {
     fprintf(stderr, "syntax: %s [REFERENCE] [NUMBER-OF-READS] [READS-SIZE]\n", argv[0]);
     return (1);
   }
  cnt_reads = atoi(argv[2]);
  read_size = atoi(argv[3]);

  ri = read_fasta(argv[1], &seq_len);
  if (seq_len < SEQ_LEN_LIMIT_LOWER)
   {
     fprintf(stderr, "Sequence size must be at least %d bp\n", SEQ_LEN_LIMIT_LOWER);
     return (1);
   }

  printf ("Sequence Length: %d\n", seq_len);

  srand(666);
  filter_unknown(ri->seq);
  gaps = (int)ceil(seq_len * GAP_RATE);
  gm   = random_gaps(seq_len, gaps, &new_seq_len);
  sm   = random_snp(ri->seq, seq_len, (int)ceil(seq_len * SNP_RATE), gm);
  printf ("Total snp: %d\n", (int)ceil(seq_len * SNP_RATE));

  printf ("New sequence length: %d\n", new_seq_len);
  cnt_gaps = 0;
  compute_reads(cnt_reads, gm, ri->seq, seq_len, read_size, gaps, sm, &cnt_gaps);
  printf ("Gaps in the reads: %lf\n", (double) cnt_gaps/cnt_reads );

  return (0);
}
