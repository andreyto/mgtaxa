//  A. L. Delcher
//
//  File:  glimmer3.cc
//
//  Last Modified:  Tue May  9 10:25:40 EDT 2006
//
//  This program finds open reading frames in the file named
//  on the command line and scores them using the probability
//  model in the file indicated by the second command-line
//  parameter.
//
//  Copyright (c) 2006 University of Maryland Center for Bioinformatics
//  & Computational Biology



#include  "glimmer3.hh"


static int  For_Edwin = 0;


// External variables

extern int  Verbose;
extern int  Global_Debug_Flag;


// Global variables

static bool  Allow_Truncated_Orfs = false;
  // If set true by -X option, then score orfs that
  // extend to the end of the sequence
static Event_Node_t  * Best_Event [6];
  // Best parse event up to the current point in each reading frame
static string  Command_Line;
  // Command, options and parameters that invoked the program
static vector <double>  Cumulative_Score [6];
  // Prefix-sum score at each position of the input sequence
  // in each reading frame, plus the independent model
  // Frames are, in order:  +1, +2, +3, -1, -2, -3, ind
static const char  * Fasta_Header;
  // Header on first line of fasta input file
static Event_Node_t  First_Event, Final_Event;
  // First and last nodes in DAG of possible parse events
static vector <Codon_t>  Fwd_Start_Pattern;
  // Bit patterns representing possible forward start codons
static vector <Codon_t>  Fwd_Stop_Pattern;
  // Bit patterns representing possible forward stop codons
static bool  GC_Frac_Set = false;
  // If true, then  Indep_GC_Frac  is set by -C option; otherwise,
  // it is determined from the input sequence data.
static int  Genbank_Xlate_Code = 0;
  // Holds the Genbank translation table number that determines
  // stop codons and codon translation.
static ICM_t  Gene_ICM;
  // The interpolated context model (ICM) of the coding
  // part of genes.
static int  Gene_ID_Ct = 0;
  // Counter used to assign ID numbers to tentative genes
static bool  Genome_Is_Circular = DEFAULT_GENOME_IS_CIRCULAR;
  // If true, input sequences are assumed to be circularly connected
  // so genes will be allowed to wrap around the end
static char  * ICM_File_Name = NULL;
  // Name of the file containing the probability model
static char  * Ignore_File_Name = NULL;
  // Name of file containing list of regions that cannot be included
  // in gene predictions
static int  Ignore_Score_Len = INT_MAX;
  // Genes at least this long do not count the independent model
  // in their score
static vector <Range_t>  Ignore_Region;
static double  Indep_GC_Frac = -1.0;
  // GC proportion used in simple independent model.
  // Set from counts of input sequences or by -C option
static ICM_t  Indep_Model (3, 2, 3);;
  // The ICM for an independent model of bases, based on GC-percentage
  // but without in-frame stop codons
static Event_Node_t  * Last_Event [6];
  // Last parse event up to the current point in each reading frame
static PWM_t  LogOdds_PWM;
  // Log odds wrt background gc-fraction of  Ribosome_PWM .
static int  Min_Gene_Len = DEFAULT_MIN_GENE_LEN;
  // Shortest (in nucleotides) gene that will be considered for scoring
static int  Max_Olap_Bases = DEFAULT_MAX_OLAP_BASES;
  // Overlaps of this many or fewer bases are allowed between adjacent
  // genes
static double  Neg_Entropy_Profile [20] = DEFAULT_NEG_ENTROPY_PROF;
  // Entropy distribution of amino-acids in non-genes
static int  Num_Start_Codons;
  // Number of different start codon patterns
static int  Num_Stop_Codons;
  // Number of different stop codon patterns
static char  * Orflist_File_Name = NULL;
  // Name of file containing a list of regions (which should be valid
  // orfs) that will be scored separately with no overlap rules
static char  * Output_Tag = NULL;
  // Prefix used for output files
static double  Pos_Entropy_Profile [20] = DEFAULT_POS_ENTROPY_PROF;
  // Entropy distribution of amino-acids in genes
static vector <Codon_t>  Rev_Start_Pattern;
  // Bit patterns representing possible reverse start codons
static vector <Codon_t>  Rev_Stop_Pattern;
  // Bit patterns representing possible reverse stop codons
static PWM_t  Ribosome_PWM;
  // Position weight matrix for the ribosome binding pattern
static int  Ribosome_Window_Size = DEFAULT_RIBOSOME_WINDOW_SIZE;
  // Width of window before starts in which to look for matches to
  //  Ribosome_PWM .
static bool  Separate_Orf_Input = false;
  // If set true by -M option then input is multifasta file
  // of orfs to be scored separately (like Orflist_Option)
static vector <Orf_Pos_t>  Orf_Pos_List;
  // List of orfs specified by the -L option to be scored separatedly
static string  Sequence;
  // The input sequence to be scored.
static int  Sequence_Ct;
  // The number of sequences in the input fasta file
static char  * Sequence_File_Name = NULL;
  // Name of the input sequence file
static int  Sequence_Len;
  // Length of genomic sequence string being processed.
static vector <const char *>  Start_Codon;
  // Sequences assumed to be start codons
static vector <double>  Start_Prob;
  // Probability of occurrence of start codons
static vector <const char *>  Stop_Codon;
  // Sequences assumed to be stop codons
static string  Tag;
  // The fasta-header lines of the sequence in  Sequence
static int  Threshold_Score = DEFAULT_THRESHOLD_SCORE;
  // Minimum score for an orf to be considered a potential gene
static bool  Use_Entropy_Profiles = false;
  // If set true (by the -E option) then show the entropy distance
  // ratio in the output.
static bool  Use_First_Start_Codon = DEFAULT_USE_FIRST_START_CODON;
  // If true, automatically use the earliest start codon in a gene;
  // otherwise, try to choose the best start codon
static bool  Use_Independent_Score = DEFAULT_USE_INDEPENDENT_SCORE;
  // If true, let the non-Markov independent model compete with
  // the periodic Markov models to score genes.
static bool  Use_PWM = false;
  // If set true (by the -b option), use the PWM matrix read in
  // to help find gene starts.



int  main
    (int argc, char * argv [])

  {
   FILE  * sequence_fp, * detail_fp, * predict_fp;
   vector <string>  seq_list, hdr_list;
   vector <Orf_t>  orf_list;
   vector <Gene_t>  gene_list;
   string  hdr, filename;
   time_t  now;
   int  i;

   try
     {
      now = time (NULL);
      cerr << "Starting at " << ctime (& now) << endl;

      Verbose = 0;

      Parse_Command_Line (argc, argv);

      if  (Ignore_File_Name != NULL)
          Get_Ignore_Regions ();

      if  (Orflist_File_Name != NULL)
          Get_Orf_Pos_List ();

      Set_Start_And_Stop_Codons ();

      if  (GC_Frac_Set)
          {
           Indep_Model . Build_Indep_WO_Stops (Indep_GC_Frac, Stop_Codon);
           Set_Ignore_Score_Len ();
          }

      filename = Output_Tag;
      filename . append (".detail");
      detail_fp = File_Open (filename, "w", __FILE__, __LINE__);

      filename = Output_Tag;
      filename . append (".predict");
      predict_fp = File_Open (filename, "w", __FILE__, __LINE__);

      sequence_fp = File_Open (Sequence_File_Name, "r", __FILE__, __LINE__);
      Gene_ICM . Read (ICM_File_Name);

      Read_Sequences (sequence_fp, seq_list, hdr_list, Sequence_Ct);
      fclose (sequence_fp);

      if  (! GC_Frac_Set)
          {
           Set_GC_Fraction (Indep_GC_Frac, seq_list);
           Indep_Model . Build_Indep_WO_Stops (Indep_GC_Frac, Stop_Codon);
           Set_Ignore_Score_Len ();
          }

      Echo_General_Settings (stderr);
      fprintf (detail_fp, "Command:  %s\n\n", Command_Line . c_str ());
      Echo_General_Settings (detail_fp);


      Prob_To_Logs (Start_Prob);
      if  (Use_PWM)
          {
           LogOdds_PWM = Ribosome_PWM;
           LogOdds_PWM . Make_Log_Odds_WRT_GC (Indep_GC_Frac);
          }

      if  (Separate_Orf_Input)
          Print_Orflist_Headings (detail_fp);

      for  (i = 0;  i < Sequence_Ct;  i ++)
        {
         if  (Separate_Orf_Input)
             {
              Score_Separate_Input (seq_list [i], hdr_list [i], i,
                   detail_fp, predict_fp);
              continue;
             }

         Sequence = seq_list [i];
         Sequence_Len = Sequence . length ();

         Fasta_Header = hdr_list [i] . c_str ();
         fprintf (detail_fp, "\n\n>%s\n", Fasta_Header);
         Echo_Specific_Settings (detail_fp, Sequence_Len);
         fprintf (predict_fp, ">%s\n", Fasta_Header);

         if  (Orflist_File_Name != NULL)
            {
             Print_Orflist_Headings (detail_fp);
             Score_Orflist (detail_fp, predict_fp);
             break;
            }

         Initialize_Terminal_Events (First_Event, Final_Event, Best_Event,
              Last_Event);

         Print_Headings (detail_fp);

         cerr << "Analyzing Sequence #" << i + 1 << endl;
         cerr << "Start Find_Orfs" << endl;
         Find_Orfs (orf_list);

         cerr << "Start Score_Orfs" << endl;
         Score_Orfs (orf_list, gene_list, detail_fp);

         if  (Verbose > 1)
             Show_Events (stdout);

         cerr << "Start Process_Events" << endl;
         Process_Events ();
         Set_Final_Event (Final_Event, Best_Event, Sequence_Len);

         cerr << "Start Trace_Back" << endl;
         Trace_Back (predict_fp, Final_Event);

         gene_list . clear ();
         orf_list . clear ();

         Clear_Events ();
        }

      fclose (detail_fp);
      fclose (predict_fp);
     }
   catch (std :: exception & e)
     {
      cerr << "** Standard Exception **" << endl;
      cerr << e << endl;
      exit (EXIT_FAILURE);
     }

   return  0;
  }



static void  Add_Events
    (const Orf_t & orf, vector <Start_t> & start_list, int id)

//  Add events for  orf  with possible coding start sites in
//   start_list  to global  Last_Event .   id  is the id number
//  of this orf (which corresponds to the numbers in the detail
//  file.

  {
   Event_Node_t  * ne;   // new event
   double  sc;
   int  fr, sub;
   int  i, n;

   fr = orf . Get_Frame ();
   n = start_list . size ();

   if  (orf . Get_Orf_Len () >= Ignore_Score_Len)
       { // artificially inflate start scores
        for  (i = 0;  i < n;  i ++)
          {
           sc = LONG_ORF_SCORE_PER_BASE * (1 + start_list [i] . j);
           if  (sc > start_list [i] . score)
               start_list [i] . score = sc;
          }
       }

   n = start_list . size ();

   if  (fr > 0)
       {
        sub = fr - 1;

        for  (i = 0;  i < n;  i ++)
          if  (1 + start_list [i] . j >= Min_Gene_Len)
              {
               ne = new Event_Node_t;
               ne -> e_type = FWD_START;
               ne -> id = id;
               ne -> pos = start_list [i] . pos + 2;
                 // event pos is last base of codon; start pos is first
               PWM_Score_Fwd_Start (start_list [i] . pos, LogOdds_PWM,
                    Ribosome_Window_Size, ne -> pwm_score, ne -> pwm_sep);
               ne -> frame = fr;
               ne -> score = start_list [i] . score;
               Add_PWM_Score (ne);
               if  (start_list [i] . which >= 0)
                   ne -> score += Start_Prob [start_list [i] . which];
                 // which will be -1 for truncated orfs with an
                 // artificial start at the start
               ne -> is_first_start = start_list [i] . first;
               ne -> truncated = start_list [i] . truncated;
               ne -> best_pred = NULL;
               ne -> frame_pred = Last_Event [sub];
               Last_Event [sub] = ne;
              }

        ne = new Event_Node_t;
        ne -> e_type = FWD_STOP;
        ne -> id = id;
        ne -> pos = orf . Get_Stop_Position () + 2;
          // event pos is last base of codon; orf pos is first
        ne -> frame = fr;
        ne -> is_first_start = false;
        ne -> truncated = false;
        ne -> score = 0.0;
        ne -> best_pred = NULL;
        ne -> frame_pred = Last_Event [sub];
        Last_Event [sub] = ne;
       }
     else
       {
        sub = 2 - fr;

        ne = new Event_Node_t;
        ne -> e_type = REV_STOP;
        ne -> id = id;
        ne -> pos = orf . Get_Stop_Position () + 2;
          // event pos is last base of codon; orf pos is first
        ne -> frame = fr;
        ne -> is_first_start = false;
        ne -> truncated = false;
        ne -> score = 0.0;
        ne -> best_pred = NULL;
        ne -> frame_pred = Last_Event [sub];
        Last_Event [sub] = ne;

        for  (i = 0;  i < n;  i ++)
          if  (1 + start_list [i] . j >= Min_Gene_Len)
              {
               ne = new Event_Node_t;
               ne -> e_type = REV_START;
               ne -> id = id;
               ne -> pos = start_list [i] . pos;
                 // both pos's are last base of codon, i.e., highest coord
               PWM_Score_Rev_Start (start_list [i] . pos, LogOdds_PWM,
                    Ribosome_Window_Size, ne -> pwm_score, ne -> pwm_sep);
               ne -> frame = fr;
               ne -> score = start_list [i] . score;
               Add_PWM_Score (ne);
               if  (start_list [i] . which >= 0)
                   ne -> score += Start_Prob [start_list [i] . which];
                 // which will be -1 for truncated orfs with an
                 // artificial start at the start
               ne -> is_first_start = start_list [i] . first;
               ne -> truncated = start_list [i] . truncated;
               ne -> best_pred = NULL;
               ne -> frame_pred = Last_Event [sub];
               Last_Event [sub] = ne;
              }
       }

   return;
  }



static void  Add_PWM_Score
    (Event_Node_t * p)

//  Add all or part of  p -> pwm_score  to  p -> score  depending
//  on the location of the PWM match.

  {
   static const int  LO_SEP = 4, HI_SEP = 10, HI_TAIL = 6;
   double  coeff;

   if  (p -> pwm_score < 0.0)
       return;

   // Use all the pwm_score if the pwm_sep is between LO_SEP and HI_SEP
   // Otherwise, use a fraction of it.
   if  (p -> pwm_sep < LO_SEP)
       coeff = double (p -> pwm_sep) / LO_SEP;
   else if  (p -> pwm_sep <= HI_SEP)
       coeff = 1.0;
   else if  (p -> pwm_sep < HI_SEP + HI_TAIL)
       coeff = double (HI_SEP + HI_TAIL - p -> pwm_sep) / HI_TAIL;
     else
       coeff = 0.0;

   if  (0.0 < coeff)
       p -> score += coeff * p -> pwm_score;

   return;
  }



static void  All_Frame_Score
    (const string & s, int len, int frame, vector <double> & af)

//  Score the first  len  characters of string  s  in all six reading
//  frames using global model  Gene_ICM .   frame  is the
//  frame position in the original genome of the first character of
//   s , where frame positions are numbered  1,2,3,1,2,3  starting
//  with the first character of the genome.  frame also has the
//  direction of the gene in the genome string.
//  **NOTE**  s  is the reverse (but not complemented) of the gene.
//  Store the results in  af  where the order of reading frames
//  is  +1,+2,+3,-1,-2,-3 .   af  is assumed to be large enough
//  to hold the results.

  {
   string  rev_compl;
   const char  * cstr = s . c_str ();

   af [0] = Gene_ICM . Score_String (cstr, len, 1);
   af [1] = Gene_ICM . Score_String (cstr, len, 2);
   af [2] = Gene_ICM . Score_String (cstr, len, 0);

   Reverse_Complement_Transfer (rev_compl, s, 0, len);

   af [3] = Gene_ICM . Score_String (rev_compl . c_str (), len, 1);
   af [4] = Gene_ICM . Score_String (rev_compl . c_str (), len, 0);
   af [5] = Gene_ICM . Score_String (rev_compl . c_str (), len, 2);

   Permute_By_Frame (af, frame);

   return;
  }



static void  Clear_Events
    (void)

//  Free memory in chains pointed to by  Last_Event .  Note that
//  the initial event is not dynamically allocated (it's the global
//  variable  First_Event ) so it is not cleared.

  {
   Event_Node_t  * p, * q;
   int  i;

   for  (i = 0;  i < 6;  i ++)
     for  (p = Last_Event [i];  p != NULL && p -> e_type != INITIAL;  p = q)
       {
        q = p -> frame_pred;
        delete p;
       }

   return;
  }



static void  Complement_Transfer
    (string & buff, const string & s, int start, int len)

//  Copy to string  buff  the substring of  s  starting at subscript
//   start  and going to the right for a length of  len .  Wraparound
//  the end of  s  if necessary.  Convert each character to its
//  Watson-Crick complement as it is copied.

  {
   int  j, n;

   n = s . length ();
   assert (start < n);
   assert (0 <= len);

   buff . resize (len);
   for  (j = 0;  j < len;  j ++, start ++)
     {
      if  (start >= n)
          start -= n;
      buff [j] = Complement (s [start]);
     }

   return;
  }



static void  Disqualify
    (Event_Node_t * p, int cutoff)

//  Set the  disqualified  bit true for nodes reachable from
//   p  by  best_pred  pointers that have  pos  values at least
//  as great as  cutoff .

  {
   Event_Node_t  * q;

   if  (p == NULL)
       return;

   for  (q = p -> best_pred;  q != NULL && cutoff <= q -> pos;  q = q -> best_pred)
     q -> disqualified = true;

   return;
  }



static void  Do_Fwd_Stop_Codon
    (int i, int frame, int prev_fwd_stop [3], int first_fwd_start [3],
     int first_fwd_stop [3], int first_base, bool hit_ignore,
     vector <Orf_t> & orf_list)

//  Create a new entry for the forward orf ending at sequence subscript  i
//  and add it to  orf_list , if it's sufficiently long.   frame  is
//  the reading frame subscript of this orf.   prev_fwd_stop  indicates
//  the location of the previous forward stop codons.   first_fwd_start
//  has the locations of the first start codon for the current forward
//  reading frames.  Set  first_fwd_stop  to this position if there
//  is no prior stop in this reading frame.   first_base  is the position
//  of the first sequence base after an ignore region, or the start of
//  the sequence if no ignore regions have been encountered, which is
//  indicated by  hit_ignore.

  {
   Orf_t  orf;
   int  gene_len, orf_len;

   if  (prev_fwd_stop [frame] == 0)
       {
        Handle_First_Forward_Stop (frame, i - 1, first_fwd_start [frame],
             first_base, gene_len, orf_len,
             Genome_Is_Circular && ! hit_ignore);
        first_fwd_stop [frame] = i - 1;
       }
     else
       {
        gene_len = i - first_fwd_start [frame] - 1;
        orf_len = i - prev_fwd_stop [frame] - 4;
       }

   if  (gene_len >= Min_Gene_Len)
       {
        orf . Set_Stop_Position (i - 1);
        orf . Set_Frame (1 + (frame + 1) % 3);
        orf . Set_Gene_Len (gene_len);
        orf . Set_Orf_Len (orf_len);
        orf_list . push_back (orf);
       }

   first_fwd_start [frame] = INT_MAX;
   prev_fwd_stop [frame] = i - 1;

   return;
  }



static void  Echo_General_Settings
    (FILE * fp)

//  Output values of global variables and parameter settings
//  to  fp .

  {
   int  i, n;

   fprintf (fp, "Sequence file = %s\n", Sequence_File_Name);
   fprintf (fp, "Number of sequences = %d\n", Sequence_Ct);
   fprintf (fp, "ICM model file = %s\n", ICM_File_Name);
   fprintf (fp, "Excluded regions file = %s\n",
        Printable (Ignore_File_Name));
   fprintf (fp, "List of orfs file = %s\n",
        Printable (Orflist_File_Name));

   fprintf (fp, "Input %s separate orfs\n",
        Separate_Orf_Input ? "is" : "is NOT");
   fprintf (fp, "Independent (non-coding) scores %s used\n",
        Use_Independent_Score ? "are" : "are NOT");
   if  (! Separate_Orf_Input)
       {
        fprintf (fp, "Circular genome = %s\n", Printable (Genome_Is_Circular));
       }
   if  (! Separate_Orf_Input && Orflist_File_Name == NULL)
       {
        fprintf (fp, "Truncated orfs = %s\n", Printable (Allow_Truncated_Orfs));
        fprintf (fp, "Minimum gene length = %d bp\n", Min_Gene_Len);
        fprintf (fp, "Maximum overlap bases = %d\n", Max_Olap_Bases);
        fprintf (fp, "Threshold score = %d\n", Threshold_Score);
        fprintf (fp, "Use first start codon = %s\n",
             Printable (Use_First_Start_Codon));
        if  (Genbank_Xlate_Code != 0)
            fprintf (fp, "Translation table = %d\n", Genbank_Xlate_Code);
        fprintf (fp, "Start codons = ");
        Print_Comma_Separated_Strings (Start_Codon, fp);
        fputc ('\n', fp);
        fprintf (fp, "Start probs = ");
        n = Start_Prob . size ();
        for  (i = 0;  i < n;  i ++)
          {
           if  (i > 0)
               fputc (',', fp);
           fprintf (fp, "%.3f", Start_Prob [i]);
          }
        fputc ('\n', fp);
        fprintf (fp, "Stop codons = ");
        Print_Comma_Separated_Strings (Stop_Codon, fp);
        fputc ('\n', fp);
       }

   fprintf (fp, "GC percentage = %.1f%%\n", 100.0 * Indep_GC_Frac);
   if  (Use_Independent_Score)
       fprintf (fp, "Ignore score on orfs longer than %s\n",
            Num_Or_Max (Ignore_Score_Len));

   return;
  }



static void  Echo_Specific_Settings
    (FILE * fp, int len)

//  Output values of variables an settings that depend on the
//  current input string, which has length  len .

  {
   fprintf (fp, "Sequence length = %d\n", len);

   return;
  }



static double  Entropy_Distance_Ratio
    (int start, int len, int fr)

//  Return the distance ratio for the entropy profile for the
//  gene starting at position  start  (in 1-based coordinates)
//  on global  Sequence with length  len  and in reading frame  fr .
//  The ratio is the distance to global  Pos_Entropy_Profile  over
//  the distance to global  Neg_Entropy_Profile .

  {
   string  buff;
   int  count [26] = {0};
   double  ep [20];
   double  pos_dist, neg_dist, ratio;
   char  aa;
   int  i;

   if  (fr > 0)
       Forward_Strand_Transfer (buff, Sequence, On_Seq_0 (start - 1), len);
     else
       Reverse_Strand_Transfer (buff, Sequence, On_Seq_0 (start - 1), len);

   for  (i = 0; i < len;  i += 3)
     {
      aa = Codon_Translation (buff . c_str () + i, Genbank_Xlate_Code);
      if  (aa != '*')
          count [aa - 'A'] ++;
     }
   Counts_To_Entropy_Profile (count, ep);

   pos_dist = neg_dist = 0.0;
   for  (i = 0;  i < 20;  i ++)
     {
      pos_dist += pow (ep [i] - Pos_Entropy_Profile [i], 2);
      neg_dist += pow (ep [i] - Neg_Entropy_Profile [i], 2);
     }

   pos_dist = sqrt (pos_dist);
   neg_dist = sqrt (neg_dist);
   if  (neg_dist == 0.0)
       {
        if  (pos_dist == 0.0)
            ratio = 1.0;
          else
            ratio = 1e3;
       }
     else
       ratio = pos_dist / neg_dist;

   return  ratio;
  }



static int  Find_Uncovered_Position
    (vector <Event_Node_t *> ep)

//  Find a position in  ep  that is not covered by any potential
//  gene, if possible.  If the first gene does not overlap the
//  last gene, then return  0 .  Also return  0  if there is
//  no uncovered position.  The position is regarded as being
//  between bases, and positions are numbered from  0  to  Sequence_Len .

  {
   int  cover_ct, zero_pos;
   int  first_pos, last_pos;
   int  i, n;

   n = ep . size ();

   if  (n <= 1)
       return  0;

   // ep is already sorted ascending by position and the initial
   // event is first in it
   first_pos = ep [1] -> pos - 3;  // between position in front of codon
   last_pos = ep [n - 1] -> pos - Sequence_Len;
     // between position after codon normalized to wrapped front position

   if  (last_pos <= first_pos)
       return  0;  // no overlap between front and back

   cover_ct = 0;
   zero_pos = ep [n - 1] -> pos;
   for  (i = 1;  i < n;  i ++)
     switch  (ep [i] -> e_type)
       {
        case  FWD_START :
          if  (ep [i] -> is_first_start)
              {
               cover_ct ++;
               if  (cover_ct == 1 && 3 <= ep [i] -> pos - zero_pos)
                   {
                    assert (zero_pos >= 1);
                    return  zero_pos;
                   }
              }
          break;

        case  FWD_STOP :
          cover_ct --;
          if  (cover_ct == 0)
              zero_pos = ep [i] -> pos;
          break;

        case  REV_START :
          if  (ep [i] -> is_first_start)
              {
               cover_ct --;
               if  (cover_ct == 0)
                   zero_pos = ep [i] -> pos;
              }
          break;

        case  REV_STOP :
          cover_ct ++;
          if  (cover_ct == 1 && 3 <= ep [i] -> pos - zero_pos)
              {
               assert (zero_pos >= 1);
               return  zero_pos;
              }
          break;

        case  INITIAL :
        case  TERMINAL :
        default :
          sprintf (Clean_Exit_Msg_Line, "ERROR:  Unexpected event type = %s",
               Print_String (ep [i] -> e_type));
          SIMPLE_THROW (Clean_Exit_Msg_Line);
       }

   return  0;
  }



static void  Find_Orfs
    (vector <Orf_t> & orf_list)

//  Put in  orf_list  all sufficiently long orfs in global
//  string  Sequence .

  {
   Orf_t  orf;
   Codon_t  codon;

   // Positions stored in these are the first (i.e., lowest-subscript)
   // base of the codon, using positions starting at 1.
   int  first_fwd_start [3] = {INT_MAX, INT_MAX, INT_MAX};
   int  last_rev_start [3] = {0};
   int  prev_fwd_stop [3] = {0}, prev_rev_stop [3] = {0};
   int  first_fwd_stop [3] = {0};
        // Used for wraparound in circular genomes
   int  ignore_start, ignore_stop;
        // indicate next beginning and ending positions of next
        // region to be ignored
   int  ignore_ct;
        // number of ignore regions
   int  ignore_sub;
        // subscript of current ignore region
   bool  hit_ignore = false;
        // indicates if any ignore region has been reached yet
   bool  ignoring = false;
        // indicates current status of ignore region
   int  first_base = 1;
        // position of the first base in the current region being
        // processed
   int  frame, gene_len;
        // frame subscripts are 0, 1, 2 for both forward and reverse
        // events.  The frame is based on the *LAST* (i.e., highest-subscript)
        // base of the codon, using positions starting at 0
   int  i, j, n;

   orf_list . clear ();
   n = Sequence_Len;

   if  (n < Min_Gene_Len)
       return;

   if  (Genome_Is_Circular)
       {
        // allow 2-base overhang to catch start and stop codons that
        // span the end of  Sequence
        n += 2;
        Sequence . push_back (Sequence [0]);
        Sequence . push_back (Sequence [1]);
       }
 
   if  (Ignore_Region . size () == 0)
       ignore_start = ignore_stop = INT_MAX;
     else
       {
        ignore_ct = Ignore_Region . size ();
        ignore_start = Ignore_Region [0] . lo;
        ignore_stop = Ignore_Region [0] . hi;
        ignore_sub = 0;
       }

   frame = 0;
   for  (i = 0;  i < n;  i ++)
     {
      // check if this position is the boundary of an ignore region
      if  (i == ignore_start)
          {
           Finish_Orfs (false, prev_rev_stop, last_rev_start, i, orf_list);
           hit_ignore = ignoring = true;
          }
      else if  (i == ignore_stop)
          {
           // reset saved positions to their initial values as if the
           // start of the genome
           for  (j = 0;  j < 3;  j ++)
             {
              first_fwd_start [j] = INT_MAX;
              last_rev_start [j] = 0;
              prev_fwd_stop [j] = 0;
              prev_rev_stop [j] = 0;
             }
           codon . Clear ();
           first_base = i + 1;
           ignoring = false;
           ignore_sub ++;
           if  (ignore_sub >= ignore_ct)
               ignore_start = ignore_stop = INT_MAX;
             else
               {
                ignore_start = Ignore_Region [ignore_sub] . lo;
                ignore_stop = Ignore_Region [ignore_sub] . hi;
               }
          }

      if  (! ignoring)
          {
           int  which, orf_stop;

           codon . Shift_In (Sequence [i]);

           if  (codon . Can_Be (Fwd_Start_Pattern, which)
                   && first_fwd_start [frame] == INT_MAX)
               first_fwd_start [frame] = i - 1;

           if  (codon . Can_Be (Rev_Start_Pattern, which))
               {
                last_rev_start [frame] = i - 1;
               }

           if  (codon . Must_Be (Fwd_Stop_Pattern, which))
               Do_Fwd_Stop_Codon (i, frame, prev_fwd_stop, first_fwd_start,
                    first_fwd_stop, first_base, hit_ignore, orf_list);

           if  (codon . Must_Be (Rev_Stop_Pattern, which))
               {
                if  (prev_rev_stop [frame] == 0)
                    Handle_First_Reverse_Stop (i - 1, last_rev_start [frame],
                         gene_len, orf_stop, hit_ignore);
                  else
                    {
                     orf_stop = prev_rev_stop [frame];
                     gene_len = last_rev_start [frame] - orf_stop;
                    }

                if  (gene_len >= Min_Gene_Len)
                    {
                     orf . Set_Stop_Position (orf_stop);
                     orf . Set_Frame (-1 - (frame + 1) % 3);
                     orf . Set_Gene_Len (gene_len);
                     orf . Set_Orf_Len (i - orf_stop - 4);
                     orf_list . push_back (orf);
                    }
                last_rev_start [frame] = 0;
                prev_rev_stop [frame] = i - 1;
               }
          }

      if  (frame == 2)
          frame = 0;
        else
          frame ++;
     }

   Finish_Orfs (Genome_Is_Circular, prev_rev_stop, last_rev_start,
        Sequence_Len, orf_list);

   if  (Genome_Is_Circular)
       Sequence . resize (Sequence_Len);
   else if  (Allow_Truncated_Orfs)
       // Treat 3 bp past the end of the sequence as stop codons
       for  ( ;  i < n + 3;  i ++)
         {
          if  (! ignoring)
              Do_Fwd_Stop_Codon (i, frame, prev_fwd_stop, first_fwd_start,
                   first_fwd_stop, first_base, hit_ignore, orf_list);

          if  (frame == 2)
              frame = 0;
            else
              frame ++;
         }

   return;
  }



static void  Find_Stops_Reverse
    (const string & s, int len, vector <bool> & has_stop)

//  Set  has_stop [i]  to true iff string  s  has a
//  stop codon in the frame corresponding to  i .
//  The order of frames is  +1,+2,+3,-1,-2,-3 .
//  Use only the first  len  characters of  s .
//   s  is the reverse (but not complemented) of the DNA strand
//  Automatically set  has_stop [6]  to  false, representing the
//  independent model "frame".

  {
   Codon_t  codon;
   int  frame_ss;    // frame subscript
   int  which;
   int  i;

   has_stop . resize (7);
   for  (i = 0;  i < 7;  i ++)
     has_stop [i] = false;

   frame_ss = 1;

   for  (i = len - 1;  i >= 0;  i --)
     {
      codon . Shift_In (s [i]);

      if  (codon . Must_Be (Fwd_Stop_Pattern, which))
          has_stop [frame_ss] = true;
      if  (codon . Must_Be (Rev_Stop_Pattern, which))
          has_stop [frame_ss + 3] = true;

      if  (frame_ss == 2)
          frame_ss = 0;
        else
          frame_ss ++;
     }

   return;
  }



static void  Finish_Orfs
    (bool use_wraparound, const int prev_rev_stop [3],
     const int last_rev_start [3], int last_position,
     vector <Orf_t> & orf_list)

//  Finish reverse-strand orfs because we've hit the end of the
//  genome (or hit an ignore region).  If  use-wraparound  is true,
//  then the orfs can wrap around the end of the (circular) genome;
//  otherwise, not.   prev_rev_stop  has the position of the last-seen
//  reverse stop codons in each frame, and  last_rev_start  has the
//  position of the last-seen reverse start codons in each frame.
//   last_position  is the last available sequence position to use.
//  Add any suitable orfs to  orf_list .

  {
   Orf_t  orf;
   int  frame, gene_len, orf_len;

   for  (frame = 0;  frame < 3;  frame ++)
     {
      Handle_Last_Reverse_Stop (frame, prev_rev_stop, last_rev_start,
           gene_len, orf_len, use_wraparound, last_position);
      if  (gene_len >= Min_Gene_Len)
          {
           orf . Set_Stop_Position (prev_rev_stop [frame]);
           orf . Set_Frame (-1 - (frame + 1) % 3);
           orf . Set_Gene_Len (gene_len);
           orf . Set_Orf_Len (orf_len);
           orf_list . push_back (orf);
          }
     }

   return;
  }



static void  Fix_Wrap
    (int & p, const int n)

//  Change position  p  so that it falls in the interval  1 .. n
//  where it should be assuming a circular coordinate scheme.

  {
   while  (p < 1)
     p += n;

   while  (p > n)
     p -= n;

   return;
  }



static int  Frame_To_Sub
    (int f)

//  Return the subscript equivalent of frame  f .

  {
   if  (f > 0)
       return  f - 1;
     else
       return  2 - f;
  }



static void  Get_Ignore_Regions
    (void)

//  Read the list of regions from the file with name in global
//   Ignore_File_Name .  Sort them and coalesce overlapping regions.
//  Put the results in global  Ignore_Region .  The format for each
//  line of input is:
//     <lo>  <hi>  <rest of line ignored>  
//  where <lo> and <hi> are integer values.  The region specified
//  is bases <lo>..<hi> inclusive, where bases are numbered starting
//  at 1.  If <hi> is less than <lo> the values are silently swapped.
//  There is no provision for circularity.  If more than one sequence
//  is read in to be searched for genes, these regions will be used
//  to screen them *ALL*, which is very likely not at all what is
//  desired.  Blank lines and lines beginning with # are skipped.

  {
   FILE  * fp;
   char  line [MAX_LINE];
   Range_t  range;
   int  i, j, n, line_ct;

   fp = File_Open (Ignore_File_Name, "r", __FILE__, __LINE__);

   line_ct = 0;
   while  (fgets (line, MAX_LINE, fp) != NULL)
     {
      char  * p;
      int  a, b;

      line_ct ++;

      // set  p  to point to the first non-blank character on the line
      for  (p = line;  * p != '\0' && isspace (* p);  p ++)
        ;
      
      if  (* p == '\0' || * p == '#')
          continue;
      else if  (sscanf (line, "%d %d", & a, & b) == 2)
          {
           if  (a < b)
               {
                range . lo = a - 1;
                  // convert to 0-based between coordinates
                range . hi = b;
               }
             else
               {
                range . lo = b - 1;
                range . hi = a;
               }
           Ignore_Region . push_back (range);
          }
        else
          {
           fprintf (stderr, "ERROR:  Following line %d in file %s is bad--skipped:\n",
                line_ct, Ignore_File_Name);
           fputs (line, stderr);
           fputc ('\n', stderr);
          }
     }

   fclose (fp);

   // sort regions by lo value
   sort (Ignore_Region . begin (), Ignore_Region . end (), Range_Cmp);

   // combine overlapping regions and move them to the front of  Ignore_Region
   n = Ignore_Region . size ();

   if  (n <= 1)
       return;

   for  (i = 0, j = 1;  j < n;  j ++)
     if  (Ignore_Region [j] . lo < Ignore_Region [i] . hi)
         {  // overlap
          if  (Ignore_Region [i] . hi < Ignore_Region [j] . hi)
              Ignore_Region [i] . hi = Ignore_Region [j] . hi;
                 // j extends i to the right
         }
       else
         {
          i ++;
          if  (i != j)
              Ignore_Region [i] = Ignore_Region [j];
                // move j region down to front of list
         }

   Ignore_Region . resize (i + 1);

   return;
  }



static void  Get_Orf_Pos_List
    (void)

//  Read the list of orfs from the file with name in global
//   Orflist_File_Name  and store them in global list
//   Orf_Pos_List .  The format for each
//  line of input is:
//     <tag>  <start>  <stop>  <dir>  <rest of line ignored>  
//  where <start> and <stop> are integer values.  The <stop> position
//  includes the ending stop codon for the orf.  The orf specified
//  is bases <start>..<stop> inclusive, where bases in the input
//  sequence are numbered starting at 1.  <dir> indicates the
//  strand of the gene for cases where it might wraparound the
//  start position of the genome sequence.
//  Blank lines and lines beginning with # are skipped.

  {
   FILE  * fp;
   char  line [MAX_LINE], t [MAX_LINE];
   Orf_Pos_t  orf;
   int  line_ct;

   fp = File_Open (Orflist_File_Name, "r", __FILE__, __LINE__);

   Orf_Pos_List . clear ();
   line_ct = 0;
   while  (fgets (line, MAX_LINE, fp) != NULL)
     {
      char  * p;
      int  a, b, d;

      line_ct ++;

      // set  p  to point to the first non-blank character on the line
      for  (p = line;  * p != '\0' && isspace (* p);  p ++)
        ;
      
      if  (* p == '\0' || * p == '#')
          continue;
      else if  (sscanf (line, "%s %d %d %d", t, & a, & b, & d) == 4)
          {
           orf . tag = strdup (t);
           orf . start = a;
           orf . stop = b;
           orf . dir = d;
           Orf_Pos_List . push_back (orf);
          }
        else
          {
           fprintf (stderr, "ERROR:  Following line %d in file %s is bad--skipped:\n",
                line_ct, Orflist_File_Name);
           fputs (line, stderr);
           fputc ('\n', stderr);
          }
     }

   fclose (fp);

   return;
  }



static void  Handle_First_Forward_Stop
     (int fr, int pos, int start_pos, int first_base, int & gene_len,
      int & orf_len, bool use_wraparound)

//  Handle the case of a forward stop codon, beginning at position
//   pos  in the global  Sequence  (counting starting at 1)  which
//  is in frame subscript  fr  (0, 1 or 2).   start_pos  is the
//  position of the first possible start codon in this frame, or else
//   INT_MAX  if none has been encountered yet.   first_base  is the
//  position of the first base in this region.  Set gene_len
//  to the length of longest possible gene for this orf.  If no gene
//  is possible (e.g., because there is no start codon), then set
//   gene_len  to  0 .  Set  orf_len  to the length of this orf.
//  If  use_wraparound  is true, allow orfs/genes to wrap around
//  through the front of the (circular) sequence.

  {
   if  (use_wraparound)
       {
        Wrap_Through_Front (fr, pos, gene_len, orf_len);
        if  (gene_len == 0 && start_pos != INT_MAX)
            gene_len = pos - start_pos;
       }
     else
       {
        // assume the orf is entirely contained in  Sequence  no
        // matter whether the odd 1 or 2 bases at the front could be
        // a stop or not
        orf_len = pos - first_base;
        orf_len -= orf_len % 3;  // round down
        if  (start_pos == INT_MAX)
            gene_len = 0;
          else
            gene_len = pos - start_pos;
        if  (Allow_Truncated_Orfs && gene_len < Min_Gene_Len)
            gene_len = orf_len;
       }

   return;
  }



static void  Handle_First_Reverse_Stop
    (int pos, int last_start, int & gene_len, int & orf_stop, bool hit_ignore)

//  Set  gene_len  to the length of the reverse-strand gene whose start
//  is at  last_start  (left base of start codon, start-at-1) and which
//  extends off the front of the sequence.  Set  orf_stop  to the first,
//  frame-correct position < 1 where the stop codon (left base) could be.
//  It doesn't matter if the 2nd or 3rd base of this stop codon placement
//  overlaps the beginning of the sequence.
//   pos  is the position (start-at-1 coords) of the right bounding stop
//  codon of this gene.   Set  gene_len  to zero and return, however,
//  if either  hit_ignore  is true or  Allow_Truncated_Orfs  is false.

  {
   if  (hit_ignore || ! Allow_Truncated_Orfs)
       {
        gene_len = 0;
        return;
       }

   orf_stop = pos % 3;
   if  (orf_stop > 0)
       orf_stop -= 3;
   gene_len = last_start - orf_stop;

   return;
  }



static void  Handle_Last_Reverse_Stop
     (int fr, const int prev_rev_stop [3], const int last_rev_start [3],
      int & gene_len, int & orf_len, bool use_wraparound, int last_position)

//  Set  orf_len  and  gene_len  to the length of the last orf, and longest
//  gene in it, resp., in reverse reading frame  fr .
//   prev_rev_stop  has the last stop position in  Sequence  in each
//  reverse reading frame, and  last_rev_start  has the corresponding
//  last start locations.    use_wraparound  indicates whether the
//  orfs are allowed to wrap around the end of the (circular) genome.
//   last_position  is the highest-numbered sequence position available

  {
   if  (prev_rev_stop [fr] == 0)
       {
        // no reverse stop in this frame at all
        gene_len = orf_len = 0;
        return;
       }

   if  (use_wraparound)
       {
        int  wrap_fr;
             // the frame at the front of the genome corresponding
             // to  fr
        wrap_fr = (3 + fr - (Sequence_Len % 3)) % 3;

        Wrap_Around_Back (wrap_fr, prev_rev_stop [fr], gene_len, orf_len);

        if  (gene_len == 0 && last_rev_start [fr] > 0)
            gene_len = last_rev_start [fr] - prev_rev_stop [fr];
       }
     else
       {
        orf_len = last_position - prev_rev_stop [fr] - 2;
             // round down to next multiple of 3
        orf_len -= orf_len % 3;

        if  (last_rev_start [fr] == 0)
            gene_len = 0;
          else
            gene_len = last_rev_start [fr] - prev_rev_stop [fr];
        if  (Allow_Truncated_Orfs && gene_len < Min_Gene_Len)
            gene_len = orf_len;
       }

   assert (orf_len % 3 == 0);
   assert (gene_len % 3 == 0);

   return;
  }



static void  Initialize_Terminal_Events
    (Event_Node_t & first_event, Event_Node_t & final_event,
     Event_Node_t * best_event [6], Event_Node_t * last_event [6])

//  Set up  first_event  and  final_event  and make all
//  entries in  best_event  and  last_event  point to
//   first_event .

  {
   int  i;

   first_event . e_type = INITIAL;
   first_event . pos = 0;
   first_event . score = 0.0;
   first_event . best_pred = NULL;
   first_event . frame_pred = NULL;

   for  (i = 0;  i < 6;  i ++)
     last_event [i] = best_event [i] = & first_event;

   final_event . e_type = TERMINAL;
   final_event . frame_pred = NULL;

   return;
  }



static void  Integerize_Scores
    (const vector <double> ds, int hi_score, const vector <bool> set_negative,
    vector <int> & is)

//  Convert the scores in  ds  to integers ranging from
//   0 .. hi_score  putting the results into  is .
//  Automatically set to  -1  entries corresponding
//  to values in  set_negative  that are true and ignore them
//  in the calculation.

  {
   vector <double>  v;
   double  min, max, sum;
   int  i, n;

   n = ds . size ();
   is . resize (n);
   v . resize (n);

   min = DBL_MAX;
   max = - DBL_MAX;
   for  (i = 0;  i < n;  i ++)
     if  (! set_negative [i])
         {
          if  (ds [i] > max)
              max = ds [i];
          if  (ds [i] < min)
              min = ds [i];
         }

   if  (min < max + MAX_LOG_DIFF)
       min = max + MAX_LOG_DIFF;
   
   sum = 0.0;
   for  (i = 0;  i < n;  i ++)
     if  (set_negative [i])
         v [i] = -1.0;
     else if  (ds [i] < min)
         v [i] = 0.0;
       else
         {
          v [i] = exp (ds [i] - min);
          sum += v [i];
         }

   for  (i = 0;  i < n;  i ++)
     if  (set_negative [i])
         is [i] = -1;
       else
         {
          is [i] = int (HI_SCORE * (v [i] / sum));
          if  (is [i] >= HI_SCORE)
              is [i] = HI_SCORE - 1;
         }

   return;
  }



static double  Olap_Score_Adjustment
    (int lo, int hi, int f1, int f2)

//  Return the larger of the frame  f1  and  frame  f2  scores
//  on the subsequence from  lo .. hi  of global  Sequence .
//   lo  and  hi  are inclusive, start at 1 coordinates.
//  Because wraparounds may have confused the frames, only the
//  sign of the frames is used.   f1  is assumed to be the
//  frame of the beginnning of the subsequence starting on
//  a codon boundary.   f2  is the corresponding frame at the
//  end of the sequence.

  {
   string  buff;
   double  s1, s2;
   int  len, fs;

   len = 1 + hi - lo;
   if  (len < 1)
       return  0.0;

   if  (lo < 1)
       lo += Sequence_Len;
   if  (lo > Sequence_Len)
       lo -= Sequence_Len;
   if  (hi < 1)
       hi += Sequence_Len;
   if  (hi > Sequence_Len)
       hi -= Sequence_Len;

   lo --;  // convert to subscript
   hi --;

   switch  (len % 3)
     {
      case  0 :
        fs = 1;
        break;
      case  1 :
        fs = 0;
        break;
      case  2 :
        fs = 2;
        break;
     }
   // fs is the frame subscript to use in scoring in the direction
   // that does not necessarily start on a codon boundary

   if  (f1 > 0)
       {
        Reverse_Transfer (buff, Sequence, hi, len);
        s1 = Gene_ICM . Score_String (buff . c_str (), len, fs)
               - Indep_Model . Score_String (buff . c_str (), len, fs);
       }
     else
       {
        Complement_Transfer (buff, Sequence, lo, len);
        s1 = Gene_ICM . Score_String (buff . c_str (), len, 1)
               - Indep_Model . Score_String (buff . c_str (), len, 1);
       }

   if  (f1 * f2 < 0)
       Reverse_Complement (buff);

   if  (f2 > 0)
       s2 = Gene_ICM . Score_String (buff . c_str (), len, 1)
              - Indep_Model . Score_String (buff . c_str (), len, 1);
     else
       s2 = Gene_ICM . Score_String (buff . c_str (), len, fs)
              - Indep_Model . Score_String (buff . c_str (), len, fs);

   return  Max (s1, s2);
  }



static int  On_Seq_0
    (int i)

//  Return the subscript equivalent to  i  on a sequence of
//  length  Sequence_Len  (with subscripts starting at 0)
//  assuming circular wraparounds.

  {
   while  (i < 0)
     i += Sequence_Len;
   while  (Sequence_Len <= i)
     i -= Sequence_Len;

   return  i;
  }



static int  On_Seq_1
    (int i)

//  Return the subscript equivalent to  i  on a sequence of
//  length  Sequence_Len  (with subscripts starting at 1)
//  assuming circular wraparounds.

  {
   while  (i < 1)
     i += Sequence_Len;
   while  (Sequence_Len < i)
     i -= Sequence_Len;

   return  i;
  }



static void  Output_Extra_Start_Info
    (FILE * fp, int i, int lo, int hi, int frame,
     vector <Start_t> & start_list)

//  Print to  fp  additional information about the start sites
//  in  start_list .   i  is the subscript of the orf, and  lo .. hi
//  are its tweeny coordinates.   frame  is the reading frame of the
//  orf.

  {
   int  stop_pos;
   int  q, r;

   if  (i == 0)
       printf (">%s\n", Fasta_Header);

   if  (frame > 0)
       stop_pos = hi + 3;
     else
       stop_pos = lo - 2;

   Fix_Wrap (stop_pos, Sequence_Len);
   printf ("# %7d  %+2d\n", stop_pos, frame);

   r = start_list . size ();
   for  (q = 0;  q < r && q < 10;  q ++)
     {
      double  score, combined_score;
      int  j, sep;

      j = start_list [q] . pos;

      if  (frame > 0)
          {
           PWM_Score_Fwd_Start (j, LogOdds_PWM, Ribosome_Window_Size, score, sep);
           combined_score = start_list [q] . score
                + Start_Prob [start_list [q] . which];
           if  (score > 0.0)
                combined_score += score;
           printf ("  %7d %c%c%c %7.3f %7.3f %7.3f %3d\n", j, Sequence [On_Seq_0 (j - 1)],
                Sequence [On_Seq_0 (j)], Sequence [On_Seq_0 (j + 1)],
                start_list [q] . score, score, combined_score, sep);
          }
        else
          {
           PWM_Score_Rev_Start (j, LogOdds_PWM, Ribosome_Window_Size, score, sep);
           combined_score = start_list [q] . score
                + Start_Prob [start_list [q] . which];
           if  (score > 0.0)
                combined_score += score;
           printf ("  %7d %c%c%c %7.3f %7.3f %7.3f %3d\n", j,
                Complement (Sequence [On_Seq_0 (j - 1)]),
                Complement (Sequence [On_Seq_0 (j - 2)]),
                Complement (Sequence [On_Seq_0 (j - 3)]),
                start_list [q] . score, score, combined_score, sep);
          }
     }

   return;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   FILE  * fp;
   char  * p, * q;
   bool  errflg = false;
   int  i, ch;

   optarg = NULL;
   Command_Line = argv [0];

#if  ALLOW_LONG_OPTIONS
   int  option_index = 0;
   static struct option  long_options [] = {
        {"start_codons", 1, 0, 'A'},
        {"rbs_pwm", 1, 0, 'b'},
        {"gc_percent", 1, 0, 'C'},
        {"entropy", 1, 0, 'E'},
        {"first_codon", 0, 0, 'f'},
        {"gene_len", 1, 0, 'g'},
        {"help", 0, 0, 'h'},
        {"ignore", 1, 0, 'g'},
        {"linear", 0, 0, 'l'},
        {"orf_coords", 1, 0, 'L'},
        {"separate_genes", 1, 0, 'M'},
        {"max_olap", 1, 0, 'o'},
        {"start_probs", 1, 0, 'P'},
        {"ignore_score_len", 1, 0, 'q'},
        {"no_indep", 0, 0, 'r'},
        {"threshold", 1, 0, 't'},
        {"extend", 0, 0, 'X'},
        {"trans_table", 1, 0, 'z'},
        {"stop_codons", 1, 0, 'Z'},
        {0, 0, 0, 0}
      };

   while  (! errflg && ((ch = getopt_long (argc, argv,
        "A:b:C:E:fg:hi:lL:Mo:P:q:rt:Xz:Z:",
        long_options, & option_index)) != EOF))
#else
   while  (! errflg && ((ch = getopt (argc, argv,
        "A:b:C:E:fg:hi:lL:Mo:P:q:rt:Xz:Z:")) != EOF))
#endif

     switch  (ch)
       {
        case  'A' :
          Command_Line . append (" -A ");
          Command_Line . append (optarg);
          Start_Codon . clear ();
          for  (p = strtok (optarg, ",");  p != NULL;  p = strtok (NULL, ","))
            {
             q = strdup (p);
             Make_Lower_Case (q);
             Start_Codon . push_back (q);
            }
          break;

        case  'b' :
          Command_Line . append (" -b ");
          Command_Line . append (optarg);
          fp = File_Open (optarg, "r", __FILE__, __LINE__);
          Ribosome_PWM . Read (fp);
          Ribosome_PWM . Counts_To_Prob ();
          Ribosome_PWM . Probs_To_Logs ();
          if  (Verbose > 1)
              Ribosome_PWM . Print (stderr);
          Use_PWM = true;
          break;

        case  'C' :
          Command_Line . append (" -C ");
          Command_Line . append (optarg);
          Indep_GC_Frac = strtod (optarg, & p) / 100.0;
          if  (p == optarg || Indep_GC_Frac < 0.0 || Indep_GC_Frac > 100.0)
              {
               fprintf (stderr, "ERROR:  Bad independent model GC fraction (-C option)\n"
                    "  value = \"%s\"", optarg);
               errflg = true;
              }
          GC_Frac_Set = true;
          break;

        case  'E' :
          Command_Line . append (" -E ");
          Command_Line . append (optarg);
          if  (strcmp (optarg, "#") != 0)
              Read_Entropy_Profiles (optarg, errflg);
          Use_Entropy_Profiles = true;
          break;

        case  'f' :
          Command_Line . append (" -f");
          Use_First_Start_Codon = true;
          break;

        case  'g' :
          Command_Line . append (" -g ");
          Command_Line . append (optarg);
          Min_Gene_Len = strtol (optarg, & p, 10);
          if  (p == optarg || Min_Gene_Len <= 0)
              {
               fprintf (stderr, "ERROR:  Bad minimum gene length (-g option)\n"
                    "  value = \"%s\"", optarg);
               errflg = true;
              }
          break;

        case  'h' :
          Command_Line . append (" -h");
          errflg = true;
          break;

        case  'i' :
          Command_Line . append (" -i ");
          Command_Line . append (optarg);
          Ignore_File_Name = optarg;
          break;

        case  'l' :
          Command_Line . append (" -l");
          Genome_Is_Circular = false;
          break;

        case  'L' :
          Command_Line . append (" -L ");
          Command_Line . append (optarg);
          Orflist_File_Name = optarg;
          break;

        case  'M' :
          Command_Line . append (" -M");
          Separate_Orf_Input = true;
          break;

        case  'o' :
          Command_Line . append (" -o ");
          Command_Line . append (optarg);
          Max_Olap_Bases = strtol (optarg, & p, 10);
          if  (p == optarg || Max_Olap_Bases < 0)
              {
               fprintf (stderr, "ERROR:  Bad max overlap bases (-o option)\n"
                    "  value = \"%s\"", optarg);
               errflg = true;
              }
          break;

        case  'P' :
          Command_Line . append (" -P ");
          Command_Line . append (optarg);
          Start_Prob . clear ();
          for  (p = strtok (optarg, ",");  p != NULL;  p = strtok (NULL, ","))
            Start_Prob . push_back (strtod (p, NULL));
          break;

        case  'q' :
          Command_Line . append (" -q ");
          Command_Line . append (optarg);
          Ignore_Score_Len = strtol (optarg, & p, 10);
          if  (p == optarg || Ignore_Score_Len < 0)
              {
               fprintf (stderr, "ERROR:  Bad ignore independent model length\n"
                    "  (-q option)  value = \"%s\"", optarg);
               errflg = true;
              }
          break;

        case  'r' :
          Command_Line . append (" -r");
          Use_Independent_Score = false;
          break;

        case  't' :
          Command_Line . append (" -t ");
          Command_Line . append (optarg);
          Threshold_Score = strtol (optarg, & p, 10);
          if  (p == optarg || Threshold_Score <= 0 || Threshold_Score >= 100)
              {
               fprintf (stderr, "ERROR:  Bad threshold score (-t option)\n"
                    "  value = \"%s\"", optarg);
               errflg = true;
              }
          break;

        case  'X' :
          Command_Line . append (" -X");
          Allow_Truncated_Orfs = true;
          Genome_Is_Circular = false;
          break;

        case  'z' :
          Command_Line . append (" -z ");
          Command_Line . append (optarg);
          Genbank_Xlate_Code = strtol (optarg, & p, 10);
          Set_Stop_Codons_By_Code (Stop_Codon, Genbank_Xlate_Code, errflg);
          break;

        case  'Z' :
          Command_Line . append (" -Z ");
          Command_Line . append (optarg);
          Stop_Codon . clear ();
          for  (p = strtok (optarg, ",");  p != NULL;  p = strtok (NULL, ","))
            {
             q = strdup (p);
             Make_Lower_Case (q);
             Stop_Codon . push_back (q);
            }
          break;

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = true;
       }

   if  (errflg)
       {
        Usage ();
        exit (EXIT_FAILURE);
       }

   if  (optind > argc - 3)
       {
        Usage ();
        exit (EXIT_FAILURE);
       }

   for  (i = optind;  i < argc;  i ++)
     {
      Command_Line . append (" ");
      Command_Line . append (argv [i]);
     }

   Sequence_File_Name = argv [optind ++];
   ICM_File_Name = argv [optind ++];
   Output_Tag = argv [optind ++];

   return;
  }



template  <class DT>
static void  Permute_By_Frame
    (vector <DT> & v, int frame)

//  Permute the first 6 entries in  v  so that they
//  represent a reverse sequence of a gene, where the first
//  base of the sequence comes from genome position with
//  frame  frame .  Positions of the genome are numbered 1,2,3,1,2,3...
//  Frame is positive for forward strand genes in the genome and negative
//  for reverse strand genes.  The input values in  v  represent
//  scores for a frame  +3  sequence.

  {
   DT  save;

   switch  (frame)
     {
      case  1 :
        save = v [0];
        v [0] = v [2];
        v [2] = v [1];
        v [1] = save;
        save = v [3];
        v [3] = v [5];
        v [5] = v [4];
        v [4] = save;
        break;
      case  2 :
        save = v [0];
        v [0] = v [1];
        v [1] = v [2];
        v [2] = save;
        save = v [3];
        v [3] = v [4];
        v [4] = v [5];
        v [5] = save;
        break;
      case  3 :
        break;
      case  -1 :
        save = v [0];
        v [0] = v [3];
        v [3] = save;
        save = v [1];
        v [1] = v [5];
        v [5] = save;
        save = v [2];
        v [2] = v [4];
        v [4] = save;
        break;
      case  -2 :
        save = v [0];
        v [0] = v [4];
        v [4] = save;
        save = v [1];
        v [1] = v [3];
        v [3] = save;
        save = v [2];
        v [2] = v [5];
        v [5] = save;
        break;
      case  -3 :
        save = v [0];
        v [0] = v [5];
        v [5] = save;
        save = v [1];
        v [1] = v [4];
        v [4] = save;
        save = v [2];
        v [2] = v [3];
        v [3] = save;
        break;
     }

   return;
  }



int  Position_To_Frame
    (int p)

//  Return the reading frame corresponding to a codon beginning in
//  position  p .  Allow  p  to be negative.  For  p = ...,-2,-1,0,1,2,3,4,...
//  frames are, respectively,  ...,1,2,3,1,2,3,1,...

  {
   if  (p >= 0)
       return  1 + ((p + 2) % 3);
     else
       return  3 - ((-1 * p) % 3);
  }



static void  Print_Comma_Separated_Strings
    (const vector <const char *> & v, FILE * fp)

//  Print the strings in  v  to  fp .  Separate them by
//  commas with no spaces.

  {
   int  i, n;

   n = v . size ();

   if  (n == 0)
       return;

   fprintf (fp, "%s", v [0]);
   for  (i = 1;  i < n;  i ++)
     fprintf (fp, ",%s", v [i]);

   return;
  }



static void  Print_Headings
    (FILE * fp)

//  Print column headings to  fp .

  {
   fputc ('\n', fp);

   fprintf (fp, "%4s %5s %17s %8s  %15s", "", "", "----- Start -----",
        "", "--- Length ----");
   if  (Use_Independent_Score)
       fprintf (fp, "  %s\n", "------------- Scores -------------");
     else
       fprintf (fp, "  %s\n", "----------- Scores ------------");
   fprintf (fp, "%4s %5s %8s %8s %8s  %7s %7s  %7s %5s %s",
        " ID ", "Frame", "of Orf", "of Gene", "Stop", "of Orf", "of Gene",
        "Raw", "InFrm", "F1 F2 F3 R1 R2 R3");
   if  (Use_Independent_Score)
       fprintf (fp, " NC");
   if  (Use_Entropy_Profiles)
       fprintf (fp, " %4s", "EDR");
   fprintf (fp, "\n");

   return;
  }



static void  Print_Orflist_Headings
    (FILE * fp)

//  Print column headings for separate orf list (-L option) to  fp .

  {
   fputc ('\n', fp);

   fprintf (fp, "%-12s %5s  %8s %8s %8s", "", "", "", "", "");
   if  (Use_Independent_Score)
       fprintf (fp, "  %s\n", "------------- Scores -------------");
     else
       fprintf (fp, "  %s\n", "----------- Scores ------------");
   fprintf (fp, "%-12s %5s  %8s %8s %8s  %7s %5s %s",
        "  ID", "Frame", "Start", "Stop", "Len", "Raw", "InFrm", "F1 F2 F3 R1 R2 R3");
   if  (Use_Independent_Score)
       fprintf (fp, " NC");
   if  (Use_Entropy_Profiles)
       fprintf (fp, " %-4s", "EDR");
   fprintf (fp, "\n");

   return;
  }



static const char  * Print_String
    (Event_t e)

//  Return a printable equivalent for  e .

  {
   switch  (e)
     {
      case  INITIAL :
        return  "Initial";
      case  FWD_START :
        return  "F_Start";
      case  FWD_STOP :
        return  "F_Stop";
      case  REV_START :
        return  "R_Start";
      case  REV_STOP :
        return  "R_Stop";
      case  TERMINAL :
        return  "Terminal";
     }
   return  "None";
  }



static void  Prob_To_Logs
    (vector <double> & v)

//  Convert the entries in  v  to their natural logarithms.
//  Add psuedo-count value for zero entries.  Normalize all
//  values in case the original values don't sum to 1.0

  {
   double  subtr;
   double  sum = 0.0, sum2 = 0.0;
   int  i, n;

   n = v . size ();
   for  (i = 0;  i < n;  i ++)
     {
      if  (v [i] < 0.0)
          {
           sprintf (Clean_Exit_Msg_Line, "ERROR:  Bad start codon probability %f\n",
                v [i]);
           Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
          }
      sum += v [i];
     }
   if  (sum == 0.0)
       {
        sprintf (Clean_Exit_Msg_Line, "ERROR:  Start codon probabilities all zero\n");
        Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
       }

   for  (i = 0;  i < n;  i ++)
     if  (v [i] == 0.0)
         {
          v [i] = sum * 1e-5;
          sum2 += v [i];
         }
   subtr = log (sum + sum2);

   for  (i = 0;  i < n;  i ++)
     v [i] = log (v [i]) - subtr;

   return;
  }



static void  Process_Events
    (void)

//  Find the best-scoring collection of genes represented by the
//  sequence of events in the global list of events pointed to by
//   Last_Event .

  {
   vector <Event_Node_t *> ep;
   Event_Node_t  * p;
   int  i, n;

   // Make  ep  point to all the events
   // Also make the initial event's position smaller than the
   // position of any other event
   for  (i = 0;  i < 6;  i ++)
     {
      int  min_pos = 0;

      for  (p = Last_Event [i];  p != NULL && p -> e_type != INITIAL ;
                p = p -> frame_pred)
        {
         ep . push_back (p);
         min_pos = Min (min_pos, p -> pos - 1);
        }
      if  (p == NULL)
          {
           sprintf (Clean_Exit_Msg_Line, "ERROR:  Missing initial event\n");
           Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
          }
      p -> pos = Min (min_pos, p -> pos);
     }
   // Add a single copy of the initial event
   ep . push_back (p);

   n = int (ep . size ());

   // Sort all events into order by their  pos  field
   sort (ep . begin (), ep . end (), Event_Pos_Cmp);

   if  (Genome_Is_Circular)
       {
        int  reference_pos;

        reference_pos = Find_Uncovered_Position (ep);
        if  (reference_pos > 0)
            Shift_Events (ep, reference_pos);
       }

   // Scan  ep  and by dynamic programming find the best predecessor
   // event for each event.  Save the best event in each frame in
   // global  Best_Event [] .

   for  (i = 0;  i < n;  i ++)
     switch  (ep [i] -> e_type)
       {
        case  INITIAL :
          Process_Initial_Event (ep [i]);
          break;
        case  FWD_START :
          Process_Fwd_Start_Event (ep [i]);
          break;
        case  FWD_STOP :
          Process_Fwd_Stop_Event (ep [i]);
          break;
        case  REV_START :
          Process_Rev_Start_Event (ep [i]);
          break;
        case  REV_STOP :
          Process_Rev_Stop_Event (ep [i]);
          break;
        default :
          sprintf (Clean_Exit_Msg_Line, "ERROR:  Unexpected event type = %d\n",
               int (ep [i] -> e_type));
          Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
       }

   return;
  }



static void  Process_Fwd_Start_Event
    (Event_Node_t * ep)

//  Process the forward-start-type event pointed to by  ep  by computing
//  the best score that can be obtained by combining it with
//  prior events.

  {
   int  i, f, mxi;

   f = Frame_To_Sub (ep -> frame);

   // Connect  ep  to the highest-scoring prior event and increment
   //  ep -> score  by that score
   mxi = 0;
   for  (i = 1;  i < 6;  i ++)
     if  (Best_Event [i] -> score > Best_Event [mxi] -> score)
         mxi = i;
   ep -> best_pred = Best_Event [mxi];
   ep -> score += Best_Event [mxi] -> score;

   // Make  ep  the last in the chain of events in this reading frame
   ep -> frame_pred = Last_Event [f];
   Last_Event [f] = ep;

   return;
  }



static void  Process_Fwd_Stop_Event
    (Event_Node_t * ep)

//  Process the forward-stop-type event pointed to by  ep  by 
//  connecting it to the best previous start codon in the same frame.
//  If that score is better than the best score in the frame, then
//  make  Best_Event  for the frame point to  ep .  Also check for
//  allowed overlaps with prior forward starts or reverse stops.

  {
   Event_Node_t  * p, * best_p;
   double  mx;
   int  i, f;

   f = Frame_To_Sub (ep -> frame);

   // Find the best preceding event and make  ep  point back to it
   mx = 0.0;
   best_p = NULL;
   for  (p = Last_Event [f];  p -> e_type == FWD_START;  p = p -> frame_pred)
     if  (p -> score > mx)
         {
          mx = p -> score;
          best_p = p;
         }
   ep -> best_pred = best_p;
   ep -> score = mx;

   // Check any events that represent genes that may overlap this one
   // by less than the allowable overlap threshold and adjust their
   // score and make them point to  ep  if it gives a better score
   if  (Best_Event [f] -> score < ep -> score)
       {
        Disqualify (best_p, 3 + ep -> pos - Max_Olap_Bases);
        Best_Event [f] = ep;
        for  (i = 0;  i < 6;  i ++)
          {
           if  (i == f)
               continue;
           for  (p = Last_Event [i];
                      p != NULL && 3 + ep -> pos - p -> pos <= Max_Olap_Bases;
                      p = p -> frame_pred)
             {
              double  score_needed;

              if  (p -> best_pred == NULL)
                  score_needed = 0.0;
                else
                  score_needed = p -> best_pred -> score;
              if  ((p -> e_type == FWD_START || p -> e_type == REV_STOP)
                        && ! p -> disqualified
                        && score_needed < ep -> score)
                  {
                   Event_Node_t  * q;
                   double  adj, diff;
                   int  lo;

                   if  (p -> e_type == FWD_START)
                       lo = p -> pos - 2;
                     else
                       lo = p -> pos + 1;
                   adj = Olap_Score_Adjustment (lo, ep -> pos - 3, p -> frame,
                             ep -> frame);
                   diff = ep -> score - p -> best_pred -> score - adj;

                   if  (diff <= 0.0)
                       continue;

                   p -> score += diff;
                   p -> best_pred = ep;
                   for  (q = Last_Event [i];  q != p;  q = q -> frame_pred)
                     if  (q -> best_pred == p)
                         q -> score += diff;
                  }
             }
          }
        Requalify (best_p, 3 + ep -> pos - Max_Olap_Bases);
       }

   // Make  ep  the last in the chain of events in this reading frame
   ep -> frame_pred = Last_Event [f];
   Last_Event [f] = ep;

   return;
  }



static void  Process_Initial_Event
    (Event_Node_t * ep)

//  Process the initial-type event pointed to by  ep  by adding
//  it to the global lists  Best_Event []  and  Last_Event [] .

  {
   int  i;

   for  (i = 0;  i < 6;  i ++)
     Best_Event [i] = Last_Event [i] = ep;

   ep -> pos = 0;
   ep -> score = 0.0;
   ep -> frame_pred = ep -> best_pred = NULL;

   return;
  }



static void  Process_Rev_Start_Event
    (Event_Node_t * ep)

//  Process the reverse-start-type event pointed to by  ep  by computing
//  the best score that can be obtained by combining it with
//  prior events.

  {
   Event_Node_t  * p;
   int  i, f;

   f = Frame_To_Sub (ep -> frame);

   // Connect  ep  to its corresponding reverse-stop event and increment
   //  ep -> score  by that score
   for  (p = Last_Event [f];  p != NULL && p -> e_type == REV_START;
              p = p -> frame_pred)
     ;
   if  (p == NULL || p -> e_type != REV_STOP)
       {
        sprintf (Clean_Exit_Msg_Line,
             "ERROR:  No reverse stop for reverse start at pos = %d\n", ep -> pos);
        Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
       }
   ep -> best_pred = p;
   ep -> score += p -> score;

   // Check any events that represent genes that may overlap this one
   // by less than the allowable overlap threshold and adjust their
   // score and make them point to  ep  if it gives a better score
   if  (Best_Event [f] -> score < ep -> score)
       {
        Disqualify (p, 3 + ep -> pos - Max_Olap_Bases);
        Best_Event [f] = ep;
        for  (i = 0;  i < 6;  i ++)
          {
           if  (i == f)
               continue;
           for  (p = Last_Event [i];
                      p != NULL && 3 + ep -> pos - p -> pos <= Max_Olap_Bases;
                      p = p -> frame_pred)
             {
              double  score_needed;

              if  (p -> best_pred == NULL)
                  score_needed = 0.0;
                else
                  score_needed = p -> best_pred -> score;
              if  ((p -> e_type == FWD_START || p -> e_type == REV_STOP)
                        && ! p -> disqualified
                        && score_needed < ep -> score)
                  {
                   Event_Node_t  * q;
                   double  adj, diff;
                   int  lo;

                   if  (p -> e_type == FWD_START)
                       lo = p -> pos - 2;
                     else
                       lo = p -> pos + 1;
                   adj = Olap_Score_Adjustment (lo, ep -> pos, p -> frame,
                             ep -> frame);
                   diff = ep -> score - p -> best_pred -> score - adj;

                   if  (diff <= 0.0)
                       continue;

                   p -> score += diff;
                   p -> best_pred = ep;
                   for  (q = Last_Event [i];  q != p;  q = q -> frame_pred)
                     if  (q -> best_pred == p)
                         q -> score += diff;
                  }
             }
          }
        Requalify (p, 3 + ep -> pos - Max_Olap_Bases);
       }

   // Make  ep  the last in the chain of events in this reading frame
   ep -> frame_pred = Last_Event [f];
   Last_Event [f] = ep;

   return;
  }



static void  Process_Rev_Stop_Event
    (Event_Node_t * ep)

//  Process the reverse-stop-type event pointed to by  ep  by computing
//  the best score that can be obtained by combining it with
//  prior events.

  {
   int  i, f, mxi;

   f = Frame_To_Sub (ep -> frame);

   // Connect  ep  to the highest-scoring prior event and increment
   //  ep -> score  by that score
   mxi = 0;
   for  (i = 1;  i < 6;  i ++)
     if  (Best_Event [i] -> score > Best_Event [mxi] -> score)
         mxi = i;
   ep -> best_pred = Best_Event [mxi];
   ep -> score = Best_Event [mxi] -> score;

   // Make  ep  the last in the chain of events in this reading frame
   ep -> frame_pred = Last_Event [f];
   Last_Event [f] = ep;

   return;
  }



static void  PWM_Score_Fwd_Start
    (int pos, const PWM_t & pwm, int window, double & score, int & separation)

//  Find the highest scoring match for  pwm
//  against the sequence in a window of length  window
//  in front of position  pos  (numbered starting at 1) in the
//  forward direction.  Set  score  to the highest score and
//  set  separation  to the number of positions between the best
//  match and  pos .

  {
   double  sc;
   int  bottom, lo, sep;
   int  j, n;

   score = 0.0;
   separation = 0;

   if  (pwm . Is_Empty ())
       return;

   n = pwm . Width ();
   bottom = pos - window - 1;

   score = - DBL_MAX;
   separation = sep = 0;
   for  (lo = pos - n - 1;  0 <= lo && bottom <= lo;  lo --, sep ++)
     {
      sc = 0.0;
      for  (j = 0;  j < n;  j ++)
        sc += pwm . Column_Score (Sequence [lo + j], j);
      if  (sc > score)
          {
           score = sc;
           separation = sep;
          }
     }

   // handle wraparound here
   if  (Genome_Is_Circular)
       for  ( ;  bottom <= lo;  lo --, sep ++)
         {
          sc = 0.0;
          for  (j = 0;  j < n;  j ++)
            {
             int  k;

             k = lo + j;
             if  (k < 0)
                 k += Sequence_Len;
             sc += pwm . Column_Score (Sequence [k], j);
            }
          if  (sc > score)
              {
               score = sc;
               separation = sep;
              }
         }

   return;
  }



static void  PWM_Score_Rev_Start
    (int pos, const PWM_t & pwm, int window, double & score, int & separation)

//  Find the highest scoring match for  pwm
//  against the sequence in a window of length  window
//  following position  pos  (numbered starting at 1) on the
//  reverse-complement strand.  Set  score  to the highest score and
//  set  separation  to the number of positions between the best
//  match and  pos .

  {
   double  sc;
   int  top, hi, sep;
   int  j, n;

   if  (pwm . Is_Empty ())
       {
        score = 0.0;
        separation = 0;
        return;
       }

   n = pwm . Width ();
   top = pos - 1 + window;

   score = - DBL_MAX;
   separation = sep = 0;
   for  (hi = pos - 1 + n;  hi < Sequence_Len && hi <= top;  hi ++, sep ++)
     {
      sc = 0.0;
      for  (j = 0;  j < n;  j ++)
        sc += pwm . Column_Score (Complement (Sequence [hi - j]), j);
      if  (sc > score)
          {
           score = sc;
           separation = sep;
          }
     }

   // handle wraparound here
   for  ( ;  hi <= top;  hi ++, sep ++)
     {
      sc = 0.0;
      for  (j = 0;  j < n;  j ++)
        {
         int  k;

         k = hi - j;
         if  (Sequence_Len <= k)
             k -= Sequence_Len;
         sc += pwm . Column_Score (Complement (Sequence [k]), j);
        }
      if  (sc > score)
          {
           score = sc;
           separation = sep;
          }
     }

   return;
  }



static void  Read_Entropy_Profiles
    (const char * fn, bool & errflg)

//  Read positive and negative entropy profiles from the
//  file name  fn .  If not successful, set  errflg  to  true .
//  Save the entropy profiles in globals  Pos_Entropy_Profile
//  and  Neg_Entropy_Profile .

  {
   FILE  * fp;
   char  line [MAX_LINE];
   int  i;

   fp = File_Open (fn, "r");
   fgets (line, MAX_LINE, fp);  // skip header line
   for  (i = 0;  i < 20;  i ++)
     if  (fscanf (fp, "%s %lf %lf\n", line, Pos_Entropy_Profile + i,
             Neg_Entropy_Profile + i) != 3)
         {
          errflg = true;
          return;
         }

   fclose (fp);

   return;
  }



static void  Read_Sequences
    (FILE * fp, vector <string> & seq_list, vector <string> & hdr_list,
     int & seq_ct)

//  Read fasta-format sequences from  fp  (which is already open),
//  convert them to lower-case, and store them in  seq_list .
//  Store the fasta header lines in  hdr_list .  Set  seq_ct  to
//  the number of sequences read.

  {
   string  seq, hdr;
   int  i, len;

   seq_list . clear ();
   hdr_list . clear ();
   seq_ct = 0;

   while  (Fasta_Read (fp, seq, hdr))
     {
      len = seq . length ();
      for  (i = 0;  i < len;  i ++)
        seq [i] = Filter (tolower (seq [i]));

      seq_list . push_back (seq);
      hdr_list . push_back (hdr);
      seq_ct ++;
     }

   return;
  }



static void  Requalify
    (Event_Node_t * p, int cutoff)

//  Set the  disqualified  bit false for nodes reachable from
//   p  by  best_pred  pointers that have  pos  values at least
//  as great as  cutoff .

  {
   Event_Node_t  * q;

   if  (p == NULL)
       return;

   for  (q = p -> best_pred;  q != NULL && cutoff <= q -> pos;  q = q -> best_pred)
     q -> disqualified = false;

   return;
  }



static void  Reverse_Complement_Transfer
    (string & buff, const string & s, int lo, int hi)

//  Copy to string  buff  the reverse complement of the substring
//  of  s  between positions  lo  and  hi  (which are
//  space-based coordinates).

  {
   int  i, j;

   assert (hi <= int (s . length ()));

   buff . resize (hi - lo);
   for  (j = 0, i = hi - 1;  i >= lo;  j ++, i --)
     buff [j] = Complement (s [i]);

   return;
  }



static void  Reverse_Transfer
    (string & buff, const string & s, int start, int len)

//  Copy to string  buff  the substring of  s  starting at subscript
//   start  and going to the left for a length of  len .  Wraparound
//  end of  s  if necessary.  Do *NOT* reverse-complement.

  {
   int  j, n;

   n = s . length ();
   assert (start < n);
   assert (0 <= len);

   buff . resize (len);
   for  (j = 0;  j < len;  j ++, start --)
     {
      buff [j] = s [start];
      if  (start <= 0)
          start += n;
     }

   return;
  }



static void  Score_Orflist
    (FILE * detail_fp, FILE * summary_fp)

//  Score the entries in global  Orf_Pos_List  using the sequence
//  in global  Sequence  sending detailed results to  detail_fp and
//  summary results to  summary_fp.

  {
   string  buff;
   vector <double>  af, score, indep_score;
   vector <int>  int_score;
   vector <bool>  has_stop;
   int  fr, frame, frame_score;
   int  lo, hi, len;
   int  i, j, m, n;

   if  (Use_Independent_Score)
       af . resize (7);
     else
       af . resize (6);

   n = Orf_Pos_List . size ();
   for  (i = 0;  i < n;  i ++)
     {
      double  gene_score;
      int  start, stop;

      start = Orf_Pos_List [i] . start;
      stop = Orf_Pos_List [i] . stop;

      if  (Orf_Pos_List [i] . dir > 0)
          {
           frame = 1 + (stop % 3);
           fr = 1 + (1 + frame) % 3;
           len = 1 + stop - start - 3;
           if  (len < 0)
               len += Sequence_Len;
           hi = stop - 3;
           if  (hi <= 0)
               hi += Sequence_Len;
           Reverse_Transfer (buff, Sequence, hi - 1, len);
          }
        else
          {
           fr = frame = - ((stop - 1) % 3) - 1;
           len = 1 + start - stop - 3;
           if  (len < 0)
               len += Sequence_Len;
           lo = stop + 2;
           if  (lo >= Sequence_Len)
               lo -= Sequence_Len;
           Complement_Transfer (buff, Sequence, lo, len);
          }

      Gene_ICM . Cumulative_Score (buff, score, 1);
      Indep_Model . Cumulative_Score (buff, indep_score, 1);
      m = score . size ();

      if  (Use_Independent_Score)
          af [6] = indep_score [m - 4];   // excludes the start codon
      All_Frame_Score (buff, m - 3, fr, af);
      Find_Stops_Reverse (buff, m - 3, has_stop);
      gene_score = 100.0 * (score [m - 4] - indep_score [m - 4]) / (m - 3);

      Permute_By_Frame (has_stop, fr);
      Integerize_Scores (af, HI_SCORE, has_stop, int_score);
      if  (frame > 0)
          frame_score = int_score [frame - 1];
        else
          frame_score = int_score [2 - frame];

      // print score details
      fprintf (detail_fp, "%-14s %+3d  %8d %8d %8d  %7.2f %5d",
           Orf_Pos_List [i] . tag, frame, start, stop, len, gene_score, frame_score);
      for  (j = 0;  j < 6;  j ++)
        if  (int_score [j] < 0)
            fprintf (detail_fp, "  -");
          else
            fprintf (detail_fp, " %2d", int_score [j]);
      if  (Use_Independent_Score)
          fprintf (detail_fp, " %2d", int_score [6]);
      if  (Use_Entropy_Profiles)
          fprintf (detail_fp, " %4.2f", Entropy_Distance_Ratio (start, m, frame));
      fputc ('\n', detail_fp);

      // print overall score
      fprintf (summary_fp, "%-14s %8d %8d %+3d %8.2f\n",
          Orf_Pos_List [i] . tag, start, stop, frame, gene_score);
     }

   return;
  }



static void  Score_Orfs
    (vector <Orf_t> & orf_list, vector <Gene_t> & gene_list, FILE * fp)

//  Compute scores for all orfs in  orf_list  using coding model
//  in global  Gene_ICM , which is assumed to have been built on reverse
//  gene strings.   Indep_Model  is the model of independent,
//  stop-codon-free sequence.  Put orfs that are candidate genes
//  onto  gene_list .  Print log information to  fp .

  {
   string  buff;
   vector <double>  af, score, indep_score;
   vector <bool>  is_start;
   vector <Start_t>  start_list;
   Start_t  start;
   char  tag [MAX_LINE];
   int  i, n, id = 0;

   if  (Use_Independent_Score)
       af . resize (7);
     else
       af . resize (6);

   gene_list . clear ();

   n = orf_list . size ();
   for  (i = 0;  i < n;  i ++)
     {
      double  first_score, best_score = - DBL_MAX;
      double  gene_score;
      vector <int>  int_score;
      vector <bool>  has_stop;
      int  first_pos = 0, best_pos = 0;
      int  first_j = 0, best_j = 0;
      double  max, max_rate;
      Codon_t  codon;
      double  s;
      bool  is_tentative_gene, orf_is_truncated = false;
      bool  first_is_truncated = false, best_is_truncated = false;
      int  which;
      int  fr, frame, max_j, orf_start, frame_score;
      int  lo, hi, len, lowest_j;
      int  j, k, m;

      frame = orf_list [i] . Get_Frame ();
      len = orf_list [i] . Get_Orf_Len ();
      if  (frame > 0)
          {
           hi = orf_list [i] . Get_Stop_Position () - 1;
           if  (hi <= 0)
               hi += Sequence_Len;
           lo = hi - len;
           Reverse_Transfer (buff, Sequence, hi - 1, len);
           fr = 1 + (1 + frame) % 3;
           orf_is_truncated = (lo < 3 && Allow_Truncated_Orfs);
           k = orf_list [i] . Get_Stop_Position () - len - 2;
          }
        else
          {
           lo = orf_list [i] . Get_Stop_Position () + 2;
           if  (lo >= Sequence_Len)
               lo -= Sequence_Len;
           hi = lo + len;
           Complement_Transfer (buff, Sequence, lo, len);
           fr = frame;
           orf_is_truncated = (Sequence_Len - hi < 3 && Allow_Truncated_Orfs);
           k = orf_list [i] . Get_Stop_Position () + len + 4;
          }
      // lo .. hi  are the between coordinates of the orf region.

      Gene_ICM . Cumulative_Score (buff, score, 1);
      Indep_Model . Cumulative_Score (buff, indep_score, 1);
      m = score . size ();

      max = 0.0;
      max_j = 0;
      is_start . resize (m, false);
      start_list . clear ();
      lowest_j = Min (3, Min_Gene_Len - 3);
      for  (j = m - 1;  j >= lowest_j;  j --)
        {
         codon . Shift_In (buff [j]);
         s = score [j] - indep_score [j];
         if  (s > max)
             {
              max = s;
              max_rate = s / (j + 1);
              max_j = j;
             }
         if  (j % 3 == 0
                 && (codon . Can_Be (Fwd_Start_Pattern, which)
                       || (first_pos == 0 && orf_is_truncated))
                 && j + 3 >= Min_Gene_Len)
             {
              double  next_s;

              next_s = score [j - 1] - indep_score [j - 1];
                // this is the score for the orf without the start
                // codon--position j is the last base of the start codon
              is_start [j + 2] = true;
              start . j = j + 2;
              start . pos = k;
                // k is the 1-based sequence coordinate of the base that
                // is 2 behind the position represented by j
              start . which = which;
              start . truncated = (which < 0);
              start . score = next_s;
              start . first = (first_pos == 0);
              start_list . push_back (start);

              if  (first_pos == 0)
                  {
                   first_score = next_s;
                   first_pos = k;
                   first_j = j + 2;
                   first_is_truncated = start . truncated;
                  }
              if  (next_s > best_score)
                  {
                   best_score = next_s;
                   best_pos = k;
                   best_j = j + 2;
                   best_is_truncated = start . truncated;
                  }
             }
         if  (frame > 0)
             k ++;
           else
             k --;
        }

      if  (Use_First_Start_Codon)
          {
           best_score = first_score;
           best_pos = first_pos;
           best_j = first_j;
           best_is_truncated = first_is_truncated;
          }

      if  (first_j + 1 < Min_Gene_Len)
          continue;

      if  (frame > 0)
          {
           k = hi + 3;
           orf_start = lo + 1;
          }
        else
          {
           k = lo - 2;
           orf_start = hi;
          }

      if  (Use_Independent_Score)
          af [6] = indep_score [best_j - 3];

//**ALD  Changed  best_j + 1  to  best_j - 2  to omit start codon
//  from score to be consistent with the independent score
      All_Frame_Score (buff, best_j - 2, fr, af);
      Find_Stops_Reverse (buff, best_j - 2, has_stop);

      Permute_By_Frame (has_stop, fr);
      Integerize_Scores (af, HI_SCORE, has_stop, int_score);
      if  (frame > 0)
          frame_score = int_score [frame - 1];
        else
          frame_score = int_score [2 - frame];

      // For now just use the score, will add more later
      is_tentative_gene
           = (best_j + 1 >= Min_Gene_Len && frame_score >= Threshold_Score);

      // If it's long enough to ignore the independent score,
      // rescue it
      if  (! is_tentative_gene && len >= Ignore_Score_Len)
          {
           best_score = first_score;
           best_pos = first_pos;
           best_j = first_j;
           is_tentative_gene = true;
          }

//**ALD  Changed to omit start codon
      gene_score = 100.0 * best_score / (best_j - 2);

      if  (For_Edwin)
          Output_Extra_Start_Info (stdout, i, lo, hi, frame, start_list);

      if  (is_tentative_gene)
          {
           sprintf (tag, "%04d", ++ Gene_ID_Ct);
           Gene_t  gene (orf_list [i]);
           gene . Set_Score (gene_score);
           gene . Set_Gene_Len (best_j + 1);
           gene_list . push_back (gene);
          }
        else
          strcpy (tag, "    ");

      if  (Genome_Is_Circular)
          {
           Fix_Wrap (orf_start, Sequence_Len);
           Fix_Wrap (best_pos, Sequence_Len);
           Fix_Wrap (k, Sequence_Len);
          }
      else if  (orf_is_truncated)
          {
           if  (frame > 0)
               {
                orf_start -= 3;
                if  (best_is_truncated)
                    best_pos -= 3;
               }
             else
               {
                orf_start += 3;
                if  (best_is_truncated)
                    best_pos += 3;
               }
          }

      fprintf (fp, "%4s %+5d %8d %8d %8d  %7d %7d  %7.2f %5d",
           tag, frame, orf_start, best_pos, k, len, best_j + 1,
           gene_score, frame_score );
      for  (j = 0;  j < 6;  j ++)
        if  (int_score [j] < 0)
            fprintf (fp, "  -");
          else
            fprintf (fp, " %2d", int_score [j]);
      if  (Use_Independent_Score)
          fprintf (fp, " %2d", int_score [6]);
      if  (Use_Entropy_Profiles)
          fprintf (fp, " %4.2f", Entropy_Distance_Ratio (best_pos,
               best_j + 1, frame));
      fputc ('\n', fp);

      if  (is_tentative_gene)
          Add_Events (orf_list [i], start_list, ++ id);
     }

   return;
  }



static void  Score_Separate_Input
    (const string & seq, const string & hdr, int seq_num, FILE * detail_fp,
     FILE * predict_fp)

//  Score the sequence  seq  with fasta header  hdr  in frame and output
//  the results to  detail_fp  and  predict_fp .

  {
   string  buff;
   vector <double>  af, score, indep_score;
   char  line [MAX_LINE], tag [MAX_LINE], * p;
   vector <int>  int_score;
   vector <bool>  has_stop;
   double  gene_score;
   int  fr, frame, frame_score;
   int  len;
   int  j, m;

   len = seq . length () - 3;  // remove stop codon
   Reverse_Transfer (buff, seq, len - 1, len);
   strcpy (line, hdr . c_str ());
   p = strtok (line, " \t\n");
   if  (p == NULL)
       sprintf (tag, "Seq%04d", seq_num);
     else
       strcpy (tag, p);

   if  (Use_Independent_Score)
       af . resize (7);
     else
       af . resize (6);
   
   frame = 1;  // assume all orfs are in correct reading frame
   fr = 3;     // shifted number for this frame

   Gene_ICM . Cumulative_Score (buff, score, 1);
   Indep_Model . Cumulative_Score (buff, indep_score, 1);
   len = m = score . size ();

   if  (Use_Independent_Score)
       af [6] = indep_score [m - 4];   // excludes the start codon
   All_Frame_Score (buff, m - 3, fr, af);
   Find_Stops_Reverse (buff, m - 3, has_stop);
   gene_score = 100.0 * (score [m - 4] - indep_score [m - 4]) / (m - 3);

   Permute_By_Frame (has_stop, fr);
   Integerize_Scores (af, HI_SCORE, has_stop, int_score);
   if  (frame > 0)
       frame_score = int_score [frame - 1];
     else
       frame_score = int_score [2 - frame];

   // print score details
   fprintf (detail_fp, "%-14s %+3d  %8d %8d %8d  %7.2f %5d",
        tag, frame, 1, len, len, gene_score, frame_score);
   for  (j = 0;  j < 6;  j ++)
     if  (int_score [j] < 0)
         fprintf (detail_fp, "  -");
       else
         fprintf (detail_fp, " %2d", int_score [j]);
   if  (Use_Independent_Score)
       fprintf (detail_fp, " %2d", int_score [6]);
   if  (Use_Entropy_Profiles)
       fprintf (detail_fp, " %4.2f", Entropy_Distance_Ratio (1, m, frame));
   fputc ('\n', detail_fp);

   // print overall score
   fprintf (predict_fp, "%-14s %8d %8d %+3d %8.2f\n",
       tag, 1, len, frame, gene_score);

   return;
  }



static void  Set_Final_Event
    (Event_Node_t & fe, Event_Node_t * best_event [6],
     int seq_len)

//  Set final event  fe , representing the end of genome,
//  and make it point back to the best event in  best_event .
//   seq_len  is the length of the entire genome sequence.

  {
   int  i;

   fe . pos = seq_len;
   fe . score = best_event [0] -> score;
   fe . best_pred = best_event [0];

   for  (i = 1;  i < 6;  i ++)
     {
      if  (best_event [i] -> score >= fe . score)
          {
           fe . score = best_event [i] -> score;
           fe . best_pred = best_event [i];
          }
     }

   return;
  }



static void  Set_GC_Fraction
    (double & gc, const vector <string> & s)

//  Set  gc  to the fraction of letters in all strings in  s  that are
//  'g' or 'c'.

  {
   int  i, j, n, m, ct = 0, total = 0;

   n = s . size ();
   for  (i = 0;  i < n;  i ++)
     {
      m = s [i] . length ();
      total += m;
      for  (j = 0;  j < m;  j ++)
        if  (s [i] [j] == 'g' || s [i] [j] == 'c')
            ct ++;
     }

   gc = double (ct) / total;

   return;
  }



static void  Set_Ignore_Score_Len
    (void)

//  Set global  Ignore_Score_Len  to the length of the longest orf
//  that would be expected to occur once at random in a million bases.
//  Assume an over-simplified model with independent stop codons.

  {

   if  (Ignore_Score_Len == INT_MAX)
       {
        double  poisson_lambda = 0.0;
        int  i, n;

        n = Stop_Codon . size ();
        for  (i = 0;  i < n;  i ++)
          {
           double  x = 1.0;
           int  j;

           for  (j = 0;  j < 3;  j ++)
             if  (Stop_Codon [i] [j] == 'c' || Stop_Codon [i] [j] == 'g')
                 x *= Indep_GC_Frac / 2.0;
               else
                 x *= (1.0 - Indep_GC_Frac) / 2.0;

           poisson_lambda += x;
          }

        assert (poisson_lambda != 0.0);
        Ignore_Score_Len
             = (long int) floor (3.0 * log (2.0 * 1000000 * poisson_lambda)
                 / poisson_lambda);
       }

   return;
  }



static void  Set_Start_And_Stop_Codons
    (void)

//  Set globals  Start_Codon  and  Stop_Codon  to the sequences
//  that are allowed to be start and stop codons for genes.

  {
   Codon_t  codon;
   int  i, n;

   if  (Start_Codon . size () == 0)
       {
        n = sizeof (DEFAULT_START_CODON) / sizeof (char *);
        for  (i = 0;  i < n;  i ++)
          Start_Codon . push_back (DEFAULT_START_CODON [i]);
        if  (Start_Prob . size () == 0)
            for  (i = 0;  i < n;  i ++)
              Start_Prob . push_back (DEFAULT_START_PROB [i]);
        else if  (Start_Codon . size () != Start_Prob . size ())
            {
             sprintf (Clean_Exit_Msg_Line,
                  "ERROR:  Different number of start codons & probs (%d & %d, resp.)\n",
                  int (Start_Codon . size ()), int (Start_Prob . size ()));
             Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
            }
       }
   else if  (Start_Prob . size () == 0)
       { // assign equal likelihood
        n = Start_Codon . size ();
        for  (i = 0;  i < n;  i ++)
          Start_Prob . push_back (1.0 / n);
       }
   else if  (Start_Codon . size () != Start_Prob . size ())
       {
        sprintf (Clean_Exit_Msg_Line,
             "ERROR:  Different number of start codons & probs (%d & %d, resp.)\n",
             int (Start_Codon . size ()), int (Start_Prob . size ()));
        Clean_Exit (Clean_Exit_Msg_Line, __FILE__, __LINE__);
       }

   if  (Stop_Codon . size () == 0)
       {
        n = sizeof (DEFAULT_STOP_CODON) / sizeof (char *);
        for  (i = 0;  i < n;  i ++)
          Stop_Codon . push_back (DEFAULT_STOP_CODON [i]);
       }

   Fwd_Start_Pattern . clear ();
   Fwd_Stop_Pattern . clear ();
   Rev_Start_Pattern . clear ();
   Rev_Stop_Pattern . clear ();

   n = Num_Start_Codons = Start_Codon . size ();
   for  (i = 0;  i < n;  i ++)
     {
      codon . Set_From (Start_Codon [i]);
      Fwd_Start_Pattern . push_back (codon);
      codon . Reverse_Complement ();
      Rev_Start_Pattern . push_back (codon);
     }

   n = Num_Stop_Codons = Stop_Codon . size ();
   for  (i = 0;  i < n;  i ++)
     {
      codon . Set_From (Stop_Codon [i]);
      Fwd_Stop_Pattern . push_back (codon);
      codon . Reverse_Complement ();
      Rev_Stop_Pattern . push_back (codon);
     }

   return;
  }



static void  Shift_Events
    (vector <Event_Node_t *> & ep, int reference_pos)

//  Change the position of all events in  ep  that are before
//   reference_pos  by adding global  Sequence_Len  to them
//  and then sort according to the new positions

  {
   Event_Node_t  * frame_last [6];
   int  f, i, n, q;

   n = ep . size ();
   if  (n <= 1)
       return;

   for  (f = 0;  f < 6;  f ++)
     frame_last [f] = Last_Event [f];

   // Find the lowest-position event in each frame after  reference_pos
   // ep [0] is the initial-state event
   for  (q = n - 1;  q > 0 && reference_pos < ep [q] -> pos;  q --)
     {
      f = Frame_To_Sub (ep [q] -> frame);
      frame_last [f] = ep [q];
     }

   // Break the chain of events in each frame to skip over events
   // before  reference_pos
   for  (f = 0;  f < 6;  f ++)
     if  (reference_pos < frame_last [f] -> pos)
         frame_last [f] -> frame_pred = ep [0];
       else
         Last_Event [f] = ep [0];

   // Add the events before  reference_pos  onto the back of the
   // frame chains after incrementing positions.
   for  (i = 1;  i <= q;  i ++)
     {
      ep [i] -> pos += Sequence_Len;
      ep [i] -> Set_Frame_From_Pos ();
      f = Frame_To_Sub (ep [i] -> frame);
      ep [i] -> frame_pred = Last_Event [f];
      Last_Event [f] = ep [i];
     }

   // Sort all events into order by their  pos  field
   sort (ep . begin (), ep . end (), Event_Pos_Cmp);

   return;
  }



static void  Show_Events
    (FILE * fp)

//  Display to  fp  the contents of the global lists of events
//  pointed to by  Last_Event .

  {
   vector <Event_Node_t *> ep;
   Event_Node_t  * p;
   int  i, n;

   for  (i = 0;  i < 6;  i ++)
     for  (p = Last_Event [i];  p != NULL;  p = p -> frame_pred)
       ep . push_back (p);

   n = int (ep . size ());

   // Sort all events into order by their  pos  field
   sort (ep . begin (), ep . end (), Event_Pos_Cmp);

   fprintf (fp, "\n%8s  %-8s  %2s  %10s\n", "Position", "Type", "Fr", "Score");
   for  (i = 0;  i < n;  i ++)
     fprintf (fp, "%8d  %-8s  %+2d  %10.2f\n", ep [i] -> pos,
          Print_String (ep [i] -> e_type), ep [i] -> frame, ep [i] -> score);

   return;
  }



static void  Trace_Back
    (FILE * fp, const Event_Node_t & final_event)

//  Trace back through the list of best events starting at
//   final_event . best_pred  and output to  fp  the corresponding
//  set of genes.

  {
   Event_Node_t  * p;
   vector <Gene_t>  gene_list;
   Gene_t  gene;
   double  prev_score;
   int  f, i, j, n, rev_start;

   for  (p = final_event . best_pred;  p -> e_type != INITIAL;  p = p -> best_pred)
     {
      switch  (p -> e_type)
        {
         case  FWD_START :
           j = gene . Get_Stop_Position ();
           gene . Set_Gene_Len (2 + j - p -> pos);
           gene . Set_Score (p -> score - p -> best_pred -> score);
           gene . Set_ID (p -> id);
           if  (p -> truncated)
               gene . Set_Status_Bit (TRUNCATED_START_FLAG);
           gene_list . push_back (gene);
           gene . Clear_Status ();
           break;
         case  FWD_STOP :
           gene . Set_Stop_Position (p -> pos - 2);
           gene . Set_Frame (1 + (p -> pos % 3));
           break;
         case  REV_START :
           rev_start = p -> pos;
           prev_score = p -> score;
           if  (p -> truncated)
               gene . Set_Status_Bit (TRUNCATED_START_FLAG);
           break;
         case  REV_STOP :
           gene . Set_Stop_Position (p -> pos - 2);
           gene . Set_Frame (- (1 + (p -> pos % 3)));
           gene . Set_Gene_Len (rev_start - p -> pos);
           gene . Set_Score (prev_score - p -> score);
           gene . Set_ID (p -> id);
           gene_list . push_back (gene);
           gene . Clear_Status ();
           break;
         default :
           printf ("Bad event type = %d\n", int (p -> e_type));
           exit (EXIT_FAILURE);
        }
     }

   n = gene_list . size ();

   // Adjust stop positions to be in the range  1 .. Sequence_Len
   // and set the frame accordingly
   for  (i = 0;  i < n;  i ++)
     {
      if  (Genome_Is_Circular)
          {
           j = On_Seq_1 (gene_list [i] . Get_Stop_Position ());
           gene_list [i] . Set_Stop_Position (j);
          }
        else
          j = gene_list [i] . Get_Stop_Position ();
      f = Position_To_Frame (j);
      if  (gene_list [i] . Get_Frame () > 0)
          gene_list [i] . Set_Frame (f);
        else
          gene_list [i] . Set_Frame (-1 * f);
     }

   sort (gene_list . begin (), gene_list . end (), By_ID);

   for  (i = 0;  i < n;  i ++)
     {
      int  start, stop;

      if  (gene_list [i] . Get_Frame () > 0)
          {
           if  (Genome_Is_Circular)
               {
                stop = On_Seq_1 (gene_list [i] . Get_Stop_Position () + 2);
                start = On_Seq_1 (stop - gene_list [i] . Get_Gene_Len () - 2);
               }
             else
               {
                stop = gene_list [i] . Get_Stop_Position () + 2;
                start = stop - gene_list [i] . Get_Gene_Len () - 2;
                if  (gene_list [i] . Get_Status_Bit (TRUNCATED_START_FLAG))
                    start -= 3;
                  // move an artificial start at the beginning of the sequence
                  // off the front to indicate the gene could extend there
               }
          }
        else
          {
           if  (Genome_Is_Circular)
               {
                stop = On_Seq_1 (gene_list [i] . Get_Stop_Position ());
                start = On_Seq_1 (stop + gene_list [i] . Get_Gene_Len () + 2);
               }
             else
               {
                stop = gene_list [i] . Get_Stop_Position ();
                start = stop + gene_list [i] . Get_Gene_Len () + 2;
                if  (gene_list [i] . Get_Status_Bit (TRUNCATED_START_FLAG))
                    start += 3;
                  // move an artificial start at the end of the sequence
                  // off the back to indicate the gene could extend there
               }
          }
      fprintf (fp, "orf%05d %8d %8d %+3d %8.2f\n",
           gene_list [i] . Get_ID (),  start, stop,
           gene_list [i] . Get_Frame (),
           100.0 * gene_list [i] . Get_Score () / gene_list [i] . Get_Gene_Len ());
     }

   return;
  }



static void  Usage
    (void)

//  Print to stderr description of options and command line for
//  this program.

  {
   fprintf (stderr,
       "USAGE:  glimmer3 [options] <sequence-file> <icm-file> <tag>\n"
       "\n"
       "Read DNA sequences in <sequence-file> and predict genes\n"
       "in them using the Interpolated Context Model in <icm-file>.\n"
       "Output details go to file <tag>.detail and predictions go to\n"
       "file <tag>.predict\n"
       "\n"
       "Options:\n"
       " -A <codon-list>\n"
       " --start_codons <codon-list>\n"
       "    Use comma-separated list of codons as start codons\n"
       "    Sample format:  -A atg,gtg\n"
       "    Use -P option to specify relative proportions of use.\n"
       "    If -P not used, then proportions will be equal\n"
       " -b <filename>\n"
       " --rbs_pwm <filename>\n"
       "    Read a position weight matrix (PWM) from <filename> to identify\n"
       "    the ribosome binding site to help choose start sites\n"
       " -C <p>\n"
       " --gc_percent <p>\n"
       "    Use <p> as GC percentage of independent model\n"
       "    Note:  <p> should be a percentage, e.g., -C 45.2\n"
       " -E <filename>\n"
       " --entropy <filename>\n"
       "    Read entropy profiles from <filename>.  Format is one header\n"
       "    line, then 20 lines of 3 columns each.  Columns are amino acid,\n"
       "    positive entropy, negative entropy.  Rows must be in order\n"
       "    by amino acid code letter\n"
       " -f\n"
       " --first_codon\n"
       "    Use first codon in orf as start codon\n"
       " -g <n>\n"
       " --gene_len <n>\n"
       "    Set minimum gene length to <n>\n"
       " -h\n"
       " --help\n"
       "    Print this message\n"
       " -i <filename>\n"
       " --ignore <filename>\n"
       "    <filename> specifies regions of bases that are off \n"
       "    limits, so that no bases within that area will be examined\n"
       " -l\n"
       " --linear\n"
       "    Assume linear rather than circular genome, i.e., no wraparound\n"
       " -L <filename>\n"
       " --orf_coords <filename>\n"
       "    Use <filename> to specify a list of orfs that should\n"
       "    be scored separately, with no overlap rules\n"
       " -M\n"
       " --separate_genes\n"
       "    <sequence-file> is a multifasta file of separate genes to\n"
       "    be scored separately, with no overlap rules\n"
       " -o <n>\n"
       " --max_olap <n>\n"
       "    Set maximum overlap length to <n>.  Overlaps this short or shorter\n"
       "    are ignored.\n"
       " -P <number-list>\n"
       " --start_probs <number-list>\n"
       "    Specify probability of different start codons (same number & order\n"
       "    as in -A option).  If no -A option, then 3 values for atg, gtg and ttg\n"
       "    in that order.  Sample format:  -P 0.6,0.35,0.05\n"
       "    If -A is specified without -P, then starts are equally likely.\n"
       " -q <n>\n"
       " --ignore_score_len <n>\n"
       "    Do not use the initial score filter on any gene <n> or more\n"
       "    base long\n"
       " -r\n"
       " --no_indep\n"
       "    Don't use independent probability score column\n"
       " -t <n>\n"
       " --threshold <n>\n"
       "    Set threshold score for calling as gene to n.  If the in-frame\n"
       "    score >= <n>, then the region is given a number and considered\n"
       "    a potential gene.\n"
       " -X\n"
       " --extend\n"
       "    Allow orfs extending off ends of sequence to be scored\n"
       " -z <n>\n"
       " --trans_table <n>\n"
       "    Use Genbank translation table number <n> for stop codons\n"
       " -Z <codon-list>\n"
       " --stop_codons <codon-list>\n"
       "    Use comma-separated list of codons as stop codons\n"
       "    Sample format:  -Z tag,tga,taa\n"
       "\n");

   return;
  }



static void  Wrap_Around_Back
    (int wfr, int pos, int & gene_len, int & orf_len)

//  Set  orf_len  to the length of the complement-strand orf that
//  wraps around the end of the sequence in global  Sequence .  The
//  stop codon for the orf is at position  pos  (first base of codon
//  numbered starting at 1).   wfr  is the frame subscript of the
//  reading frame to use at the beginning of  Sequence  (i.e., it
//  allows for  Sequence_Len  not being a multiple of 3).  The
//  maximum possible orf length is  Sequence_Len - 3  rounded down
//  to the nearest multiple of 3.  Set  gene_len  to the longest
//  possible gene in that orf, looking only for starts that are completely
//  contained in the start of  Sequence .  If no starts are found,
//  set  gene_len  to  0  (even though there may be starts between
//   pos  and the end of  Sequence ).

  {
   Codon_t  codon;
   int  start_at, check_len, frame, orf_add, which;
   int  i;

   assert (pos > 0);
   check_len = pos - 1;

   start_at = -1;
   orf_add = 0;
     // this is the number of extra bases at the front of the sequence
     // to add to the orf at the back
   frame = 0;
   for  (i = 0;  i < check_len;  i ++)
     {
      codon . Shift_In (Sequence [i]);

      if  (frame == wfr)
          {
           if  (codon . Must_Be (Rev_Stop_Pattern, which))
               {
                orf_add = i - 2;
                break;
               }
             else
               orf_add = i + 1;
          }
      if  (frame == wfr && codon . Can_Be (Rev_Start_Pattern, which))
          start_at = i + 1;

      if  (frame == 2)
          frame = 0;
        else
          frame ++;
     }

   orf_len = orf_add + Sequence_Len - pos - 2;
   orf_len -= orf_len % 3;
   if  (start_at == -1)
       gene_len = 0;
     else
       gene_len = start_at + Sequence_Len - pos - 2;
   
   return;
  }



static void  Wrap_Through_Front
    (int fr, int pos, int & gene_len, int & orf_len)

//  Set  orf_len  to the length of the orf with forward frame subscript
//   fr  with stop codon at position  pos  that wraps around and begins
//  at the end of the sequence in global  Sequence .  Set  gene_len
//  to the longest possible gene in that orf.  Start looking at the
//  beginning of  Sequence  and assume there are no stops between
//  there and  pos .  If no starts are found, set  gene_len  to  0
//  (even though there may be starts between  0  and  pos in  Sequence ).

  {
   Codon_t  codon;
   int  start_at, check_len, which;
   int  i, j, s;

   assert (pos > 0);
   start_at = -1;
   s = (pos - 1) % 3;
   check_len = Sequence_Len + s - pos - 4;

   // Loop back to at most original stop codon.  Do not allow the
   // orf to overlap that stop codon.
   for  (i = 0;  i < check_len;  i += 3)
     {
      for  (j = 0;  j < 3;  j ++)
        {
         s --;
         if  (s < 0)
             s += Sequence_Len;
         codon . Reverse_Shift_In (Sequence [s]);
        }

      if  (codon . Must_Be (Fwd_Stop_Pattern, which))
          break;
      if  (codon . Can_Be (Fwd_Start_Pattern, which))
          start_at = i + 3;

     }

   orf_len = i + 3 * ((pos - 1) / 3);
   if  (start_at == -1)
       gene_len = 0;
     else
       gene_len = start_at + 3 * ((pos - 1) / 3);
   
   return;
  }



void  Event_Node_t :: Set_Frame_From_Pos
    (void)

// Set the  frame  field of this node to the frame corresponding
// to the value in the  pos  field  but retaining the sign of
// the  frame  field.

  {
   int  f;

   assert (pos > 2);

   f = 1 + (pos % 3);
   if  (frame > 0)
       frame = f;
     else
       frame = -1 * f;

   return;
  }



