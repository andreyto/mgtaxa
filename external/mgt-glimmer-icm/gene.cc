//  A. L. Delcher
//
//  File:  gene.cc
//
//  Last Modified:  23 October 2003
//
//  DNA- and gene-related routines.


#include "delcher.hh"
#include "gene.hh"


static const char  COMPLEMENT_TABLE []
    = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn"
      " nnnnnnnnn*nn-.nnnnnnnnnnnnnnnnn"
      "nTVGHNNCDNNMNKNNNNYSANBWNRNnnnn_"
      "ntvghnncdnnmnknnnnysanbwnrnnnnnn";

static const char  CONVERSION_STRING [] = "acgtn";

static const int  FINAL_STATE = 6;

static const int  Transition [FINAL_STATE] [5]
  //    a  c  g  t  n
    = {
       {1, 1, 1, 1, 1},    // state 0 transitions
       {2, 2, 2, 3, 3},    // state 1 transitions
       {0, 0, 0, 0, 0},    // state 2 transitions
       {4, 0, 5, 0, 4},    // state 3 transitions
       {6, 1, 6, 1, 6},    // state 4 transitions
       {6, 1, 1, 1, 6}     // state 5 transitions
      };
  // Encodes FSA to recognize in-frame stop codons



bool  Codon_t :: Can_Be
    (const vector <Codon_t> & a, int & which)

//  Return  true  iff this codon could match any of the codons
//  in  a .  "could match" means this codon could be a string
//  that equals a string that an entry of  a  could be.
//  Set  which  to the subscript of the first matching entry in  a ,
//  or else -1 if there is no match.

  {
   unsigned int  x;
   int  i, n;

   n = a . size ();
   for  (i = 0;  i < n;  i ++)
     {
      x = data & a [i] . data;
      if  ((x & 0xf00) && (x & 0xf0) && (x & 0x0f))
          {
           which = i;
           return  true;
          }
     }

   which = -1;
   return  false;
  }



bool  Codon_t :: Must_Be
    (const vector <Codon_t> & a, int & which)

//  Return  true  iff this codon must match one of the codons
//  in  a .  "must match" means that any string this codon
//  could be equals a string that an entry of  a  could be.
//  Set  which  to the subscript of the first matching entry in  a .
//  or else -1 if there is no match.

  {
   int  i, n;

   n = a . size ();
   for  (i = 0;  i < n;  i ++)
     if  ((data & a [i] . data) == data
             && (data & 0xf00) && (data & 0xf0) && (data & 0x0f))
         {
          which = i;
          return  true;
         }

   which = -1;
   return  false;
  }



void  Codon_t :: Reverse_Complement
    (void)

//  Set this codon to the reverse complement of
//  the value in it.  E.g., "atg" changes to "cat"

  {
   unsigned int  x = 0x0;
   int  i;

   for  (i = 0;  i < 12;  i ++)
     {
      x = (x << 1) | (data & 0x1);
      data >>= 1;
     }

   data = x;
   return;
  }



void  Codon_t :: Reverse_Shift_In
    (char ch)

//  Add  ch  onto the left of this codon, shifting the rightmost
//  character off the right end.

  {
   data = (data & reverse_shift_mask) >> 4;
   data |= (Ch_Mask (ch) << 8);

   return;
  }



void  Codon_t :: Set_From
    (const char * s)

//  Set this codon to the equivalent of the characters in string  s .

  {
   int  i;

   Clear ();
   for  (i = 0;  i < 3 && * s != '\0';  i ++)
     Shift_In (* s ++);

   return;
  }



void  Codon_t :: Shift_In
    (char ch)

//  Add  ch  onto the right of this codon, shifting the leftmost
//  character off the left end.

  {
   data = (data & shift_mask) << 4;
   data |= Ch_Mask (ch);

   return;
  }



double  PWM_t :: Column_Score
    (char ch, int j)  const

//  Return the entry for character  ch  in column subscript  j .
//  If  ch  is not a valid nucleotide, return  0 .

  {
   int  i;

   i = Nucleotide_To_Subscript (ch);
   if  (i < 0)
       return  0.0;
     else
       return  col [j] . p [i];
  }



void  PWM_t :: Counts_To_Prob
    (void)

//  Convert the counts in this PWM to probabilities by dividing
//  each entry by the sum of its column.  Convert zero probabilities
//  to a small positive value.

  {
   const double  ZERO_EQUIV = 1e-6;
   int  width;
   int  i, j;

   width = col . size ();

   for  (j = 0;  j < width;  j ++)
     {
      double  sum = 0.0;
      int  zero_count = 0;

      for  (i = 0;  i < 4;  i ++)
        {
         sum += col [j] . p [i];
         if  (col [j] . p [i] == 0.0)
             zero_count ++;
        }

      if  (sum > 0.0)
          for  (i = 0;  i < 4;  i ++)
            {
             col [j] . p [i] /= sum;
             if  (col [j] . p [i] == 0)
                 col [j] . p [i] = ZERO_EQUIV;
               else
                 col [j] . p [i] /= (1.0 + zero_count * ZERO_EQUIV);
            }
     }

   return;
  }



void  PWM_t :: Make_Log_Odds_WRT_GC
    (double gc_frac)

//  Convert the probabilities in this PWM to log odds
//  by subtracting the log of the base probabilities implied by
//  a GC portion of  gc_frac .

  {
   double  at_log, gc_log;
   int  width;
   int  j;

   if  (gc_frac <= 0.0)
       SIMPLE_THROW ("ERROR:  Non-positive gc-fraction");

   gc_log = log (0.5 * gc_frac);
   at_log = log (0.5 * (1.0 - gc_frac));

   width = col . size ();

   for  (j = 0;  j < width;  j ++)
     {
      col [j] . p [0] -= at_log;
      col [j] . p [1] -= gc_log;
      col [j] . p [2] -= gc_log;
      col [j] . p [3] -= at_log;
     }

   return;
  }



void  PWM_t :: Print
    (FILE * fp)

//  Print the contents of this PWM to  fp .

  {
   char  * tag = "acgt";
   int  width;
   int  i, j;

   width = col . size ();

   fprintf (fp, "PWM:\n");
   for  (i = 0;  i < 4;  i ++)
     {
      fprintf (fp, " %c", tag [i]);
      for  (j = 0;  j < width;  j ++)
        fprintf (fp, " %12.5e", col [j] . p [i]);
      fputc ('\n', fp);
     }

   return;
  }



void  PWM_t :: Probs_To_Logs
    (void)

//  Convert the probabilities in this PWM to natural logarithms.

  {
   int  width;
   int  i, j;

   width = col . size ();

   for  (j = 0;  j < width;  j ++)
     for  (i = 0;  i < 4;  i ++)
       if  (col [j] . p [i] <= 0.0)
           SIMPLE_THROW ("ERROR:  Log of non-positive value");
         else
           col [j] . p [i] = log (col [j] . p [i]);

   return;
  }



bool  PWM_t :: Read
    (FILE * fp)

//  Set this PWM to values read in from  fp , which must
//  already be open.  Return  true  if the values are read
//  successfully, and  false  otherwise.

  {
   char  tag [1000];
   int  width;
   int  i, j;

   fscanf (fp, "%d", & width);
   if  (width <= 0)
       {
        fprintf (stderr, "ERROR:  Bad width = %d in PWM\n", width);
        return  false;
       }

   col  . resize (width);

   for  (i = 0;  i < 4;  i ++)
     {
      double  x;

      fscanf (fp, "%s", tag);   // skip tag in first column
      for  (j = 0;  j < width;  j ++)
        {
         fscanf (fp, "%lf", & x);
         col [j] . p [i] = x;
        }
     }

   return  true;
  }



PWM_t &  PWM_t :: operator =
    (const PWM_t & src)

//  Assign this PWM the value in  src .

  {
   int  width;
   int  i, j;

   if  (this != & src)
       {
        width = src . col . size ();

        col . clear ();
        col . resize (width);
        for  (j = 0;  j < width;  j ++)
          for  (i = 0;  i < 4;  i ++)
            col [j] . p [i] = src . col [j] . p [i];
       }

   return  (* this);
  }



unsigned int  Gene_t :: Get_Status_Bit
    (unsigned int u)  const

//  Return  0  if the status bit(s) matching pattern  u  are
//  all zero; otherwise, return  1 .

  {
   if  ((status & u) == 0)
       return  0;
     else
       return  1;
  }



bool  By_ID
    (const Gene_t & a, const Gene_t & b)

//  Return true iff  a 's  id  field is less than  b 's.

  {
   return  (a . Get_ID () < b . Get_ID ());
  }



unsigned  Ch_Mask
    (char Ch)

/* Returns a bit mask representing character  Ch . */

  {
   switch  (tolower (Ch))
     {
      case  'a' :
        return  0x1;
      case  'c' :
        return  0x2;
      case  'g' :
        return  0x4;
      case  't' :
        return  0x8;
      case  'r' :     // a or g
        return  0x5;
      case  'y' :     // c or t
        return  0xA;
      case  's' :     // c or g
        return  0x6;
      case  'w' :     // a or t
        return  0x9;
      case  'm' :     // a or c
        return  0x3;
      case  'k' :     // g or t
        return  0xC;
      case  'b' :     // c, g or t
        return  0xE;
      case  'd' :     // a, g or t
        return  0xD;
      case  'h' :     // a, c or t
        return  0xB;
      case  'v' :     // a, c or g
        return  0x7;
      case  'n' :     // anything
        return  0xF;
      default :       // nothing
        return  0x0;
     }
  }



int  Char_Sub
    (char ch)

//  Return a subscript corresponding to character  ch .

  {
   char  * p;

   p = strchr (CONVERSION_STRING, tolower (ch));
   if  (p == NULL)
       return  4;

   return  p - CONVERSION_STRING;
  }



char  Codon_Translation
    (const char * c, int transl_tabl)

//  Return the character code for the amino acid that
//  the triplet starting at  c  translates to using
//  NCBI translation table  t .  Return  'X'  if the
//  triplet has non-acgt characters and return '*' for
//  a stop codon.

  {

   int  i, j, sub = 0;

   for  (i = 0;  i < 3;  i ++)
     {
      j = Nucleotide_To_Subscript (c [i]);
      if  (j < 0)
          return  'X';
      sub = 4 * sub + j;
     }

   switch  (transl_tabl)
     {
      case  0 :  // unspecified code--assume standard
      case  1 :  // The Standard Code
      case  11 : // The Bacterial and Plant Plastid Code
        return  CODON_XLATE_TABLE_1 [sub];
      case  2 :  // The Vertebrate Mitochondrial Code
        return  CODON_XLATE_TABLE_2 [sub];
      case  3 :  // The Yeast Mitochondrial Code
        return  CODON_XLATE_TABLE_3 [sub];
      case  4 :  // The Mold, Protozoan, and Coelenterate Mitochondrial Code
                 //   and the Mycoplasma/Spiroplasma Code
        return  CODON_XLATE_TABLE_4 [sub];
      case  5 :  // The Invertebrate Mitochondrial Code
        return  CODON_XLATE_TABLE_5 [sub];
      case  6 :  // The Ciliate, Dasycladacean and Hexamita Nuclear Code
        return  CODON_XLATE_TABLE_6 [sub];
      case  9 :  // The Echinoderm and Flatworm Mitochondrial Code
        return  CODON_XLATE_TABLE_9 [sub];
      case  10 :  // The Euplotid Nuclear Code
        return  CODON_XLATE_TABLE_10 [sub];
      case  12 :  // The Alternative Yeast Nuclear Code
        return  CODON_XLATE_TABLE_12 [sub];
      case  13 :  // The Ascidian Mitochondrial Code
        return  CODON_XLATE_TABLE_13 [sub];
      case  14 :  // The Alternative Flatworm Mitochondrial Code
        return  CODON_XLATE_TABLE_14 [sub];
      case  15 :  // Blepharisma Nuclear Code
        return  CODON_XLATE_TABLE_15 [sub];
      case  16 :  // Chlorophycean Mitochondrial Code
        return  CODON_XLATE_TABLE_16 [sub];
      case  21 :  // Trematode Mitochondrial Code
        return  CODON_XLATE_TABLE_21 [sub];
      case  22 :  // Scenedesmus obliquus mitochondrial Code
        return  CODON_XLATE_TABLE_22 [sub];
      case  23 :  // Thraustochytrium Mitochondrial Code
        return  CODON_XLATE_TABLE_23 [sub];

      default :
        sprintf (Clean_Exit_Msg_Line,
             "ERROR:  Bad translation table = %d", transl_tabl);
        SIMPLE_THROW (Clean_Exit_Msg_Line);
     }
  }



char  Complement
    (char ch)

// Returns the DNA complement of  ch

  {
   return  COMPLEMENT_TABLE [unsigned (ch)];
  }



void  Counts_To_Entropy_Profile
    (int count [26], double ep [20])

//  Convert the amino-acid counts in  count  to their
//  entropy profile in  ep .

  {
   double  sum;
   int  i, j;

   sum = 0.0;
   for  (i = 0;  i < 26;  i ++)
     if  (IS_AMINO [i])
         sum += count [i];

   if  (sum == 0.0)
       {
        for  (j = 0;  j < 20;  j ++)
          ep [j] = 0.0;
        return;
       }

   for  (i = j = 0;  i < 26;  i ++)
     if  (IS_AMINO [i])
         ep [j ++] = count [i] / sum;

   sum = 0.0;
   for  (j = 0;  j < 20;  j ++)
     {
      if  (ep [j] <= 0.0)
          ep [j] = 0.0;
        else
          ep [j] = -1.0 * ep [j] * log (ep [j]);
      sum += ep [j];
     }

   for  (j = 0;  j < 20;  j ++)
     ep [j] /= sum;

   return;
  }



int  Filter
    (char Ch)

//  Return a single  a, c, g or t  for  Ch .

  {
   switch  (tolower (Ch))
     {
      case  'a' :
      case  'c' :
      case  'g' :
      case  't' :
        return  Ch;
      case  'r' :     // a or g
        return  'g';
      case  'y' :     // c or t
        return  'c';
      case  's' :     // c or g
        return  'c';
      case  'w' :     // a or t
        return  't';
      case  'm' :     // a or c
        return  'c';
      case  'k' :     // g or t
        return  't';
      case  'b' :     // c, g or t
        return  'c';
      case  'd' :     // a, g or t
        return  'g';
      case  'h' :     // a, c or t
        return  'c';
      case  'v' :     // a, c or g
        return  'c';
      default :       // anything
        return  'c';
    }
  }



void  Find_Stop_Codons
    (const char * X, int T, int Stop [])

//  Set  Stop [0 .. 6]  TRUE  or  FALSE   according to whether
//  X [1 .. T] has a stop codon in the corresponding reading frame.
//  Stop [6]  is always set  FALSE .

  {
   unsigned  Codon;
   int  i;

   for  (i = 0;  i < 7;  i ++)
     Stop [i] = 0;

   if  (T < 3)
       return;

   Codon = Ch_Mask (X [1]) << 4 | Ch_Mask (X [2]);

   for  (i = 3;  i <= T;  i ++)
     {
      Codon = (Codon & SHIFT_MASK) << 4;
      Codon |= Ch_Mask (X [i]);
      
      if  (Is_Forward_Stop (Codon))
          Stop [i % 3] = TRUE;
      if  (Is_Reverse_Stop (Codon))
          Stop [3 + i % 3] = TRUE;
     }

   return;
  }



int  First_In_Frame_Stop
    (char * s, int frame)

//  Return the subscript of the first base of the first
//  in-frame stop codon in string  s  whose first base
//  is in frame  frame .  Return the length of  s  if
//  there is no in-frame stop codon.

  {
   int  i, state;

   state = frame;
   for  (i = 0;  s [i] != '\0' && state < FINAL_STATE;  i ++)
     state = Transition [state] [Char_Sub (s [i])];

   if  (state == FINAL_STATE)
       return  i - 3;

   return  i;
  }



void  Forward_Strand_Transfer
    (string & t, const string & s, int start, int len)

//  Copy the sequence starting at subscript  start  on  s
//  with length  len  to string  t .  Wrap circularly around
//  the end of  s  if necessary.

  {
   int  i, n;

   t . resize (len);
   n = s . length ();
   assert (0 <= start && start < n);

   for  (i = 0;  i < len;  i ++)
     {
      t [i] = s [start];
      start ++;
      if  (start >= n)
          start = 0;
     }

   return;
  }



int  Is_Forward_Start
    (unsigned Codon)

//  Return  TRUE  iff bit pattern  Codon  represents a start codon in the
//  forward direction.

  {
   return  (
            (Codon & ATG_MASK) == Codon
//               || (Codon & CTG_MASK) == Codon
               || (Codon & GTG_MASK) == Codon
               || (Codon & TTG_MASK) == Codon
           );
  }



int  Is_Forward_Stop
    (unsigned Codon)

//  Return  TRUE  iff bit pattern  Codon  represents a stop codon in the
//  forward direction.

  {
   return  (
            (Codon & TAA_MASK) == Codon
               || (Codon & TAG_MASK) == Codon
               || (Codon & TGA_MASK) == Codon
           );
  }



int  Is_Reverse_Start
    (unsigned Codon)

//  Return  TRUE  iff bit pattern  Codon  represents a start codon in the
//  reverse direction.

  {
   return  (
            (Codon & CAT_MASK) == Codon
//               || (Codon & CAG_MASK) == Codon
               || (Codon & CAC_MASK) == Codon
               || (Codon & CAA_MASK) == Codon
           );
  }



int  Is_Reverse_Stop
    (unsigned Codon)

//  Return  TRUE  iff bit pattern  Codon  represents a stop codon in the
//  reverse direction.

  {
   return  (
            (Codon & TTA_MASK) == Codon
               || (Codon & CTA_MASK) == Codon
               || (Codon & TCA_MASK) == Codon
           );
  }



int  Is_Start
    (const char * S)

/* Return  TRUE  iff  S  is a start codon. */

  {
   return  (
            strncmp (S, "atg", 3) == 0
//               || strncmp (S, "ctg", 3) == 0
               || strncmp (S, "gtg", 3) == 0
               || strncmp (S, "ttg", 3) == 0
           );
  }



int  Is_Stop
    (char * S)

/* Return  TRUE  iff  S  is a stop codon. */

  {
   return  (
            strncmp (S, "taa", 3) == 0
               || strncmp (S, "tag", 3) == 0
               || strncmp (S, "tga", 3) == 0
           );
  }



int  Nucleotide_To_Subscript
    (char ch)

//  Return the subscript that corresponds to nucleotide  ch .
//  Return  -1  if  ch  is not a, c, g or t.

  {
   switch  (tolower (ch))
     {
      case  'a' :
        return  0;
      case  'c' :
        return  1;
      case  'g' :
        return  2;
      case  't' :
        return  3;
      default :
        return -1;
     }
  }



int  Read_String
    (FILE * fp, char * & T, long int & Size, char Name [],
     int Partial)

/* Read next string from  fp  (assuming FASTA format) into  T [1 ..]
*  which has  Size  characters.  Allocate extra memory if needed
*  and adjust  Size  accordingly.  Return  TRUE  if successful,  FALSE
*  otherwise (e.g., EOF).  Partial indicates if first line has
*  numbers indicating a subrange of characters to read. */

  {
   char  * P, Line [MAX_LINE];
   long int  Len, Lo, Hi;
   int  Ch, Ct;

   while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
     ;

   if  (Ch == EOF)
       return  FALSE;

   fgets (Line, MAX_LINE, fp);
   Len = strlen (Line);
   assert (Len > 0 && Line [Len - 1] == '\n');
   P = strtok (Line, " \t\n");
   if  (P != NULL)
       strcpy (Name, P);
     else
       Name [0] = '\0';
   Lo = 0;  Hi = LONG_MAX;
   if  (Partial)
       {
        P = strtok (NULL, " \t\n");
        if  (P != NULL)
            {
             Lo = strtol (P, NULL, 10);
             P = strtok (NULL, " \t\n");
             if  (P != NULL)
                 Hi = strtol (P, NULL, 10);
            }
        assert (Lo <= Hi);
       }

   Ct = 0;
   T [0] = '\0';
   Len = 1;
   while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
     {
      if  (isspace (Ch))
          continue;

      Ct ++;
      if  (Ct < Lo || Ct > Hi)
          continue;

      if  (Len >= Size)
          {
           Size += INCR_SIZE;
           T = (char *) Safe_realloc (T, Size);
          }
      Ch = tolower (Ch);
      switch  (Ch)
        {
         case  'a' :
         case  'c' :
         case  'g' :
         case  't' :
         case  's' :
         case  'w' :
         case  'r' :
         case  'y' :
         case  'm' :
         case  'k' :
         case  'b' :
         case  'd' :
         case  'h' :
         case  'v' :
         case  'n' :
           break;
         default :
           fprintf (stderr, "Unexpected character `%c\' in string %s\n",
                                 Ch, Name);
           Ch = 'n';
        }
      T [Len ++] = char (Ch);
     }

   T [Len] = '\0';
   if  (Ch == '>')
       ungetc (Ch, fp);

   return  TRUE;
  }



void  Reverse_Complement
    (char * s)

//  Set string  s  to its DNA Watson-Crick reverse complement

  {
   int  i, j, n;

   n = strlen (s);
   for  (i = 0, j = n - 1;  i < j;  i ++, j --)
     {
      char  ch;

      ch = s [j];
      s [j] = Complement (s [i]);
      s [i] = Complement (ch);
     }

   if  (i == j)
       s [i] = Complement (s [i]);

   return;
  }



void  Reverse_Complement
    (string & s)

//  Set string  s  to its DNA Watson-Crick reverse complement

  {
   int  i, j, n;

   n = s . length ();
   for  (i = 0, j = n - 1;  i < j;  i ++, j --)
     {
      char  ch;

      ch = s [j];
      s [j] = Complement (s [i]);
      s [i] = Complement (ch);
     }

   if  (i == j)
       s [i] = Complement (s [i]);

   return;
  }



void  Reverse_Strand_Transfer
    (string & t, const string & s, int start, int len)

//  Copy the reverse-complement sequence starting at subscript
//   start  on  s   with length  len  to string  t .  Wrap circularly
//  around the end of  s  if necessary.

  {
   int  i, n;

   t . resize (len);
   n = s . length ();
   assert (0 <= start && start < n);

   for  (i = 0;  i < len;  i ++)
     {
      t [i] = Complement (s [start]);
      start --;
      if  (start < 0)
          start = n - 1;
     }

   return;
  }



void  Set_Stop_Codons_By_Code
    (vector <const char *> & stop_codon, int code, bool & errflg)

//  Put stop codon values in  stop_codon  according
//  to the Genbank translation table code  code .  If code is
//  not recognized, then set  errflg  true  and leave  stop_codon
//  empty.

  {
   stop_codon . clear ();
   switch  (code)
     {
      case  1 :  // Standard
      case  11 :  // Bacterial
      case  12 :  // Alternative yeast
        stop_codon . push_back ("taa");
        stop_codon . push_back ("tag");
        stop_codon . push_back ("tga");
        break;
      case  2 :  // Vertebrate mitochondrial
        stop_codon . push_back ("taa");
        stop_codon . push_back ("tag");
        stop_codon . push_back ("aga");
        stop_codon . push_back ("agg");
        break;
      case  3 :  // Yeast mitochondrial
      case  4 :  // Mold mitochondrial
      case  5 :  // Invertebrate mitochondrial
      case  9 :  // Echinoderm mitochondrial
      case  10 :  // Euplotid nuclear
      case  13 :  // Ascidian mitochondrial
      case  21 :  // Trematode mitochondrial
        stop_codon . push_back ("taa");
        stop_codon . push_back ("tag");
        break;
      case  6 :  // Ciliate nuclear
        stop_codon . push_back ("tga");
        break;
      case  14 :  // Flatworm mitochondrial
        stop_codon . push_back ("tag");
        break;
      case  15 :  // Blepharisma mitochondrial
      case  16 :  // Chlorophycean mitochondrial
        stop_codon . push_back ("taa");
        stop_codon . push_back ("tga");
        break;
      case  22 :  // Scenedesmus obliquus mitochondrial
        stop_codon . push_back ("taa");
        stop_codon . push_back ("tga");
        stop_codon . push_back ("tca");
        break;
      case  23 :  // Thraustochytrium aureum mitochondrial
        stop_codon . push_back ("taa");
        stop_codon . push_back ("tag");
        stop_codon . push_back ("tga");
        stop_codon . push_back ("tta");
        break;
      default :
        fprintf (stderr, "ERROR:  Unknown translation-table number = %d\n",
             code);
        errflg = true;
     }

   return;
  }



