//    Programmer:  Arthur L. Delcher
//          File:  simple-score.cc
//  Last Updated:  Fri Jan 13 10:20:08 EST 2006
//
//  Compute scores for each sequence in an input multi-fasta
//  file (read from stdin) using both positive and negative
//  ICMs, in files named on the command line.


#include  "simple-score.hh"


static char  * Pos_Model_Path;
  // Name of file containing the positive model
static char  * Neg_Model_Path;
  // Name of file containing the negative model
static bool  Include_Per_Base = false;
  // If true also output the length of each string
  // and the per-base net score
static bool  Include_Total = false;
  // If true also output the total positive model
  // and negative model scores
static bool  Use_Null_Neg_Model = false;
  // If true then negative model is ignored and value
  // from it is automatically zero


//**ALD  Gets rid of make undefined reference error
int  Unused = Filter ('a');



int  main
    (int argc, char * argv [])

  {
   ICM_t  pos_model, neg_model;
   char  * string = NULL, * tag = NULL;
   long int  string_size = 0, tag_size = 0;
   int  string_num = 0;

   Parse_Command_Line (argc, argv);

   pos_model . Read (Pos_Model_Path);
   fprintf (stderr, "Positive Model = %s\n", Pos_Model_Path);
   fprintf (stderr, "  len = %d  depth = %d  periodicity = %d\n",
        pos_model . Get_Model_Len (),
        pos_model . Get_Model_Depth (),
        pos_model . Get_Periodicity ());
   if  (Use_Null_Neg_Model)
       fprintf (stderr, "No negative model\n");
     else
       {
        neg_model . Read (Neg_Model_Path);
        fprintf (stderr, "Negative Model = %s\n", Neg_Model_Path);
        fprintf (stderr, "  len = %d  depth = %d  periodicity = %d\n",
             neg_model . Get_Model_Len (),
             neg_model . Get_Model_Depth (),
             neg_model . Get_Periodicity ());
       }

   while  (Read_String (stdin, string, string_size, tag, tag_size))
     {
      char  * token;
      double  pos_score, neg_score;
      int  len;

      string_num ++;
      len = strlen (string);
      token = strtok (tag, " \t\n");

      pos_score = pos_model . Score_String (string, len, 1);
      if  (! Use_Null_Neg_Model)
          neg_score = neg_model . Score_String (string, len, 1);
        else
          neg_score = 0.0;

      printf ("%-20s ", token);
      if  (Include_Total)
          printf (" %11.4f %11.4f", pos_score, neg_score);
      printf (" %11.4f", pos_score - neg_score);
      if  (Include_Per_Base)
          {
           double  pb;

           pb = (len > 0 ? (pos_score - neg_score) / len : 0.0);
           printf (" %9d %8.4f", len, pb);
          }
      putchar ('\n');
     }

   fprintf (stderr, "Done scoring, exiting.\n");
   return  0;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   bool  errflg = false;
   int  ch, option_index = 0;
   static struct option  long_options [] = {
        {"per_base", 1, 0, 'B'},
        {"help", 0, 0, 'h'},
        {"no_neg_model", 1, 0, 'N'},
        {"total", 1, 0, 'T'},
        {0, 0, 0, 0}
      };

   optarg = NULL;

   while  (! errflg && ((ch = getopt_long (argc, argv,
        "BhNT", long_options, & option_index)) != EOF))
     switch  (ch)
       {
        case  'B' :
          Include_Per_Base = true;
          break;

        case  'h' :
          errflg = true;
          break;

        case  'N' :
          Use_Null_Neg_Model = true;
          break;

        case  'T' :
          Include_Total = true;
          break;

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = TRUE;
       }

   if  (errflg || (Use_Null_Neg_Model && optind > argc - 1)
          || (! Use_Null_Neg_Model && optind != argc - 2))
       {
        Usage ();
        exit (EXIT_FAILURE);
       }

   Pos_Model_Path = argv [optind ++];
   if  (! Use_Null_Neg_Model)
       Neg_Model_Path = argv [optind ++];

   return;
  }



int  Read_String
    (FILE * fp, char * & s, long int & s_size, char * & tag,
     long int & tag_size)

//  Read next string from  fp  (assuming FASTA format) into  s [0 .. ]
//  which has  s_size  characters.  Allocate extra memory if needed
//  and adjust  s_size  accordingly.  Return  TRUE  if successful,  FALSE
//  otherwise (e.g., EOF).  Put FASTA header line into  tag [0 .. ]
//  (and adjust  tag_size  if needed).

  {
   int  ch, ct;

   while  ((ch = fgetc (fp)) != EOF && ch != '>')
     ;

   if  (ch == EOF)
       return  FALSE;

   ct = 0;
   while  ((ch = fgetc (fp)) != EOF && ch != '\n' && isspace (ch))
     ;
   if  (ch == EOF)
       return  FALSE;
   if  (ch != '\n' && ! isspace (ch))
       ungetc (ch, fp);
   while  ((ch = fgetc (fp)) != EOF && ch != '\n')
     {
      if  (ct >= tag_size - 1)
          {
           tag_size += INCR_SIZE;
           tag = (char *) Safe_realloc (tag, tag_size);
          }
      tag [ct ++] = char (ch);
     }
   tag [ct ++] = '\0';

   ct = 0;
   while  ((ch = fgetc (fp)) != EOF && ch != '>')
     {
      if  (isspace (ch))
          continue;
      //AT: we skip everything that is not 'actg', otherwise
      //Filter() will create runs of 'c' from runs of 'n' later on,
      //which is a very bad thing for CA scaffolds. Our solution
      //will create chimerich k-mers and shift frame, but for taxonomic
      //signatures it should be acceptable (assuming gaps are uncorellated
      //with their flanking k-mers).
      //The cleaner solution would be to split string and add scores in the
      //main() loop above.
      ch = tolower(ch);
      if  (ch != Filter(ch))
          continue;

      if  (ct >= s_size - 1)
          {
           s_size += INCR_SIZE;
           s = (char *) Safe_realloc (s, s_size);
          }
      s [ct ++] = char (ch);
     }
   s [ct ++] = '\0';

   if  (ch == '>')
       ungetc (ch, fp);

   return  TRUE;
  }



static void  Usage
    (void)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
       "USAGE:  simple-score [options] <pos-model> <neg-model> < input-file\n"
       "\n"
       "Read sequences from  stdin  and score each using the ICM's in\n"
       "<pos-model> and <neg-model> .  Output scores to  stdout\n"
       "\n"
       "Options:\n"
       " -B\n"
       " --per_base\n"
       "    Also output the length of each string and the net score per\n"
       "    character as the last two fields\n"
       " -h\n"
       " --help\n"
       "    Print this message\n"
       " -N\n"
       " --no_neg_model\n"
       "    No negative model, i.e., negative model score is constant zero\n"
       "    Omit <neg-model> parameter in this case\n"
       " -T\n"
       " --total\n"
       "    Output the total positive and negative scores for each\n"
       "    string before the net score, which is the difference\n"
       "    between them\n"
       "\n");

   return;
  }



