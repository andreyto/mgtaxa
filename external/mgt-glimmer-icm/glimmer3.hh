//  A. L. Delcher
//
//  File:  glimmer3.hh
//
//  Last Modified:  Tue May  9 10:25:40 EDT 2006
//
//  Declarations for  Glimmer3



#ifndef  __GLIMMER3_HH_INCLUDED
#define  __GLIMMER3_HH_INCLUDED


#include  "delcher.hh"
#include  "fasta.hh"
#include  "gene.hh"
#include  "icm.hh"


// Default values of global variables

static const bool  DEFAULT_GENOME_IS_CIRCULAR = true;
static const int  DEFAULT_MIN_GENE_LEN = 100;
static const int  DEFAULT_MAX_OLAP_BASES = 30;
static const int  DEFAULT_RIBOSOME_WINDOW_SIZE = 20;
static const double  DEFAULT_START_PROB []
     = {0.60, 0.30, 0.10};
static const int  DEFAULT_THRESHOLD_SCORE = 30;
static const int  DEFAULT_USE_FIRST_START_CODON = false;
static const int  DEFAULT_USE_INDEPENDENT_SCORE = true;
static const int  HI_SCORE = 100;
  // the highest possible ICM score for an orf
static const double  LONG_ORF_SCORE_PER_BASE = 0.03;
  // artificially good score value for sufficiently long orfs
  //**ALD Should maybe change to a lower value like 0.01 ??


enum  Event_t
  {INITIAL, FWD_START, FWD_STOP, REV_START, REV_STOP, TERMINAL};


struct  Event_Node_t
  {
   int  id : 24;
   int  frame : 3;
   unsigned  is_first_start : 1;
   unsigned  disqualified : 1;
   unsigned  truncated : 1;
   Event_t  e_type;
   int  pos, pwm_sep;
     // pos is the last base of the codon, numbered starting at 1
   double  score, pwm_score;
   Event_Node_t  * frame_pred;
   Event_Node_t  * best_pred;

   Event_Node_t  ()   // default constructor
     { is_first_start = disqualified = truncated = 0; }

   void  Set_Frame_From_Pos
       (void);
  };


static bool  Event_Pos_Cmp
    (Event_Node_t * const & a, Event_Node_t * const & b)
  { return  (a -> pos < b -> pos); }


struct  Orf_Pos_t
  {
   int  start, stop, dir;
   char  * tag;
  };


struct  Range_t
  {
   int  lo, hi;
  };


static bool  Range_Cmp
    (const Range_t & a, const Range_t & b)
  { return  (a . lo < b . lo); }


struct  Position_t
  {
   int  lo, hi, max_prev;
  };


struct  Start_t
  {
   int  j, pos;
   double  score, rate;
   int  which : 8;
   unsigned  truncated : 1;  
   bool  first;
  };



static void  Add_Events
    (const Orf_t & orf, vector <Start_t> & start_list, int id);
static void  Add_PWM_Score
    (Event_Node_t * p);
static void  All_Frame_Score
    (const string & s, int offset, int frame, vector <double> & af);
static void  Clear_Events
    (void);
static void  Complement_Transfer
    (string & buff, const string & s, int lo, int hi);
static void  Disqualify
    (Event_Node_t * p, int cutoff);
static void  Do_Fwd_Stop_Codon
    (int i, int frame, int prev_fwd_stop [3], int first_fwd_start [3],
     int first_fwd_stop [3], int first_base, bool hit_ignore,
     vector <Orf_t> & orf_list);
static void  Echo_General_Settings
    (FILE * fp);
static void  Echo_Specific_Settings
    (FILE * fp, int len);
static double  Entropy_Distance_Ratio
    (int start, int len, int fr);
static int  Find_Uncovered_Position
    (vector <Event_Node_t *> ep);
static void  Find_Orfs
    (vector <Orf_t> & orf_list);
static void  Find_Stops_Reverse
    (const string & s, int len, vector <bool> & has_stop);
static void  Finish_Orfs
    (bool use_wraparound, const int prev_rev_stop [3],
     const int last_rev_start [3], int last_position,
     vector <Orf_t> & orf_list);
static void  Fix_Wrap
    (int & p, const int n);
static int  Frame_To_Sub
    (int f);
static void  Get_Ignore_Regions
    (void);
static void  Get_Orf_Pos_List
    (void);
static void  Handle_First_Forward_Stop
    (int fr, int pos, int start_pos, int first_base, int & gene_len,
     int & orf_len, bool use_wraparound);
static void  Handle_First_Reverse_Stop
    (int pos, int last_start, int & gene_len, int & orf_stop, bool hit_ignore);
static void  Handle_Last_Reverse_Stop
    (int fr, const int prev_rev_stop [3], const int last_rev_start [3],
     int & gene_len, int & orf_len, bool use_wraparound, int last_position);
static void  Initialize_Terminal_Events
    (Event_Node_t & first_event, Event_Node_t & final_event,
     Event_Node_t * best_event [6], Event_Node_t * last_event [6]);
static void  Integerize_Scores
    (const vector <double> ds, int hi_score, const vector <bool> set_zero,
    vector <int> & is);
static double  Olap_Score_Adjustment
    (int lo, int hi, int f1, int f2);
static int  On_Seq_0
    (int i);
static int  On_Seq_1
    (int i);
static void  Output_Extra_Start_Info
    (FILE * fp, int i, int lo, int hi, int frame,
     vector <Start_t> & start_list);
static void  Parse_Command_Line
    (int argc, char * argv []);
template  <class DT>
static void  Permute_By_Frame
    (vector <DT> & v, int frame);
int  Position_To_Frame
    (int p);
static void  Print_Comma_Separated_Strings
    (const vector <const char *> & v, FILE * fp);
static void  Print_Headings
    (FILE * fp);
static void  Print_Orflist_Headings
    (FILE * fp);
static const char  * Print_String
    (Event_t e);
static void  Prob_To_Logs
    (vector <double> & v);
static void  Process_Events
    (void);
static void  Process_Fwd_Start_Event
    (Event_Node_t * ep);
static void  Process_Fwd_Stop_Event
    (Event_Node_t * ep);
static void  Process_Initial_Event
    (Event_Node_t * ep);
static void  Process_Rev_Start_Event
    (Event_Node_t * ep);
static void  Process_Rev_Stop_Event
    (Event_Node_t * ep);
static void  PWM_Score_Fwd_Start
    (int pos, const PWM_t & pwm, int window, double & score, int & separation);
static void  PWM_Score_Rev_Start
    (int pos, const PWM_t & pwm, int window, double & score, int & separation);
static void  Read_Entropy_Profiles
    (const char * fn, bool & errflg);
static void  Read_Sequences
    (FILE * fp, vector <string> & seq_list, vector <string> & hdr_list,
     int & seq_ct);
static void  Requalify
    (Event_Node_t * p, int cutoff);
static void  Reverse_Complement_Transfer
    (string & buff, const string & s, int lo, int hi);
static void  Reverse_Transfer
    (string & buff, const string & s, int start, int len);
static void  Score_Orflist
    (FILE * detail_fp, FILE * summary_fp);
static void  Score_Orfs
    (vector <Orf_t> & orf_list, vector <Gene_t> & gene_list, FILE * fp);
static void  Score_Separate_Input
    (const string & seq, const string & hdr, int seq_num, FILE * detail_fp,
     FILE * predict_fp);
static void  Set_Final_Event
    (Event_Node_t & fe, Event_Node_t * best_event [6],
     int seq_len);
static void  Set_GC_Fraction
    (double & gc, const vector <string> & s);
static void  Set_Ignore_Score_Len
    (void);
static void  Set_Start_And_Stop_Codons
    (void);
static void  Shift_Events
    (vector <Event_Node_t *> & ep, int reference_pos);
static void  Show_Events
    (FILE * fp);
static void  Trace_Back
    (FILE * fp, const Event_Node_t & final_event);
static void  Usage
    (void);
static void  Wrap_Around_Back
    (int wfr, int pos, int & gene_len, int & orf_len);
static void  Wrap_Through_Front
    (int fr, int pos, int & gene_len, int & orf_len);

#endif
