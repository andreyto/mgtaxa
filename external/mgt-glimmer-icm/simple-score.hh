//    Programmer:  Arthur L. Delcher
//          File:  simple-score.hh
//  Last Updated:  Thu Nov  3 08:03:26 EST 2005
//
//  Declarations for  simple-score.cc


#ifndef _SIMPLE_SCORE_HH
#define _SIMPLE_SCORE_HH


#include  "icm.hh"


static void  Parse_Command_Line
    (int argc, char * argv []);
static int  Read_String
    (FILE * fp, char * & s, long int & s_size, char * & tag,
     long int & tag_size);
static void  Usage
    (void);


#endif
