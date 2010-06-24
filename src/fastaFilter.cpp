//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the MGTAXA package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include <string>
#include <cstdlib>
#include <cstdio>
#include <assert.h>

using namespace std;

string degenC = "nN"; //degenarate symbols to compress

int maxDegen = 20; //any run of consequitive degenerate 
//symbols will be compressed to this number

int maxLine = 32768; // max output line length for residue codes

int nLine = 0; // output line length counter for residue codes

bool isDegen[256];

void setDegen(const string& degenC) {
    for ( int i = 0; i < 256; i++ ) isDegen[i] = false;
    for ( int i = 0; i < degenC.size(); i++ ) {
        int ind = degenC[i];
        assert(ind < 256 && ind >= 0);
        isDegen[ind] = true;
    }
}

/*
inline bool isDegen(int c) {
    
    //for(int i = 0; i < degenC.size(); i++) {
    //    if( c == degenC[i] ) return true;
    //}
    
    if ( c == 'N' || c == 'n' ) return true;
    return false;
}
*/

inline void putCode(int c) {
    putchar_unlocked(c);
    if ( ++nLine >= maxLine ) {
        putchar_unlocked('\n');
        nLine = 0;
    }
}

void usage() {
    const char *msg = "Clip consequtive runs of degenerate symbols in FASTA input.\n"
                      "Read FASTA from standard input and write FASTA to standard output.\n"
                      "Equal line length is restored after deleting the degenerates.\n"
                      "Arguments:\n"
                      "-d string of symbols to consider degenerates, case-sensitive [Nn]\n"
                      "-n maximum number of consequtive degenerates to keep [20]\n"
                      "-l output line length (excluding fasta deflines)[32768]\n"
                      "Example (with -n 5): 'actgNNNnnnnnnnactg' will become 'actgNNNnnactg'\n";
    fprintf(stderr,msg);
}

int main(int argc, char**argv) {

    int iArg = 1;
    for( ; iArg < argc; iArg++ ) {
        string arg = argv[iArg];
        if ( arg == "-d" ) {
            degenC = string(argv[++iArg]);
        }
        else if ( arg == "-n" ) {
            maxDegen = atoi(argv[++iArg]);
        }
        else if ( arg == "-l" ) {
            maxLine = atoi(argv[++iArg]);
        }
        else {
            usage();
            return 1;
        }
    }

    if ( iArg < argc ) {
        usage();
        return 1;
    }
    
    setDegen(degenC);

    nLine = 0;
    int nDegen = 0;
    int c = EOF;
    while( (c = getchar_unlocked()) != EOF ) {
        if ( c == '>' ) {
            if ( nLine > 0 ) {
                putchar_unlocked('\n');
                nLine = 0;
            }
            putchar_unlocked(c);
            do {
                c = getchar_unlocked();
                if ( c == EOF ) goto eof;
                putchar_unlocked(c);
                if ( c == '\n' ) break;
            } while (true);
        }
        else {
            if ( isDegen[c] ) {
                if ( ++nDegen <= maxDegen ) putCode(c);
            }
            else {
                if ( c != '\n' ) {
                    nDegen = 0;
                    putCode(c);
                }
            }
        }
    }
eof:
   
return 0;
}
