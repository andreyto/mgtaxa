#include <fstream>
#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <ext/hash_map>
#include <map>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <utility>

using namespace std;
using namespace __gnu_cxx;

char sym(int x)
{
    if (x==0)
    {
        return 'A';
    }
    else if (x==1)
    {
        return 'C';
    }
    else if (x==2)
    {
        return 'G';
    }
    else if (x==3)
    {
        return 'T';
    }
    else
    {
        return 'X';
    }
}

void compute_map(int * kmap, int size, int k)
{
    vector<int> mult(k);
    mult[0]=1;
    for (int i=1;i<k;i++)
        mult[i]=4*mult[i-1];

    vector<int> val(k);
    for(int i=0;i<k;i++)
    {
        val[i]=0;
    }

    int num=1;

    for(int i=0;i<size;i++)
    {
        if (kmap[i]==-1)
        {
            //cerr << "Column " << num << " ";
            int rem=i;
            for(int j=k-1;j>=0;j--)
            {
                val[j]=rem/mult[j];
                rem=rem%mult[j];
                //cerr << sym(val[j]);
            }
            //cerr << " ";
            // complement
            for(int j=0;j<k;j++)
            {
                val[j]=3-val[j];
                //cerr << sym(val[j]);
            }

            //cerr << endl;
            num++;

            int irc=0;
            for(int j=0;j<k;j++)
            {
                irc+=(val[j]*mult[k-1-j]);
                val[j]=0;
            }

            // check and set map
            if (irc>=size)
            {
                cerr << "Error mapping: i " << i << " irc " << irc << endl;
                exit(0);
            }
            kmap[i]=irc;
            if (i!=irc)
            {
                kmap[irc]=i;
            }
        }
    }
}

int fun (char x)
{
    if (x == 'A' || x == 'a')
    {
        return 0;
    }
    else if (x == 'C' || x == 'c')
    {
        return 1;
    }
    else if (x == 'G' || x == 'g')
    {
        return 2;
    }
    else if (x == 'T' || x == 't')
    {
        return 3;
    }
    else
    {
        return -1;
    }
}

int rfun (char x)
{
    if (x == 'A' || x == 'a')
    {
        return 3;
    }
    else if (x == 'C' || x == 'c')
    {
        return 2;
    }
    else if (x == 'G' || x == 'g')
    {
        return 1;
    }
    else if (x == 'T' || x == 't')
    {
        return 0;
    }
    else
    {
        return -1;
    }
}

void update(int * vec, vector<int> & win, int & total_sum, int size, int k)
{
    bool flag=true;
    int index=0;
    for (int i=0;i<k;i++)
    {
        if (win[i]<0)
        {
            flag=false;
            break;
        }
        else
        {
            index+=win[i];
        }
    }
    if (index>=size)
    {
        cerr << "Error: index exceeds size " << index << " " << size << endl;
        exit(0);
    }
    if (flag)
    {
        vec[index]++;
        total_sum++;
    }
}

void proc_count(string & seq, int * vec, int start, int end, string & name, ostream & outs, int * kmap, int size, int k)
{

    int total_sum=0;

    vector<int> win(k);
    // forward strand
    // initialization
    string::iterator sit=seq.begin();
    int j=0;
    while(j<k && sit!=seq.end())
    {
        win[j]=((int) pow (4.0,k-1-j))*fun(*sit);
        j++;
        sit++;
    }

    if (j==k)
    {
        update(vec,win,total_sum,size,k);
    }

    while(sit!=seq.end())
    {
        for(j=0;j<k-1;j++)
        {
            win[j]=4*win[j+1];
        }
        win[k-1]=fun(*sit);
        sit++;
        update(vec,win,total_sum,size,k);
    }

    // print output

    if (total_sum>0)
    {
        outs << name << " " << start << " " << end << " " << k << " " << total_sum;
        //    cout << name << " " << start << " " << end << " " << k << " " << total_sum << endl;
        for (int i=0;i<size;i++)
        {
            if (i==kmap[i])
            {
                outs << " " << 0.0001*(int(1000000.0*vec[i]/total_sum));
                vec[i]=0; // re-initializing
            }
            else if(i<kmap[i])
            {
                int sum=vec[i]+vec[kmap[i]];
                outs << " " << 0.0001*(int(1000000.0*sum/total_sum));
                vec[i]=0;
                vec[kmap[i]]=0;
            }
        }

        outs << endl;
    }

}

bool get_char(FILE *inp, char& ch)
{
///
	if( feof(inp) ) { return false; }
	ch = fgetc(inp);
	return ch != EOF;
///
}

void process_fasta (FILE *ins, ostream & outs, int * kmap, int size, int k, int chunk)
{

    int * vec;
    vec = new int [size];
    for(int i=0;i<size;i++)
    {
        vec[i]=0;
    }

    string str;
    char ch;
    string name;
    string seq;
    int start=0;
    int end=0;
    bool rflag=get_char(ins,ch);
    while(rflag)
    {
        if (ch == '>')
        {
            seq.clear();
            name.clear();
            while(ch != ' ' && ch != '\n' && rflag)
            {
                rflag=get_char(ins,ch);
                if (ch!=' ' && ch!='\n')
                    name += ch;
            }
            start=0;
            end=0;

            //      cout << name << endl;

            // remove rest of line
            while (ch != '\n' && rflag)
            {
                rflag=get_char(ins,ch);
            }

            // reading sequence
            bool sflag=true;
            if (rflag)
            {
                if (chunk<=0)
                {
                    while(sflag && rflag)
                    {
                        rflag=get_char(ins,ch);
                        if (ch=='>')
                        {
                            sflag=false;
                        }
                        else
                        {
                            if (ch!='\n' && ch!=' ')
                                seq += ch;
                        }
                    }
                    proc_count(seq,vec,1,seq.size(),name,outs,kmap,size,k);
                }
                else
                {
                    while(sflag && rflag)
                    {
                        start=end+1;
                        rflag=get_char(ins,ch);
                        if (ch=='>' || !rflag)
                        {
                            sflag=false;
                            if (seq.size()==chunk)
                            {
                                end=seq.size()+start-1;
                                proc_count(seq,vec,start,end,name,outs,kmap,size,k);
                            }
                        }
                        else
                        {
                            if (seq.size()==chunk)
                            {
                                end=seq.size()+start-1;
                                proc_count(seq,vec,start,end,name,outs,kmap,size,k);
                                seq.clear();
                            }

                            if (ch!='\n' && ch!=' ')
                                seq += ch;
                        }
                    }
                }
            }
        }
        else
        {
            rflag=get_char(ins,ch);
        }
    }

    delete[] vec;
}

int main (int argc, char* argv[])
{
    int k=6; // default k-mer size
    int chunk=0; // default chunk size, <=0 values mean that full sequence is considered
    int c;

    while((c = getopt(argc, argv,"k:l:"))!=-1)
    {
        switch(c)
        {
        case 'k':
            sscanf(optarg,"%d",&k);
            break;
        case 'l':
            sscanf(optarg,"%d",&chunk);
            break;
        default:
            break;
        }
    }


    if (k<=0 || k>7)
    {
        cerr << "kmer size is set to 6\n";
        k=6;
    }

    if (chunk<=0)
    {
        cerr << "Analysis for full sequence lengths\n";
    }

    time_t rawtime;
    struct tm * cur_time;

    time ( &rawtime );
    cur_time = localtime ( &rawtime );

    cerr << "Program started " << asctime(cur_time) << endl;

    cerr << "-k " << k << endl
    << "-l " << chunk << endl;


    int * kmap; // array storing map of kmer to its reverse complement
    // A:0,C:1,G:2,T:3
    // for sequences that contain non-bases, the kmers containing the non-bases
    // are not included in the count

    int size = (int) pow(4.0,k);

    kmap = new int [size];

    for(long int i=0;i<size;i++)
        kmap[i]=-1;

    compute_map(kmap,size,k);
    cerr << "Compute map finished size: " << size << endl;
	
	const int inBufferSize = 1024*1024;
	char inBuffer[inBufferSize];
	setbuffer(stdin,inBuffer,inBufferSize);
    //istream& inp = cin;
    ostream& outs = cout;
    //ifstream inp("test.fas");
    //ofstream outs("test.kmer");
    process_fasta(stdin,outs,kmap,size,k,chunk);

    delete[] kmap;

    time ( &rawtime );
    cur_time = localtime ( &rawtime );

    cerr << "Program finished " << asctime(cur_time) << endl;
}

