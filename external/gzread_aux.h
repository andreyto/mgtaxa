#ifndef GZ_H
#define GZ_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <zlib.h>
#include <string.h>

namespace std{

class INBUF {
public:
  INBUF(string file);
  bool get_line(string & str);
  bool get_char(char & ch);
  ~INBUF(){gzclose(fp);delete[] buf;}
private:
  char * buf;
  gzFile fp;
  static const int LEN=65536;
  int CUR;
  int pos;
  bool proceed;
};

INBUF::INBUF(string file){
  fp=gzopen(file.c_str(),"rb");
  buf=new char [LEN];
  CUR==-1;
  pos=-1;
  proceed=true;
}

bool INBUF::get_line(string & str){
  str.clear();
  if (!proceed)
    return false;


  //  char * temp=new char [2000]; // max len of line
  bool flag=true;
  int i;
  int tind=0;
  while(flag){
    if (CUR==-1||pos==-1){
      CUR=gzread(fp,buf,LEN);
      pos=0;
    }

    for(i=pos;i<CUR;i++){
      if (buf[i] != '\n'){
	//	temp[tind]=buf[i];
	//	tind++;
	str.push_back(buf[i]);
      }else{
	flag=false;
	//	temp[tind]='\0';
	break;
      }
    }
    pos=i+1;
    if (pos>=CUR){
      if (!gzeof(fp)){
	pos=-1;
      }else{
	proceed=false;
      }
    }
  }
  //  str=string(temp);
  //  delete[] temp;
  return true;
}

bool INBUF::get_char(char & ch){
  ch=' ';

  if (!proceed)
    return false;
  
  if (CUR==-1||pos==-1){
    CUR=gzread(fp,buf,LEN);
    pos=0;
  }
  
  ch=buf[pos];
  pos++;
  
  if (pos>=CUR){
    if (!gzeof(fp)){
      pos=-1;
    }else{
      proceed=false;
    }
  }
  return true;
}
}
#endif
