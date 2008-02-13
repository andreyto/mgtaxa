#ifndef HM_H
#define HM_H

namespace __gnu_cxx{

#define HSIZE 1000000 // initial hash size 
#define HBSIZE 1000000 // initial hash size 

  class eqstr
    {
    public:
      bool operator()(std::string const & s1, std::string const & s2) const
	{
	  return strcmp(s1.c_str(), s2.c_str()) == 0;
	}
    };
  
  class hashstr
    {
    public:
      size_t operator()(std::string const &str) const
	{
	  return hash<char const *>()(str.c_str());
	}
    };
  
  class TBOOL
    {
    public:
      TBOOL(){len=0;one=false;two=false;three=false;}
      void set(int ulen, bool uone,bool utwo,bool uthree){
	len=ulen;
	one=uone;
	two=utwo;
	three=uthree;
      }
      int get_len(){return len;}
      bool first(){return one;}
      bool second(){return two;}
      bool third(){return three;}
      TBOOL(const TBOOL & tmp){
	len=tmp.len;
	one=tmp.one;
	two=tmp.two;
	three=tmp.three;
      }
      ~TBOOL(){}
    private:
      int len;
      bool one;
      bool two;
      bool three;
    };

  typedef hash_map<std::string,int,hashstr, eqstr> SIHMAP;
  typedef hash_map<std::string,double,hashstr, eqstr> SDHMAP;
  typedef hash_map<std::string,std::string,hashstr, eqstr> SSHMAP;
  typedef hash_map<std::string,bool,hashstr, eqstr> SBHMAP;  
  typedef hash_map<int,bool> IBHMAP;
  typedef hash_map<int,int> IIHMAP;
  typedef hash_map<std::string,TBOOL,hashstr, eqstr> STHMAP;
  typedef hash_map<std::string,pair<int,int>,hashstr, eqstr> SPHMAP;
  typedef hash_map<std::string,int*,hashstr, eqstr> SIPHMAP;
  typedef hash_map<std::string,std::list<pair<int,std::string> >,hashstr, eqstr> SLHMAP;

  bool present(SBHMAP & shash,std::string & s){
    bool flag=false;
    SBHMAP::iterator it=shash.find(s);
    if (it!=shash.end()){
      flag=true;
    }
    return flag;
  }

  bool present(STHMAP & shash,std::string & s){
    bool flag=false;
    STHMAP::iterator it=shash.find(s);
    if (it!=shash.end()){
      flag=true;
    }
    return flag;
  }

  void get_indices(SIHMAP & sn,std::vector<std::list<int> > & G, std::vector<std::string> & ns,std::string & query,std::string  & subject,int & counter,int & qlen, int & slen, int & qind,int & sind){
    SIHMAP::iterator pq=sn.find(query);
    SIHMAP::iterator ps=sn.find(subject);
  
    if(pq!=sn.end()){
      qind=pq->second;
    }else{
      qind=counter;
      sn[query]=qind;
      ns.push_back(query);
      counter++;
      std::list<int> tl;
      G.push_back(tl);
    }
    if(ps!=sn.end()){
      sind=ps->second;
    }else{
      sind=counter;
      sn[subject]=sind;
      ns.push_back(subject);
      counter++;
      std::list<int> tl;
      G.push_back(tl);
    }
  }
}

#endif
