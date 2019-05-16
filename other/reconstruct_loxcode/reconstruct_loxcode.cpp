//compile g++ -O4 -I ../edlib/include -I ../../seqan/include  main_reconstruct_shortest_path.cpp -std=c++1y ../edlib/src/edlib.cpp
//or g++ -O4 -I edlib/include  main_reconstruct_loxcode.cpp -std=c++1y edlib/src/edlib.cpp

#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <assert.h>
#include "edlib.h"

using namespace std;

template<typename T>
void pop_front(std::vector<T>& vec)
{
    assert(!vec.empty());
    vec.erase(vec.begin());
}

void readFASTA(ifstream &F , std::vector<string> & lines)
{
  getline(F,lines[0]);getline(F,lines[1]);getline(F,lines[2]);getline(F,lines[3]);
}

//barcode elements;
std::map<string, string> ele_R1={{"GCTCGAATTTGCAC","start"},{"ACTCCGCA","1"},{"TCCAGAATTTGTAT","2"},{"ACATCCAC","3"},{"AAAGGAATTTCTCC","4"},{"ATTTCCTC","5"},{"GCCCGAATTTTTTC","6"},{"GCTACTGG","7"},{"ATGAGAATTTATGG","8"},{"AACTAGAA","9"},{"TGCAGAATTTCCTC","10"},{"CGACACTT","11"},{"AACGGAATTTTCAA","12"},{"CGTGTTTG","13"},{"CAAACACG","-13"},{"TTGAAAATTCCGTT","-12"},{"AAGTGTCG","-11"},{"GAGGAAATTCTGCA","-10"},{"TTCTAGTT","-9"},{"CCATAAATTCTCAT","-8"},{"CCAGTAGC","-7"},{"GAAAAAATTCGGGC","-6"},{"GAGGAAAT","-5"},{"GGAGAAATTCCTTT","-4"},{"GTGGATGT","-3"},{"ATACAAATTCTGGA","-2"},{"TGCGGAGT","-1"},{"ACACGAATTCATCC","end"}};
std::map<string, string> ele_R2={{"GGATGAATTCGTGT","end"},{"ACTCCGCA","-1"},{"TCCAGAATTTGTAT","-2"},{"ACATCCAC","-3"},{"AAAGGAATTTCTCC","-4"},{"ATTTCCTC","-5"},{"GCCCGAATTTTTTC","-6"},{"GCTACTGG","-7"},{"ATGAGAATTTATGG","-8"},{"AACTAGAA","-9"},{"TGCAGAATTTCCTC","-10"},{"CGACACTT","-11"},{"AACGGAATTTTCAA","-12"},{"CGTGTTTG","-13"},{"CAAACACG","13"},{"TTGAAAATTCCGTT","12"},{"AAGTGTCG","11"},{"GAGGAAATTCTGCA","10"},{"TTCTAGTT","9"},{"CCATAAATTCTCAT","8"},{"CCAGTAGC","7"},{"GAAAAAATTCGGGC","6"},{"GAGGAAAT","5"},{"GGAGAAATTCCTTT","4"},{"GTGGATGT","3"},{"ATACAAATTCTGGA","2"},{"TGCGGAGT","1"},{"GTGCAAATTCGAGC","start"}};
//////////////

std::vector<string> Consensus(vector<string> a, vector<string> b)
/**
 * Used if R1, R2 overlap
 */
{
    int offset=0; bool do_align=false;
    int max_overlap=0,max_offset=0;
    int asize=a.size(),bsize=b.size();
    
    for(offset=-bsize; offset<asize; offset++) //check all possible starting positions of b
    {
        int match=0,not_match=0;
        
        for(int k=0; k<b.size();k++)
        {
            if(k+offset<0) ;
            else if(k+offset>=asize) ;
            else
            {
                if(a[k+offset]==b[k] && b[k]!="?") match++;
                if(a[k+offset]!=b[k] && a[k+offset]!="?" && b[k]!="?") not_match++;
            }
        }
      if(match>max_overlap && not_match==0) {max_overlap=match; max_offset=offset;}
    }
    
    std::vector<string> consensus,empty;
    if(max_overlap==0)
    {
        //add special case of full length
        if(a.size()==8 && b.size()==6 && count(a.begin(),a.end(),"?")==0 && count(b.begin(),b.end(),"?")==0)
        {
            vector<string> full_length(a);
            full_length.push_back("?");
            full_length.insert(full_length.end(), b.begin(), b.end());
            return full_length;
        }
        else return empty;
    }
    
    
    for(int i=0; i<max_offset; i++) consensus.push_back(a[i]);
    for(int k=0; k<b.size();k++)
    {
      if(k+max_offset<a.size())
        {
          if(a[k+max_offset]==b[k]) {consensus.push_back(b[k]); continue;}
          if(a[k+max_offset]!="?" && b[k]=="?"){consensus.push_back(a[k+max_offset]); continue;}
          if(a[k+max_offset]=="?" && b[k]!="?"){consensus.push_back(b[k]); continue;}
        }
      else consensus.push_back(b[k]);
    }
    
    if(max_offset<0) { for(int i=max_offset;i<0;i++) pop_front(consensus); }
    
    
    if(count(consensus.begin(),consensus.end(),"start")==1 &&
       count(consensus.begin(),consensus.end(),"end")==1   &&
       (count(consensus.begin(),consensus.end(),"?")==0)) 
    {
       return consensus;
    }
    else return empty;
}

int main(int argc, char* argv[])
{
    string prefix=argv[1];
    string dir=argv[2];
   
    string nR1 = dir+prefix+"_R1_001.fastq";
    string nR2 = dir+prefix+"_R2_001.fastq";
    
    int counter=0;

    ofstream fs("saturation_"+prefix+".txt");    
    ifstream fileR1(nR1);ifstream fileR2(nR2);
  
    map< std::vector<string> , int > CODES;
  
    double keep=0; //keep track of reads that are useful
    
    if(!fileR1.is_open() || !fileR2.is_open()){ cout<<"file not found"<<endl; return 0;}
  
  for (int i=0; ;i++)
    {
      if(fileR1.eof() || fileR2.eof()) break;
      
      std::vector<string> lines(4);
      readFASTA(fileR1, lines);
      string R1=lines[1];
      readFASTA(fileR2, lines);
      string R2=lines[1];
        
     if(R1.length()<290 || R2.length()<244) continue; //this depends on sequencing run
     string start="GCTCGAATTTGCAC",end="GGATGAATTCGTGT";
     int start_loc=0,end_loc=0;
     
      ///using edit ditance to find start (in read 1) and end (in read 2)
      EdlibAlignResult result;
      result = edlibAlign(&start[0u],start.size(),&R1[0u],R1.size(), edlibNewAlignConfig(1, EDLIB_MODE_HW,  EDLIB_TASK_LOC, NULL, 0));
      if(result.numLocations!=1) continue;
      start_loc=result.startLocations[0];
      
      result = edlibAlign(&end[0u],end.size(),&R2[0u],R2.size(), edlibNewAlignConfig(1, EDLIB_MODE_HW,  EDLIB_TASK_LOC, NULL, 0));
      if(result.numLocations!=1) continue;
      end_loc=result.startLocations[0];
      
      edlibFreeAlignResult(result);
      
      std::size_t loc=-1;
      
      map<int,int> locs={{0,14},{48,8},{90,14},{138,8},{180,14},{228,8},{270,14},/*{318,8}*/};
      
      bool discard=false, discard_R1=false, discard_R2=false;
      vector<string> loxcode_R1, loxcode_R2;
      
      for(auto l : locs)
        {
            if(l.first+start_loc<0 || l.first+start_loc>R1.length()) {discard=true; break;}
          string el=R1.substr(l.first+start_loc,l.second);
          auto it = ele_R1.find(el); // exact
          if(it!=ele_R1.end()) {loxcode_R1.push_back(it->second);  continue;}
          
          else // can't find exact, relax and look +/- 2bp
          {
            if(l.first+start_loc-2>0 && l.first+start_loc-2<R1.length()) el=R1.substr(l.first+start_loc-2,l.second+4);
            else el=R1.substr(l.first+start_loc,l.second);
            
            map<int,string> pos;
            for(auto e : ele_R1) pos[(edlibAlign(&e.first[0u], e.first.size(),&el[0u], el.size(), edlibNewAlignConfig(1, EDLIB_MODE_HW,  EDLIB_TASK_DISTANCE, NULL, 0)).editDistance)]=e.second;
            
            if(pos.count(0)==1 && pos.find(0)->second!="5" && pos.find(0)->second!="-5") // 5 and -5 very close to some other sequence 
              {
                  loxcode_R1.push_back(pos.find(0)->second); //if(pos.find(0)->second=="end") break;
              }
            else {loxcode_R1.push_back("?"); }
          }
        }
      
      map<int,int> locs_R2={{0,14},{48,8},{90,14},{138,8},{180,14},{228,8}};
      
      for(auto l : locs_R2)
      {
        if(l.first+end_loc<0 || l.first+end_loc>R2.length()) {discard=true; break;}
        string el=R2.substr(l.first+end_loc,l.second);
        auto it = ele_R2.find(el);
        if(it!=ele_R2.end()) {loxcode_R2.push_back(it->second); continue;}
        
        else
        {
          if(l.first+end_loc-2>0 && l.first+end_loc-2<R2.length()) el=R2.substr(l.first+end_loc-2,l.second+4);
          else 
            el=R2.substr(l.first+end_loc,l.second);
          
          map<int,string> pos;
          for(auto e : ele_R2) pos[(edlibAlign( &e.first[0u], e.first.size(),&el[0u], el.size(), edlibNewAlignConfig(1, EDLIB_MODE_HW,  EDLIB_TASK_DISTANCE, NULL, 0)).editDistance)]=e.second;

          if(pos.count(0)==1 && pos.find(0)->second!="5" && pos.find(0)->second!="-5")
            {
                loxcode_R2.push_back(pos.find(0)->second);
            }
          else {loxcode_R2.push_back("?");}
          
        }
        
      }
        
      
        {
            std::vector<string>r1=loxcode_R1;
            std::vector<string>r2=loxcode_R2;
	    for(auto i : r1) cout << i << " "; cout << endl;
	    for(auto i : r2) cout << i << " "; cout << endl;
	    cout << "********" << endl;
	    std::reverse(r2.begin(),r2.end());
            std::vector<string>r3=Consensus(r1,r2);
            if(r3.size()>0) // consensus gave us empty - we know read wasn't sensible bc
            {
            keep++;
            if(CODES.find(r3)==CODES.end()) CODES[r3]=1;
            else CODES[r3]++;
            }
            fs<<i<<" "<<CODES.size()<<" "<<keep/i<<endl;
        }
        
     }
    
  ofstream f("loxcode_counts_"+prefix+".csv");
  for(auto c : CODES)
  {
      f<<c.second<<","; for(auto i : c.first) f<<i<<" "; f<<", "<<c.first.size()-2<<endl;
  }
  
    
    fileR1.close();fileR2.close();f.close(); fs.close();
 
    
    return 0;
}
      
        
    
       
   
   
