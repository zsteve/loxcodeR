/**
 * LoxcodeR
 */

#include <Rcpp.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <assert.h>
#include "edlib/edlib.h"

using namespace Rcpp;
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

//' Decode FASTQ
//'
//' Recover loxcodes from raw Illumina FASTQ output
//' @param r Paths of R1, 2 respectively
//' @param meta User-defined data-frame for sample metadata
//' @param min_read_length min read length for R1, R2 filter respectively
//' @return S4 loxcode_sample object with decoded results
//' @export
// [[Rcpp::export]]
Rcpp::S4 decode(std::vector<std::string> r, Rcpp::DataFrame meta,
                int min_r1_len = 354, int min_r2_len = 244){
  ifstream fileR1(r[0]); ifstream fileR2(r[1]); // input files
  int counter=0;
  /*
   * code_readout.first = cassette sequence
   * code_readout.second = line numbers in fastq files corresponding to each
   */
  map< std::vector<string>, std::vector<int> > code_readout;
  vector<int> saturation; // track saturation

  int keep=0; //keep track of reads that are useful

  if(!fileR1.is_open() || !fileR2.is_open()){
    // file not found
    Rcpp::stop("no such file");
  }

  for(int i=0; ;i++){
    // i is the read counter (starts from zero)
    if(fileR1.eof() || fileR2.eof()) break;

    std::vector<string> lines(4);
    readFASTA(fileR1, lines);
    string R1=lines[1];
    readFASTA(fileR2, lines);
    string R2=lines[1];

    if(R1.length()<min_r1_len || R2.length()<min_r2_len) continue; //this depends on sequencing run

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

    map<int,int> locs={{0,14},{48,8},{90,14},{138,8},{180,14},{228,8},{270,14},{318,8}};

    bool discard=false, discard_R1=false, discard_R2=false;
    vector<string> loxcode_R1, loxcode_R2;

    for(auto l : locs){
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

    for(auto l : locs_R2){
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
       std::vector<string>r2=loxcode_R2;  std::reverse(r2.begin(),r2.end());
       std::vector<string>r3=Consensus(r1,r2);
       if(r3.size()>0) // consensus gave us empty - we know read wasn't sensible bc
       {
         keep++;
         if(code_readout.find(r3)==code_readout.end()){
            code_readout[r3].reserve(10); // reserve for 10 cassettes, if we need more then reallocate
         }
         code_readout[r3].push_back(i);
       }
       saturation.push_back(keep);
    }

  }

  std::vector<int> output_code_sizes; output_code_sizes.reserve(code_readout.size());
  std::vector<string> output_code_readout; output_code_readout.reserve(code_readout.size());
  std::vector<int> output_code_counts; output_code_counts.reserve(code_readout.size());
  std::vector<std::vector<int> > output_code_readids; output_code_readids.reserve(code_readout.size());
  for(auto c : code_readout){
      output_code_readout.push_back("");
      for(int i = 1; i < c.first.size()-1; ++i){ // we suppress start and end
              output_code_readout.back() += c.first[i];
	      // remove trailing ' ' - very important when converting back to integer form 
              if(i < c.first.size()-2) output_code_readout.back() += " ";
      }
      output_code_counts.push_back(c.second.size());
      output_code_sizes.push_back(c.first.size() - 2);
      output_code_readids.push_back(c.second);
  }

  Rcpp::DataFrame output_df = Rcpp::DataFrame::create(Named("count") = wrap(output_code_counts),
                                                      Named("code") = wrap(output_code_readout),
                                                      Named("size") = wrap(output_code_sizes),
						      Named("stringsAsFactors") = false);
  Rcpp::S4 decode_output("decode_output");
  decode_output.slot("data") = output_df;
  decode_output.slot("read_ids") = output_code_readids;
  decode_output.slot("saturation") = saturation;
  Rcpp::S4 output("loxcode_sample");

  output.slot("decode") = decode_output;
  output.slot("meta") = meta;
  output.slot("files") = r;

  fileR1.close();fileR2.close();

  return output;
}
