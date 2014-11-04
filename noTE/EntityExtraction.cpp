#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cctype>
#include <ctime>
#include <algorithm>
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <cmath>
#include <cstdlib>

#include "Util.h"
#include "HeapBinary.h"

#define _multi
//#define _output
using namespace tr1;
using namespace std;

vector<int> sum;
vector<unsigned long long> dpt;
vector<int> listCount;
vector<int> optimalIndex;
vector<int> greedyIndex;
vector<int> randomIndex;
vector<int> endIndex;
vector<int> beginIndex;
vector<pair<unsigned int,unsigned int> >indexTransfer;
int lastEndIndex;
int thisEndIndex;
bool isFirstTime;
string algorithm;
string type;

#define lmax 200
#define all 200

InvertList *invertList_source=new InvertList();

int Q(int begin,int end)
{
  cout<<"enter the function Q"<<endl;
  int sum=0;
  for(int k=begin;k<=end;k++)
  {
    sum+=invertList_source->at(k).second->size();
  }
  return sum*log(end-begin+1);
}

void EditDistanceEntry(int argc, char** argv);
void EditSimilarityEntry(int argc, char** argv);
void JaccardEntry(int argc, char** argv);
void EditDistanceEntry(int argc, char** argv, int limit);//ED
void EditSimilarityEntry(int argc, char** argv, int limit);//EDS
void EditSimilarityVerify(int argc, char** argv, int limit);//EDS
void JaccardEntry(int argc, char** argv, int limit);//JS
void CosineTokenEntry(int argc, char** argv);
void CosineEntry(int argc, char** argv, int limit);//CS
void DiceTokenEntry(int argc, char** argv);
void DiceEntry(int argc, char** argv, int limit);//DS
void printHelp(char* argv)
{
  std::cerr << "Usage:" << std::endl;
  std::cerr << "  " << argv << " <Approximate Type> [<q-gram>] <Threshold>";
  std::cerr << " <Dictionary File> <Document File>" << endl;
  std::cerr << " [Scale]" << std::endl;
  std::cerr << "Parameters:" << std::endl;
  std::cerr << "  <Approximate Type> :" << std::endl;
  std::cerr << "      Approximate Function, select from: EditDistance/EditSimilarity/Jaccard/Cosine/Dice." << std::endl;
  std::cerr << "  [<q-gram>] :" << std::endl;
  std::cerr << "      value of q, only for EditDistance/EditSimilarity." << std::endl;
  std::cerr << "  <Threshold> :" << std::endl;
  std::cerr << "      Approximate Function threshold, integer for EditDistance, float in (0,1] for others." << std::endl;
  std::cerr << "  <Dictionary File> :" << std::endl;
  std::cerr << "      Dictionary file path and name." << endl;
  std::cerr << "  <Document File> :" << std::endl;
  std::cerr << "      Document file path and name." << endl;
  std::cerr << "  [Scale] :" << std::endl;
  std::cerr << "      Test Scale or not." << endl;
  std::cerr << "  <Algorithm> :" << std::endl;
  std::cerr << "      Alternative Algorithm, select from: Greedy." << endl;
  std::cerr << "  <Param> :" << std::endl;
  std::cerr << "      Param for Greedy." << endl;
}

int main(int argc, char **argv)
{
  if(argc < 5)
  {
    std::cerr << "Invalid program arguments!!" << std::endl;
    printHelp(argv[0]);
    std::cerr << "Program exits with -1" << std::endl;
    return -1;
  }
  string function(argv[1]);
  cerr << "# Approximate Type : " << function << endl;
  if(function == "EditDistance")
  {
    if(argc > 6 && string(argv[6]) == "Scale")
    {
      for(int i = 5; i >= 1; i--)
      {
        entityLengthMap.clear();
        keywords.clear();
        gramListMap.clear();
        canda = candb = candc = candd = 0;
        maxLen = -1;
        minLen = MAX;
        MAXKEY = new Array<int>(1);
        MINKEY = new Array<int>(1);
        EditDistanceEntry(argc, argv, 20000 * i);
      }
    }
    else if(argc > 6 && string(argv[6]) == "Best")
      EditDistanceEntry(argc, argv, 0x7FFFFFFF);
    else
      EditDistanceEntry(argc, argv);
  }
  else if(function == "Jaccard")
  {
    if(argc > 6 && string(argv[6]) == "Scale")
    {
      for(int i = 5; i >= 1; i--)
      {
        entityLengthMap.clear();
        keywords.clear();
        gramListMap.clear();
        canda = candb = candc = candd = 0;
        maxLen = -1;
        minLen = MAX;
        MAXKEY = new Array<int>(1);
        MINKEY = new Array<int>(1);
        JaccardEntry(argc, argv, 20000 * i);
      }
    }
    else if(argc > 6 && string(argv[6]) == "Best")
      JaccardEntry(argc, argv, 0x7FFFFFFF);
    else
      JaccardEntry(argc, argv);
  }
  else if(function == "EditSimilarity")
  {
    if(argc > 6 && string(argv[6]) == "Scale")
    {
      for(int i = 5; i >= 1; i--)
      {
        entityLengthMap.clear();
        keywords.clear();
        gramListMap.clear();
        canda = candb = candc = candd = 0;
        maxLen = -1;
        minLen = MAX;
        MAXKEY = new Array<int>(1);
        MINKEY = new Array<int>(1);
        EditSimilarityEntry(argc, argv, 20000 * i);
      }
    }
    else if(argc > 6 && string(argv[6]) == "Best")
    {
      if(string(argv[3]) == "0.8")
        EditSimilarityVerify(argc, argv, 0X7FFFFFFF);
      else
        EditSimilarityEntry(argc, argv, 0x7FFFFFFF);
    }
    else
      EditSimilarityEntry(argc, argv);
  }
  else if(function == "Cosine")
  {
    if(argc > 6 && string(argv[6]) == "Scale")
    {
      for(int i = 5; i >= 1; i--)
      {
        entityLengthMap.clear();
        keywords.clear();
        gramListMap.clear();
        canda = candb = candc = candd = 0;
        maxLen = -1;
        minLen = MAX;
        MAXKEY = new Array<int>(1);
        MINKEY = new Array<int>(1);
        CosineEntry(argc, argv, 20000 * i);
      }
    }
    else
      CosineEntry(argc, argv, 0x7FFFFFFF);
  }
  else if(function == "Dice")
  {
    if(argc > 6 && string(argv[6]) == "Scale")
    {
      for(int i = 5; i >= 1; i--)
      {
        entityLengthMap.clear();
        keywords.clear();
        gramListMap.clear();
        canda = candb = candc = candd = 0;
        maxLen = -1;
        minLen = MAX;
        MAXKEY = new Array<int>(1);
        MINKEY = new Array<int>(1);
        DiceEntry(argc, argv, 20000 * i);
      }
    }
    else
      DiceEntry(argc, argv, 0x7FFFFFFF);
  }
  else {
    cerr << "# Wrong Approximate Type !! Programme exit"<< endl;
    printHelp(argv[0]);
    exit(1);
  }
  return 0;
}

//**************************************************************************************************************//
//                           Greedy,Optimal,Random,Average                                                                                                               //
//**************************************************************************************************************//

void getListCount(InvertList* invertList)
{
    listCount.clear();
	for(int i = 0;i < invertList->size();i++)
	{
	    int c = (invertList->at(i).second->size() > 0) ? invertList->at(i).second->size() - 1 : 0;
	    listCount.push_back(c);
	}
}

int getTE()
{
    int Elength = 0;
    for(int i = 0;i < keywords.size();i++)
    {
        if(keywords[i].size() > Elength)
            Elength = keywords[i].size();
    }
    if(type == "ED")
        return Elength - q + 1 + tau;
    else if(type == "JS")
        return Elength / det;
    else if(type == "CS")
        return (Elength - q + 1) / (det * det);
    else if(type == "DS")
        return (Elength - q + 1) * (2 - det) / det;
    else if(type == "ES")
        return Elength / det - q + 1;
}

void initialSum(int n)
{
    sum.clear();
	sum.push_back(listCount[0]);
	for(int i = 1;i < n;i++)
	{
		int s = sum[i - 1] + listCount[i];
		sum.push_back(s);
	}
}

void initialDpt(int n)
{
    dpt.clear();
    for(int i = 0;i < n;i++)
        dpt.push_back(0);
}

void initialOptimalIndex(int n)
{
    optimalIndex.clear();
    for(int i = 0;i < n;i++)
        optimalIndex.push_back(0);
}

void initialGreedyIndex(int n)
{
    greedyIndex.clear();
    for(int i = 0;i < n;i++)
        greedyIndex.push_back(0);
}

void initialRandomIndex(int n)
{
    randomIndex.clear();
    for(int i = 0;i < n;i++)
        randomIndex.push_back(0);
}


int optimal(int n,int param)
{
    unsigned long long min;
    optimalIndex[n - 1] = 0;
    unsigned long long temp;
    int k,TE = getTE();
    cout << "TE:" << TE << endl;

    dpt[0] = 0;
    for(int i = 1;i < TE;i++)
    {
        unsigned long long s = 0;

        //s = sum[i] * (log10(indexTransfer[i].second+1) / log2 + TE) + indexTransfer[i].second + 1;
        dpt[i] = 0;
    }

    for(int i = TE;i < n;i++)
    {
        //min = 44444444444444444;
        min = sum[i] * (log10(indexTransfer[i].second+1) / log2 + TE) + indexTransfer[i].second+1 + dpt[0];
        for(int j = TE - 1;j < i;j++)
        {
            int s = 0;
            k = j - TE + 1;
//            if(k < 0)
//                k = 0;
//            if(k == 0)
//            {
//                s = (sum[i - 1] - sum[0]) * (log10(i - k) / log2 + 2 * maxLen * tau) + i + dpt[j];
//            }
//            else
            if(k>0)
            {
                s = (sum[i] - sum[k-1]) * (log10(indexTransfer[i].second-indexTransfer[k-1].first+1) / log2 + TE) + indexTransfer[i].second-indexTransfer[k-1].first+1 + dpt[j];
            }
            else
            {
                s = sum[i] * (log10(indexTransfer[i].second+1) / log2 + TE) + indexTransfer[i].second+1 + dpt[j];
            }
            if(s < min)
            {
                cout << "i: " << i << " j: " << j << " s: " << s << " min: " << min << endl;
                min = s;
                optimalIndex[i] = j;
            }

        }
        dpt[i] = min;
    }
    return dpt[n];
}

void getEndIndex()
{
    endIndex.clear();
    endIndex.push_back(optimalIndex.size() - 1);
    for(int i = optimalIndex[optimalIndex.size() - 1];i >= 0;)
    {
        if(i == 0)
            break;
        endIndex.push_back(i);
        i = optimalIndex[i];
    }
}

void getBeginIndex()
{
    beginIndex.clear();
    int te = getTE();
    for(int i = optimalIndex[optimalIndex.size() - 1];i >= 0;)
    {
        if(i - te <= 0)
            break;
        beginIndex.push_back(i - te);
        i = optimalIndex[i];
    }
    beginIndex.push_back(0);
}

void greedy(int n,double param)
{
    beginIndex.clear();
    endIndex.clear();

    int te = getTE();
    cout  << "te: " << te << endl;
    unsigned long long s = sum[n-1] * (log10(indexTransfer[n-1].second+1) / log2 + te) + indexTransfer[n-1].second+1;
    int begin = 0;
    for(int i = te + 1;i <= n;i++)
    {
        unsigned long long ts1,ts2;
        if(begin>0)
        {
            ts1 = (sum[i-1] - sum[begin-1]) * (log10(indexTransfer[i-1].second-indexTransfer[begin-1].first+1) / log2 + te) + indexTransfer[i-1].second-indexTransfer[begin-1].first+1;
            ts2 = (sum[n-1] - sum[i - te - 1]) * (log10(indexTransfer[n-1].second-indexTransfer[i-te-1].first+1) / log2 + te) + indexTransfer[n-1].second-indexTransfer[i-te-1].first+1;
        }
        else
        {
            ts1 = sum[i-1] * (log10(indexTransfer[i-1].second+1) / log2 + te) + indexTransfer[i-1].second+1;
            ts2 = (sum[n-1] - sum[i - te - 1]) * (log10(indexTransfer[n-1].second-indexTransfer[i-te-1].first+1) / log2 + te) + indexTransfer[n-1].second-indexTransfer[i-te-1].first+1;
        }
        if(ts1+ts2 < s)
        {
            cout << "ts1: " << ts1 << " ts2: " << ts2 << " s: " << s << " i: " << i << endl;

            if(i<n)
            {
                beginIndex.push_back(begin);
                endIndex.push_back(i-1);
            }
            begin = i - te;

            if(begin>0)
            {
                s = (sum[n-1] - sum[begin-1]) * (log10(indexTransfer[n-1].second-indexTransfer[begin-1].first+1) / log2 + te) + indexTransfer[n-1].second-indexTransfer[begin-1].first+1;
            }
            else
            {
                s = sum[n-1] * (log10(indexTransfer[n-1].second+1) / log2 + te) + indexTransfer[n-1].second+1;
            }

        }
    }
    if(begin<n)
    {
        beginIndex.push_back(begin);
        endIndex.push_back(n-1);
    }
}

void getGreedyBeginIndex()
{
//    beginIndex.clear();
//    beginIndex.push_back(0);
//    int te = getTE();
//    for(int i = greedyIndex[0];i < greedyIndex.size();)
//    {
//        beginIndex.push_back(i - te);
//        i = greedyIndex[i - te];
//        if(i == -1)
//            break;
//    }
}

void getGreedyEndIndex()
{
//    endIndex.clear();
//    int te = getTE();
//    for(int i = greedyIndex[0];i < greedyIndex.size();)
//    {
//        endIndex.push_back(i);
//        i = greedyIndex[i - te];
//        if(i == -1)
//            break;
//    }
//    endIndex.push_back(greedyIndex.size() - 1);
}

void random(int n,int param)
{
    srand(time(0));
    int te = getTE();
    int begin = 0;
    vector<int> ranindex;
    for(int i = 0;i < param;i++)
    {
        ranindex.push_back(rand() % n + te);
    }
    sort(ranindex.begin(),ranindex.end());

    for(int i = 0;i < ranindex.size();i++)
    {
        randomIndex[begin] = ranindex[i];
        begin = ranindex[i] - te;
    }
    randomIndex[begin] = -1;
}

void getRandomBeginIndex()
{
    beginIndex.clear();
    beginIndex.push_back(0);
    int te = getTE();
    for(int i = randomIndex[0];i < randomIndex.size();)
    {
        beginIndex.push_back(i - te);
        i = randomIndex[i - te];
        if(i == -1)
            break;
    }
}

void getRandomEndIndex()
{
    endIndex.clear();
    int te = getTE();
    for(int i = randomIndex[0];i < randomIndex.size();)
    {
        endIndex.push_back(i);
        i = randomIndex[i - te];
        if(i == -1)
            break;
    }
    endIndex.push_back(randomIndex.size() - 1);
}



void initialBeginIndex()
{
    beginIndex.clear();
    int te = getTE();
    int k = lmax - te;
    for(int i = optimalIndex.size() - 1 - lmax;i >= 0;)
    {
        if(i <= 0)
            break;
        beginIndex.push_back(i);
        i -= k;
    }
    beginIndex.push_back(0);
}

void initialEndIndex()
{
    endIndex.clear();
    int te = getTE();
    int k = lmax - te;
    for(int i = optimalIndex.size() - 1;i >= 0;)
    {
        if(i == 0)
            break;
        endIndex.push_back(i);
        i -= k;
    }
}


//**************************************************************************************************************//
//                              EDIT DISTANCE                                                                                                               //
//**************************************************************************************************************//

void EditDistanceEntry(int argc, char** argv)
{
  q = atoi(argv[2]);
  tau = atoi(argv[3]);
  entity = string(argv[4]);
  document = string(argv[5]);
  int limit = MAX;
  MAXKEY->append(MAX);

  vector<string> documents;
  int realResult = 0;
  timeval tb, te;
  double tt = 0;
  timeval sbegin,send;

  cerr <<"# Threshold: " << tau << endl;
  cerr <<"# Loading Entities..."<< endl;
  ifstream input(entity.c_str(), ios::in);
  string keyword;
  int lim = 0;
  while(getline(input,keyword)&& lim < limit)
  {
    lim++;
    for(unsigned int j=0; j<keyword.length(); j++)
      keyword[j] = tolower(keyword[j]);
    keywords.push_back(keyword);
  }
  input.close();
  cerr << "# Entities Number: "<<keywords.size()<<endl;

  cerr << "# Loading Documents:"<<endl;
  ifstream querys(document.c_str(), ios::in);
  string query;
  while(getline(querys,query))
  {
    if(query.length() - q < 0)  continue;
    for(unsigned int j=0; j<query.length(); j++)
      query[j] = tolower(query[j]);
    documents.push_back(query);
  }
  querys.close();
  cerr << "# Documents Number:"<<documents.size()<<endl;

  cerr << "# Create Gram List:" <<endl;
  for(unsigned int id = 0; id < keywords.size(); id++)
    appendGramList(keywords[id], id);
  addTail();
  cerr <<"# Minimum Length: "<<minLen+q-1 << endl;
  cerr <<"# Maximum Length: "<<maxLen+q-1 << endl;
  cerr <<"# Gram List Number: "<<gramListMap.size()<<endl;

  if(minLen <= tau * q) {
    cerr << "!!q is too short, choose a bigger q, Program will exit" << endl;
    exit(1);
  }

  HeapBinary hbi(EditDistance);
  InvertList *invertList = new InvertList();
  vector<int> *index = new vector<int>(invertList->size(), 0);
  if(string(argv[6]) == "Batch")
  {
    gettimeofday(&sbegin,0);
    int prevPos;
    for(unsigned int i = 0; i < documents.size(); i++)
    {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      prevPos = 0;
      for(int j = 0; j < (int)documents[i].length() - q + 1; j++) {
        GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
        if(it != gramListMap.end())
        {
          if (j - prevPos > tau * q + 1)
          {
            if((int)invertList->size() >= minLen - tau * q)
            {
              index->assign(invertList->size(), 0);
              hbi.EntityExtract(invertList,index,documents[i],i);
              invertList->clear();
            }
          }
          invertList->push_back(make_pair(j, it->second));
          prevPos = j;
        }
      }

      int n=documents[i].length()-q+1;
      int state[n];
      state[1]=invertList_source->at(1).second->size();
      int pos[n];
      pos[1]=1;
      int lower=minLen+q-1-tau;
      int upper=maxLen+q-1+tau;
      for(int k=2;k<=n;k++)
      {
        int postemp;
        int minimum=MAX;
        for(int m=k-1;m>=1;m--)
        {
          int m1=m-upper;
          int temp=state[m]+Q(m1,k);
          if(temp<minimum)
          {
            minimum=temp;
            postemp=m1;
          }
        }
        state[k]=minimum;
        pos[k]=postemp;
      }
      int begin_pos=pos[n];
      int end_pos=n;
      cerr<<"going enter while"<<endl;
      while(begin_pos>1)
      {
        cerr<<"enter while"<<endl;
        cerr<<end_pos<<"  "<<begin_pos<<endl;
        for(int h=begin_pos;h<=end_pos;h++)
        {
          invertList->at(h-begin_pos).first=invertList_source->at(h).first;
          for(int indextemp=1;indextemp<=invertList_source->at(h).second->size();indextemp++)
          {
            invertList->at(h-begin_pos).second->at(indextemp)=invertList_source->at(h).second->at(indextemp);
          }
        }
        if((int)invertList->size() >= minLen - tau * q)
        {
          index->assign(invertList->size(), 0);
          hbi.EntityExtract(invertList,index,documents[i],i);
          invertList->clear();
        }
        end_pos=begin_pos+upper;
        begin_pos=pos[end_pos];

        if(begin_pos<1)
          begin_pos=1;
      }
      cout<<end_pos<<"  "<<begin_pos<<endl;
      for(int h=begin_pos;h<=end_pos;h++)
      {
        invertList->at(h-begin_pos).first=invertList_source->at(h).first;
        for(int indextemp=1;indextemp<=invertList_source->at(h).second->size();indextemp++)
        {
          invertList->at(h-begin_pos).second->at(indextemp)=invertList_source->at(h).second->at(indextemp);
        }
      }
      if((int)invertList->size() >= minLen - tau * q)
      {
        index->assign(invertList->size(), 0);
        hbi.EntityExtract(invertList,index,documents[i],i);
        invertList->clear();
      }
    }

    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(edth_imp(keywords[result[i].ent_id].c_str(),  documents[result[i].doc_id].substr(result[i].pos, result[i].len).c_str(), entityLengthMap[result[i].ent_id]+q-1, result[i].len, tau))
      {
        //cout << endl << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
        //cout << "##" << keywords[result[i].ent_id] << "$$" << endl;
        //cout << "##" << documents[result[i].doc_id].substr(result[i].pos,result[i].len) << "$$" << endl;
        realResult++;
      }
    }
    cerr<<"bengkui"<<endl;
    gettimeofday(&te, 0);
    tt = 0;
    cerr<<".................................."<<endl;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With Binary====="<<endl;
    cerr << "#   Candidate Num : " << canda << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With Binary====="<<endl;

    result.clear();
    realResult = 0;
  }
  if(string(argv[6]) == "Binary")
  {
    gettimeofday(&sbegin,0);
    for(unsigned int i = 0; i < documents.size(); i++)
    {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      for(unsigned int j = 0; j < documents[i].length() - q + 1; j++) {
        GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
        if(it != gramListMap.end())
          invertList->push_back(make_pair(j, it->second));
      }
      if(invertList->size() == 0)
        continue;
      index->assign(invertList->size(), 0);
      hbi.EntityExtract(invertList,index,documents[i],i);
      invertList->clear();
    }
    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(edth_imp(keywords[result[i].ent_id].c_str(),  documents[result[i].doc_id].substr(result[i].pos, result[i].len).c_str(), entityLengthMap[result[i].ent_id]+q-1, result[i].len, tau))
      {
        //cout << endl << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
        //cout << "##" << keywords[result[i].ent_id] << "$$" << endl;
        //cout << "##" << documents[result[i].doc_id].substr(result[i].pos,result[i].len) << "$$" << endl;
        realResult++;
      }
    }
    gettimeofday(&te, 0);
    tt = 0;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With Binary====="<<endl;
    cerr << "#   Candidate Num : " << canda << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With Binary====="<<endl;

    result.clear();
    realResult = 0;
  }
  else if(string(argv[6]) == "Bucket")
  {
    gettimeofday(&sbegin,0);
    for(unsigned int i = 0; i < documents.size(); i++)
    {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      for(unsigned int j = 0; j < documents[i].length() - q + 1; j++) {
        GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
        if(it != gramListMap.end())
          invertList->push_back(make_pair(j, it->second));
      }
      if(invertList->size() == 0)
        continue;
      index->assign(invertList->size(), 0);
      hbi.EntityExtractWithBucketAndLazy(invertList,index,documents[i],i);
      invertList->clear();
    }
    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(edth_imp(keywords[result[i].ent_id].c_str(),  documents[result[i].doc_id].substr(result[i].pos, result[i].len).c_str(), entityLengthMap[result[i].ent_id]+q-1, result[i].len, tau))
      {
        //cout << endl << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
        //cout << "##" << keywords[result[i].ent_id] << "$$" << endl;
        //cout << "##" << documents[result[i].doc_id].substr(result[i].pos,result[i].len) << "$$" << endl;
        realResult++;
      }
    }
    gettimeofday(&te, 0);
    tt = 0;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With Bucket and Lazy====="<<endl;
    cerr << "#   Candidate Num : " << candb << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With Bucket and Lazy====="<<endl;

    result.clear();
    realResult = 0;
  }
  else if(string(argv[6]) == "Lazy")
  {
    gettimeofday(&sbegin,0);
    for(unsigned int i = 0; i < documents.size(); i++)
    {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      for(unsigned int j = 0; j < documents[i].length() - q + 1; j++) {
        GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
        if(it != gramListMap.end())
          invertList->push_back(make_pair(j, it->second));
      }
      if(invertList->size() == 0)
        continue;
      index->assign(invertList->size(), 0);
      hbi.EntityExtractWithLazyUpdate(invertList,index,documents[i],i);
      invertList->clear();
    }
    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(edth_imp(keywords[result[i].ent_id].c_str(),  documents[result[i].doc_id].substr(result[i].pos, result[i].len).c_str(), entityLengthMap[result[i].ent_id]+q-1, result[i].len, tau))
      {
        //cout << endl << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
        //cout << "##" << keywords[result[i].ent_id] << "$$" << endl;
        //cout << "##" << documents[result[i].doc_id].substr(result[i].pos,result[i].len) << "$$" << endl;
        realResult++;
      }
    }
    gettimeofday(&te, 0);
    tt = 0;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With Lazy Update====="<<endl;
    cerr << "#   Candidate Num : " << candc << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With Lazy Update====="<<endl;

    result.clear();
    realResult = 0;
  }
  else if(string(argv[6]) == "Single")
  {
    gettimeofday(&sbegin,0);
    for(unsigned int i = 0; i < documents.size(); i++)
    {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      for(unsigned int j = 0; j < documents[i].length() - q + 1; j++) {
        GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
        if(it != gramListMap.end())
          invertList->push_back(make_pair(j, it->second));
      }
      if(invertList->size() == 0)
        continue;
      index->assign(invertList->size(), 0);
      hbi.EntityExtractWithNoPruning(invertList,index,documents[i],i);
      invertList->clear();
    }
    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(edth_imp(keywords[result[i].ent_id].c_str(),  documents[result[i].doc_id].substr(result[i].pos, result[i].len).c_str(), entityLengthMap[result[i].ent_id]+q-1, result[i].len, tau))
      {
        //cout << endl << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
        //cout << "##" << keywords[result[i].ent_id] << "$$" << endl;
        //cout << "##" << documents[result[i].doc_id].substr(result[i].pos,result[i].len) << "$$" << endl;
        realResult++;
      }
    }
    gettimeofday(&te, 0);
    tt = 0;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With No Pruning====="<<endl;
    cerr << "#   Candidate Num : " << candd << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With No Pruning====="<<endl;

    result.clear();
    realResult = 0;
  }
  else if(string(argv[6]) == "Multi")
  {
#ifdef _multi
    gettimeofday(&sbegin,0);

    int lower = minLen + q - 1 - tau;
    int upper = maxLen + q - 1 + tau;

    for(unsigned int i = 0; i < documents.size(); i++) {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      for(int len = lower; len <= upper; len++) {
        for(unsigned int st = 0; st + len < documents[i].length() + 1; st++) {
          for(unsigned int j = st; j < st + len - q + 1; j++) {
            GramListMap::iterator it = gramListMap.find(documents[i].substr(j, q));
            if(it != gramListMap.end())
              invertList->push_back(make_pair(j - st, it->second));
          }
          if(invertList->size() == 0)
            continue;
          index->assign(invertList->size(), 0);
          hbi.EntityExtractWithMultHeap(invertList,index,documents[i], i, st, len);
          invertList->clear();
        }
      }
    }
    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(edth_imp(keywords[result[i].ent_id].c_str(),  documents[result[i].doc_id].substr(result[i].pos, result[i].len).c_str(), entityLengthMap[result[i].ent_id]+q-1, result[i].len, tau))
      {
        //cout << endl << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
        //cout << "##" << keywords[result[i].ent_id] << "$$" << endl;
        //cout << "##" << documents[result[i].doc_id].substr(result[i].pos,result[i].len) << "$$" << endl;
        realResult++;
      }
    }
    gettimeofday(&te, 0);
    tt = 0;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With Multi Heap====="<<endl;
    cerr << "#   Candidate Num : " << "" << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With Multi Heap====="<<endl;

    result.clear();
    realResult = 0;
#endif
  }
  delete index;
  delete invertList;
  GramListMap::iterator gbegin = gramListMap.begin();
  GramListMap::iterator gend = gramListMap.end();
  while(gbegin!=gend)
  {
    delete gbegin->second;
    gbegin++;
  }
  delete MINKEY;
  delete MAXKEY;
}

//**************************************************************************************************************//
//                              EDIT DISTANCE SIMILARITY                                                                                           //
//**************************************************************************************************************//

void EditSimilarityEntry(int argc, char** argv)
{
  q = atoi(argv[2]);
  det = atof(argv[3]);
  entity = string(argv[4]);
  document = string(argv[5]);
  int limit = MAX;
  MAXKEY->append(MAX);

  vector<string> documents;
  int realResult = 0;
  timeval tb, te;
  double tt = 0;
  timeval sbegin,send;

  cerr <<"# Threshold: " << det << endl;
  cerr<<"# Loading Entities:"<<endl;
  ifstream input(entity.c_str(), ios::in);
  string keyword;
  int lim = 0;
  while(getline(input,keyword)&& lim < limit)
  {
    lim++;
    for(unsigned int j=0; j<keyword.length(); j++)
      keyword[j] = tolower(keyword[j]);
    keywords.push_back(keyword);
  }
  input.close();
  cerr << "# Entities Number: "<<keywords.size()<<endl;

  cerr << "# Loading Documents:"<<endl;
  ifstream querys(document.c_str(), ios::in);
  string query;
  while(getline(querys,query))
  {
    if(query.length() - q < 0)  continue;
    for(unsigned int j=0; j<query.length(); j++)
      query[j] = tolower(query[j]);
    documents.push_back(query);
  }
  querys.close();
  cerr << "# Documents Number:"<<documents.size()<<endl;

  cerr << "# Create Gram List:" <<endl;
  for(unsigned int id = 0; id < keywords.size(); id++)
    appendGramList(keywords[id], id);
  addTail();
  cerr <<"# Minimum Length: "<<minLen+q-1 << endl;
  cerr <<"# Maximum Length: "<<maxLen+q-1 << endl;
  cerr <<"# Gram List Number: "<<gramListMap.size()<<endl;

  if(minLen <= (minLen + q - 1) * (1 - det) * q) {
    cerr << "!!q is too short, choose a bigger q, Program will exit" << endl;
    exit(1);
  }

  HeapBinary hbi(EditDistanceSimilarity);
  InvertList *invertList = new InvertList();
  vector<int> *index = new vector<int>(invertList->size(), 0);
  if(string(argv[6]) == "Binary")
  {
    gettimeofday(&sbegin,0);
    for(unsigned int i = 0; i < documents.size(); i++)
    {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      for(unsigned int j = 0; j < documents[i].length() - q + 1; j++) {
        GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
        if(it != gramListMap.end())
          invertList->push_back(make_pair(j, it->second));
      }
      if(invertList->size() == 0)
        continue;
      index->assign(invertList->size(), 0);
      hbi.edsEntityExtract(invertList,index,documents[i],i);
      invertList->clear();
    }
    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(eds_imp(keywords[result[i].ent_id].c_str(),  documents[result[i].doc_id].substr(result[i].pos, result[i].len).c_str(), entityLengthMap[result[i].ent_id]+q-1, result[i].len, det))
      {
        //cout << endl << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
        //cout << "##" << keywords[result[i].ent_id] << "$$" << endl;
        //cout << "##" << documents[result[i].doc_id].substr(result[i].pos,result[i].len) << "$$" << endl;
        realResult++;
      }
    }
    gettimeofday(&te, 0);
    tt = 0;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With Binary====="<<endl;
    cerr << "#   Candidate Num : " << canda << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With Binary====="<<endl;

    result.clear();
    realResult = 0;
  }
  else if(string(argv[6]) == "Bucket")
  {

    gettimeofday(&sbegin,0);
    for(unsigned int i = 0; i < documents.size(); i++)
    {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      for(unsigned int j = 0; j < documents[i].length() - q + 1; j++) {
        GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
        if(it != gramListMap.end())
          invertList->push_back(make_pair(j, it->second));
      }
      if(invertList->size() == 0)
        continue;
      index->assign(invertList->size(), 0);
      hbi.edsEntityExtractWithBucketAndLazy(invertList,index,documents[i],i);
      invertList->clear();
    }
    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(eds_imp(keywords[result[i].ent_id].c_str(),  documents[result[i].doc_id].substr(result[i].pos, result[i].len).c_str(), entityLengthMap[result[i].ent_id]+q-1, result[i].len, det))
      {
        //cout << endl << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
        //cout << "##" << keywords[result[i].ent_id] << "$$" << endl;
        //cout << "##" << documents[result[i].doc_id].substr(result[i].pos,result[i].len) << "$$" << endl;
        realResult++;
      }
    }
    gettimeofday(&te, 0);
    tt = 0;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With Bucket and Lazy====="<<endl;
    cerr << "#   Candidate Num : " << candb << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With Bucket and Lazy====="<<endl;

    result.clear();
    realResult = 0;
  }
  else if(string(argv[6]) == "Lazy")
  {

    gettimeofday(&sbegin,0);
    for(unsigned int i = 0; i < documents.size(); i++)
    {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      for(unsigned int j = 0; j < documents[i].length() - q + 1; j++) {
        GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
        if(it != gramListMap.end())
          invertList->push_back(make_pair(j, it->second));
      }
      if(invertList->size() == 0)
        continue;
      index->assign(invertList->size(), 0);
      hbi.edsEntityExtractWithLazyUpdate(invertList,index,documents[i],i);
      invertList->clear();
    }
    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(eds_imp(keywords[result[i].ent_id].c_str(),  documents[result[i].doc_id].substr(result[i].pos, result[i].len).c_str(), entityLengthMap[result[i].ent_id]+q-1, result[i].len, det))
      {
        //cout << endl << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
        //cout << "##" << keywords[result[i].ent_id] << "$$" << endl;
        //cout << "##" << documents[result[i].doc_id].substr(result[i].pos,result[i].len) << "$$" << endl;
        realResult++;
      }
    }
    gettimeofday(&te, 0);
    tt = 0;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With Lazy Update====="<<endl;
    cerr << "#   Candidate Num : " << candc << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With Lazy Update====="<<endl;

    result.clear();
    realResult = 0;
  }
  else if(string(argv[6]) == "Single")
  {
    gettimeofday(&sbegin,0);
    for(unsigned int i = 0; i < documents.size(); i++)
    {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      for(unsigned int j = 0; j < documents[i].length() - q + 1; j++) {
        GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
        if(it != gramListMap.end())
          invertList->push_back(make_pair(j, it->second));
      }
      if(invertList->size() == 0)
        continue;
      index->assign(invertList->size(), 0);
      hbi.edsEntityExtractWithNoPruning(invertList,index,documents[i],i);
      invertList->clear();
    }
    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(eds_imp(keywords[result[i].ent_id].c_str(),  documents[result[i].doc_id].substr(result[i].pos, result[i].len).c_str(), entityLengthMap[result[i].ent_id]+q-1, result[i].len, det))
      {
        //cout << endl << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
        //cout << "##" << keywords[result[i].ent_id] << "$$" << endl;
        //cout << "##" << documents[result[i].doc_id].substr(result[i].pos,result[i].len) << "$$" << endl;
        realResult++;
      }
    }
    gettimeofday(&te, 0);
    tt = 0;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With No Pruning====="<<endl;
    cerr << "#   Candidate Num : " << candd << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With No Pruning====="<<endl;

    result.clear();
    realResult = 0;
  }
  else if(string(argv[6]) == "Multi")
  {
#ifdef _multi
    if(det > 0.87) {
      gettimeofday(&sbegin,0);

      int lower = int(ceil((minLen + q - 1 ) * det));
      int upper = int(floor((maxLen + q - 1) / det));

      for(unsigned int i = 0; i < documents.size(); i++) {
#ifdef _output
        cerr << i << "/" << documents.size() << "\r\r";
#endif
        for(int len = lower; len <= upper; len++) {
          for(int st = 0; st + len < (int)documents[i].length() + 1; st++) {
            for(int j = st; j < st + len - q + 1; j++) {
              GramListMap::iterator it = gramListMap.find(documents[i].substr(j, q));
              if(it != gramListMap.end())
                invertList->push_back(make_pair(j - st, it->second));
            }
            if(invertList->size() == 0)
              continue;
            index->assign(invertList->size(), 0);
            hbi.edsEntityExtractWithMultHeap(invertList,index,documents[i], i, st, len);
            invertList->clear();
          }
        }
      }
      gettimeofday(&send,0);

      gettimeofday(&tb,0);
      for(unsigned int i = 0; i < result.size(); i++)
      {
        if(eds_imp(keywords[result[i].ent_id].c_str(),  documents[result[i].doc_id].substr(result[i].pos, result[i].len).c_str(), entityLengthMap[result[i].ent_id]+q-1, result[i].len, det))
        {
          //cout << endl << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
          //cout << "##" << keywords[result[i].ent_id] << "$$" << endl;
          //cout << "##" << documents[result[i].doc_id].substr(result[i].pos,result[i].len) << "$$" << endl;
          realResult++;
        }
      }
      gettimeofday(&te, 0);
      tt = 0;
      tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
      cerr << "=====Search With Multi Heap====="<<endl;
      cerr << "#   Candidate Num : " << "" << endl;
      cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
      cerr << "#     Verify Time : " << tt << "s" <<endl;
      cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
      cerr << "#  To Be Verified : " << result.size() << endl;
      cerr << "# Real Result Num : " << realResult << endl;
      cerr << "=====Search With Multi Heap====="<<endl;

      result.clear();
      realResult = 0;
    }else{
      cerr << "=====Search With Multi Heap====="<<endl;
      cerr << "# This will take too long time. " << endl;
      cerr << "# It cannot finish in 72000s. " << endl;
      cerr << "# So just let's assume it is 72000s" <<endl;
      cerr << "#      Total Time : " << 72000 << "s"<< endl;
      cerr << "#  To Be Verified : " << endl;
      cerr << "# Real Result Num : " << endl;
      cerr << "=====Search With Multi Heap====="<<endl;
    }
#endif
  }
  delete index;
  delete invertList;
  GramListMap::iterator gbegin = gramListMap.begin();
  GramListMap::iterator gend = gramListMap.end();
  while(gbegin!=gend)
  {
    delete gbegin->second;
    gbegin++;
  }
  delete MINKEY;
  delete MAXKEY;
}

//**************************************************************************************************************//
//                              JACCARD  SIMILARITY                                                                                                    //
//**************************************************************************************************************//

void JaccardEntry(int argc, char** argv)
{
  //q = atoi(argv[2]);
  det = atof(argv[2]);
  entity = string(argv[3]);
  document = string(argv[4]);
  int limit = MAX;
  MAXKEY->append(MAX);

  vector<string> documents;
  int realResult = 0;
  timeval tb, te;
  double tt = 0;
  timeval sbegin,send;

  cerr <<"# Threshold: " << det << endl;
  cerr<<"# Loading Entities:"<<endl;
  ifstream input(entity.c_str(), ios::in);
  string keyword;
  int lim = 0;
  while(getline(input,keyword)&& lim < limit)
  {
    lim++;
    for(unsigned int j=0; j<keyword.length(); j++)
      keyword[j] = tolower(keyword[j]);
    keywords.push_back(keyword);
  }
  input.close();
  cerr << "# Entities Number: "<<keywords.size()<<endl;

  cerr << "# Loading Documents:"<<endl;
  ifstream querys(document.c_str(), ios::in);
  string query;
  while(getline(querys,query))
  {
    for(unsigned int j=0; j<query.length(); j++)
      query[j] = tolower(query[j]);
    documents.push_back(query);
  }
  querys.close();
  cerr << "# Documents Number:"<<documents.size()<<endl;

  cerr << "# Create Gram List:" <<endl;
  for(unsigned int id = 0; id < keywords.size(); id++)
    appendTokenList(keywords[id], id);
  addTail();
  cerr <<"# Minimum Length: "<< minLen << endl;
  cerr <<"# Maximum Length: "<< maxLen << endl;
  cerr <<"# Gram List Number: "<< gramListMap.size() <<endl;

  vector<vector<string>*> tokendoc;
  for(unsigned int i = 0; i < documents.size(); i++) {
    vector<string> *tokens = new vector<string>();
    strToTokens(documents[i], tokens, delim);
    tokendoc.push_back(tokens);
  }

  HeapBinary hbi(Jaccard);
  InvertList *invertList = new InvertList();
  vector<int> *index = new vector<int>(invertList->size(), 0);
  gettimeofday(&sbegin,0);
  if(string(argv[5]) == "Binary")
  {
    for(unsigned int i = 0; i < documents.size(); i++)
    {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      for(unsigned int j=0; j<tokendoc[i]->size(); j++) {
        GramListMap::iterator it = gramListMap.find(tokendoc[i]->at(j));
        if(it!=gramListMap.end())
          invertList->push_back(make_pair(j,it->second));
      }
      if(invertList->size() == 0)
        continue;
      index->assign(invertList->size(), 0);
      hbi.jacEntityExtract(invertList,index,tokendoc[i],i);
      invertList->clear();
    }
    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(jaccardToken(keywords[result[i].ent_id], tokendoc[result[i].doc_id], result[i].pos, result[i].len) >= det)
      {
        //cerr << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
        realResult++;
      }
    }
    gettimeofday(&te, 0);
    tt = 0;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With Binary====="<<endl;
    cerr << "#   Candidate Num : " << canda << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With Binary====="<<endl;

    result.clear();
    realResult = 0;
  }
  else if(string(argv[5]) == "Bucket")
  {

    gettimeofday(&sbegin,0);
    for(unsigned int i = 0; i < documents.size(); i++)
    {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      for(unsigned int j=0; j<tokendoc[i]->size(); j++) {
        GramListMap::iterator it = gramListMap.find(tokendoc[i]->at(j));
        if(it!=gramListMap.end())
          invertList->push_back(make_pair(j,it->second));
      }
      if(invertList->size() == 0)
        continue;
      index->assign(invertList->size(), 0);
      hbi.jacEntityExtractWithBucketAndLazy(invertList,index,tokendoc[i],i);
      invertList->clear();
    }
    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(jaccardToken(keywords[result[i].ent_id], tokendoc[result[i].doc_id], result[i].pos, result[i].len) >= det)
        realResult++;
    }
    gettimeofday(&te, 0);
    tt = 0;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With Bucket and Lazy====="<<endl;
    cerr << "#   Candidate Num : " << candb << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With Bucket and Lazy====="<<endl;

    result.clear();
    realResult = 0;
  }
  else if(string(argv[5]) == "Lazy")
  {

    gettimeofday(&sbegin,0);
    for(unsigned int i = 0; i < documents.size(); i++)
    {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      for(unsigned int j=0; j<tokendoc[i]->size(); j++) {
        GramListMap::iterator it = gramListMap.find(tokendoc[i]->at(j));
        if(it!=gramListMap.end())
          invertList->push_back(make_pair(j,it->second));
      }
      if(invertList->size() == 0)
        continue;
      index->assign(invertList->size(), 0);
      hbi.jacEntityExtractWithLazyUpdate(invertList,index,tokendoc[i],i);
      invertList->clear();
    }
    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(jaccardToken(keywords[result[i].ent_id], tokendoc[result[i].doc_id], result[i].pos, result[i].len) >= det)
        realResult++;
    }
    gettimeofday(&te, 0);
    tt = 0;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With Lazy Update====="<<endl;
    cerr << "#   Candidate Num : " << candc << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With Lazy Update====="<<endl;

    result.clear();
    realResult = 0;
  }
  else if(string(argv[5]) == "Single")
  {

    gettimeofday(&sbegin,0);
    for(unsigned int i = 0; i < documents.size(); i++)
    {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      for(unsigned int j=0; j<tokendoc[i]->size(); j++) {
        GramListMap::iterator it = gramListMap.find(tokendoc[i]->at(j));
        if(it!=gramListMap.end())
          invertList->push_back(make_pair(j,it->second));
      }
      if(invertList->size() == 0)
        continue;
      index->assign(invertList->size(), 0);
      hbi.jacEntityExtractWithNoPruning(invertList,index,tokendoc[i],i);
      invertList->clear();
    }
    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(jaccardToken(keywords[result[i].ent_id], tokendoc[result[i].doc_id], result[i].pos, result[i].len) >= det)
        realResult++;
    }
    gettimeofday(&te, 0);
    tt = 0;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With No Pruning====="<<endl;
    cerr << "#   Candidate Num : " << candd << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With No Pruning====="<<endl;

    result.clear();
    realResult = 0;
  }
  else if(string(argv[5]) == "Multi")
  {
#ifdef _multi
    gettimeofday(&sbegin,0);

    int lower = int(ceil(minLen * det));
    int upper = int(floor(maxLen / det));

    for(unsigned int i = 0; i < documents.size(); i++) {
#ifdef _output
      cerr << i << "/" << documents.size() << "\r\r";
#endif
      for(int len = lower; len <= upper; len++) {
        for(int st = 0; st + len < (int)tokendoc[i]->size() + 1; st++) {
          for(int j = st; j < st + len; j++) {
            GramListMap::iterator it = gramListMap.find(tokendoc[i]->at(j));
            if(it != gramListMap.end())
              invertList->push_back(make_pair(j - st, it->second));
          }
          if(invertList->size() == 0)
            continue;
          index->assign(invertList->size(), 0);
          hbi.jacEntityExtractWithMultHeap(invertList,index,tokendoc[i], i, st, len);
          invertList->clear();
        }
      }
    }
    gettimeofday(&send,0);

    gettimeofday(&tb,0);
    for(unsigned int i = 0; i < result.size(); i++)
    {
      if(jaccardToken(keywords[result[i].ent_id], tokendoc[result[i].doc_id], result[i].pos, result[i].len) >= det)
        realResult++;
    }
    gettimeofday(&te, 0);
    tt = 0;
    tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
    cerr << "=====Search With Multi Heap====="<<endl;
    cerr << "#   Candidate Num : " << "" << endl;
    cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
    cerr << "#     Verify Time : " << tt << "s" <<endl;
    cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
    cerr << "#  To Be Verified : " << result.size() << endl;
    cerr << "# Real Result Num : " << realResult << endl;
    cerr << "=====Search With Multi Heap====="<<endl;

    result.clear();
    realResult = 0;
#endif
  }
  vector< vector<string>* >::iterator ibegin = tokendoc.begin();
  vector< vector<string>* >::iterator iend = tokendoc.end();
  while(ibegin!=iend)
  {
    delete *ibegin;
    ibegin++;
  }
  delete index;
  delete invertList;
  GramListMap::iterator gbegin = gramListMap.begin();
  GramListMap::iterator gend = gramListMap.end();
  while(gbegin!=gend)
  {
    delete gbegin->second;
    gbegin++;
  }
  delete MINKEY;
  delete MAXKEY;
}

//**************************************************************************************************************//
//                              COSINE  SIMILARITY TOKEN                                                                                          //
//**************************************************************************************************************//

void CosineTokenEntry(int argc, char** argv)
{
  //q = atoi(argv[2]);
  det = atof(argv[2]);
  entity = string(argv[3]);
  document = string(argv[4]);
  int limit = MAX;
  if(argc>5)
    limit = atoi(argv[5]);
  MAXKEY->append(MAX);

  vector<string> documents;
  int realResult = 0;
  timeval tb, te;
  double tt = 0;
  timeval sbegin,send;

  cerr<<"# Loading Entities:"<<endl;
  ifstream input(entity.c_str(), ios::in);
  string keyword;
  int lim = 0;
  while(getline(input,keyword)&& lim < limit)
  {
    lim++;
    for(unsigned int j=0; j<keyword.length(); j++)
      keyword[j] = tolower(keyword[j]);
    keywords.push_back(keyword);
  }
  input.close();
  cerr << "# Entities Number: "<<keywords.size()<<endl;

  cerr << "# Loading Documents:"<<endl;
  ifstream querys(document.c_str(), ios::in);
  string query;
  while(getline(querys,query))
  {
    for(unsigned int j=0; j<query.length(); j++)
      query[j] = tolower(query[j]);
    documents.push_back(query);
  }
  querys.close();
  cerr << "# Documents Number:"<<documents.size()<<endl;

  cerr << "# Create Gram List:" <<endl;
  for(unsigned int id = 0; id < keywords.size(); id++)
    appendTokenList(keywords[id], id);
  addTail();
  cerr <<"# Minimum Length: "<< minLen << endl;
  cerr <<"# Maximum Length: "<< maxLen << endl;
  cerr <<"# Gram List Number: "<< gramListMap.size() <<endl;

  vector<vector<string>*> tokendoc;
  for(unsigned int i = 0; i < documents.size(); i++) {
    vector<string> *tokens = new vector<string>();
    strToTokens(documents[i], tokens, delim);
    tokendoc.push_back(tokens);
  }

  gettimeofday(&sbegin,0);
  HeapBinary hbi(Cosine);
  InvertList *invertList = new InvertList();
  vector<int> *index = new vector<int>(invertList->size(), 0);
  for(unsigned int i = 0; i < documents.size(); i++)
  {
#ifdef _output
    cerr << i << "/" << documents.size() << "\r\r";
#endif
    for(unsigned int j=0; j<tokendoc[i]->size(); j++) {
      GramListMap::iterator it = gramListMap.find(tokendoc[i]->at(j));
      if(it!=gramListMap.end())
        invertList->push_back(make_pair(j,it->second));
    }
    if(invertList->size() == 0)
      continue;
    index->assign(invertList->size(), 0);
    hbi.cosEntityExtract(invertList,index,tokendoc[i],i);
    invertList->clear();
  }
  gettimeofday(&send,0);

  gettimeofday(&tb,0);
  for(unsigned int i = 0; i < result.size(); i++)
  {
    if(cosineToken(keywords[result[i].ent_id], tokendoc[result[i].doc_id], result[i].pos, result[i].len) >= det)
      realResult++;
  }
  gettimeofday(&te, 0);
  tt = 0;
  tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
  cerr << "=====Search With Binary====="<<endl;
  cerr << "#   Candidate Num : " << "" << endl;
  cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
  cerr << "#     Verify Time : " << tt << "s" <<endl;
  cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
  cerr << "#  To Be Verified : " << result.size() << endl;
  cerr << "# Real Result Num : " << realResult << endl;
  cerr << "=====Search With Binary====="<<endl;

  result.clear();
  realResult = 0;

  vector< vector<string>* >::iterator ibegin = tokendoc.begin();
  vector< vector<string>* >::iterator iend = tokendoc.end();
  while(ibegin!=iend)
  {
    delete *ibegin;
    ibegin++;
  }
  delete index;
  delete invertList;
  GramListMap::iterator gbegin = gramListMap.begin();
  GramListMap::iterator gend = gramListMap.end();
  while(gbegin!=gend)
  {
    delete gbegin->second;
    gbegin++;
  }
  delete MINKEY;
  delete MAXKEY;
}

//**************************************************************************************************************//
//                              DICE  SIMILARITY TOKEN                                                                                               //
//**************************************************************************************************************//

void DiceTokenEntry(int argc, char** argv)
{
  //q = atoi(argv[2]);
  det = atof(argv[2]);
  entity = string(argv[3]);
  document = string(argv[4]);
  int limit = MAX;
  if(argc>5)
    limit = atoi(argv[5]);
  MAXKEY->append(MAX);

  vector<string> documents;
  int realResult = 0;
  timeval tb, te;
  double tt = 0;
  timeval sbegin,send;

  cerr<<"# Loading Entities:"<<endl;
  ifstream input(entity.c_str(), ios::in);
  string keyword;
  while(getline(input,keyword))
  {
    for(unsigned int j=0; j<keyword.length(); j++)
      keyword[j] = tolower(keyword[j]);
    keywords.push_back(keyword);
  }
  input.close();
  cerr << "# Entities Number: "<<keywords.size()<<endl;

  cerr << "# Loading Documents:"<<endl;
  ifstream querys(document.c_str(), ios::in);
  string query;
  int lim = 0;
  while(getline(querys,query)&& lim < limit)
  {
    lim++;
    for(unsigned int j=0; j<query.length(); j++)
      query[j] = tolower(query[j]);
    documents.push_back(query);
  }
  querys.close();
  cerr << "# Documents Number:"<<documents.size()<<endl;

  cerr << "# Create Gram List:" <<endl;
  for(unsigned int id = 0; id < keywords.size(); id++)
    appendTokenList(keywords[id], id);
  addTail();
  cerr <<"# Minimum Length: "<< minLen << endl;
  cerr <<"# Maximum Length: "<< maxLen << endl;
  cerr <<"# Gram List Number: "<< gramListMap.size() <<endl;

  vector<vector<string>*> tokendoc;
  for(unsigned int i = 0; i < documents.size(); i++) {
    vector<string> *tokens = new vector<string>();
    strToTokens(documents[i], tokens, delim);
    tokendoc.push_back(tokens);
  }

  gettimeofday(&sbegin,0);
  HeapBinary hbi(Dice);
  InvertList *invertList = new InvertList();
  vector<int> *index = new vector<int>(invertList->size(), 0);
  for(unsigned int i = 0; i < documents.size(); i++)
  {
#ifdef _output
    cerr << i << "/" << documents.size() << "\r\r";
#endif
    for(unsigned int j=0; j<tokendoc[i]->size(); j++) {
      GramListMap::iterator it = gramListMap.find(tokendoc[i]->at(j));
      if(it!=gramListMap.end())
        invertList->push_back(make_pair(j,it->second));
    }
    if(invertList->size() == 0)
      continue;
    index->assign(invertList->size(), 0);
    hbi.diceEntityExtract(invertList,index,tokendoc[i],i);
    invertList->clear();
  }
  gettimeofday(&send,0);

  gettimeofday(&tb,0);
  for(unsigned int i = 0; i < result.size(); i++)
  {
    if(diceToken(keywords[result[i].ent_id], tokendoc[result[i].doc_id], result[i].pos, result[i].len) >= det)
      realResult++;
  }
  gettimeofday(&te, 0);
  tt = 0;
  tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
  cerr << "=====Search With Binary====="<<endl;
  cerr << "#   Candidate Num : " << "" << endl;
  cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC << "s" << endl;
  cerr << "#     Verify Time : " << tt << "s" <<endl;
  cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
  cerr << "#  To Be Verified : " << result.size() << endl;
  cerr << "# Real Result Num : " << realResult << endl;
  cerr << "=====Search With Binary====="<<endl;

  result.clear();
  realResult = 0;

  vector< vector<string>* >::iterator ibegin = tokendoc.begin();
  vector< vector<string>* >::iterator iend = tokendoc.end();
  while(ibegin!=iend)
  {
    delete *ibegin;
    ibegin++;
  }
  delete index;
  delete invertList;
  GramListMap::iterator gbegin = gramListMap.begin();
  GramListMap::iterator gend = gramListMap.end();
  while(gbegin!=gend)
  {
    delete gbegin->second;
    gbegin++;
  }
  delete MINKEY;
  delete MAXKEY;
}

//**************************************************************************************************************//
//                              EDIT DISTANCE SCALE                                                                                                    //
//**************************************************************************************************************//

void EditDistanceEntry(int argc, char** argv, int limit)
{
  q = atoi(argv[2]);
  tau = atoi(argv[3]);
  entity = string(argv[4]);
  document = string(argv[5]);
  algorithm = string(argv[7]);
  MAXKEY->append(MAX);

  type = "ED";

  vector<string> documents;
  int realResult = 0;
  timeval tb, te;
  double tt = 0;
  double pp = 0;
  timeval sbegin,send;
  timeval pb, pe;

  cerr <<"# Threshold: " << tau << endl;
  cerr<<"# Loading Entities:"<<endl;
  ifstream input(entity.c_str(), ios::in);
  string keyword;
  int lim = 0;
  while(getline(input,keyword)&& lim < limit)
  {
    lim++;
    for(unsigned int j=0; j<keyword.length(); j++)
      keyword[j] = tolower(keyword[j]);
    keywords.push_back(keyword);
  }
  input.close();
  cerr << "# Entities Number: "<<keywords.size()<<endl;

  cerr << "# Loading Documents:"<<endl;
  ifstream querys(document.c_str(), ios::in);
  string query;
  while(getline(querys,query))
  {
    if(query.length() - q < 0)  continue;
    for(unsigned int j=0; j<query.length(); j++)
      query[j] = tolower(query[j]);
    documents.push_back(query);
  }
  querys.close();
  cerr << "# Documents Number:"<<documents.size()<<endl;

  cerr << "# Create Gram List:" <<endl;
  for(unsigned int id = 0; id < keywords.size(); id++)
    appendGramList(keywords[id], id);
  addTail();
  cerr <<"# Minimum Length: "<<minLen+q-1 << endl;
  cerr <<"# Maximum Length: "<<maxLen+q-1 << endl;
  cerr <<"# Gram List Number: "<<gramListMap.size()<<endl;

  if(minLen <= tau * q) {
    cerr << "!!q is too short, choose a bigger q, Program will exit" << endl;
    exit(1);
  }

  gettimeofday(&sbegin,0);
  HeapBinary hbi(EditDistance);
  InvertList *invertList = new InvertList();
  vector<int> *index = new vector<int>(invertList->size(), 0);

  int size = 0;
  for(unsigned int i = 0; i < documents.size(); i++)
  {
#ifdef _output
    cerr << i << "/" << documents.size() << "\r\r";
#endif

    gettimeofday(&pb,0);

    lastEndIndex = documents[i].size();
    isFirstTime = true;

    indexTransfer.clear();
    unsigned int b = -1,e = -1;
    bool f = false,isFirst = true;
    for(unsigned int j = 0; j < documents[i].length() - q + 1; j++) {
      GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
      if(it != gramListMap.end())
      {
        invertList->push_back(make_pair(j, it->second));
        if(!f)
        {
            ++b;
            ++e;
            indexTransfer.push_back(make_pair(b,e));
        }
        else
        {
            ++e;
            f = false;
            indexTransfer.push_back(make_pair(b,e));
        }
        isFirst = false;
      }
      else
      {
        pair<int,Array<int>*> empty = make_pair(-1, new Array<int>());
        invertList->push_back(empty);
        if(isFirst)
        {
            indexTransfer.push_back(make_pair(0,0));
        }
        else if(!f)
        {
            ++b;
            f = true;
            indexTransfer.push_back(make_pair(b,e));
        }
        else
        {
            indexTransfer.push_back(make_pair(b,e));
        }

      }
    }
    if(invertList->size() == 0)
      continue;


    getListCount(invertList);
    cout << "invertList size:  " << invertList->size() << endl;
    cout << "count size:  " << listCount.size() << endl;
    cout << "documents[" << i << "].size - q + 1:  " << documents[i].size() - q + 1 << endl;

    for(int i = 0;i < listCount.size();i++)
    {
        //int t = (listCount[i]>0)?(listCount[i]-1):0;
        cout << listCount[i] << " ";
    }
    cout << endl;
    for(int i = 0;i < indexTransfer.size();i++)
    {
        cout << indexTransfer[i].first << " ";
    }
    cout << endl;
    for(int i = 0;i < indexTransfer.size();i++)
    {
        cout << indexTransfer[i].second << " ";
    }
    cout << endl;

    initialSum(documents[i].size() - q + 1);

    if(algorithm == "Optimal")
    {
        int param = atoi(argv[8]);
        initialOptimalIndex(documents[i].size() - q + 1);
        initialDpt(documents[i].size() - q + 1);
        optimal(documents[i].size() - q + 1,param);

        for(int i = 0;i < sum.size();i++)
        {
            cout << sum[i] << " ";
        }
        cout << endl;

        for(int i = 0;i < dpt.size();i++)
        {
            cout << dpt[i] << " ";
        }
        cout << endl;

        for(int i = 0;i < optimalIndex.size();i++)
        {
            cout << optimalIndex[i] << " ";
        }
        cout << endl;

        getBeginIndex();
        getEndIndex();
    }
    else if(algorithm == "Greedy")
    {
        double param = atof(argv[8]);
        initialGreedyIndex(documents[i].size() - q + 1);
        greedy(documents[i].size() - q + 1,param);

        for(int i = 0;i < sum.size();i++)
        {
            cout << sum[i] << " ";
        }
        cout << endl;

        for(int i = 0;i < greedyIndex.size();i++)
        {
            cout << greedyIndex[i] << " ";
        }
        cout << endl;

        getGreedyBeginIndex();
        getGreedyEndIndex();
    }
    else if(algorithm == "Random")
    {
        int param = atoi(argv[8]);
        initialRandomIndex(documents[i].size() - q + 1);
        random(documents[i].size() - q + 1,param);
        getRandomBeginIndex();
        getRandomEndIndex();
    }
//    //average
//    initialBeginIndex();
//    initialEndIndex();
//    algorithm = "Average";



    for(int j = beginIndex.size() - 1;j >= 0;j--)
    {
        cout << beginIndex[j] << " ";
    }
    cout << endl;
    for(int j = endIndex.size() - 1;j >= 0;j--)
    {
        cout << endIndex[j] << " ";
    }
    cout << endl;



    gettimeofday(&pe,0);
    pp += (pe.tv_sec-pb.tv_sec+(pe.tv_usec-pb.tv_usec)*1.0/CLOCKS_PER_SEC);
//    InvertList *invertList_part = new InvertList();
//    vector<int> *index_part = new vector<int>(invertList_part->size(),0);
    for(int j = 0;j < beginIndex.size();j++)
//    int m = lmax - getTE();
//    for(int j = optimalIndex.size() - 1 - lmax;j >= 0;j -= m)
    {
//        for(int k = beginIndex[j];k <= endIndex[j];k++)
//        {
//            if(invertList->at(k).first != -1)
//            {
//                invertList_part->push_back(invertList->at(k));
//            }
//
//        }
        if(invertList->size() == 0)
            continue;
        index->assign(invertList->size(),0);
        if(algorithm == "Optimal")
        {
            if(j == beginIndex.size() - 1)
                lastEndIndex == 0;
            else
                lastEndIndex = endIndex[j+1];
        }
        else
        {
            if(j == 0)
                lastEndIndex = 0;
            else
                lastEndIndex = endIndex[j-1];
        }

        if(beginIndex.size() == 1)
            lastEndIndex = 0;

        thisEndIndex = endIndex[j];
//        cout << "index.size:" << index->size() << endl;
        //cout << beginIndex[j] << " " << endIndex[j] << " " << lastEndIndex << endl;
        hbi.EntityExtractScale(invertList,index,documents[i],i,beginIndex[j],endIndex[j]);
        isFirstTime = false;
//        hbi.EntityExtractScale(invertList,index,documents[i],i,j,j+lmax);
        //invertList->clear();
    }
    size += beginIndex.size();
    //index->assign(invertList->size(), 0);
    //hbi.EntityExtractScale(invertList,index,documents[i],i);
    invertList->clear();
  }
  gettimeofday(&send,0);

  gettimeofday(&tb,0);
  for(unsigned int i = 0; i < result.size(); i++)
  {
    if(edth_imp(keywords[result[i].ent_id].c_str(),  documents[result[i].doc_id].substr(result[i].pos, result[i].len).c_str(), entityLengthMap[result[i].ent_id]+q-1, result[i].len, tau))
    {
      cout << endl << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
      cout << "##" << keywords[result[i].ent_id] << "$$" << endl;
      cout << "##" << documents[result[i].doc_id].substr(result[i].pos,result[i].len) << "$$" << endl;
      realResult++;
    }
  }
//    unordered_set<Cand>::iterator it = result2.begin();
//  while(it != result2.end())
//  {
//    if(edth_imp(keywords[it->ent_id].c_str(),  documents[it->doc_id].substr(it->pos, it->len).c_str(), entityLengthMap[it->ent_id]+q-1, it->len, tau))
//    {
//      cout << endl << it->ent_id << " " << it->doc_id << " " << it->pos << " " << it->len << endl;
//      cout << "##" << keywords[it->ent_id] << "$$" << endl;
//      cout << "##" << documents[it->doc_id].substr(it->pos,it->len) << "$$" << endl;
//      realResult++;
//    }
//    ++it;
//  }
  gettimeofday(&te, 0);
  tt = 0;
  tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
  cerr << "=====Search With Binary====="<<endl;
  cerr << "#       Algorithm : " << algorithm << endl;
  cerr << "#  Partition Time : " << pp << "s" << endl;
  cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC - pp << "s" << endl;
  cerr << "#     Verify Time : " << tt << "s" <<endl;
  cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
  cerr << "#  To Be Verified : " << result.size() << endl;
  cerr << "# Real Result Num : " << realResult << endl;
  cerr << "#    Average size : " << size / documents.size() << endl;
  cerr << "=====Search With Binary====="<<endl;

FILE *fp = fopen("record_ED.csv","a");
  //fprintf(fp, "\n");
  fprintf(fp, "#%d\n",q );
  fprintf(fp, "#%d\n" , tau );
  fprintf(fp, "#param %.2f\n" , atof(argv[8]) );
  fprintf(fp, "#PartitionTime %5.4f\n" , pp);
  fprintf(fp, "#FilterTime %5.4f\n" , send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC - pp);
  fprintf(fp, "#VerifyTime %5.4f\n" , tt );
  fprintf(fp, "#TotalTime %5.4f\n" , tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC);
  //fprintf(fp, "#ToBeVerified %d\n" , result.size() );
  //fprintf(fp, "#RealResultNum %d\n" , realResult );
  //fprintf(fp, "#AverageSize %d\n" , size / 100 );
  fprintf(fp, "\n");
  fclose(fp);

  result.clear();
  realResult = 0;

  delete index;
  delete invertList;
  GramListMap::iterator gbegin = gramListMap.begin();
  GramListMap::iterator gend = gramListMap.end();
  while(gbegin!=gend)
  {
    delete gbegin->second;
    gbegin++;
  }
  delete MINKEY;
  delete MAXKEY;

}


//**************************************************************************************************************//
//                              EDIT DISTANCE SIMILARITY  SCALE                                                                               //
//**************************************************************************************************************//

void EditSimilarityEntry(int argc, char** argv, int limit)
{
  q = atoi(argv[2]);
  det = atof(argv[3]);
  entity = string(argv[4]);
  document = string(argv[5]);
  algorithm = string(argv[7]);
  MAXKEY->append(MAX);

  type = "ES";

  vector<string> documents;
  int realResult = 0;
  timeval tb, te;
  double tt = 0;
  timeval sbegin,send;
  timeval pb,pe;
  double pp = 0;

  cerr <<"# Threshold: " << det << endl;
  cerr<<"# Loading Entities:"<<endl;
  ifstream input(entity.c_str(), ios::in);
  string keyword;
  int lim = 0;
  while(getline(input,keyword)&& lim < limit)
  {
    lim++;
    for(unsigned int j=0; j<keyword.length(); j++)
      keyword[j] = tolower(keyword[j]);
    keywords.push_back(keyword);
  }
  input.close();
  cerr << "# Entities Number: "<<keywords.size()<<endl;

  cerr << "# Loading Documents:"<<endl;
  ifstream querys(document.c_str(), ios::in);
  string query;
  while(getline(querys,query))
  {
    if(query.length() - q < 0)  continue;
    for(unsigned int j=0; j<query.length(); j++)
      query[j] = tolower(query[j]);
    documents.push_back(query);
  }
  querys.close();
  cerr << "# Documents Number:"<<documents.size()<<endl;

  cerr << "# Create Gram List:" <<endl;
  for(unsigned int id = 0; id < keywords.size(); id++)
    appendGramList(keywords[id], id);
  addTail();
  cerr <<"# Minimum Length: "<<minLen+q-1 << endl;
  cerr <<"# Maximum Length: "<<maxLen+q-1 << endl;
  cerr <<"# Gram List Number: "<<gramListMap.size()<<endl;

  if(minLen <= (minLen + q - 1) * (1 - det) * q) {
    cerr << "!!q is too short, choose a bigger q, Program will exit" << endl;
    exit(1);
  }

  gettimeofday(&sbegin,0);
  HeapBinary hbi(EditDistanceSimilarity);
  InvertList *invertList = new InvertList();
  vector<int> *index = new vector<int>(invertList->size(), 0);

  int size = 0;
  for(unsigned int i = 0; i < documents.size(); i++)
  {
#ifdef _output
    cerr << i << "/" << documents.size() << "\r\r";
#endif

    gettimeofday(&pb,0);
    lastEndIndex = documents[i].size();

    indexTransfer.clear();
    unsigned int b = 0,e = 0;
    bool f = false,isFirst = true;
    for(unsigned int j = 0; j < documents[i].length() - q + 1; j++) {
      GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
      if(it != gramListMap.end())
      {
        invertList->push_back(make_pair(j, it->second));
        if(!f)
        {
            ++b;
            ++e;
            indexTransfer.push_back(make_pair(b,e));
        }
        else
        {
            ++e;
            f = false;
            indexTransfer.push_back(make_pair(b,e));
        }
        isFirst = false;
      }
      else
      {
        pair<int,Array<int>*> empty = make_pair(-1, new Array<int>());
        invertList->push_back(empty);
        if(isFirst)
        {
            indexTransfer.push_back(make_pair(0,0));
        }
        else if(!f)
        {
            ++b;
            f = true;
            indexTransfer.push_back(make_pair(b,e));
        }
        else
        {
            indexTransfer.push_back(make_pair(b,e));
        }

      }
    }
    if(invertList->size() == 0)
      continue;


    getListCount(invertList);
    cout << "invertList size:  " << invertList->size() << endl;
    cout << "count size:  " << listCount.size() << endl;
    cout << "documents[" << i << "].size - q + 1:  " << documents[i].size() - q + 1 << endl;

    for(int i = 0;i < listCount.size();i++)
    {
        //int t = (listCount[i]>0)?(listCount[i]-1):0;
        cout << listCount[i] << " ";
    }
    cout << endl;

    initialSum(documents[i].size() - q + 1);

    if(algorithm == "Optimal")
    {
        int param = atoi(argv[8]);
        initialOptimalIndex(documents[i].size() - q + 1);
        initialDpt(documents[i].size() - q + 1);
        optimal(documents[i].size() - q + 1,param);

        for(int i = 0;i < sum.size();i++)
        {
            cout << sum[i] << " ";
        }
        cout << endl;

        for(int i = 0;i < dpt.size();i++)
        {
            cout << dpt[i] << " ";
        }
        cout << endl;

        for(int i = 0;i < optimalIndex.size();i++)
        {
            cout << optimalIndex[i] << " ";
        }
        cout << endl;

        getBeginIndex();
        getEndIndex();
    }
    else if(algorithm == "Greedy")
    {
        double param = atof(argv[8]);
        initialGreedyIndex(documents[i].size() - q + 1);
        greedy(documents[i].size() - q + 1,param);
        getGreedyBeginIndex();
        getGreedyEndIndex();
    }
    else if(algorithm == "Random")
    {
        int param = atoi(argv[8]);
        initialRandomIndex(documents[i].size() - q + 1);
        random(documents[i].size() - q + 1,param);
        getRandomBeginIndex();
        getRandomEndIndex();
    }

    for(int j = beginIndex.size() - 1;j >= 0;j--)
    {
        cout << beginIndex[j] << " ";
    }
    cout << endl;
    for(int j = endIndex.size() - 1;j >= 0;j--)
    {
        cout << endIndex[j] << " ";
    }
    cout << endl;

    gettimeofday(&pe,0);
    pp += (pe.tv_sec-pb.tv_sec+(pe.tv_usec-pb.tv_usec)*1.0/CLOCKS_PER_SEC);

    for(int j = 0;j < beginIndex.size();j++)
    {
        if(invertList->size() == 0)
            continue;
        index->assign(invertList->size(),0);
        if(algorithm == "Optimal")
        {
            if(j == beginIndex.size() - 1)
                lastEndIndex == 0;
            else
                lastEndIndex = endIndex[j+1];
        }
        else
        {
            if(j == 0)
                lastEndIndex = 0;
            else
                lastEndIndex = endIndex[j-1];
        }

        if(beginIndex.size() == 1)
            lastEndIndex = 0;

        thisEndIndex = endIndex[j];
//        cout << "index.size:" << index->size() << endl;
        hbi.edsEntityExtract(invertList,index,documents[i],i,beginIndex[j],endIndex[j]);
        isFirstTime = false;
    }

    size += beginIndex.size();
    invertList->clear();
  }
  gettimeofday(&send,0);

  gettimeofday(&tb,0);
  for(unsigned int i = 0; i < result.size(); i++)
  {
    if(eds_imp(keywords[result[i].ent_id].c_str(),  documents[result[i].doc_id].substr(result[i].pos, result[i].len).c_str(), entityLengthMap[result[i].ent_id]+q-1, result[i].len, det))
    {
      //cout << endl << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
      //cout << "##" << keywords[result[i].ent_id] << "$$" << endl;
      //cout << "##" << documents[result[i].doc_id].substr(result[i].pos,result[i].len) << "$$" << endl;
      realResult++;
    }
  }
  gettimeofday(&te, 0);
  tt = 0;
  tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
  cerr << "=====Search With Binary====="<<endl;
  cerr << "#       Algorithm : " << algorithm << endl;
  cerr << "#  Partition Time : " << pp << "s" << endl;
  cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC - pp << "s" << endl;
  cerr << "#     Verify Time : " << tt << "s" <<endl;
  cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
  cerr << "#  To Be Verified : " << result.size() << endl;
  cerr << "# Real Result Num : " << realResult << endl;
  cerr << "#    Average size : " << size / documents.size() << endl;
  cerr << "=====Search With Binary====="<<endl;

  FILE *fp = fopen("record_EDS.csv","a");
  //fprintf(fp, "\n");
  fprintf(fp, "#%d\n",q );
  fprintf(fp, "#%f\n" , det );
  fprintf(fp, "#param %.2f\n" , atof(argv[8]) );
  fprintf(fp, "#PartitionTime %5.4f\n" , pp);
  fprintf(fp, "#FilterTime %5.4f\n" , send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC - pp);
  fprintf(fp, "#VerifyTime %5.4f\n" , tt );
  fprintf(fp, "#TotalTime %5.4f\n" , send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC);
  //fprintf(fp, "#ToBeVerified %d\n" , result.size() );
  //fprintf(fp, "#RealResultNum %d\n" , realResult );
  //fprintf(fp, "#AverageSize %d\n" , size / 100 );
  fprintf(fp, "\n");
  fclose(fp);

  result.clear();
  realResult = 0;

  delete index;
  delete invertList;
  GramListMap::iterator gbegin = gramListMap.begin();
  GramListMap::iterator gend = gramListMap.end();
  while(gbegin!=gend)
  {
    delete gbegin->second;
    gbegin++;
  }
  delete MINKEY;
  delete MAXKEY;
}

void EditSimilarityVerify(int argc, char** argv, int limit)
{
  q = atoi(argv[2]);
  det = atof(argv[3]);
  entity = string(argv[4]);
  document = string(argv[5]);
  algorithm = string(argv[7]);
  MAXKEY->append(MAX);

  type = "ES";

  vector<string> documents;
  int realResult = 0;
  timeval sbegin,send;
  timeval pb,pe;
  double pp;

  cerr <<"# Threshold: " << det << endl;
  cerr<<"# Loading Entities:"<<endl;
  ifstream input(entity.c_str(), ios::in);
  string keyword;
  int lim = 0;
  while(getline(input,keyword)&& lim < limit)
  {
    lim++;
    for(unsigned int j=0; j<keyword.length(); j++)
      keyword[j] = tolower(keyword[j]);
    keywords.push_back(keyword);
  }
  input.close();
  cerr << "# Entities Number: "<<keywords.size()<<endl;

  cerr << "# Loading Documents:"<<endl;
  ifstream querys(document.c_str(), ios::in);
  string query;
  while(getline(querys,query))
  {
    if(query.length() - q < 0)  continue;
    for(unsigned int j=0; j<query.length(); j++)
      query[j] = tolower(query[j]);
    documents.push_back(query);
  }
  querys.close();
  cerr << "# Documents Number:"<<documents.size()<<endl;

  cerr << "# Create Gram List:" <<endl;
  for(unsigned int id = 0; id < keywords.size(); id++)
    appendGramList(keywords[id], id);
  addTail();
  cerr <<"# Minimum Length: "<<minLen+q-1 << endl;
  cerr <<"# Maximum Length: "<<maxLen+q-1 << endl;
  cerr <<"# Gram List Number: "<<gramListMap.size()<<endl;

  if(minLen <= (minLen + q - 1) * (1 - det) * q) {
    cerr << "!!q is too short, choose a bigger q, Program will exit" << endl;
    exit(1);
  }

  gettimeofday(&sbegin,0);
  HeapBinary hbi(EditDistanceSimilarity);
  InvertList *invertList = new InvertList();
  vector<int> *index = new vector<int>(invertList->size(), 0);

  int size = 0;
  for(unsigned int i = 0; i < documents.size(); i++)
  {
#ifdef _output
    cerr << i << "/" << documents.size() << "\r\r";
#endif

    gettimeofday(&pb,0);
    lastEndIndex = documents[i].size();

    indexTransfer.clear();
    unsigned int b = 0,e = 0;
    bool f = false,isFirst = true;
    for(unsigned int j = 0; j < documents[i].length() - q + 1; j++) {
      GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
      if(it != gramListMap.end())
      {
        invertList->push_back(make_pair(j, it->second));
        if(!f)
        {
            ++b;
            ++e;
            indexTransfer.push_back(make_pair(b,e));
        }
        else
        {
            ++e;
            f = false;
            indexTransfer.push_back(make_pair(b,e));
        }
        isFirst = false;
      }
      else
      {
        pair<int,Array<int>*> empty = make_pair(-1, new Array<int>());
        invertList->push_back(empty);
        if(isFirst)
        {
            indexTransfer.push_back(make_pair(0,0));
        }
        else if(!f)
        {
            ++b;
            f = true;
            indexTransfer.push_back(make_pair(b,e));
        }
        else
        {
            indexTransfer.push_back(make_pair(b,e));
        }

      }
    }
    if(invertList->size() == 0)
      continue;


    getListCount(invertList);
    cout << "invertList size:  " << invertList->size() << endl;
    cout << "count size:  " << listCount.size() << endl;
    cout << "documents[" << i << "].size - q + 1:  " << documents[i].size() - q + 1 << endl;

    for(int i = 0;i < listCount.size();i++)
    {
        //int t = (listCount[i]>0)?(listCount[i]-1):0;
        cout << listCount[i] << " ";
    }
    cout << endl;

    initialSum(documents[i].size() - q + 1);

    if(algorithm == "Optimal")
    {
        int param = atoi(argv[8]);
        initialOptimalIndex(documents[i].size() - q + 1);
        initialDpt(documents[i].size() - q + 1);
        optimal(documents[i].size() - q + 1,param);

        for(int i = 0;i < sum.size();i++)
        {
            cout << sum[i] << " ";
        }
        cout << endl;

        for(int i = 0;i < dpt.size();i++)
        {
            cout << dpt[i] << " ";
        }
        cout << endl;

        for(int i = 0;i < optimalIndex.size();i++)
        {
            cout << optimalIndex[i] << " ";
        }
        cout << endl;

        getBeginIndex();
        getEndIndex();
    }
    else if(algorithm == "Greedy")
    {
        double param = atof(argv[8]);
        initialGreedyIndex(documents[i].size() - q + 1);
        greedy(documents[i].size() - q + 1,param);
        getGreedyBeginIndex();
        getGreedyEndIndex();
    }
    else if(algorithm == "Random")
    {
        int param = atoi(argv[8]);
        initialRandomIndex(documents[i].size() - q + 1);
        random(documents[i].size() - q + 1,param);
        getRandomBeginIndex();
        getRandomEndIndex();
    }

    for(int j = beginIndex.size() - 1;j >= 0;j--)
    {
        cout << beginIndex[j] << " ";
    }
    cout << endl;
    for(int j = endIndex.size() - 1;j >= 0;j--)
    {
        cout << endIndex[j] << " ";
    }
    cout << endl;

    gettimeofday(&pe,0);
    pp += (pe.tv_sec-pb.tv_sec+(pe.tv_usec-pb.tv_usec)*1.0/CLOCKS_PER_SEC);

    for(int j = 0;j < beginIndex.size();j++)
    {
        if(invertList->size() == 0)
            continue;
        index->assign(invertList->size(),0);
        if(algorithm == "Optimal")
        {
            if(j == beginIndex.size() - 1)
                lastEndIndex == 0;
            else
                lastEndIndex = endIndex[j+1];
        }
        else
        {
            if(j == 0)
                lastEndIndex = 0;
            else
                lastEndIndex = endIndex[j-1];
        }
        if(beginIndex.size() == 1)
            lastEndIndex = 0;

        thisEndIndex = endIndex[j];
//        cout << "index.size:" << index->size() << endl;
        hbi.edsEntityExtractVerify(invertList,index,documents[i],i,beginIndex[j],endIndex[j]);
        isFirstTime = false;
    }
    size += beginIndex.size();

    invertList->clear();
  }
  gettimeofday(&send,0);

  cerr << "=====Search With Binary====="<<endl;
  cerr << "#  Partition Time : " << pp << "s" << endl;
  cerr << "#       Algorithm : " << algorithm << endl;
  cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC - pp<<"s"<< endl;
  cerr << "#      Total Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
  cerr << "#  To Be Verified : " << canda << endl;
  cerr << "# Real Result Num : " << result.size() << endl;
  cerr << "#    Average size : " << size / documents.size() << endl;
  cerr << "=====Search With Binary====="<<endl;

  FILE *fp = fopen("record_EDS.csv","a");
  //fprintf(fp, "\n");
  fprintf(fp, "#%d\n",q );
  fprintf(fp, "#%f\n" , det );
  fprintf(fp, "#param %.2f\n" , atof(argv[8]) );
  fprintf(fp, "#PartitionTime %5.4f\n" , pp);
  fprintf(fp, "#FilterTime %5.4f\n" , send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC - pp);
  //fprintf(fp, "#     Verify Time : %.6fs\n" , tt );
  fprintf(fp, "#TotalTime %5.4f\n" , send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC);
  //fprintf(fp, "#ToBeVerified %d\n" , result.size() );
  //fprintf(fp, "#RealResultNum %d\n" , realResult );
  //fprintf(fp, "#AverageSize %d\n" , size / 100 );
  fprintf(fp, "\n");
  fclose(fp);

  result.clear();
  realResult = 0;

  delete index;
  delete invertList;
  GramListMap::iterator gbegin = gramListMap.begin();
  GramListMap::iterator gend = gramListMap.end();
  while(gbegin!=gend)
  {
    delete gbegin->second;
    gbegin++;
  }
  delete MINKEY;
  delete MAXKEY;
}

//**************************************************************************************************************//
//                              JACCARD  SIMILARITY  SCALE                                                                                        //
//**************************************************************************************************************//

void JaccardEntry(int argc, char** argv, int limit)
{
  q = atoi(argv[2]);
  det = atof(argv[3]);
  entity = string(argv[4]);
  document = string(argv[5]);
  algorithm = string(argv[7]);
  MAXKEY->append(MAX);

  type = "JS";

  vector<string> documents;
  int realResult = 0;
  timeval tb, te;
  double tt = 0;
  timeval sbegin,send;
  timeval pb,pe;
  double pp = 0;
  timeval bb,be;
  double bt = 0;

  cerr <<"# Threshold: " << det << endl;
  cerr<<"# Loading Entities:"<<endl;
  ifstream input(entity.c_str(), ios::in);
  string keyword;
  int lim = 0;
  while(getline(input,keyword)&& lim < limit)
  {
    lim++;
    for(unsigned int j=0; j<keyword.length(); j++)
      keyword[j] = tolower(keyword[j]);
    keywords.push_back(keyword);
  }
  input.close();
  cerr << "# Entities Number: "<<keywords.size()<<endl;

  cerr << "# Loading Documents:"<<endl;
  ifstream querys(document.c_str(), ios::in);
  string query;
  while(getline(querys,query))
  {
    for(unsigned int j=0; j<query.length(); j++)
      query[j] = tolower(query[j]);
    documents.push_back(query);
  }
  querys.close();
  cerr << "# Documents Number:"<<documents.size()<<endl;

  cerr << "# Create Gram List:" <<endl;
  for(unsigned int id = 0; id < keywords.size(); id++)
    appendGramList(keywords[id], id);
    //appendTokenList(keywords[id], id);
  addTail();
  cerr <<"# Minimum Length: "<< minLen << endl;
  cerr <<"# Maximum Length: "<< maxLen << endl;
  cerr <<"# Gram List Number: "<< gramListMap.size() <<endl;

//  vector<vector<string>*> tokendoc;
//  for(unsigned int i = 0; i < documents.size(); i++) {
//    vector<string> *tokens = new vector<string>();
//    strToTokens(documents[i], tokens, delim);
//    tokendoc.push_back(tokens);
//  }

  gettimeofday(&sbegin,0);
  HeapBinary hbi(Jaccard);
  InvertList *invertList = new InvertList();
  vector<int> *index = new vector<int>(invertList->size(), 0);

  int size = 0;
  for(unsigned int i = 0; i < documents.size(); i++)
  {
#ifdef _output
    cerr << i << "/" << documents.size() << "\r\r";
#endif

    gettimeofday(&pb,0);
    lastEndIndex = documents[i].size();

    indexTransfer.clear();
    unsigned int b = 0,e = 0;
    bool f = false,isFirst = true;
    for(unsigned int j = 0; j < documents[i].length() - q + 1; j++) {
      GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
      if(it != gramListMap.end())
      {
        invertList->push_back(make_pair(j, it->second));
        if(!f)
        {
            ++b;
            ++e;
            indexTransfer.push_back(make_pair(b,e));
        }
        else
        {
            ++e;
            f = false;
            indexTransfer.push_back(make_pair(b,e));
        }
        isFirst = false;
      }
      else
      {
        pair<int,Array<int>*> empty = make_pair(-1, new Array<int>());
        invertList->push_back(empty);
        if(isFirst)
        {
            indexTransfer.push_back(make_pair(0,0));
        }
        else if(!f)
        {
            ++b;
            f = true;
            indexTransfer.push_back(make_pair(b,e));
        }
        else
        {
            indexTransfer.push_back(make_pair(b,e));
        }

      }
    }
    if(invertList->size() == 0)
      continue;


    getListCount(invertList);
    cout << "invertList size:  " << invertList->size() << endl;
    cout << "count size:  " << listCount.size() << endl;
    cout << "documents[" << i << "].size - q + 1:  " << documents[i].size() - q + 1 << endl;

    for(int i = 0;i < listCount.size();i++)
    {
        //int t = (listCount[i]>0)?(listCount[i]-1):0;
        cout << listCount[i] << " ";
    }
    cout << endl;

    initialSum(documents[i].size() - q + 1);

    if(algorithm == "Optimal")
    {
        int param = atoi(argv[8]);
        initialOptimalIndex(documents[i].size() - q + 1);
        initialDpt(documents[i].size() - q + 1);
        optimal(documents[i].size() - q + 1,param);

        for(int i = 0;i < sum.size();i++)
        {
            cout << sum[i] << " ";
        }
        cout << endl;

        for(int i = 0;i < dpt.size();i++)
        {
            cout << dpt[i] << " ";
        }
        cout << endl;

        for(int i = 0;i < optimalIndex.size();i++)
        {
            cout << optimalIndex[i] << " ";
        }
        cout << endl;

        getBeginIndex();
        getEndIndex();
    }
    else if(algorithm == "Greedy")
    {
        double param = atof(argv[8]);
        initialGreedyIndex(documents[i].size() - q + 1);
        greedy(documents[i].size() - q + 1,param);
        getGreedyBeginIndex();
        getGreedyEndIndex();
    }
    else if(algorithm == "Random")
    {
        int param = atoi(argv[8]);
        initialRandomIndex(documents[i].size() - q + 1);
        random(documents[i].size() - q + 1,param);
        getRandomBeginIndex();
        getRandomEndIndex();
    }


//    cout << "beginIndex.size: " << beginIndex.size() << endl;
//    for(int j = beginIndex.size() - 1;j >= 0;j--)
//    {
//        cout << beginIndex[j] << " ";
//    }
//    cout << endl;
//    cout << "endIndex.size: " << endIndex.size() << endl;
//    for(int j = endIndex.size() - 1;j >= 0;j--)
//    {
//        cout << endIndex[j] << " ";
//    }
//    cout << endl;

    gettimeofday(&pe,0);
    pp += (pe.tv_sec-pb.tv_sec+(pe.tv_usec-pb.tv_usec)*1.0/CLOCKS_PER_SEC);


    for(int j = 0;j < beginIndex.size();j++)
    {
        //int j = 0;
        index->assign(invertList->size(),0);
        if(algorithm == "Optimal")
        {
            if(j == beginIndex.size() - 1)
                lastEndIndex == 0;
            else
                lastEndIndex = endIndex[j+1];
        }
        else
        {
            if(j == 0)
                lastEndIndex = 0;
            else
                lastEndIndex = endIndex[j-1];
        }
        if(beginIndex.size() == 1)
            lastEndIndex = 0;

        thisEndIndex = endIndex[j];
        gettimeofday(&bb,0);
        hbi.jacEntityExtract(invertList,index,documents[i],i,beginIndex[j],endIndex[j]);
        gettimeofday(&be,0);

        bt += (be.tv_sec-bb.tv_sec+(be.tv_usec-bb.tv_usec)*1.0/CLOCKS_PER_SEC);
//        cout << "index.size:" << index->size() << endl;
//        hbi.EntityExtractScale(invertList,index,documents[i],i,beginIndex[j],endIndex[j]);
    }
    size += beginIndex.size();
    cout << "ok" << endl;

    invertList->clear();
  }
  gettimeofday(&send,0);

  gettimeofday(&tb,0);
//  for(unsigned int i = 0; i < result.size(); i++)
//  {
//    if(jaccardToken(keywords[result[i].ent_id], tokendoc[result[i].doc_id], result[i].pos, result[i].len) >= det)
//    {
//      //cerr << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
//      realResult++;
//    }
//  }
  gettimeofday(&te, 0);
  tt = 0;
  tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
  cerr << "=====Search With Binary====="<<endl;
  cerr << "#       Algorithm : " << algorithm << endl;
  cerr << "#  Partition Time : " << pp << "s" << endl;
  cerr << "#     Binary Time : " << bt << "s" << endl;
  cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC - pp << "s" << endl;
//  cerr << "#     Verify Time : " << tt << "s" <<endl;
  cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
//  cerr << "#  To Be Verified : " << result.size() << endl;
  cerr << "# Real Result Num : " << result.size() << endl;
  cerr << "#    Average size : " << size / documents.size() << endl;
  cerr << "=====Search With Binary====="<<endl;

  FILE *fp = fopen("record_JAC.csv","a");
  //fprintf(fp, "\n");
  fprintf(fp, "#%d\n",q );
  fprintf(fp, "#%f\n" , det );
  fprintf(fp, "#param %.2f\n" , atof(argv[7]) );
  fprintf(fp, "#PartitionTime %5.4f\n" , pp);
  fprintf(fp, "#FilterTime %5.4f\n" , send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC - pp);
  fprintf(fp, "#VerifyTime %5.4f\n" , tt );
  fprintf(fp, "#TotalTime %5.4f\n" , tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC);
  //fprintf(fp, "#ToBeVerified %d\n" , result.size() );
  //fprintf(fp, "#RealResultNum %d\n" , realResult );
  //fprintf(fp, "#AverageSize %d\n" , size / 100 );
  fprintf(fp, "\n");
  fclose(fp);

  result.clear();
  realResult = 0;

//  vector< vector<string>* >::iterator ibegin = tokendoc.begin();
//  vector< vector<string>* >::iterator iend = tokendoc.end();
//  while(ibegin!=iend)
//  {
//    delete *ibegin;
//    ibegin++;
//  }
  delete index;
  delete invertList;
  GramListMap::iterator gbegin = gramListMap.begin();
  GramListMap::iterator gend = gramListMap.end();
  while(gbegin!=gend)
  {
    delete gbegin->second;
    gbegin++;
  }
  delete MINKEY;
  delete MAXKEY;
}

//**************************************************************************************************************//
//                              DICE  SIMILARITY SCALE                                                                                                //
//**************************************************************************************************************//

void DiceEntry(int argc, char** argv, int limit)
{
  q = atoi(argv[2]);
  det = atof(argv[3]);
  entity = string(argv[4]);
  document = string(argv[5]);
  algorithm = string(argv[7]);
  MAXKEY->append(MAX);

  type = "DS";

  vector<string> documents;
  int realResult = 0;
  timeval tb, te;
  double tt = 0;
  timeval sbegin,send;
  timeval pb,pe;
  double pp = 0;

  cerr <<"# Threshold: " << det << endl;
  cerr<<"# Loading Entities:"<<endl;
  ifstream input(entity.c_str(), ios::in);
  string keyword;
  int lim = 0;
  while(getline(input,keyword)&& lim < limit)
  {
    lim++;
    for(unsigned int j=0; j<keyword.length(); j++)
      keyword[j] = tolower(keyword[j]);
    keywords.push_back(keyword);
  }
  input.close();
  cerr << "# Entities Number: "<<keywords.size()<<endl;

  cerr << "# Loading Documents:"<<endl;
  ifstream querys(document.c_str(), ios::in);
  string query;
  while(getline(querys,query))
  {
    if(query.length() - q < 0)  continue;
    for(unsigned int j=0; j<query.length(); j++)
      query[j] = tolower(query[j]);
    documents.push_back(query);
  }

  querys.close();

  cerr << "# Documents Number:"<<documents.size()<<endl;
  cerr << "# Create Gram List:" <<endl;

  for(unsigned int id = 0; id < keywords.size(); id++)
    appendGramList(keywords[id], id);
  addTail();

  cerr <<"# Minimum Length: "<<minLen+q-1 << endl;
  cerr <<"# Maximum Length: "<<maxLen+q-1 << endl;
  cerr <<"# Gram List Number: "<<gramListMap.size()<<endl;

  gettimeofday(&sbegin,0);
  HeapBinary hbi(Dice);

  InvertList *invertList = new InvertList();
  vector<int> *index = new vector<int>(invertList->size(), 0);

  int size = 0;
  for(unsigned int i = 0; i < documents.size(); i++)
  {
#ifdef _output
    cerr << i << "/" << documents.size() << "\r\r";
#endif

    gettimeofday(&pb,0);
    lastEndIndex = documents[i].size();

    indexTransfer.clear();
    unsigned int b = 0,e = 0;
    bool f = false,isFirst = true;
    for(unsigned int j = 0; j < documents[i].length() - q + 1; j++) {
      GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
      if(it != gramListMap.end())
      {
        invertList->push_back(make_pair(j, it->second));
        if(!f)
        {
            ++b;
            ++e;
            indexTransfer.push_back(make_pair(b,e));
        }
        else
        {
            ++e;
            f = false;
            indexTransfer.push_back(make_pair(b,e));
        }
        isFirst = false;
      }
      else
      {
        pair<int,Array<int>*> empty = make_pair(-1, new Array<int>());
        invertList->push_back(empty);
        if(isFirst)
        {
            indexTransfer.push_back(make_pair(0,0));
        }
        else if(!f)
        {
            ++b;
            f = true;
            indexTransfer.push_back(make_pair(b,e));
        }
        else
        {
            indexTransfer.push_back(make_pair(b,e));
        }

      }
    }
    if(invertList->size() == 0)
      continue;


    getListCount(invertList);
    cout << "invertList size:  " << invertList->size() << endl;
    cout << "count size:  " << listCount.size() << endl;
    cout << "documents[" << i << "].size - q + 1:  " << documents[i].size() - q + 1 << endl;

    for(int i = 0;i < listCount.size();i++)
    {
        //int t = (listCount[i]>0)?(listCount[i]-1):0;
        cout << listCount[i] << " ";
    }
    cout << endl;

    initialSum(documents[i].size() - q + 1);

    if(algorithm == "Optimal")
    {
        int param = atoi(argv[8]);
        initialOptimalIndex(documents[i].size() - q + 1);
        initialDpt(documents[i].size() - q + 1);
        optimal(documents[i].size() - q + 1,param);

        for(int i = 0;i < sum.size();i++)
        {
            cout << sum[i] << " ";
        }
        cout << endl;

        for(int i = 0;i < dpt.size();i++)
        {
            cout << dpt[i] << " ";
        }
        cout << endl;

        for(int i = 0;i < optimalIndex.size();i++)
        {
            cout << optimalIndex[i] << " ";
        }
        cout << endl;

        getBeginIndex();
        getEndIndex();
    }
    else if(algorithm == "Greedy")
    {
        double param = atof(argv[8]);
        initialGreedyIndex(documents[i].size() - q + 1);
        greedy(documents[i].size() - q + 1,param);
        getGreedyBeginIndex();
        getGreedyEndIndex();
    }
    else if(algorithm == "Random")
    {
        int param = atoi(argv[8]);
        initialRandomIndex(documents[i].size() - q + 1);
        random(documents[i].size() - q + 1,param);
        getRandomBeginIndex();
        getRandomEndIndex();
    }

    for(int j = beginIndex.size() - 1;j >= 0;j--)
    {
        cout << beginIndex[j] << " ";
    }
    cout << endl;
    for(int j = endIndex.size() - 1;j >= 0;j--)
    {
        cout << endIndex[j] << " ";
    }
    cout << endl;

    gettimeofday(&pe,0);
    pp += (pe.tv_sec-pb.tv_sec+(pe.tv_usec-pb.tv_usec)*1.0/CLOCKS_PER_SEC);
    for(int j = 0;j < beginIndex.size();j++)
    {
        if(invertList->size() == 0)
            continue;
        index->assign(invertList->size(),0);
        if(algorithm == "Optimal")
        {
            if(j == beginIndex.size() - 1)
                lastEndIndex == 0;
            else
                lastEndIndex = endIndex[j+1];
        }
        else
        {
            if(j == 0)
                lastEndIndex = 0;
            else
                lastEndIndex = endIndex[j-1];
        }
        if(beginIndex.size() == 1)
            lastEndIndex = 0;

        thisEndIndex = endIndex[j];
//        cout << "index.size:" << index->size() << endl;
        hbi.diceEntityExtract(invertList,index,documents[i],i,beginIndex[j],endIndex[j]);
        isFirstTime = false;
    }
    size += beginIndex.size();

    invertList->clear();
  }
  gettimeofday(&send,0);

  gettimeofday(&tb,0);
  /*for(unsigned int i = 0; i < result.size(); i++)
    {
    realResult++;
    }*/
  gettimeofday(&te, 0);
  tt = 0;
  tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
  cerr << "=====Search With Binary====="<<endl;
  cerr << "#       Algorithm : " << algorithm << endl;
  cerr << "#  Partition Time : " << pp << "s" << endl;
  cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC - pp << "s" << endl;
  //cerr << "#     Verify Time : " << tt << "s" <<endl;
  cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
  //cerr << "#  To Be Verified : " << result.size() << endl;
  cerr << "# Real Result Num : " << result.size() << endl;
  cerr << "#    Average size : " << size / 100 << endl;
  cerr << "=====Search With Binary====="<<endl;

  result.clear();
  realResult = 0;

  delete index;
  delete invertList;

  GramListMap::iterator gbegin = gramListMap.begin();
  GramListMap::iterator gend = gramListMap.end();
  while(gbegin!=gend)
  {
    delete gbegin->second;
    gbegin++;
  }
  delete MINKEY;
  delete MAXKEY;
}

//**************************************************************************************************************//
//                              COSINE  SIMILARITY SCALE                                                                                          //
//**************************************************************************************************************//

void CosineEntry(int argc, char** argv, int limit)
{
  q = atoi(argv[2]);
  det = atof(argv[3]);
  entity = string(argv[4]);
  document = string(argv[5]);
  algorithm = string(argv[7]);
  MAXKEY->append(MAX);

    type = "CS";

  vector<string> documents;
  int realResult = 0;
  timeval tb, te;
  double tt = 0;
  timeval sbegin,send;
  timeval pb,pe;
  double pp = 0;

  cerr <<"# Threshold: " << det << endl;
  cerr<<"# Loading Entities:"<<endl;
  ifstream input(entity.c_str(), ios::in);
  string keyword;
  int lim = 0;
  while(getline(input,keyword)&& lim < limit)
  {
    lim++;
    for(unsigned int j=0; j<keyword.length(); j++)
      keyword[j] = tolower(keyword[j]);
    keywords.push_back(keyword);
  }
  input.close();
  cerr << "# Entities Number: "<<keywords.size()<<endl;

  cerr << "# Loading Documents:"<<endl;
  ifstream querys(document.c_str(), ios::in);
  string query;
  while(getline(querys,query))
  {
    if(query.length() - q < 0)  continue;
    for(unsigned int j=0; j<query.length(); j++)
      query[j] = tolower(query[j]);
    documents.push_back(query);
  }

  querys.close();

  cerr << "# Documents Number:"<<documents.size()<<endl;
  cerr << "# Create Gram List:" <<endl;

  for(unsigned int id = 0; id < keywords.size(); id++)
    appendGramList(keywords[id], id);
  addTail();

  cerr <<"# Minimum Length: "<<minLen+q-1 << endl;
  cerr <<"# Maximum Length: "<<maxLen+q-1 << endl;
  cerr <<"# Gram List Number: "<<gramListMap.size()<<endl;

  gettimeofday(&sbegin,0);
  HeapBinary hbi(Cosine);

  InvertList *invertList = new InvertList();
  vector<int> *index = new vector<int>(invertList->size(), 0);

  int size = 0;
  for(unsigned int i = 0; i < documents.size(); i++)
  {
#ifdef _output
    cerr << i << "/" << documents.size() << "\r\r";
#endif

    gettimeofday(&pb,0);
    lastEndIndex = documents[i].size();

    indexTransfer.clear();
    unsigned int b = 0,e = 0;
    bool f = false,isFirst = true;
    for(unsigned int j = 0; j < documents[i].length() - q + 1; j++) {
      GramListMap::iterator it = gramListMap.find(documents[i].substr(j,q));
      if(it != gramListMap.end())
      {
        invertList->push_back(make_pair(j, it->second));
        if(!f)
        {
            ++b;
            ++e;
            indexTransfer.push_back(make_pair(b,e));
        }
        else
        {
            ++e;
            f = false;
            indexTransfer.push_back(make_pair(b,e));
        }
        isFirst = false;
      }
      else
      {
        pair<int,Array<int>*> empty = make_pair(-1, new Array<int>());
        invertList->push_back(empty);
        if(isFirst)
        {
            indexTransfer.push_back(make_pair(0,0));
        }
        else if(!f)
        {
            ++b;
            f = true;
            indexTransfer.push_back(make_pair(b,e));
        }
        else
        {
            indexTransfer.push_back(make_pair(b,e));
        }

      }
    }
    if(invertList->size() == 0)
      continue;


    getListCount(invertList);
    cout << "invertList size:  " << invertList->size() << endl;
    cout << "count size:  " << listCount.size() << endl;
    cout << "documents[" << i << "].size - q + 1:  " << documents[i].size() - q + 1 << endl;

    for(int i = 0;i < listCount.size();i++)
    {
        //int t = (listCount[i]>0)?(listCount[i]-1):0;
        cout << listCount[i] << " ";
    }
    cout << endl;

    initialSum(documents[i].size() - q + 1);

    if(algorithm == "Optimal")
    {
        int param = atoi(argv[8]);
        initialOptimalIndex(documents[i].size() - q + 1);
        initialDpt(documents[i].size() - q + 1);
        optimal(documents[i].size() - q + 1,param);

        for(int i = 0;i < sum.size();i++)
        {
            cout << sum[i] << " ";
        }
        cout << endl;

        for(int i = 0;i < dpt.size();i++)
        {
            cout << dpt[i] << " ";
        }
        cout << endl;

        for(int i = 0;i < optimalIndex.size();i++)
        {
            cout << optimalIndex[i] << " ";
        }
        cout << endl;

        getBeginIndex();
        getEndIndex();
    }
    else if(algorithm == "Greedy")
    {
        double param = atof(argv[8]);
        initialGreedyIndex(documents[i].size() - q + 1);
        greedy(documents[i].size() - q + 1,param);
        getGreedyBeginIndex();
        getGreedyEndIndex();
    }
    else if(algorithm == "Random")
    {
        int param = atoi(argv[8]);
        initialRandomIndex(documents[i].size() - q + 1);
        random(documents[i].size() - q + 1,param);
        getRandomBeginIndex();
        getRandomEndIndex();
    }

    for(int j = beginIndex.size() - 1;j >= 0;j--)
    {
        cout << beginIndex[j] << " ";
    }
    cout << endl;
    for(int j = endIndex.size() - 1;j >= 0;j--)
    {
        cout << endIndex[j] << " ";
    }
    cout << endl;

    gettimeofday(&pe,0);
    pp += (pe.tv_sec-pb.tv_sec+(pe.tv_usec-pb.tv_usec)*1.0/CLOCKS_PER_SEC);

    for(int j = 0;j < beginIndex.size();j++)
    {
        if(invertList->size() == 0)
            continue;
        index->assign(invertList->size(),0);
        if(algorithm == "Optimal")
        {
            if(j == beginIndex.size() - 1)
                lastEndIndex == 0;
            else
                lastEndIndex = endIndex[j+1];
        }
        else
        {
            if(j == 0)
                lastEndIndex = 0;
            else
                lastEndIndex = endIndex[j-1];
        }
        if(beginIndex.size() == 1)
            lastEndIndex = 0;

        thisEndIndex = endIndex[j];
//        cout << "index.size:" << index->size() << endl;
        hbi.cosEntityExtract(invertList,index,documents[i],i,beginIndex[j],endIndex[j]);
        isFirstTime = false;
//        hbi.EntityExtractScale(invertList,index,documents[i],i,j,j+lmax);
        //invertList->clear();
    }
    size += beginIndex.size();

    invertList->clear();
  }
  gettimeofday(&send,0);

  gettimeofday(&tb,0);
  for(unsigned int i = 0; i < result.size(); i++)
  {
    //cout << result[i].ent_id << " " << result[i].doc_id << " " << result[i].pos << " " << result[i].len << endl;
    realResult++;
  }
  gettimeofday(&te, 0);
  tt = 0;
  tt += (te.tv_sec-tb.tv_sec+(te.tv_usec-tb.tv_usec)*1.0/CLOCKS_PER_SEC);
  cerr << "=====Search With Binary====="<<endl;
  cerr << "#       Algorithm : " << algorithm << endl;
  cerr << "#  Partition Time : " << pp << "s" << endl;
  cerr << "#     Filter Time : " << send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC - pp << "s" << endl;
  //cerr << "#     Verify Time : " << tt << "s" <<endl;
  cerr << "#      Total Time : " << tt+send.tv_sec-sbegin.tv_sec+(send.tv_usec-sbegin.tv_usec)*1.0/CLOCKS_PER_SEC<<"s"<< endl;
  //cerr << "#  To Be Verified : " << result.size() << endl;
  cerr << "# Real Result Num : " << result.size() << endl;
  cerr << "#    Average size : " << size / 100 << endl;
  cerr << "=====Search With Binary====="<<endl;

  result.clear();
  realResult = 0;

  delete index;
  delete invertList;

  GramListMap::iterator gbegin = gramListMap.begin();
  GramListMap::iterator gend = gramListMap.end();
  while(gbegin!=gend)
  {
    delete gbegin->second;
    gbegin++;
  }
  delete MINKEY;
  delete MAXKEY;
}
