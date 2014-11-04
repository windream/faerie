#ifndef _UTIL_H_
#define _UTIL_H_
#include <sys/time.h>
#include <ctime>
#include <string>
#include <set>
#ifndef _WIN32
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#else
#include <unordered_map>
#include <unordered_set>
#endif
#include "array.h"

#define MAX 0x7FFFFFFF             //inf big
#define uint unsigned

#define log2 0.301029995663981198017467022509663365781307220458984375
#define pow_size 100000


enum ApproximateType{
	EditDistance,
    EditDistanceSimilarity,
	Jaccard,
    Cosine,
    Dice
};
typedef unordered_map<string,Array<int>*> GramListMap;
typedef vector< pair<int,Array<int>*> > InvertList;
typedef struct cand{
	int ent_id;
	int doc_id;
	int pos;
	int len;

	bool operator==(const cand& c2) const
    {
        return((this->ent_id==c2.ent_id)&&(this->doc_id==c2.doc_id)&&(this->pos==c2.pos)&&(this->len==c2.len));
    }
}Cand;

class CandHash {
    public:
        std::size_t operator()(const Cand &c) const
        {
            return c.ent_id * pow(pow_size,3) + c.doc_id * pow(pow_size,2) + c.pos * pow_size + c.len;
        }
};

extern int maxLen;
extern int minLen;
extern long long canda;
extern long long candb;
extern long long candc;
extern long long candd;
extern string entity;
extern string document;
extern string delim;
extern int q;
extern double det;
extern int tau;
extern GramListMap gramListMap;
extern vector<string> keywords;
extern vector<int> entityLengthMap;
extern vector<Cand> result;
extern unordered_set<Cand, CandHash> result2;
extern Array<int>* MAXKEY;
extern Array<int>* MINKEY;




void strToTokens(const string &s, vector<string> *res, const string &delims);
void strToTokens(const string &s, unordered_set<string> &res, const string &delims);
void strToTokens(const string &s, set<string> &res, const string &delims);
inline void appendEntityLengthMap(int length);
void appendGramList(string keyword, int id);
void appendTokenList(string keyword, int id);
void addTail();
bool edth_imp(const char* x,const char* y, int lx, int ly, int D);
bool eds_imp(const char* x,const char* y, int lx, int ly, double d);
double jaccardToken(const string &s1, const vector<string> *s2, int begin, int len);
double cosineToken(const string &s1, const vector<string> *s2, int begin, int len);
double diceToken(const string &s1, const vector<string> *s2, int begin, int len);
#endif
