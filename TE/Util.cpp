#include "Util.h"

long long canda = 0;
long long candb = 0;
long long candc = 0;
long long candd = 0;
int maxLen = -1;
int minLen = MAX;
string entity;
string document;
string delim = " \t";
int q;
double det;
int tau;
GramListMap gramListMap;
vector<string> keywords;
vector<int> entityLengthMap;
unordered_set<Cand, CandHash> result2;
vector<Cand> result;
Array<int>* MAXKEY = new Array<int>(1);
Array<int>* MINKEY = new Array<int>(1);

extern string type;

void strToTokens(const string &s, vector<string> *res, const string &delims) {
        string::size_type begIdx, endIdx;
        begIdx = s.find_first_not_of(delims);
        while (begIdx != string::npos) {
                endIdx = s.find_first_of(delims, begIdx);
                if (endIdx == string::npos)
                    endIdx = s.length();
                res->push_back(s.substr(begIdx, endIdx - begIdx));
                begIdx = s.find_first_not_of(delims, endIdx);
        }
}
void strToTokens(const string &s, unordered_set<string> &res, const string &delims) {
        string::size_type begIdx, endIdx;
        begIdx = s.find_first_not_of(delims);
        while (begIdx != string::npos) {
                endIdx = s.find_first_of(delims, begIdx);
                if (endIdx == string::npos)
                        endIdx = s.length();
                res.insert(s.substr(begIdx, endIdx - begIdx));
                begIdx = s.find_first_not_of(delims, endIdx);
        }
}
void strToTokens(const string &s, set<string> &res, const string &delims) {
        string::size_type begIdx, endIdx;
        begIdx = s.find_first_not_of(delims);
        while (begIdx != string::npos) {
                endIdx = s.find_first_of(delims, begIdx);
                if (endIdx == string::npos)
                        endIdx = s.length();
                res.insert(s.substr(begIdx, endIdx - begIdx));
                begIdx = s.find_first_not_of(delims, endIdx);
        }
}
inline void appendEntityLengthMap(int length)
{
	if(maxLen < length)
		maxLen = length;
	if(minLen > length)
		minLen = length;
	entityLengthMap.push_back(length);
}

void appendGramList(string keyword, int id)
{
	appendEntityLengthMap(keyword.length()-q+1);
	int len = keyword.length() - q + 1;
	unordered_set<string> grams;
	for(int j=0; j<len; j++)
		grams.insert(keyword.substr(j, q));
	unordered_set<string>::const_iterator begin = grams.begin();
	unordered_set<string>::const_iterator end = grams.end();
	while(begin!=end)
	{
		GramListMap::iterator it = gramListMap.find(*begin);
		if(it==gramListMap.end())
		{
			Array<int>* newArray = new Array<int>();
			newArray->append(id);
			gramListMap.insert(GramListMap::value_type(*begin,newArray));
		}
		else
			it->second->append(id);
		++begin;
	}
}

void appendTokenList(string keyword, int id)
{
        unordered_set<string> grams;
        strToTokens(keyword, grams, delim);
        appendEntityLengthMap(grams.size());
        unordered_set<string>::const_iterator begin = grams.begin();
        unordered_set<string>::const_iterator end = grams.end();
        while(begin!=end)
        {
                GramListMap::iterator it = gramListMap.find(*begin);
                if(it==gramListMap.end())
                {
                        Array<int>* newArray = new Array<int>();
                        newArray->append(id);
                        gramListMap.insert(GramListMap::value_type(*begin,newArray));
                }
                else
                        it->second->append(id);
                ++begin;
        }
}

void addTail()
{
	GramListMap::iterator begin = gramListMap.begin();
	GramListMap::iterator end = gramListMap.end();
	int i = 0;
	while(begin!=end)
	{
		Array<int>::iterator b = (*begin).second->begin();
		Array<int>::iterator e = (*begin).second->end();
		unsigned int max = 0;

		int j = 0;
		while(b != e)
		{
		    if(keywords[*b].size() > max)
                max = keywords[*b].size();
            ++b;
            //cout << "i:" << i++ << " " << "j:" << j << endl;
		}
		if(type == "ED")
            max = max - q + 1 + tau;
        else if(type == "JS")
            max = (max - q + 1) / det;
        else if(type == "CS")
            max = (max - q + 1) / (det * det);
        else if(type == "DS")
            max = (max - q + 1) * (2 - det) / det;
        else if(type == "ES")
            max = max / det - q + 1;
		(*begin).second->sette(max);

//		cout << (*begin).second->gette() << endl;

		(*begin).second->append(MAX);
		++begin;
	}
}

int min(int a, int b, int c)
{
        if (a <= b && a <= c) return a;
        if (b <= a && b <= c) return b;
        return c;
}

bool edth_imp(const char* x,const char* y, int lx, int ly, int D) {
    int i, j;
    int* mat[2]; bool row, valid, ans;
    for (i = 0; i < 2; i++)
        mat[i] = new int[ly + 1];
    for (i = 0; i <= D; i++)
        mat[0][ly - i] = i;
    for (i = 1, row = 1; i <= lx; i++, row = !row) {
        valid = 0;
        if (i <= D) mat[row][ly] = i;
        for (j = (i - D >= 1 ? i - D : 1); j <= (i + D <= ly ? i + D : ly); j++) {
                if (x[lx - i] == y[ly - j])
                    mat[row][ly - j] = mat[!row][ly - j + 1];
                else
                    mat[row][ly - j] = min(j - 1 >= i - D ? mat[row][ly - j + 1] : D, mat[!row][ly - j + 1], j + 1 <= i + D ? mat[! row][ly - j] : D) + 1;
                if (mat[row][ly - j] <= D) valid = 1;
            }
            if (!valid) {
                for (i = 0; i < 2; i++)
                    delete [] mat[i];
                return false;
            }
    }
    ans = mat[!row][0] <= D;
    for (i = 0; i < 2; i++)
            delete [] mat[i];
    return ans;
}

bool eds_imp(const char* x,const char* y, int lx, int ly, double d) {
    int D = int(floor((1-d)*max(lx, ly)));
    int i, j;
    int* mat[2]; bool row, valid, ans;
    for (i = 0; i < 2; i++)
        mat[i] = new int[ly + 1];
    for (i = 0; i <= D; i++)
        mat[0][ly - i] = i;
    for (i = 1, row = 1; i <= lx; i++, row = !row) {
        valid = 0;
        if (i <= D) mat[row][ly] = i;
        for (j = (i - D >= 1 ? i - D : 1); j <= (i + D <= ly ? i + D : ly); j++) {
                if (x[lx - i] == y[ly - j])
                    mat[row][ly - j] = mat[!row][ly - j + 1];
                else
                    mat[row][ly - j] = min(j - 1 >= i - D ? mat[row][ly - j + 1] : D, mat[!row][ly - j + 1], j + 1 <= i + D ? mat[! row][ly - j] : D) + 1;
                if (mat[row][ly - j] <= D) valid = 1;
            }
            if (!valid) {
                for (i = 0; i < 2; i++)
                    delete [] mat[i];
                return false;
            }
    }
    ans = mat[!row][0] <= D;
    for (i = 0; i < 2; i++)
            delete [] mat[i];
    return ans;
}

double jaccardToken(const string &s1, const vector<string> *s2, int begin, int len)
{
    set<string> tokenSet1;
    strToTokens(s1, tokenSet1, delim);
    unsigned news1Size = tokenSet1.size();
    unsigned news2Size = len;
    tokenSet1.insert(s2->begin()+begin,s2->begin()+begin+len);
    unsigned unionSize = tokenSet1.size();
    unsigned intersectionSize = news1Size+news2Size-unionSize;
    return intersectionSize*1.0/unionSize;
}

double cosineToken(const string &s1, const vector<string> *s2, int begin, int len)
{
    set<string> tokenSet1;
    strToTokens(s1, tokenSet1, delim);
    unsigned news1Size = tokenSet1.size();
    unsigned news2Size = len;
    tokenSet1.insert(s2->begin()+begin,s2->begin()+begin+len);
    unsigned unionSize = tokenSet1.size();
    unsigned intersectionSize = news1Size+news2Size-unionSize;
    return intersectionSize*1.0/sqrt(news1Size*news2Size*1.0);
}

double diceToken(const string &s1, const vector<string> *s2, int begin, int len)
{
    set<string> tokenSet1;
    strToTokens(s1, tokenSet1, delim);
    unsigned news1Size = tokenSet1.size();
    unsigned news2Size = len;
    tokenSet1.insert(s2->begin()+begin,s2->begin()+begin+len);
    unsigned unionSize = tokenSet1.size();
    unsigned intersectionSize = news1Size+news2Size-unionSize;
    return intersectionSize*2.0/(news1Size + news2Size);
}
