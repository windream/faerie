#ifndef _HEAPBINARY_H_
#define _HEAPBINARY_H_
#include "Util.h"
#include <algorithm>
class HeapBinary
{
private:
	enum ApproximateType type;
	int Tl;
	int Le;
	int Te;
	vector<int>indexContainer;
public:
	HeapBinary(enum ApproximateType type);
	/********Edit Distance**********/
	void EntityExtract(InvertList* invertList, vector<int>* index, const string &query,const int id);
	void EntityExtractScale(InvertList* invertList, vector<int>* index, const string &query,const int id,const int begin,const int end);
	void EntityExtractWithNoPruning(InvertList* invertList, vector<int>* index, const string &query,const int id);
	void EntityExtractWithLazyUpdate(InvertList* invertList, vector<int>* index, const string &query,const int id);
	void EntityExtractWithBucketAndLazy(InvertList* invertList, vector<int>* index, const string &query,const int id);
	void EntityExtractWithMultHeap(InvertList* invertList, vector<int>* index, const string &query,const int id, int st, int len);
	/*******Edit Similarity********/
    void edsEntityExtract(InvertList* invertList, vector<int>* index, const string &query,const int id,const int begin,const int end);
    void edsEntityExtract(InvertList* invertList, vector<int>* index, const string &query,const int id);
    void edsEntityExtractVerify(InvertList* invertList, vector<int>* index, const string &query,const int id,const int begin,const int end);
    void edsEntityExtractVerify(InvertList* invertList, vector<int>* index, const string &query,const int id);
	void edsEntityExtractWithNoPruning(InvertList* invertList, vector<int>* index, const string &query,const int id);
	void edsEntityExtractWithLazyUpdate(InvertList* invertList, vector<int>* index, const string &query,const int id);
	void edsEntityExtractWithBucketAndLazy(InvertList* invertList, vector<int>* index, const string &query,const int id);
	void edsEntityExtractWithMultHeap(InvertList* invertList,vector<int>* index, const string &query, const int id, int st, int len);
	/************Jaccard***********/
    void jacEntityExtract(InvertList* invertList, vector<int>* index, const vector<string>*query,const int id,const int begin,const int end);
    void jacEntityExtract(InvertList* invertList, vector<int>* index, const string &query,const int id,const int begin,const int end);
    void jacEntityExtract(InvertList* invertList, vector<int>* index, const vector<string>*query,const int id);
	void jacEntityExtractWithNoPruning(InvertList* invertList, vector<int>* index, const vector<string>*query,const int id);
	void jacEntityExtractWithLazyUpdate(InvertList* invertList, vector<int>* index, const vector<string>*query,const int id);
	void jacEntityExtractWithBucketAndLazy(InvertList* invertList, vector<int>* index, const vector<string>*query,const int id);
	void jacEntityExtractWithMultHeap(InvertList* invertList,vector<int>* index, const vector<string>*query, const int id, int st, int len);
	/************Cosine ***********/
	void cosEntityExtract(InvertList* invertList, vector<int>* index, const vector<string>*query,const int id);
	void cosEntityExtract(InvertList* invertList, vector<int>* index, const string &query,const int id,const int begin,const int end);
	void cosEntityExtract(InvertList* invertList, vector<int>* index, const string &query,const int id);
	/************Dice ***********/
	void diceEntityExtract(InvertList* invertList, vector<int>* index, const vector<string>*query,const int id);
	void diceEntityExtract(InvertList* invertList, vector<int>* index, const string &query,const int id,const int begin,const int end);
	void diceEntityExtract(InvertList* invertList, vector<int>* index, const string &query,const int id);
	/*******Span*******/
	void BinarySpan(int &i, int lower, int upper, int id, int minValue,const string &query);
	void BinarySpanScale(int &i, int lower, int upper, int id, int minValue,const string &query);
	void edsBinarySpan(int &i, int lower, int upper, int id, int minValue,const string &query);
	void edsBinarySpanVerify(int &i, int lower, int upper, int id, int minValue,const string &query);
	void jacBinarySpan(int &i, int lower, int upper, int id, int minValue,const vector<string>*query);
	void jacBinarySpan(int &i, int lower, int upper, int id, int minValue,const string &query);
	void cosBinarySpan(int &i, int lower, int upper, int id, int minValue,const vector<string>*query);
	void cosBinarySpan(int &i, int lower, int upper, int id, int minValue,const string &query);
	void diceBinarySpan(int &i, int lower, int upper, int id, int minValue,const vector<string>*query);
	void diceBinarySpan(int &i, int lower, int upper, int id, int minValue,const string &query);
 	/*******Set/Get Arguments***/
	int GetEditDistanceT(int s, int e) { return  max(s, e) - tau * q; }
	int GetEditSimilarityT(int s, int e) { return int(ceil(max(s, e) - (max(s, e) + q - 1) * (1 - det) * q)); }
	int GetJaccardT(int s, int e) { return int(ceil((e + s) * det / (1 + det))); }
	int GetCosineT(int s, int e) { return int(ceil(sqrt(e * s) * det)); }
	int GetDiceT(int s, int e) { return int(ceil((s + e) * det / 2)); }
	void SetEditDistance(int e){
		Tl = e - tau * q;
		Te = e + tau;
		Le = e - tau;
	}
	void SetEditSimilarity(int e){
		Tl = int(ceil(e - (e + q - 1) * (1 - det) * q));
		Te = int(floor((e + q - 1) / det - (q - 1)));
		Le = int(ceil((e + q - 1) * det - (q - 1)));
	}
	void SetJaccard(int e){
		Tl = int(ceil(e * det));
		Te = int(floor(e / det));
		Le = int(ceil(e * det));
	}
	void SetCosine(int e){
		Tl = int(ceil(e * det * det));
		Te = int(floor(e / det / det));
		Le = int(ceil( e * det * det));
	}
	void SetDice(int e){
		Tl = int(ceil(e * det / (2 - det)));
		Te = int(floor(e * (2 - det) / det));
		Le = int(ceil(e * det / (2 - det)));
	}
};
#endif
