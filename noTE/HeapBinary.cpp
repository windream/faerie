#include "HeapBinary.h"
#include <memory.h>

extern int lastEndIndex;
extern int thisEndIndex;
extern bool isFirstTime;
HeapBinary::HeapBinary(enum ApproximateType approximateType)
{
	type = approximateType;
}

inline void adjustStableHeap(InvertList* invertList, vector<int>* index, vector<int>* hp, int pos, int n)
{
	int j;
	int temp = hp->at(pos);
	j = 2 * pos;
//	cout << "pos: " << pos << endl;
//	cout << "j: " << j << endl;
//	cout << "temp: " << temp << endl;
	while(j <= n)
	{
//	  cout << "inner j:" << j << endl;
		if(j < n)
		{
			if( invertList->at(hp->at(j)).second->at(index->at(hp->at(j)))
			  > invertList->at(hp->at(j+1)).second->at(index->at(hp->at(j+1)))
			 || ( invertList->at(hp->at(j)).second->at(index->at(hp->at(j)))
			 == invertList->at(hp->at(j+1)).second->at(index->at(hp->at(j+1)))
			 && hp->at(j)>hp->at(j+1) ) )
				j++;
		}

		if( invertList->at(temp).second->at(index->at(temp))
		  < invertList->at(hp->at(j)).second->at(index->at(hp->at(j)))
		 || ( invertList->at(temp).second->at(index->at(temp))
  		 == invertList->at(hp->at(j)).second->at(index->at(hp->at(j)))
		 && temp < hp->at(j) ) )
			break;

		hp->at(j / 2) = hp->at(j);
		j = 2 * j;
	}
	hp->at(j / 2) = temp;
}

inline void createStableHeap(InvertList* invertList,vector<int>* index,vector<int>* hp,const int begin,const int end)
{
	int i;
  hp->push_back(-1);
	for(i = begin; i <= end; i++)
	{
	  if(invertList->at(i).first != -1)
	  {
	    hp->push_back(i);
	  }
	}
  int n = hp->size() - 1;
  //index->assign(n - 1,0);
//  cout << "inner  hp.size:" << hp->size() << endl;
	for(i = n/2; i >= 1; i--)
		adjustStableHeap(invertList,index,hp,i,n);
}

inline void createStableHeap(InvertList* invertList,vector<int>* index,vector<int>* hp)
{
	int i;
  int n = hp->size();
	for(i = 0; i <= n; i++)
	{
	  hp->push_back(i - 1);
	}

  //index->assign(n - 1,0);
	for(i = n/2; i >= 1; i--)
		adjustStableHeap(invertList,index,hp,i,n);
}

void HeapBinary::BinarySpan(int &i, int lower, int upper, int id, int minValue,const string &query)
{
	//cout << "Span: " << i << endl;
	int mid = -1;
	int X = Te + indexContainer[i] - 1;
	while(lower <= upper)
	{
		mid = (lower + upper + 1) / 2;
		if (indexContainer[mid] > X)
			upper = mid - 1;
		else
			lower = mid + 1;
	}
	mid = upper;

	for(int j = i + Tl - 1; j <= mid; j++)
	{
		for(int len = Le; len <= Te; len++)
		{
			if(j - i + 1 < len - tau*q)
				break;

			int left = max(indexContainer[j] - len + 1,
				       i == 0 ? 0 : indexContainer[i - 1] + 1 );

			int right = min(indexContainer[i] + len - 1,
					j == (int)indexContainer.size() - 1 ?
					int(query.length()) - q : indexContainer[j + 1] - 1 );

			for(int cur = left; cur + len - 1 <= right; cur++)
			{
				canda++;
				int c = 0;
				int suffix,prefix;
				//string str =  query.substr(cur ,len + q - 1);
				for(prefix = 0; prefix < q - 1; prefix++)
				{
					if(keywords[minValue][prefix] == query[cur + prefix])
						c++;
					else
						break;
				}
				for(suffix = 1; suffix < q; suffix++)
				{
					if(keywords[minValue][keywords[minValue].length()-suffix]==query[cur+len+q-1-suffix])
					    c++;
					else
					    break;
				}
				if( j - i + 1 + c >= GetEditDistanceT(len + 2 * q - 2, entityLengthMap[minValue] + 2 * q - 2))
				{
					 Cand res1 = {minValue,id, cur ,len + q - 1};
					 result.push_back(res1);
				}
			}
		}
	}

	//cout << "End Span: " << i+1 << endl;
	for(i = i + 1; i <= mid - Tl + 1;)
		BinarySpan(i, mid, min(i + Te - 1,int(indexContainer.size()) - 1), id, minValue, query);
}

void HeapBinary::BinarySpanScale(int &i, int lower, int upper, int id, int minValue,const string &query)
{
	//cout << "Span: " << i << endl;
	int mid = -1;
	int X = Te + indexContainer[i] - 1;
	while(lower <= upper)
	{
		mid = (lower + upper + 1) / 2;
		if (indexContainer[mid] > X)
			upper = mid - 1;
		else
			lower = mid + 1;
	}
	mid = upper;

	for(int j = i + Tl - 1; j <= mid; j++)
	{
		for(int len = Le; len <= Te; len++)
		{
			if(j - i + 1 < len - tau*q)
				break;

			int left = max(indexContainer[j] - len + 1,
				       i == 0 ? 0 : indexContainer[i - 1] + 1 );

			int right = min(indexContainer[i] + len - 1,
					j == (int)indexContainer.size() - 1 ?
					int(query.length()) - q : indexContainer[j + 1] - 1 );

			for(int cur = left; cur + len - 1 <= right; cur++)
			{
				int c = 0;
				int suffix,prefix;
				//string str =  query.substr(cur ,len + q - 1);
				for(prefix = 0; prefix < q - 1; prefix++)
				{
					if(keywords[minValue][prefix] == query[cur + prefix])
						c++;
					else
						break;
				}
				for(suffix = 1; suffix < q; suffix++)
				{
					if(keywords[minValue][keywords[minValue].length()-suffix]==query[cur+len+q-1-suffix])
					    c++;
					else
					    break;
				}
				if( j - i + 1 + c >= GetEditDistanceT(len + 2 * q - 2, entityLengthMap[minValue] + 2 * q - 2))
				{
					 Cand res1 = {minValue,id, cur ,len + q - 1};
					 //cout << cur << " " << len+q-1 << " " << lastEndIndex << " " << thisEndIndex << endl;
//					 if(isFirstTime || (!isFirstTime && cur+len+q-1>lastEndIndex))
					 if(cur+len+q-1>lastEndIndex&&cur+len+q-1<=thisEndIndex)
                        result.push_back(res1);
                    //cout << "lastEndIndex:"  << lastEndIndex << endl;
//                     result2.insert(res1);
				}
			}
		}
	}

	//cout << "End Span: " << i+1 << endl;
	for(i = i + 1; i <= mid - Tl + 1;)
		BinarySpanScale(i, mid, min(i + Te - 1,int(indexContainer.size()) - 1), id, minValue, query);
}

void HeapBinary::edsBinarySpan(int &i, int lower, int upper, int id, int minValue,const string &query)
{
	//cout << "Span: " << i << endl;
	int mid = -1;
	int X = Te + indexContainer[i] - 1;
	while(lower <= upper)
	{
		mid = (lower + upper + 1) / 2;
		if (indexContainer[mid] > X)
			upper = mid - 1;
		else
			lower = mid + 1;
	}
	mid = upper;

	for(int j = i + Tl - 1; j <= mid; j++)
	{
		for(int len = Le; len <= Te; len++)
		{
			if(j - i + 1 < int(ceil(len - (len + q - 1) * (1 - det) * q)))
				break;

			int left = max(indexContainer[j] - len + 1,
				       i == 0 ? 0 : indexContainer[i - 1] + 1 );

			int right = min(indexContainer[i] + len - 1,
					j == (int)indexContainer.size() - 1 ?
					int(query.length()) - q : indexContainer[j + 1] - 1 );

			for(int cur = left; cur + len - 1 <= right; cur++)
			{
				canda++;
				int c = 0;
				int suffix,prefix;
				//string str =  query.substr(cur ,len + q - 1);
				for(prefix = 0; prefix < q - 1; prefix++)
				{
					if(keywords[minValue][prefix] == query[cur + prefix])
						c++;
					else
						break;
				}
				for(suffix = 1; suffix < q; suffix++)
				{
					if(keywords[minValue][keywords[minValue].length()-suffix]==query[cur+len+q-1-suffix])
					    c++;
					else
					    break;
				}
				if( j - i + 1 + c >= GetEditSimilarityT(len + 2 * q - 2, entityLengthMap[minValue] + 2 * q - 2))
				{
					 Cand res1 = {minValue,id, cur ,len + q - 1};
					 if(cur+len+q-1>lastEndIndex&&cur+len+q-1<=thisEndIndex)
                        result.push_back(res1);
				}
			}
		}
	}

	//cout << "End Span: " << i+1 << endl;
	for(i = i + 1; i <= mid - Tl + 1;)
		edsBinarySpan(i, mid, min(i + Te - 1,int(indexContainer.size()) - 1), id, minValue, query);
}

void HeapBinary::edsBinarySpanVerify(int &i, int lower, int upper, int id, int minValue,const string &query)
{
	//cout << "Span: " << i << endl;
	int mid = -1;
	int X = Te + indexContainer[i] - 1;
	while(lower <= upper)
	{
		mid = (lower + upper + 1) / 2;
		if (indexContainer[mid] > X)
			upper = mid - 1;
		else
			lower = mid + 1;
	}
	mid = upper;

	for(int j = i + Tl - 1; j <= mid; j++)
	{
		for(int len = Le; len <= Te; len++)
		{
			if(j - i + 1 < int(ceil(len - (len + q - 1) * (1 - det) * q)))
				break;

			int left = max(indexContainer[j] - len + 1,
				       i == 0 ? 0 : indexContainer[i - 1] + 1 );

			int right = min(indexContainer[i] + len - 1,
					j == (int)indexContainer.size() - 1 ?
					int(query.length()) - q : indexContainer[j + 1] - 1 );

			for(int cur = left; cur + len - 1 <= right; cur++)
			{
				canda++;
				int c = 0;
				int suffix,prefix;
				//string str =  query.substr(cur ,len + q - 1);
				for(prefix = 0; prefix < q - 1; prefix++)
				{
					if(keywords[minValue][prefix] == query[cur + prefix])
						c++;
					else
						break;
				}
				for(suffix = 1; suffix < q; suffix++)
				{
					if(keywords[minValue][keywords[minValue].length()-suffix]==query[cur+len+q-1-suffix])
					    c++;
					else
					    break;
				}
				if( j - i + 1 + c >= GetEditSimilarityT(len + 2 * q - 2, entityLengthMap[minValue] + 2 * q - 2))
				{
                    canda++;
                    if(eds_imp(keywords[minValue].c_str(), query.substr(cur, len+q-1).c_str(), entityLengthMap[minValue]+q-1, len+q-1, det))
                    {
                         Cand res1 = {minValue,id, cur ,len + q - 1};
                         if(cur+len+q-1>lastEndIndex&&cur+len+q-1<=thisEndIndex)
                            result.push_back(res1);
                    }
				}
			}
		}
	}

	//cout << "End Span: " << i+1 << endl;
	for(i = i + 1; i <= mid - Tl + 1;)
		edsBinarySpanVerify(i, mid, min(i + Te - 1,int(indexContainer.size()) - 1), id, minValue, query);
}

void HeapBinary::jacBinarySpan(int &i, int lower, int upper, int id, int minValue,const vector<string>*query)
{
	//cout << "Span: " << i << endl;
	int mid = -1;
	int X = Te + indexContainer[i] - 1;
	while(lower <= upper)
	{
		mid = (lower + upper + 1) / 2;
		if (indexContainer[mid] > X)
			upper = mid - 1;
		else
			lower = mid + 1;
	}
	mid = upper;

	for(int j = i + Tl - 1; j <= mid; j++)
	{
		for(int len = Le; len <= Te; len++)
		{
			if(j - i + 1 < GetJaccardT(len, entityLengthMap[minValue]))
				break;

			int left = max(indexContainer[j] - len + 1,
				       i == 0 ? 0 : indexContainer[i - 1] + 1 );

			int right = min(indexContainer[i] + len - 1,
					j == (int)indexContainer.size() - 1 ?
					int(query->size()) - 1 : indexContainer[j + 1] - 1 );

			for(int cur = left; cur + len - 1 <= right; cur++)
			{
				canda++;
                Cand res1 = {minValue,id, cur ,len};
                if(cur+len>lastEndIndex&&cur+len<=thisEndIndex)
                    result.push_back(res1);
			}
		}
	}

	//cout << "End Span: " << i+1 << endl;
	for(i = i + 1; i <= mid - Tl + 1;)
		jacBinarySpan(i, mid, min(i + Te - 1,int(indexContainer.size()) - 1), id, minValue, query);
}

void HeapBinary::jacBinarySpan(int &i, int lower, int upper, int id, int minValue,const string &query)
{
	cout << "Span: " << i << endl;
	int mid = -1;
	int X = Te + indexContainer[i] - 1;
	while(lower <= upper)
	{
		mid = (lower + upper + 1) / 2;
		if (indexContainer[mid] > X)
			upper = mid - 1;
		else
			lower = mid + 1;
	}
	mid = upper;

	for(int j = i + Tl - 1; j <= mid; j++)
	{
		for(int len = Le; len <= Te; len++)
		{
			if(j - i + 1 < GetJaccardT(len, entityLengthMap[minValue]))
				break;

			int left = max(indexContainer[j] - len + 1,
				       i == 0 ? 0 : indexContainer[i - 1] + 1 );

			int right = min(indexContainer[i] + len - 1,
					j == (int)indexContainer.size() - 1 ?
					(int)query.length() - q: indexContainer[j + 1] - 1 );

			for(int cur = left; cur + len - 1 <= right; cur++)
			{
				canda++;
                Cand res1 = {minValue,id, cur ,len};
                if(cur+len>lastEndIndex&&cur+len<=thisEndIndex)
                    result.push_back(res1);
			}
		}
	}

	//cout << "End Span: " << i+1 << endl;
	for(i = i + 1; i <= mid - Tl + 1;)
		jacBinarySpan(i, mid, min(i + Te - 1,int(indexContainer.size()) - 1), id, minValue, query);
}

void HeapBinary::cosBinarySpan(int &i, int lower, int upper, int id, int minValue,const vector<string>*query)
{
	//cout << "Span: " << i << endl;
	int mid = -1;
	int X = Te + indexContainer[i] - 1;
	while(lower <= upper)
	{
		mid = (lower + upper + 1) / 2;
		if (indexContainer[mid] > X)
			upper = mid - 1;
		else
			lower = mid + 1;
	}
	mid = upper;

	for(int j = i + Tl - 1; j <= mid; j++)
	{
		for(int len = Le; len <= Te; len++)
		{
			if(j - i + 1 < GetCosineT(len, entityLengthMap[minValue]))
				break;

			int left = max(indexContainer[j] - len + 1,
				       i == 0 ? 0 : indexContainer[i - 1] + 1 );

			int right = min(indexContainer[i] + len - 1,
					j == (int)indexContainer.size() - 1 ?
					int(query->size()) - 1 : indexContainer[j + 1] - 1 );

			for(int cur = left; cur + len - 1 <= right; cur++)
			{
                 Cand res1 = {minValue,id, cur ,len};
                 result.push_back(res1);
			}
		}
	}

	//cout << "End Span: " << i+1 << endl;
	for(i = i + 1; i <= mid - Tl + 1;)
		cosBinarySpan(i, mid, min(i + Te - 1,int(indexContainer.size()) - 1), id, minValue, query);
}

void HeapBinary::cosBinarySpan(int &i, int lower, int upper, int id, int minValue,const string &query)
{
	//cout << "Span: " << i << endl;
	int mid = -1;
	int X = Te + indexContainer[i] - 1;
	while(lower <= upper)
	{
		mid = (lower + upper + 1) / 2;
		if (indexContainer[mid] > X)
			upper = mid - 1;
		else
			lower = mid + 1;
	}
	mid = upper;

	for(int j = i + Tl - 1; j <= mid; j++)
	{
		for(int len = Le; len <= Te; len++)
		{
			if(j - i + 1 < GetCosineT(len, entityLengthMap[minValue]))
				break;

			int left = max(indexContainer[j] - len + 1,
				       i == 0 ? 0 : indexContainer[i - 1] + 1 );

			int right = min(indexContainer[i] + len - 1,
					j == (int)indexContainer.size() - 1 ?
					(int)query.length() - q: indexContainer[j + 1] - 1 );

			for(int cur = left; cur + len - 1 <= right; cur++)
			{
                 Cand res1 = {minValue,id, cur ,len + q - 1};
                 if(cur+len+q-1>lastEndIndex&&cur+len+q-1<=thisEndIndex)
                    result.push_back(res1);
			}
		}
	}

	//cout << "End Span: " << i+1 << endl;
	for(i = i + 1; i <= mid - Tl + 1;)
		cosBinarySpan(i, mid, min(i + Te - 1,int(indexContainer.size()) - 1), id, minValue, query);
}

void HeapBinary::diceBinarySpan(int &i, int lower, int upper, int id, int minValue,const vector<string>*query)
{
	//cout << "Span: " << i << endl;
	int mid = -1;
	int X = Te + indexContainer[i] - 1;
	while(lower <= upper)
	{
		mid = (lower + upper + 1) / 2;
		if (indexContainer[mid] > X)
			upper = mid - 1;
		else
			lower = mid + 1;
	}
	mid = upper;

	for(int j = i + Tl - 1; j <= mid; j++)
	{
		for(int len = Le; len <= Te; len++)
		{
			if(j - i + 1 < GetDiceT(len, entityLengthMap[minValue]))
				break;

			int left = max(indexContainer[j] - len + 1,
				       i == 0 ? 0 : indexContainer[i - 1] + 1 );

			int right = min(indexContainer[i] + len - 1,
					j == (int)indexContainer.size() - 1 ?
					int(query->size()) - 1 : indexContainer[j + 1] - 1 );

			for(int cur = left; cur + len - 1 <= right; cur++)
			{
                 Cand res1 = {minValue,id, cur ,len};
                 result.push_back(res1);
			}
		}
	}

	//cout << "End Span: " << i+1 << endl;
	for(i = i + 1; i <= mid - Tl + 1;)
		diceBinarySpan(i, mid, min(i + Te - 1,int(indexContainer.size()) - 1), id, minValue, query);
}

void HeapBinary::diceBinarySpan(int &i, int lower, int upper, int id, int minValue,const string &query)
{
	//cout << "Span: " << i << endl;
	int mid = -1;
	int X = Te + indexContainer[i] - 1;
	while(lower <= upper)
	{
		mid = (lower + upper + 1) / 2;
		if (indexContainer[mid] > X)
			upper = mid - 1;
		else
			lower = mid + 1;
	}
	mid = upper;

	for(int j = i + Tl - 1; j <= mid; j++)
	{
		for(int len = Le; len <= Te; len++)
		{
			if(j - i + 1 < GetDiceT(len, entityLengthMap[minValue]))
				break;

			int left = max(indexContainer[j] - len + 1,
				       i == 0 ? 0 : indexContainer[i - 1] + 1 );

			int right = min(indexContainer[i] + len - 1,
					j == (int)indexContainer.size() - 1 ?
					int(query.length()) - q : indexContainer[j + 1] - 1 );

			for(int cur = left; cur + len - 1 <= right; cur++)
			{
                 Cand res1 = {minValue,id, cur ,len + q - 1};
                 if(cur+len+q-1>lastEndIndex&&cur+len+q-1<=thisEndIndex)
                    result.push_back(res1);
			}
		}
	}

	//cout << "End Span: " << i+1 << endl;
	for(i = i + 1; i <= mid - Tl + 1;)
		diceBinarySpan(i, mid, min(i + Te - 1,int(indexContainer.size()) - 1), id, minValue, query);
}

/*******Edit Distance*********/

void HeapBinary::EntityExtract(InvertList* invertList, vector<int>* index, const string &query, const int id)
{
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

  cout << "MAX: " << MAX << endl;
  cout << "Tl: " << Tl << endl;
	while(minValue != MAX)
	{
		do{
		  cout << "minValue: " << minValue << endl;
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetEditDistance(length);

		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			cout << "indexContainer.size:" << indexContainer.size() << endl;
			int j, i = 0;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
					BinarySpan(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::EntityExtractScale(InvertList* invertList, vector<int>* index, const string &query, const int id, const int begin, const int end)
{
	indexContainer.clear();
	//int K = index->size();

	vector<int>* hp = new vector<int>();
	//hp->reserve(K+1);
	createStableHeap(invertList,index,hp,begin,end);
//	cout << "hp.size:" << hp->size() << endl;
  //index->assign(index->size(),0);
  int K = hp->size() - 1;
    if(K==0)
        return;
	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

//  cout << "MAX: " << MAX << endl;
//  cout << "Tl: " << Tl << endl;
  int i = 0;
	while(minValue != MAX)
	{
//	  cout << "第" << ++i << "次" << endl;
		do{
		  //cout << "minValue: " << minValue << endl;
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetEditDistance(length);

//    cout << "indexContainer.size()" << indexContainer.size() << endl;
		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			int j, i = 0;
//			cout << "indexContainer.size() - Tl + 1" << (int)indexContainer.size() - Tl + 1 << endl;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
					BinarySpanScale(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::EntityExtractWithNoPruning(InvertList* invertList, vector<int>* index, const string &query, const int id)
{
    vector<int> count(query.length() - q + 1, 0);
    vector<int> posRec;
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetEditDistance(length);

	    for(int offset = Le; offset <= Te; offset++) {
            for(unsigned lp = 0; lp < indexContainer.size(); lp++){
                for(int i = indexContainer[lp] - offset + 1 > 0 ? indexContainer[lp] - offset + 1 : 0;i <= (indexContainer[lp] < (int)query.length() - q + 1 - offset ? indexContainer[lp] : (int)query.length() - q + 1 - offset); i++)
                {
                    count[i]++;
                    if(count[i]==1)
                        posRec.push_back(i);
                }
            }
            for(int i = 0; i < (int)posRec.size(); i++) {
                if(count[posRec[i]] >= GetEditDistanceT(offset, length)) {
                    int c = 0;
                    int prefix, suffix;
                    for(prefix = 0; prefix < q - 1; prefix++) {
                        if(keywords[minValue][prefix] == query[posRec[i] + prefix])
                            c++;
                        else
                            break;
                    }
                    for(suffix = 1; suffix < q; suffix++) {
                        if(keywords[minValue][keywords[minValue].length() - suffix] == query[posRec[i] + offset + q - 1 - suffix])
                            c++;
                        else
                            break;
                    }
                    if(c + count[posRec[i]] >= GetEditDistanceT(offset + 2 * q - 2, length + 2 * q - 2)) {
                        Cand res = {minValue, id, posRec[i], offset + q - 1};
                        result.push_back(res);
                    }
                }
                candd++;
                count[posRec[i]] = 0;
            }
            posRec.clear();
        }
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::EntityExtractWithLazyUpdate(InvertList* invertList, vector<int>* index, const string &query, const int id)
{
    vector<int> count(query.length() - q + 1, 0);
    vector<int> posRec;
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetEditDistance(length);

		if((int)indexContainer.size() >= Tl) {
            for(int offset = Le; offset <= Te; offset++) {
                for(unsigned lp = 0; lp < indexContainer.size(); lp++) {
                    for(int i = indexContainer[lp] - offset + 1 > 0 ? indexContainer[lp] - offset + 1 : 0; i <= (indexContainer[lp] < (int)query.length() - q + 1 - offset ? indexContainer[lp] : (int)query.length() - q + 1 - offset); i++)
                    {
                        count[i]++;
                        if(count[i]==1)
                            posRec.push_back(i);
                    }
                }
                for(int i = 0; i < (int)posRec.size(); i++) {
                    if(count[posRec[i]] >= GetEditDistanceT(offset, length)) {
                        int c = 0;
                        int prefix, suffix;
                        for(prefix = 0; prefix < q - 1; prefix++) {
                            if(keywords[minValue][prefix] == query[posRec[i] + prefix])
                                c++;
                            else
                                break;
                        }
                        for(suffix = 1; suffix < q; suffix++) {
                            if(keywords[minValue][keywords[minValue].length() - suffix] == query[posRec[i] + offset + q - 1 - suffix])
                                c++;
                            else
                                break;
                        }
                        if(c + count[posRec[i]] >= GetEditDistanceT(offset + 2 * q - 2, length + 2 * q - 2)) {
                            Cand res = {minValue, id, posRec[i], offset + q - 1};
                            result.push_back(res);
                        }
                    }
                    candc++;
                    count[posRec[i]] = 0;
                }
                posRec.clear();
            }
        }
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::EntityExtractWithBucketAndLazy(InvertList* invertList, vector<int>* index, const string &query, const int id)
{
    vector<int> count(query.length() - q + 1, 0);
    vector<int> posRec;
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));
    int nowIndex;
	while(minValue != MAX)
	{
		do{
            nowIndex = minIndex;
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)) && minIndex - nowIndex - 1 <= q * tau);

		int length = entityLengthMap[minValue];
		SetEditDistance(length);

		if((int)indexContainer.size() >= Tl) {
            for(int offset = Le; offset <= Te; offset++) {
                for(unsigned lp = 0; lp < indexContainer.size(); lp++) {
                    for(int i = indexContainer[lp] - offset + 1 > 0 ? indexContainer[lp] - offset + 1 : 0; i <= (indexContainer[lp] < (int)query.length() - q + 1 - offset ? indexContainer[lp] : (int)query.length() - q + 1 - offset); i++)
                    {
                        count[i]++;
                        if(count[i]==1)
                            posRec.push_back(i);
                    }
                }
                for(int i = 0; i < (int)posRec.size(); i++) {
                    if(count[posRec[i]] >= GetEditDistanceT(offset, length)) {
                        int c = 0;
                        int prefix, suffix;
                        for(prefix = 0; prefix < q - 1; prefix++) {
                            if(keywords[minValue][prefix] == query[posRec[i] + prefix])
                                c++;
                            else
                                break;
                        }
                        for(suffix = 1; suffix < q; suffix++) {
                            if(keywords[minValue][keywords[minValue].length() - suffix] == query[posRec[i] + offset + q - 1 - suffix])
                                c++;
                            else
                                break;
                        }
                        if(c + count[posRec[i]] >= GetEditDistanceT(offset + 2 * q - 2, length + 2 * q - 2)) {
                            Cand res = {minValue, id, posRec[i], offset + q - 1};
                            result.push_back(res);
                        }
                    }
                    candb++;
                    count[posRec[i]] = 0;
                }
                posRec.clear();
            }
        }
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::EntityExtractWithMultHeap(InvertList* invertList, vector<int>* index, const string &query, const int id, int st, int len)
{
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetEditDistance(length);

        if(len - q + 1 >= Le && len - q + 1 <= Te)
        {
            if((int)indexContainer.size() >= GetEditDistanceT(len - q + 1, length)) {
                int c = 0;
                int prefix, suffix;
                for(prefix = 0; prefix < q - 1; prefix++) {
                    if(keywords[minValue][prefix] == query[st + prefix])
                        c++;
                    else
                        break;
                }
                for(suffix = 1; suffix < q; suffix++) {
                    if(keywords[minValue][keywords[minValue].length() - suffix] == query[st + len - suffix])
                        c++;
                    else
                        break;
                }
                if(c + (int)indexContainer.size() >= GetEditDistanceT(len + q - 1, length + 2 * q - 2)) {
                    Cand res = {minValue, id, st, len};
                    result.push_back(res);
                }
            }
        }
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

/*******Edit Distance Similarity*********/

void HeapBinary::edsEntityExtract(InvertList* invertList, vector<int>* index, const string &query, const int id, const int begin, const int end)
{
	indexContainer.clear();
	vector<int>* hp = new vector<int>();
	createStableHeap(invertList,index,hp,begin,end);

	int K = hp->size() - 1;
//    cout << "hp.size: " << hp->size() << endl;
    if(K == 0)
        return;

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetEditSimilarity(length);

		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			int j, i = 0;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
					edsBinarySpan(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::edsEntityExtract(InvertList* invertList, vector<int>* index, const string &query, const int id)
{
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetEditSimilarity(length);

		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			int j, i = 0;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
					edsBinarySpan(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::edsEntityExtractVerify(InvertList* invertList, vector<int>* index, const string &query, const int id)
{
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetEditSimilarity(length);

		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			int j, i = 0;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
					edsBinarySpanVerify(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::edsEntityExtractVerify(InvertList* invertList, vector<int>* index, const string &query, const int id, const int begin, const int end)
{
	indexContainer.clear();

	vector<int>* hp = new vector<int>();
	createStableHeap(invertList,index,hp,begin,end);
	int K = hp->size() - 1;

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetEditSimilarity(length);

		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			int j, i = 0;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
					edsBinarySpanVerify(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::edsEntityExtractWithNoPruning(InvertList* invertList, vector<int>* index, const string &query, const int id)
{
    vector<int> posRec;
    vector<int> count(query.length() - q + 1, 0);
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetEditSimilarity(length);

	    for(int offset = Le; offset <= Te; offset++) {
            for(unsigned lp = 0; lp < indexContainer.size(); lp++) {
                for(int i = indexContainer[lp] - offset + 1 > 0 ? indexContainer[lp] - offset + 1 : 0; i <= (indexContainer[lp] < (int)query.length() - q + 1 - offset ? indexContainer[lp] : (int)query.length() - q + 1 - offset); i++)
                {
                    count[i]++;
                    if(count[i] == 1)
                        posRec.push_back(i);
                }
            }
            for(int i = 0; i < (int)posRec.size(); i++) {
                if(count[posRec[i]] >= GetEditSimilarityT(offset, length)) {
                    int c = 0;
                    int prefix, suffix;
                    for(prefix = 0; prefix < q - 1; prefix++) {
                        if(keywords[minValue][prefix] == query[posRec[i] + prefix])
                            c++;
                        else
                            break;
                    }
                    for(suffix = 1; suffix < q; suffix++) {
                        if(keywords[minValue][keywords[minValue].length() - suffix] == query[posRec[i] + offset + q - 1 - suffix])
                            c++;
                        else
                            break;
                    }
                    if(c + count[posRec[i]] >= GetEditSimilarityT(offset + 2 * q - 2, length + 2 * q - 2)) {
                        Cand res = {minValue, id, posRec[i], offset + q - 1};
                        result.push_back(res);
                    }
                }
                candd++;
                count[posRec[i]] = 0;
            }
            posRec.clear();
        }
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::edsEntityExtractWithLazyUpdate(InvertList* invertList, vector<int>* index, const string &query, const int id)
{
    vector<int> count(query.length() - q + 1, 0);
    vector<int> posRec;
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetEditSimilarity(length);

		if((int)indexContainer.size() >= Tl) {
            for(int offset = Le; offset <= Te; offset++) {
                for(unsigned lp = 0; lp < indexContainer.size(); lp++) {
                    for(int i = indexContainer[lp] - offset + 1 > 0 ? indexContainer[lp] - offset + 1 : 0; i <= (indexContainer[lp] < (int)query.length() - q + 1 - offset ? indexContainer[lp] : (int)query.length() - q + 1 - offset); i++)
                    {
                        count[i]++;
                        if(count[i] == 1)
                            posRec.push_back(i);
                    }
                }
                for(int i = 0; i < (int)posRec.size(); i++) {
                    if(count[posRec[i]] >= GetEditSimilarityT(offset, length)) {
                        int c = 0;
                        int prefix, suffix;
                        for(prefix = 0; prefix < q - 1; prefix++) {
                            if(keywords[minValue][prefix] == query[posRec[i] + prefix])
                                c++;
                            else
                                break;
                        }
                        for(suffix = 1; suffix < q; suffix++) {
                            if(keywords[minValue][keywords[minValue].length() - suffix] == query[posRec[i] + offset + q - 1 - suffix])
                                c++;
                            else
                                break;
                        }
                        if(c + count[posRec[i]] >= GetEditSimilarityT(offset + 2 * q - 2, length + 2 * q - 2)) {
                            Cand res = {minValue, id, posRec[i], offset + q - 1};
                            result.push_back(res);
                        }
                    }
                    candc++;
                    count[posRec[i]] = 0;
                }
                posRec.clear();
            }
        }
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::edsEntityExtractWithBucketAndLazy(InvertList* invertList, vector<int>* index, const string &query, const int id)
{
    vector<int> count(query.length() - q + 1, 0);
    vector<int> posRec;
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));
    int nowIndex;
	while(minValue != MAX)
	{
		do{
            nowIndex = minIndex;
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)) && minIndex - nowIndex - 1 <= int(ceil((entityLengthMap[minValue] + q - 1) * (1 - det) * q / det)));

		int length = entityLengthMap[minValue];
		SetEditSimilarity(length);

		if((int)indexContainer.size() >= Tl) {
            for(int offset = Le; offset <= Te; offset++) {
                for(unsigned lp = 0; lp < indexContainer.size(); lp++) {
                    for(int i = indexContainer[lp] - offset + 1 > 0 ? indexContainer[lp] - offset + 1 : 0; i <= (indexContainer[lp] < (int)query.length() - q + 1 - offset ? indexContainer[lp] : (int)query.length() - q + 1 - offset); i++)
                    {
                        count[i]++;
                        if(count[i] == 1)
                            posRec.push_back(i);
                    }
                }
                for(int i = 0; i < (int)posRec.size(); i++) {
                    if(count[posRec[i]] >= GetEditSimilarityT(offset, length)) {
                        int c = 0;
                        int prefix, suffix;
                        for(prefix = 0; prefix < q - 1; prefix++) {
                            if(keywords[minValue][prefix] == query[posRec[i] + prefix])
                                c++;
                            else
                                break;
                        }
                        for(suffix = 1; suffix < q; suffix++) {
                            if(keywords[minValue][keywords[minValue].length() - suffix] == query[posRec[i] + offset + q - 1 - suffix])
                                c++;
                            else
                                break;
                        }
                        if(c + count[posRec[i]] >= GetEditSimilarityT(offset + 2 * q - 2, length + 2 * q - 2)) {
                            Cand res = {minValue, id, posRec[i], offset + q - 1};
                            result.push_back(res);
                        }
                    }
                    candb++;
                    count[posRec[i]] = 0;
                }
                posRec.clear();
            }
        }
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::edsEntityExtractWithMultHeap(InvertList* invertList, vector<int>* index, const string &query, const int id, int st, int len)
{
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetEditSimilarity(length);

		if(len - q + 1 >= Le && len - q + 1 <= Te)
		{
			if((int)indexContainer.size() >= GetEditSimilarityT(len - q + 1, length)) {
				int c = 0;
				int prefix, suffix;
				for(prefix = 0; prefix < q - 1; prefix++) {
					if(keywords[minValue][prefix] == query[st + prefix])
						c++;
					else
						break;
				}
				for(suffix = 1; suffix < q; suffix++) {
					if(keywords[minValue][keywords[minValue].length() - suffix] == query[st + len - suffix])
						c++;
					else
						break;
				}
				if(c + (int)indexContainer.size() >= GetEditSimilarityT(len + q - 1, length + 2 * q - 2)) {
					Cand res = {minValue, id, st, len};
					result.push_back(res);
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

/**********Jaccard Similariy**********/

void HeapBinary::jacEntityExtract(InvertList* invertList, vector<int>* index, const vector<string>* query, const int id, const int begin, const int end)
{
	indexContainer.clear();
	//int K = index->size();

	vector<int>* hp = new vector<int>();
	//hp->reserve(K+1);
	createStableHeap(invertList,index,hp,begin,end);

    int K = hp->size() - 1;
    if(K==0) return;
//    cout << "hp.size: " << hp->size() << endl;
	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetJaccard(length);

		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			int j, i = 0;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
					jacBinarySpan(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::jacEntityExtract(InvertList* invertList, vector<int>* index, const string &query, const int id, const int begin, const int end)
{
	indexContainer.clear();
	//int K = index->size();

	vector<int>* hp = new vector<int>();
	//hp->reserve(K+1);
	createStableHeap(invertList,index,hp,begin,end);

    int K = hp->size() - 1;
    if(K==0) return;
    //cout << "hp.size: " << hp->size() << endl;
	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{

		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetJaccard(length);

//        cout << "indexContainer.size: " << indexContainer.size() << endl;
//        cout << "Tl: " << Tl << endl;
		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			int j, i = 0;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
				    cout << "binary" << endl;
					jacBinarySpan(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::jacEntityExtract(InvertList* invertList, vector<int>* index, const vector<string>* query, const int id)
{
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetJaccard(length);

		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			int j, i = 0;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
					jacBinarySpan(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::jacEntityExtractWithNoPruning(InvertList* invertList, vector<int>* index, const vector<string>*query,const int id)
{
    vector<int> count(query->size(), 0);
    vector<int> posRec;
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetJaccard(length);

	    for(int offset = Le; offset <= Te; offset++) {
            for(unsigned lp = 0; lp < indexContainer.size(); lp++) {
                for(int i = indexContainer[lp] - offset + 1 > 0 ? indexContainer[lp] - offset + 1 : 0;i <= (indexContainer[lp] < (int)query->size() - offset ? indexContainer[lp] : (int)query->size() - offset);i++)
                {
                    count[i]++;
                    if(count[i] == 1)
                        posRec.push_back(i);
                }
            }
            for(int i = 0; i < (int)posRec.size(); i++) {
                candd++;
                if(count[posRec[i]] >= GetJaccardT(offset, length)) {
                    Cand res = {minValue, id, posRec[i], offset};
                    result.push_back(res);
                }
                count[posRec[i]] = 0;
            }
            posRec.clear();
        }
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::jacEntityExtractWithLazyUpdate(InvertList* invertList, vector<int>* index, const vector<string>*query,const int id)
{
    vector<int> count(query->size(), 0);
    vector<int> posRec;
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetJaccard(length);

		if((int)indexContainer.size() >= Tl) {
            for(int offset = Le; offset <= Te; offset++) {
                for(unsigned lp = 0; lp < indexContainer.size(); lp++) {
                    for(int i = indexContainer[lp] - offset + 1 > 0 ? indexContainer[lp] - offset + 1 : 0;i <= (indexContainer[lp] < (int)query->size() - offset ? indexContainer[lp] : (int)query->size() - offset);i++)
                    {
                        count[i]++;
                        if(count[i] == 1)
                            posRec.push_back(i);
                    }
                }
                for(int i = 0; i < (int)posRec.size(); i++) {
                    candc++;
                    if(count[posRec[i]] >= GetJaccardT(offset, length)) {
                        Cand res = {minValue, id, posRec[i], offset};
                        result.push_back(res);
                    }
                    count[posRec[i]] = 0;
                }
                posRec.clear();
            }
        }
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::jacEntityExtractWithBucketAndLazy(InvertList* invertList, vector<int>* index, const vector<string>*query,const int id)
{
    vector<int> count(query->size(), 0);
    vector<int> posRec;
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));
    int nowIndex;
	while(minValue != MAX)
	{
		do{
            nowIndex = minIndex;
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)) && minIndex - nowIndex - 1 <= (Te - Tl));

		int length = entityLengthMap[minValue];
		SetJaccard(length);

		if((int)indexContainer.size() >= Tl) {
            for(int offset = Le; offset <= Te; offset++) {
                for(unsigned lp = 0; lp < indexContainer.size(); lp++) {
                    for(int i = indexContainer[lp] - offset + 1 > 0 ? indexContainer[lp] - offset + 1 : 0;i <= (indexContainer[lp] < (int)query->size() - offset ? indexContainer[lp] : (int)query->size() - offset);i++)
                    {
                        count[i]++;
                        if(count[i] == 1)
                            posRec.push_back(i);
                    }
                }
                for(int i = 0; i < (int)posRec.size(); i++) {
                    candb++;
                    if(count[posRec[i]] >= GetJaccardT(offset, length)) {
                            Cand res = {minValue, id, posRec[i], offset};
                            result.push_back(res);
                    }
                    count[posRec[i]] = 0;
                }
                posRec.clear();
            }
        }
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::jacEntityExtractWithMultHeap(InvertList* invertList,vector<int>* index, const vector<string>*query, const int id, int st, int len)
{
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetJaccard(length);
		if(len >= Le && len <= Te)
		{
			if((int)indexContainer.size() >= GetJaccardT(len, length))
			{
				Cand res = {minValue, id, st, len};
				result.push_back(res);
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

/**********Cosine Similariy**********/

void HeapBinary::cosEntityExtract(InvertList* invertList, vector<int>* index, const vector<string>* query, const int id)
{
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetCosine(length);

		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			int j, i = 0;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
					cosBinarySpan(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::cosEntityExtract(InvertList* invertList, vector<int>* index, const string &query, const int id, const int begin, const int end)
{
	indexContainer.clear();
	vector<int>* hp = new vector<int>();
	createStableHeap(invertList,index,hp,begin,end);
	int K = hp->size() - 1;
//    cout << "hp.size: " << hp->size() << endl;
    if(K==0)
        return;

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetCosine(length);

		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			int j, i = 0;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
					cosBinarySpan(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::cosEntityExtract(InvertList* invertList, vector<int>* index, const string &query, const int id)
{
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetCosine(length);

		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			int j, i = 0;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
					cosBinarySpan(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

/**********Dice Similariy**********/

void HeapBinary::diceEntityExtract(InvertList* invertList, vector<int>* index, const vector<string>* query, const int id)
{
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetDice(length);

		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			int j, i = 0;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
				    cout << "binary" << endl;
					diceBinarySpan(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::diceEntityExtract(InvertList* invertList, vector<int>* index, const string &query, const int id, const int begin, const int end)
{
	indexContainer.clear();
	vector<int>* hp = new vector<int>();
	createStableHeap(invertList,index,hp,begin,end);
	int K = hp->size() - 1;

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetDice(length);

		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			int j, i = 0;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
				    cout << "binary" << endl;
					diceBinarySpan(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}

void HeapBinary::diceEntityExtract(InvertList* invertList, vector<int>* index, const string &query, const int id)
{
	indexContainer.clear();
	int K = index->size();

	vector<int>* hp = new vector<int>();
	hp->reserve(K+1);
	createStableHeap(invertList,index,hp);

	int minIndex = hp->at(1);
	int minValue = invertList->at(minIndex).second->at(index->at(minIndex));

	while(minValue != MAX)
	{
		do{
			indexContainer.push_back(invertList->at(minIndex).first);
			index->at(minIndex)++;
			adjustStableHeap(invertList,index,hp,1,K);
			minIndex = hp->at(1);
		}while(minValue == invertList->at(minIndex).second->at(index->at(minIndex)));

		int length = entityLengthMap[minValue];
		SetDice(length);

		if((int)indexContainer.size() >= Tl)
		{
			//cout << "BEGIN:" << Tl << endl;
			int j, i = 0;
			while (i < (int)indexContainer.size() - Tl + 1)
			{
				j = i + Tl - 1;
				if (indexContainer[j] - indexContainer[i] + 1 <= Te)
				{
					diceBinarySpan(i, j, min(i + Te - 1, int(indexContainer.size()) - 1), id, minValue, query);
				}
				else
				{
					//cout << "Shift: " << i << endl;
					//Binary Shift
					int lower = i;
					int upper = j;
					int mid = -1;
					while (lower <= upper)
					{
						mid = (lower + upper + 1)/2;
						if ((indexContainer[j] + (mid - i)) - indexContainer[mid] + 1 > Te)
							lower = mid + 1;
						else
							upper = mid - 1;
					}
					i = lower;
					//cout << "End Shift: " << i << endl;
				}
			}
		}
		minValue = invertList->at(minIndex).second->at(index->at(minIndex));
		indexContainer.clear();
	}
	delete hp;
}
