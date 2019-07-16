/* ハンガリー法の計算（C++）*/
/* calculations of Hungarian method (C++) */

#define _GLIBCXX_DEBUG
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "mex.h"
#include "matrix.h"
#include <cstdlib>

using namespace std;

#define REP(i,n) for(int i=0;i<(int)n;++i)
#define FOR(i,c) for(__typeof((c).begin())i=(c).begin();i!=(c).end();++i)
#define ALL(c) (c).begin(), (c).end()

typedef double weight;

void hungarian(const double *b,double *permutation,const double *n_data) {
  vector<double> *a;
  int n, p, q;
  double n_ = *n_data;
  n = int(n_);
  int inf = 100;
  a = new vector<double>;
  for(int i = 0;i<n*n;i++) a->push_back(b[i]);
  vector<double> fx(n, inf), fy(n, 0);
  vector<int> x(n, -1), y(n, -1);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
		fx[i] = max(fx[i], a->at(n*j+i));
  for (int i = 0; i < n; ) {
    vector<int> t(n, -1), s(n+1, i);
    for (p = q = 0; p <= q && x[i] < 0; ++p)
      for (int k = s[p], j = 0; j < n && x[i] < 0; ++j)
        if (abs(fx[k] + fy[j] - a->at(n*j+k))<0.000000001 && t[j] < 0) {
          s[++q] = y[j], t[j] = k;
          if (s[q] < 0)
            for (p = j; p >= 0; j = p)
              y[j] = k = t[j], p = x[k], x[k] = j;
        }
    if (x[i] < 0) {
      weight d = inf;
      for (int k = 0; k <= q; ++k)
        for (int j = 0; j < n; ++j)
          if (t[j] < 0) d = min(d, fx[s[k]] + fy[j] - a->at(n*j+s[k]));
      for (int j = 0; j < n; ++j) fy[j] += (t[j] < 0 ? 0 : d);
      for (int k = 0; k <= q; ++k) fx[s[k]] -= d;
    } else ++i;
  }
  /*
  weight ret = 0;
  for (int i = 0; i < n; i++) {
	  ret += a->at(n*x[i]+i);
  }*/
  for(int i = 0;i<n;++i) permutation[n*x[i]+i] = 1;
  delete a;
}

void mexFunction( int Nreturned, mxArray *returned[], int Noperand, const mxArray *operand[] ){
  double *a;
  double *permutation;
  double *n_data;
  double *permutation_out;
  int rows;

/* aの行の数を保存 */
/* store the number of rows of a */
  rows = mxGetM(operand[0]);

  *returned = mxCreateDoubleMatrix(rows * rows,1,mxREAL);
  
/* Matlab側の変数のアドレスをC側の変数にコピーする */
/* copy the address of Matlab-variables to C-variables */
  a = mxGetPr(operand[0]);
  permutation = mxGetPr(operand[1]);  
  n_data = mxGetPr(operand[2]);
  
  permutation_out = mxGetPr(returned[0]);
  
  hungarian(a,permutation,n_data);
  for (int i = 0; i < rows * rows; i++)
  {
      permutation_out[i] = permutation[i];
  }
}

