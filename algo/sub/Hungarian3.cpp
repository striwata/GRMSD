/* ?ｽn?ｽ?ｽ?ｽK?ｽ?ｽ?ｽ[?ｽ@?ｽﾌ計?ｽZ?ｽiC++?ｽj?ｽi?ｽT?ｽC?ｽY?ｽﾌ異なゑｿｽ鼾?ｿｽj*/

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

void hungarian3(const double *b,double *permutation,const double *n_data1,const double *n_data2) {
  vector<double> *a;
  int n1,n2, p, q;
  double n1_ = *n_data1;
  double n2_ = *n_data2;
  n1 = int(n1_);
  n2 = int(n2_);
  int inf = 100;
  a = new vector<double>;
  for(int i = 0;i<n1*n2;i++) a->push_back(b[i]);
  vector<double> fx(n2, inf), fy(n1, 0);
  vector<int> x(n2, -1), y(n1, -1);
  for (int i = 0; i < n2; ++i)
    for (int j = 0; j < n1; ++j)
		fx[i] = max(fx[i], a->at(n2*j+i));
  for (int i = 0; i < n2; ) {
    vector<int> t(n1, -1), s(n2+1, i);
    for (p = q = 0; p <= q && x[i] < 0; ++p)
      for (int k = s[p], j = 0; j < n1 && x[i] < 0; ++j)
        if (abs(fx[k] + fy[j] - a->at(n2*j+k))<0.000000001 && t[j] < 0) {
          s[++q] = y[j], t[j] = k;
          if (s[q] < 0)
            for (p = j; p >= 0; j = p)
              y[j] = k = t[j], p = x[k], x[k] = j;
        }
    if (x[i] < 0) {
      weight d = inf;
      for (int k = 0; k <= q; ++k)
        for (int j = 0; j < n1; ++j)
          if (t[j] < 0) d = min(d, fx[s[k]] + fy[j] - a->at(n2*j+s[k]));
      for (int j = 0; j < n1; ++j) fy[j] += (t[j] < 0 ? 0 : d);
      for (int k = 0; k <= q; ++k) fx[s[k]] -= d;
    } else ++i;
  }
  /*
  weight ret = 0;
  for (int i = 0; i < n; i++) {
	  ret += a->at(n*x[i]+i);
  }*/
  for(int i = 0;i<n2;++i) permutation[n2*x[i]+i] = 1;
  delete a;
}

void mexFunction( int Nreturned, mxArray *returned[], int Noperand, const mxArray *operand[] ){
  double *a;
  double *permutation;
  double *n_data1;
  double *n_data2;
  double *permutation_out;
  int rows,cols;

  /* a?ｽﾌ行?ｽﾌ撰ｿｽ?ｽ?ｽﾛ托ｿｽ */  
  rows = mxGetM(operand[0]);
  cols = mxGetN(operand[0]);

  *returned = mxCreateDoubleMatrix(rows * cols,1,mxREAL);
  
  /* Matlab?ｽ?ｽ?ｽﾌ変撰ｿｽ?ｽﾌア?ｽh?ｽ?ｽ?ｽX?ｽ?ｽC?ｽ?ｽ?ｽﾌ変撰ｿｽ?ｽﾉコ?ｽs?ｽ[?ｽ?ｽ?ｽ?ｽ */
  a = mxGetPr(operand[0]);
  permutation = mxGetPr(operand[1]);  
  n_data1 = mxGetPr(operand[2]);
  n_data2 = mxGetPr(operand[3]);
  
  permutation_out = mxGetPr(returned[0]);
  
  hungarian3(a,permutation,n_data1,n_data2);
  for (int i = 0; i < rows * cols; i++)
  {
      permutation_out[i] = permutation[i];
  }
}
