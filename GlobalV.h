#ifndef GLB_HEADER
#define GLB_HEADER

#include "TMinuit.h"
#include "TMath.h"
#include "TFile.h"
#include <iostream>
#include <sstream>
#include <string>


using namespace std;

extern double Bin_LLinclu[20];
extern double Bin_LTinclu[20];
extern double Bin_TLinclu[20];
extern double Bin_TTinclu[20];
extern double BinError_LLinclu[20];
extern double BinError_LTinclu[20];
extern double BinError_TLinclu[20];
extern double BinError_TTinclu[20];

extern double Bin_LL1[20];
extern double Bin_LL2[20];
extern double Bin_LL3[20];
extern double Bin_LL4[20];

extern double Bin_LT1[20];
extern double Bin_LT2[20];
extern double Bin_LT3[20];
extern double Bin_LT4[20];

extern double Bin_TL1[20];
extern double Bin_TL2[20];
extern double Bin_TL3[20];
extern double Bin_TL4[20];

extern double Bin_TT1[20];
extern double Bin_TT2[20];
extern double Bin_TT3[20];
extern double Bin_TT4[20];

extern double Alpha[5];
extern double FitError[5];

extern string ToString(double num);

#endif
