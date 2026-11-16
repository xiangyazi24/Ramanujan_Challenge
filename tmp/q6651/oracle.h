#ifndef EXACT_D_ORACLE_H
#define EXACT_D_ORACLE_H
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
enum { ED_EXACT=1, ED_SHORT=2, ED_NONDET=3 };
#define ED_PREFIX (1u<<22)
static unsigned char *ed_factors[9], *ed_al;
static unsigned ed_boundary_difference_count[9];
static uint64_t ed_count(int m){uint64_t a=1,b=2;for(int j=0;j<m;j++){uint64_t c=a+b;a=b;b=c;}return a;}
static unsigned char ed_alpha_symbol(uint64_t k){
 uint64_t doubles=k&(k>>1);unsigned a=2;
 if(doubles){unsigned t=(unsigned)__builtin_ctzll(doubles);uint64_t tail=k&((1ULL<<(t+2))-1);int L=64-__builtin_clzll(tail);uint64_t rank=3;
  for(int l=2;l<L;l++)rank+=ed_count(l>2?l-3:0);
  int m=L-2;for(int p=m-1;p>=0;p--)if((tail>>p)&1)rank+=ed_count(p);
  unsigned bl=k?64-__builtin_clzll(k):1;if(rank<=bl)a=(k>>rank)&1;
 }
 return (unsigned char)(2*a+(__builtin_popcountll(k)&1));
}
static void ed_oracle_init(int max_len){
 if(max_len<1||max_len>8){fprintf(stderr,"oracle length outside 1..8\n");exit(2);}if(ed_al)return;
 ed_al=malloc(ED_PREFIX+8);if(!ed_al)exit(2);for(unsigned i=0;i<ED_PREFIX+8;i++)ed_al[i]=ed_alpha_symbol(i);
 unsigned size=1;for(int l=1;l<=8;l++){size*=6;ed_factors[l]=calloc(size,1);if(!ed_factors[l])exit(2);}
 for(unsigned at=0;at<ED_PREFIX;at++){unsigned code=0;for(unsigned l=1;l<=8 && at+l<=ED_PREFIX;l++){code=code*6+ed_al[at+l-1];ed_factors[l][code]=1;}}
 for(unsigned l=1;l<=8;l++)for(unsigned at=ED_PREFIX-l+1;at<ED_PREFIX;at++){unsigned code=0;for(unsigned j=0;j<l;j++)code=code*6+ed_al[at+j];if(!ed_factors[l][code])ed_boundary_difference_count[l]++;}
}
static int ed_illegal(const unsigned char*w,int len,int policy){
 static const unsigned char triples[11][3]={{5,0,0},{5,0,1},{4,4,5},{4,5,0},{0,5,0},{4,0,5},{5,0,4},{5,1,4},{0,3,5},{5,4,1},{4,1,1}};
 if(len<1||len>8){fprintf(stderr,"oracle query length outside 1..8\n");exit(2);}
 for(int i=0;i<len;i++)if(w[i]>=6){fprintf(stderr,"oracle symbol outside 0..5\n");exit(2);}
 if(policy==ED_SHORT){for(int i=0;i+2<len;i++){if((w[i]&1)==(w[i+1]&1)&&(w[i]&1)==(w[i+2]&1))return 1;for(int j=0;j<11;j++)if(w[i]==triples[j][0]&&w[i+1]==triples[j][1]&&w[i+2]==triples[j][2])return 1;}for(int i=0;i+3<len;i++)if((w[i]==4&&w[i+1]==5&&w[i+2]==4&&w[i+3]==5)||(w[i]==5&&w[i+1]==4&&w[i+2]==5&&w[i+3]==5))return 1;return 0;}
 if(policy!=ED_EXACT||!ed_al)exit(2);unsigned code=0;for(int i=0;i<len;i++)code=code*6+w[i];return !ed_factors[len][code];
}
static void exactd_oracle_init(int n){ed_oracle_init(n);}
static unsigned exactd_alpha(unsigned k){return ed_alpha_symbol(k);}
static int exactd_illegal(const unsigned char*w,int len,int short_mode){return ed_illegal(w,len,short_mode?ED_SHORT:ED_EXACT);}
#endif
