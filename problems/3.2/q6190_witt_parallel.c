/* Q6190 scratch-only parallel first-Witt/companion summary census. */
#include <inttypes.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned __int128 u128;

typedef struct { uint32_t r, eta, comp; } Zero;
typedef struct {
  uint64_t primes, total_zeros, reflection_pairs, consecutive_pairs, midpoint_zeros;
  uint64_t g[3], lambda[3], bock[3];
  uint32_t maxz, maxp;
} Stats;

static uint64_t modpow64(uint64_t a, uint64_t e, uint64_t m) {
  uint64_t out = 1 % m;
  while (e) { if (e & 1) out = (uint64_t)((u128)out*a % m); a = (uint64_t)((u128)a*a % m); e >>= 1; }
  return out;
}
static int color(uint32_t a, uint32_t p) {
  if (a == 0) return 0;
  uint32_t q = (uint32_t)modpow64(a, (p-1)/2, p);
  if (q == 1) return 1;
  if (q == p-1) return 2;
  abort();
}
static uint32_t subp(uint32_t a, uint32_t b, uint32_t p) { return a >= b ? a-b : a+p-b; }
static uint32_t addp(uint32_t a, uint32_t b, uint32_t p) { uint64_t s=(uint64_t)a+b; return (uint32_t)(s>=p?s-p:s); }
static uint64_t polyP(uint64_t n, uint64_t mod) { return (uint64_t)((((u128)34*n+51)*n+27)*n+5)%mod; }
static uint64_t pow6(uint64_t n, uint64_t mod) {
  uint64_t n2=(uint64_t)((u128)n*n%mod), n3=(uint64_t)((u128)n2*n%mod);
  return (uint64_t)((u128)n3*n3%mod);
}
static uint64_t nextv(uint64_t Pn,uint64_t cur,uint64_t n6,uint64_t prev,uint64_t mod) {
  uint64_t a=(uint64_t)((u128)Pn*cur%mod), b=(uint64_t)((u128)n6*prev%mod);
  return a>=b?a-b:a+mod-b;
}
static uint8_t *sieve(uint32_t limit) {
  uint8_t *a=calloc((size_t)limit+1,1); if(!a) abort(); memset(a,1,(size_t)limit+1); a[0]=a[1]=0;
  for(uint32_t d=2;(uint64_t)d*d<=limit;++d) if(a[d]) for(uint32_t m=d*d;m<=limit;m+=d) a[m]=0;
  return a;
}
static void merge(Stats *a, const Stats *b) {
  a->primes+=b->primes; a->total_zeros+=b->total_zeros; a->reflection_pairs+=b->reflection_pairs;
  a->consecutive_pairs+=b->consecutive_pairs; a->midpoint_zeros+=b->midpoint_zeros;
  for(int i=0;i<3;++i){a->g[i]+=b->g[i];a->lambda[i]+=b->lambda[i];a->bock[i]+=b->bock[i];}
  if(b->maxz>a->maxz || (b->maxz==a->maxz && b->maxz && b->maxp<a->maxp)){a->maxz=b->maxz;a->maxp=b->maxp;}
}
int main(int argc,char **argv){
  uint32_t limit=1000000; if(argc>1) limit=(uint32_t)strtoul(argv[1],NULL,10);
  if(limit<5 || limit>1000000){fputs("limit out of range\n",stderr);return 2;}
  uint8_t *isprime=sieve(limit); Stats total={0}; int threads=0;
#pragma omp parallel
  {
    Stats s={0};
#pragma omp single
    threads=omp_get_num_threads();
#pragma omp for schedule(dynamic,8)
    for(int64_t pp=5;pp<=(int64_t)limit;++pp){
      uint32_t p=(uint32_t)pp; if(!isprime[p]) continue; ++s.primes;
      uint64_t p2=(uint64_t)p*p, bprev=1%p2,bcur=5%p2; uint32_t cprev=0,ccur=1%p,fact=1;
      Zero z[256]; uint32_t nz=0;
      if(bcur%p==0){ z[nz++]=(Zero){1,(uint32_t)((bcur/p)%p),ccur}; }
      for(uint32_t n=1;n<p-1;++n){
        uint32_t r=n+1; uint64_t P2=polyP(n,p2),N62=pow6(n,p2),bn=nextv(P2,bcur,N62,bprev,p2);
        uint32_t P1=(uint32_t)(P2%p),N61=(uint32_t)(N62%p),cn=(uint32_t)nextv(P1,ccur,N61,cprev,p);
        fact=(uint32_t)((uint64_t)fact*r%p);
        if(bn%p==0){
          if(nz>=256) abort(); uint32_t f2=(uint32_t)((uint64_t)fact*fact%p),f3=(uint32_t)((uint64_t)f2*fact%p);
          uint32_t invf3=(uint32_t)modpow64(f3,p-2,p), en=(uint32_t)((bn/p)%p);
          z[nz++]=(Zero){r,(uint32_t)((uint64_t)en*invf3%p),(uint32_t)((uint64_t)cn*invf3%p)};
        }
        bprev=bcur;bcur=bn;cprev=ccur;ccur=cn;
      }
      s.total_zeros+=nz; if(nz>s.maxz || (nz==s.maxz && nz && p<s.maxp)){s.maxz=nz;s.maxp=p;}
      for(uint32_t i=0;i<nz;++i){
        uint32_t j=nz-1-i; if(z[i].r+z[j].r!=p-1) abort();
        if(z[i].r==(p-1)/2) ++s.midpoint_zeros;
      }
      for(uint32_t i=0;i<nz/2;++i){uint32_t j=nz-1-i,g=subp(z[i].eta,z[j].eta,p);++s.g[color(g,p)];++s.reflection_pairs;}
      for(uint32_t i=0;i+1<nz;++i){
        Zero x=z[i],y=z[i+1]; if(x.comp==0||y.comp==0) abort();
        uint32_t inv=(uint32_t)modpow64(x.comp,p-2,p),lam=(uint32_t)((uint64_t)y.comp*inv%p); lam=lam?p-lam:0;
        uint32_t B=addp(y.eta,(uint32_t)((uint64_t)lam*x.eta%p),p);
        ++s.lambda[color(lam,p)];++s.bock[color(B,p)];++s.consecutive_pairs;
      }
    }
#pragma omp critical
    merge(&total,&s);
  }
  free(isprime);
  printf("Q6190_WITT_PARALLEL limit=%u threads=%d primes=%"PRIu64" total_zeros=%"PRIu64" maxZ=%u maxZ_p=%u\n",limit,threads,total.primes,total.total_zeros,total.maxz,total.maxp);
  printf("midpoint_zeros=%"PRIu64"\n",total.midpoint_zeros);
  printf("reflection_pairs=%"PRIu64" g_color_zero=%"PRIu64" g_color_square=%"PRIu64" g_color_nonsquare=%"PRIu64"\n",total.reflection_pairs,total.g[0],total.g[1],total.g[2]);
  printf("consecutive_zero_pairs=%"PRIu64" lambda_zero=%"PRIu64" lambda_square=%"PRIu64" lambda_nonsquare=%"PRIu64"\n",total.consecutive_pairs,total.lambda[0],total.lambda[1],total.lambda[2]);
  printf("bockstein_zero=%"PRIu64" bockstein_square=%"PRIu64" bockstein_nonsquare=%"PRIu64"\n",total.bock[0],total.bock[1],total.bock[2]);
  return 0;
}
