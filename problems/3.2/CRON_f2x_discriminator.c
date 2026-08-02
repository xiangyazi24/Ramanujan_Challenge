/* CRON_f2x_discriminator.c — GARQI-1' numerical discriminator (Q6444 appendix AD action item 1).
 *
 * For dyadic windows (X, 2X], X = 2^13..2^19, and every prime p with the Apery
 * first block reaching the window (p in (X/2, 2X]):
 *   N_p(X) = #{ r in Z_p : X < p + r <= 2X },  Z_p = { 0<=r<=p-2 : b_r = 0 mod p }
 *   S_X    = sum_p N_p(X)
 *   F2_X   = sum_p N_p(X) (N_p(X)-1)          (dyadic two-zero factorial moment)
 *   F2r_X  = reflection-forced part: pairs (r, p-1-r), r != p-1-r, both in window.
 * GARQI-1' predicts F2_X <= C X/log X, i.e. F2_X * log X / X bounded.
 * Globals over all primes 7 <= p < 2^20: mean |Z_p|, E|Z_p|^2 (Poisson-pair predicts ~3),
 * freq(|Z_p|=0) vs e^{-1/2} = 0.60653, odd-|Z_p| primes (expect exactly {11, 3137}).
 *
 * Ground-truth gates: Z_13 = {}, Z_17 = {3,13}, Z_29 = {}. Abort if violated.
 * cc -O2 -pthread -o CRON_f2x_discriminator CRON_f2x_discriminator.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

#define PMAX (1u<<20)
#define KMIN 13
#define KMAX 19
#define NW (KMAX-KMIN+1)
#define MAXZ 64

typedef unsigned long long u64;
typedef unsigned int u32;

static u32 *primes; static u32 nprimes;

static void sieve(void){
    static unsigned char comp[PMAX+1];
    primes = malloc(sizeof(u32)*90000);
    for (u32 i=2;i<=PMAX;i++) if(!comp[i]){
        if (i>=7) primes[nprimes++]=i;
        for (u64 j=(u64)i*i;j<=PMAX;j+=i) comp[j]=1;
    }
}

typedef struct {
    double S[NW], F2[NW], F2r[NW]; u64 npr[NW];
    double sumZ, sumZ2; u64 nzero, nprime;
    u32 oddp[64]; int nodd;
    int tid, nthreads;
    u32 *inv, *row_unused;
} W;

static void scan_prime(u32 p, W*w){
    u32 *inv = w->inv;
    inv[1]=1;
    for (u32 i=2;i<p;i++) inv[i] = (u32)(p - (u64)(p/i)*inv[p%i]%p);
    /* recurrence, collect zeros */
    u32 z[MAXZ]; int nz=0;
    u64 bm1=1%p, b=5%p;             /* b_0, b_1 */
    if (bm1==0){ z[nz++]=0; }
    if (p==5){ /* b_1=0 mod 5 */ }
    if (b==0 && nz<MAXZ) z[nz++]=1;
    for (u32 n=1;n<=p-3;n++){        /* produce b_{n+1}, r up to p-2 */
        u64 nn=n%p, n2=nn*nn%p, n3=n2*nn%p;
        u64 poly=(34*n3%p + 51*n2%p + 27*nn%p + 5)%p;
        u64 iv=inv[(n+1)%p]; u64 iv3=iv*iv%p*iv%p;
        u64 bn1=( (poly*b%p + p - n3*bm1%p) % p ) * iv3 % p;
        bm1=b; b=bn1;
        if (b==0 && nz<MAXZ) z[nz++]=n+1;
    }
    /* gates */
    if (p==13 && nz!=0){fprintf(stderr,"GATE FAIL p=13\n");exit(2);}
    if (p==17 && !(nz==2&&z[0]==3&&z[1]==13)){fprintf(stderr,"GATE FAIL p=17\n");exit(2);}
    if (p==29 && nz!=0){fprintf(stderr,"GATE FAIL p=29\n");exit(2);}
    w->nprime++; w->sumZ+=nz; w->sumZ2+=(double)nz*nz; if(nz==0)w->nzero++;
    if (nz&1){ if(w->nodd<64) w->oddp[w->nodd]=p; w->nodd++; }
    /* windows */
    for (int k=KMIN;k<=KMAX;k++){
        u64 X=1ull<<k;
        if (p>2*X) continue; if ((u64)2*p<X+2) continue;
        long long rlo = (long long)X+1-(long long)p; if (rlo<0) rlo=0;
        long long rhi = (long long)2*X-(long long)p; if (rhi>(long long)p-2) rhi=p-2;
        if (rhi<rlo) continue;
        int N=0;
        for (int i=0;i<nz;i++) if (z[i]>=rlo && z[i]<=rhi) N++;
        int idx=k-KMIN;
        w->npr[idx]++; w->S[idx]+=N; w->F2[idx]+=(double)N*(N-1);
        /* reflection-forced pairs: r and p-1-r both zeros (guaranteed) and both in window */
        for (int i=0;i<nz;i++){
            u32 r=z[i], rm=p-1-r;
            if (rm==r) continue;
            if (r<rm) { /* count each mirror pair once, then x2 for ordered convention */
                if (r>=rlo&&r<=rhi&&rm>=(u64)rlo&&rm<=(u64)rhi) w->F2r[idx]+=2.0;
            }
        }
    }
}

static void* worker(void*arg){
    W*w=(W*)arg;
    w->inv=malloc(sizeof(u32)*(PMAX+2));
    time_t t0=time(0), last=t0;
    for (u32 i=w->tid;i<nprimes;i+=w->nthreads){
        scan_prime(primes[i],w);
        if (w->tid==0 && time(0)-last>=5){
            last=time(0);
            fprintf(stderr,"[progress] %.1f%% (p=%u) elapsed %lds\n",
                100.0*i/nprimes, primes[i], last-t0);
        }
    }
    free(w->inv); return 0;
}

int main(int argc,char**argv){
    int NT = argc>1?atoi(argv[1]):8;
    sieve();
    fprintf(stderr,"primes 7..%u: %u, threads=%d\n",PMAX,nprimes,NT);
    W *ws=calloc(NT,sizeof(W)); pthread_t th[64];
    for (int t=0;t<NT;t++){ws[t].tid=t;ws[t].nthreads=NT;pthread_create(&th[t],0,worker,&ws[t]);}
    for (int t=0;t<NT;t++) pthread_join(th[t],0);
    W g; memset(&g,0,sizeof g);
    for (int t=0;t<NT;t++){
        g.nprime+=ws[t].nprime; g.sumZ+=ws[t].sumZ; g.sumZ2+=ws[t].sumZ2; g.nzero+=ws[t].nzero;
        for(int i=0;i<ws[t].nodd&&g.nodd<64;i++) g.oddp[g.nodd++]=ws[t].oddp[i];
        for(int k=0;k<NW;k++){g.S[k]+=ws[t].S[k];g.F2[k]+=ws[t].F2[k];g.F2r[k]+=ws[t].F2r[k];g.npr[k]+=ws[t].npr[k];}
    }
    printf("== GARQI-1' discriminator, primes 7..%u (%llu primes) ==\n",PMAX,(u64)g.nprime);
    printf("global: mean|Z|=%.4f  E|Z|^2=%.4f (Poisson-pair predicts 3)  freq(|Z|=0)=%.5f (e^-1/2=0.60653)\n",
        g.sumZ/g.nprime, g.sumZ2/g.nprime, (double)g.nzero/g.nprime);
    printf("odd-|Z| primes (%d): ",g.nodd); for(int i=0;i<g.nodd;i++)printf("%u ",g.oddp[i]); printf("\n\n");
    printf("%-3s %10s %8s %12s %12s %12s | %10s %10s %10s\n",
        "k","X","#primes","S_X","F2_X","F2r_X","S*lnX/X","F2*lnX/X","F2r/F2");
    for (int k=KMIN;k<=KMAX;k++){
        int i=k-KMIN; double X=(double)(1ull<<k), lx=log(X);
        printf("%-3d %10llu %8llu %12.0f %12.0f %12.0f | %10.4f %10.4f %10.4f\n",
            k,(u64)1<<k,g.npr[i],g.S[i],g.F2[i],g.F2r[i],
            g.S[i]*lx/X, g.F2[i]*lx/X, g.F2[i]>0?g.F2r[i]/g.F2[i]:0.0);
    }
    return 0;
}
