/* Q6651 public CI compile interface only. The delivered engine includes the
 * unchanged private literal parent instead of this stub. */
#include <stdint.h>
enum Addr { AD_0=0, AD_1=1, AD_HASH=2 };
enum Trail { WT_TODO=0, WT_D0=1, WT_D1=2, WT_F0=3, WT_F1=4 };
enum NBit { NB_X=0, NB_0=1, NB_1=2, NB_END=3 };
enum CBit { CB_X=0, CB_0=1, CB_1=2, CB_END=3 };
enum Phase { PH_IDLE=0, PH_CLEAR=1, PH_SCAN=2, PH_INIT=3, PH_VALID=4,
 PH_EQUAL=5, PH_READ=6, PH_ADV=7, PH_INC=8, PH_WRITE=9, PH_DONE=10,
 PH_ABORT=11, PH_SINGLE0=12, PH_SINGLE1=13, PH_BAD=14, PH_RSVD=15,
 PH_COUNT=16 };
enum Op { OP_SCAN=0, OP_INIT=1, OP_VALID=2, OP_EQUAL=3, OP_READ=4,
 OP_ADV=5, OP_INC=6, OP_WRITE=7, OP_WRITE_FINAL=8, OP_CLEAR=9,
 OP_COUNT=10 };
enum HDir { HD_NONE=0, HD_ABORT=1, HD_OUT=2, HD_RET=3 };
typedef struct { uint8_t dir,op,acc; } Head;
typedef struct Cell {
 uint8_t L,B,H,S,e,mu,w,rw,fresh;
 uint8_t ad,wt,nb,cb,rb,fu;
 uint8_t phase,par,final,result;
 Head h;
 uint8_t repair_pending,repair_prev;
} Cell;
enum { CR_NONE=6 };
static int K=2,cutauth_control,last_gap_defect;
static Head h_none(void){Head h={HD_NONE,0,0};return h;}
static int chain_anchor(const Cell*s){return s->rw&&!s->w;}
static int old_follower(const Cell*s){return !s->L&&!s->rw;}
static int auth_frontier(const Cell*l,const Cell*m){return old_follower(m)&&chain_anchor(l);}
static int effective_cut(const Cell*l,const Cell*m){return m->L||(l->L&&!l->rw)||auth_frontier(l,m)||l->repair_pending||m->repair_pending;}
static void cutauth_after_birth(Cell*a,Cell*b){(void)a;b->L=1;b->B=0;b->H=1;b->S=0;b->e=(uint8_t)((b->e+1)%K);b->ad=AD_HASH;b->mu=0;b->rw=b->w=b->repair_pending=0;b->repair_prev=CR_NONE;b->phase=PH_CLEAR;b->wt=WT_D0;b->h=h_none();}
static int parent_pair_step(Cell*a,Cell*b,int coin,int*changed){(void)a;(void)b;(void)coin;last_gap_defect=0;*changed=0;return 0;}
