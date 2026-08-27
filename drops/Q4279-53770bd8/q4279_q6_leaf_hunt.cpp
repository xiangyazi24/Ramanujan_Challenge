// Q4279: exact candidate-first hunt for an actual selected q6 Apéry leaf.
//
// Build (standard library only):
//   g++ -O3 -march=native -std=c++20 -pthread -Wall -Wextra -Wpedantic \
//       q4279_q6_leaf_hunt.cpp -o q4279_q6_leaf_hunt
//
// The scan computes Apéry zero sets with the division-free cleared recurrence
//
//   B_(n+1) = (34n^3+51n^2+27n+5) B_n - n^6 B_(n-1) (mod p),
//   B_n=(n!)^3 b_n,
//
// only through (p-1)/2 and uses the proved reflection symmetry.  It then
// projects ordered lower-zero pairs by the exact inverse formulae
//
//   plus:  q=(7p-r+a)/6,      t=2r-a-p,
//   minus: q=p+(r+a+1)/7,    t=(r-6a-6)/7,
//
// applies the raw minus-first tie-break, and only then tests t in Z_q.
// Every reported nonexistence statement is explicitly finite-range.

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i64 = std::int64_t;
using i32 = std::int32_t;

struct Options {
    u32 qmin_exclusive = 1000000;
    u32 qmax = 1250000;
    u32 threads = 0;
    u32 block_width = 250000;
    std::string out_dir = "q4279-scan";
};

Options parse_options(int argc, char **argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        auto need = [&](const char *name) {
            if (++i >= argc) throw std::runtime_error(std::string("missing value after ") + name);
            return std::string(argv[i]);
        };
        if (arg == "--qmin-exclusive") o.qmin_exclusive = static_cast<u32>(std::stoul(need("--qmin-exclusive")));
        else if (arg == "--qmax") o.qmax = static_cast<u32>(std::stoul(need("--qmax")));
        else if (arg == "--threads") o.threads = static_cast<u32>(std::stoul(need("--threads")));
        else if (arg == "--block-width") o.block_width = static_cast<u32>(std::stoul(need("--block-width")));
        else if (arg == "--out") o.out_dir = need("--out");
        else if (arg == "--help" || arg == "-h") {
            std::cout << "usage: q4279_q6_leaf_hunt --qmin-exclusive N --qmax N "
                         "[--threads N] [--block-width N] [--out DIR]\n";
            std::exit(0);
        } else throw std::runtime_error("unknown option: " + arg);
    }
    if (o.qmin_exclusive < 1000) throw std::runtime_error("qmin-exclusive must be >=1000");
    if (o.qmax <= o.qmin_exclusive) throw std::runtime_error("qmax must exceed qmin-exclusive");
    if (o.qmax > 20000000U) throw std::runtime_error("certified release cap is qmax=20,000,000");
    if (o.block_width == 0) throw std::runtime_error("block-width must be positive");
    if (o.threads == 0) o.threads = std::max(1U, std::thread::hardware_concurrency());
    return o;
}

std::vector<u32> primes_upto(u32 limit) {
    std::vector<std::uint8_t> composite(static_cast<std::size_t>(limit / 2) + 1, 0);
    const u32 root = static_cast<u32>(std::sqrt(static_cast<long double>(limit)));
    for (u32 p = 3; p <= root; p += 2) if (!composite[p >> 1]) {
        for (u64 n = u64(p) * p; n <= limit; n += 2ULL * p)
            composite[static_cast<std::size_t>(n >> 1)] = 1;
    }
    std::vector<u32> out;
    out.reserve(static_cast<std::size_t>(1.1L * limit / std::max<long double>(2, std::log(limit))));
    if (limit >= 2) out.push_back(2);
    for (u32 n = 3; n <= limit; n += 2) if (!composite[n >> 1]) out.push_back(n);
    return out;
}

// Montgomery arithmetic for odd p<2^30.  Values stay in [0,p).
struct Mont32 {
    u32 p, nprime;
    explicit Mont32(u32 modulus) : p(modulus), nprime(0) {
        if ((p & 1U) == 0 || p >= (1U << 30)) throw std::runtime_error("bad Montgomery modulus");
        u32 inv = p;
        for (int i = 0; i < 5; ++i) inv *= (2U - p * inv);
        if (u32(u64(inv) * p) != 1U) throw std::runtime_error("Montgomery inverse failed");
        nprime = 0U - inv;
    }
    u32 reduce(u64 x) const {
        const u32 m = u32(x) * nprime;
        u64 u = (x + u64(m) * p) >> 32;
        if (u >= p) u -= p;
        return static_cast<u32>(u);
    }
    u32 to_mont(u32 x) const { return static_cast<u32>((u64(x % p) << 32) % p); }
    u32 from_mont(u32 x) const { return reduce(x); }
    u32 mul(u32 a, u32 b) const { return reduce(u64(a) * b); }
    u32 add(u32 a, u32 b) const { u32 s = a + b; if (s >= p) s -= p; return s; }
    u32 sub(u32 a, u32 b) const { return a >= b ? a - b : a + p - b; }
};

u64 P_exact(u64 n) { return (2*n+1)*(17*n*n+17*n+5); }

u64 pow_mod(u64 a, u64 e, u64 p) {
    u64 r = 1;
    while (e) {
        if (e & 1U) r = static_cast<u64>((__uint128_t(r) * a) % p);
        a = static_cast<u64>((__uint128_t(a) * a) % p);
        e >>= 1U;
    }
    return r;
}

std::vector<u32> zeros_mont(u32 p, u64 *steps = nullptr) {
    if (p <= 5) return {};
    Mont32 m{p}; // braces: no most-vexing parse.
    u32 y0 = m.to_mont(1), y1 = m.to_mont(5);

    // Forward differences at n=1 for P(n).
    u32 pv=m.to_mont(117), p1=m.to_mont(418), p2=m.to_mont(510), p3=m.to_mont(204);
    // Forward differences at n=1 for n^6.
    u32 nv=m.to_mont(1), n1=m.to_mont(63), n2=m.to_mont(602), n3=m.to_mont(2100);
    u32 n4=m.to_mont(3360), n5=m.to_mont(2520), n6=m.to_mont(720);

    const u32 center=(p-1)/2;
    std::vector<u32> low;
    low.reserve(8);
    for (u32 n=1; n<center; ++n) {
        const u32 y2=m.sub(m.mul(pv,y1),m.mul(nv,y0));
        if (y2==0) low.push_back(n+1);
        y0=y1; y1=y2;
        pv=m.add(pv,p1); p1=m.add(p1,p2); p2=m.add(p2,p3);
        nv=m.add(nv,n1); n1=m.add(n1,n2); n2=m.add(n2,n3);
        n3=m.add(n3,n4); n4=m.add(n4,n5); n5=m.add(n5,n6);
    }
    if (steps) *steps = center ? center-1 : 0;
    std::vector<u32> z;
    z.reserve(low.size()*2);
    for (u32 r:low) {
        z.push_back(r);
        const u32 s=p-1-r;
        if (s!=r) z.push_back(s);
    }
    std::sort(z.begin(),z.end());
    z.erase(std::unique(z.begin(),z.end()),z.end());
    return z;
}

std::vector<u32> zeros_full(u32 p) {
    std::vector<u32> z;
    u64 b0=1%p,b1=5%p;
    if (b0==0) z.push_back(0);
    if (b1==0) z.push_back(1);
    for (u32 n=1;n<=p-2;++n) {
        const u64 nn=n;
        const u64 pp=P_exact(nn)%p;
        const u64 n3=nn*nn%p*nn%p;
        const u64 rhs=(pp*b1+p-n3*b0%p)%p;
        u64 den=u64(n+1)%p; den=den*den%p*u64(n+1)%p;
        const u64 b2=static_cast<u64>((__uint128_t(rhs)*pow_mod(den,p-2,p))%p);
        if (b2==0) z.push_back(n+1);
        b0=b1;b1=b2;
    }
    return z;
}

bool has_zero(const std::vector<u32>&z,i64 r) {
    return r>=0 && r<=std::numeric_limits<u32>::max() &&
           std::binary_search(z.begin(),z.end(),static_cast<u32>(r));
}

struct RawHit { u32 q,t,p,r,a; char sign; };

bool raw_less(const RawHit&x,const RawHit&y) {
    return std::tie(x.q,x.t,x.p,x.sign,x.r,x.a)<std::tie(y.q,y.t,y.p,y.sign,y.r,y.a);
}

struct ProjectionStats {
    u64 p_seen=0,p_eligible=0,ordered_pairs=0,parity_pruned=0,gap_pruned=0;
    u64 plus_integral=0,plus_prime=0,minus_integral=0,minus_prime=0;
};

std::vector<RawHit> project_hits(const std::vector<u32>&primes,
                                 const std::vector<i32>&index,
                                 const std::vector<std::vector<u32>>&zeros,
                                 u32 pmin,u32 qlo,u32 qmax,
                                 ProjectionStats*st=nullptr) {
    std::vector<RawHit> h;
    auto begin=std::lower_bound(primes.begin(),primes.end(),pmin);
    for (auto it=begin;it!=primes.end();++it) {
        const std::size_t pi=static_cast<std::size_t>(it-primes.begin());
        const u32 p=*it; const auto&zp=zeros[pi];
        if(st)++st->p_seen;
        if(zp.size()<3)continue;
        if(st)++st->p_eligible;
        for(u32 r:zp)for(u32 a:zp) {
            if(r==a)continue;
            if(st)++st->ordered_pairs;
            if(((r^a)&1U)==0){if(st)++st->parity_pruned;continue;}
            const u32 gap=r>a?r-a:a-r;
            if(gap<3){if(st)++st->gap_pruned;continue;}

            const i64 pn=7LL*p-i64(r)+i64(a);
            if(pn>0 && pn%6==0) {
                if(st)++st->plus_integral;
                const i64 q=pn/6,t=2LL*r-i64(a)-p;
                if(q>=qlo&&q<=qmax&&q>p&&t>=0&&t<q&&index[static_cast<std::size_t>(q)]>=0) {
                    const i64 d=q-p;
                    if(i64(r)!=t+6*d||i64(a)!=t+13*d-q||!(q-13*d<=t&&t<q-7*d))
                        throw std::runtime_error("plus inverse self-check failed");
                    h.push_back({static_cast<u32>(q),static_cast<u32>(t),p,r,a,'+'});
                    if(st)++st->plus_prime;
                }
            }
            const i64 mn=i64(r)+i64(a)+1;
            if(mn%7==0) {
                if(st)++st->minus_integral;
                const i64 d=mn/7,q=i64(p)+d,t=d-i64(a)-1;
                if(q>=qlo&&q<=qmax&&q>p&&t>=0&&t<q&&index[static_cast<std::size_t>(q)]>=0) {
                    if(i64(r)!=t+6*d||i64(a)!=d-t-1||!(0<=t&&t<=d-1&&7*d+t<q))
                        throw std::runtime_error("minus inverse self-check failed");
                    h.push_back({static_cast<u32>(q),static_cast<u32>(t),p,r,a,'-'});
                    if(st)++st->minus_prime;
                }
            }
        }
    }
    std::sort(h.begin(),h.end(),raw_less);
    h.erase(std::unique(h.begin(),h.end(),[](const RawHit&x,const RawHit&y){
        return x.q==y.q&&x.t==y.t&&x.p==y.p&&x.r==y.r&&x.a==y.a&&x.sign==y.sign;
    }),h.end());
    return h;
}

struct Occurrence {
    u32 q=0,t=0,p=0,r=0; i64 ap=-1,am=-1; bool rawp=false,rawm=false;
    char assigned()const{return rawm?'-':'+';}
};

std::vector<Occurrence> merge_hits(const std::vector<RawHit>&h) {
    std::vector<Occurrence> out;
    for(std::size_t i=0;i<h.size();) {
        std::size_t j=i+1;while(j<h.size()&&h[j].q==h[i].q&&h[j].t==h[i].t&&h[j].p==h[i].p)++j;
        Occurrence o;o.q=h[i].q;o.t=h[i].t;o.p=h[i].p;
        bool have_r=false;
        for(std::size_t k=i;k<j;++k) {
            if(have_r&&o.r!=h[k].r)throw std::runtime_error("primary residue merge mismatch");
            o.r=h[k].r;have_r=true;
            if(h[k].sign=='+'){o.rawp=true;o.ap=h[k].a;}else{o.rawm=true;o.am=h[k].a;}
        }
        out.push_back(o);i=j;
    }
    std::sort(out.begin(),out.end(),[](const Occurrence&x,const Occurrence&y){
        return std::tie(x.q,x.t,x.p)<std::tie(y.q,y.t,y.p);
    });
    return out;
}

std::vector<RawHit> direct_hits(const std::vector<u32>&primes,
                                const std::vector<i32>&index,
                                const std::vector<std::vector<u32>>&zeros,u32 qmax) {
    std::unordered_map<u64,std::vector<u32>> row6;
    for(std::size_t i=0;i<primes.size();++i)if(primes[i]>=7)
        for(u32 r:zeros[i])row6[6ULL*primes[i]+r].push_back(primes[i]);
    std::vector<RawHit> out;
    for(u32 q:primes)if(q>=17&&q<=qmax)for(u32 t=0;t<q;++t) {
        const u64 n=6ULL*q+t;auto it=row6.find(n);if(it==row6.end())continue;
        for(u32 p:it->second) {
            if(p>=q||n>=7ULL*p)continue;const i32 pi=index[p];if(pi<0)continue;
            const u32 r=static_cast<u32>(n-6ULL*p);
            const i64 am=i64(q)-1-t-p;
            if(has_zero(zeros[static_cast<std::size_t>(pi)],am))out.push_back({q,t,p,r,static_cast<u32>(am),'-'});
            const i64 ap=12LL*q+t-13LL*p;
            if(has_zero(zeros[static_cast<std::size_t>(pi)],ap))out.push_back({q,t,p,r,static_cast<u32>(ap),'+'});
        }
    }
    std::sort(out.begin(),out.end(),raw_less);return out;
}

struct Validation {
    u64 mont_products=0,mont_fail=0,rec_checks=0,rec_fail=0;
    u64 inverse_hits=0,direct_hits=0,candidate_fail=0,states5000=0,leaves5000=0;
};

Validation validate() {
    Validation v;
    const std::array<u32,6> ps={7,11,101,1009,10007,999983};
    u64 s=0x9e3779b97f4a7c15ULL;
    auto rng=[&](){s^=s<<7;s^=s>>9;s^=s<<8;return s;};
    for(u32 p:ps){Mont32 m{p};for(int k=0;k<2000;++k){
        u32 a=static_cast<u32>(rng()%p),b=static_cast<u32>(rng()%p);
        u32 got=m.from_mont(m.mul(m.to_mont(a),m.to_mont(b))),want=static_cast<u32>(u64(a)*b%p);
        ++v.mont_products;if(got!=want)++v.mont_fail;
    }}
    if(v.mont_fail)throw std::runtime_error("Montgomery validation failed");
    const u32 lim=5000;auto primes=primes_upto(lim);
    std::vector<i32> index(lim+1,-1);std::vector<std::vector<u32>> zeros(primes.size());
    for(std::size_t i=0;i<primes.size();++i){index[primes[i]]=static_cast<i32>(i);zeros[i]=zeros_mont(primes[i]);
        if(primes[i]>=7&&primes[i]<=500){++v.rec_checks;if(zeros[i]!=zeros_full(primes[i]))++v.rec_fail;}
        if(primes[i]>=17)v.states5000+=zeros[i].size();
    }
    if(v.rec_fail||v.states5000!=605)throw std::runtime_error("recurrence/state validation failed");
    ProjectionStats st;auto inv=project_hits(primes,index,zeros,7,17,lim,&st);auto dir=direct_hits(primes,index,zeros,lim);
    v.inverse_hits=inv.size();v.direct_hits=dir.size();
    if(inv.size()!=dir.size())++v.candidate_fail;else for(std::size_t i=0;i<inv.size();++i)
        if(std::tie(inv[i].q,inv[i].t,inv[i].p,inv[i].sign,inv[i].r,inv[i].a)!=
           std::tie(dir[i].q,dir[i].t,dir[i].p,dir[i].sign,dir[i].r,dir[i].a)){++v.candidate_fail;break;}
    if(v.candidate_fail)throw std::runtime_error("candidate-first/direct validation failed");
    for(const auto&x:inv){const i32 qi=index[x.q];if(qi>=0&&has_zero(zeros[static_cast<std::size_t>(qi)],x.t))++v.leaves5000;}
    if(v.leaves5000)throw std::runtime_error("unexpected selected leaf below 5000");
    return v;
}

struct Nearest { i64 linear=-1,cyclic=-1,zero=-1; };
Nearest nearest(u32 q,u32 t,const std::vector<u32>&z) {
    Nearest a;
    for(u32 r:z){i64 d=std::llabs(i64(t)-r),c=std::min<i64>(d,i64(q)-d);
        if(a.linear<0||d<a.linear)a.linear=d;
        if(a.cyclic<0||c<a.cyclic||(c==a.cyclic&&r<a.zero)){a.cyclic=c;a.zero=r;}
    }return a;
}

double quantile(std::vector<double>v,double x){
    if(v.empty())return std::numeric_limits<double>::quiet_NaN();std::sort(v.begin(),v.end());
    double p=x*(v.size()-1);std::size_t i=static_cast<std::size_t>(std::floor(p)),j=static_cast<std::size_t>(std::ceil(p));
    return i==j?v[i]:v[i]+(p-i)*(v[j]-v[i]);
}

struct Block {
    u64 lo=0,hi=0,primes=0,qz=0,qa=0,qboth=0,totalz=0,totala=0,occ=0,actual=0,maxz=0,maxa=0;
    long double heuristic=0;std::vector<double>norm,scaled;
};

void add(Block&b,u32 q,std::size_t z,std::size_t a,std::size_t o,std::size_t actual,i64 gap){
    ++b.primes;b.totalz+=z;b.totala+=a;b.occ+=o;b.actual+=actual;b.maxz=std::max<u64>(b.maxz,z);b.maxa=std::max<u64>(b.maxa,a);
    if(z)++b.qz;if(a)++b.qa;if(z&&a){++b.qboth;b.heuristic+=static_cast<long double>(z)*a/q;
        if(gap>=0){b.norm.push_back(double(gap)/q);b.scaled.push_back(2.0*z*gap/q);}}
}

void write_blocks(const fs::path&path,const std::vector<Block>&v){
    std::ofstream o(path);o<<std::setprecision(17)
      <<"lo_exclusive,hi_inclusive,prime_count,q_with_z,q_with_a,q_with_both,total_z,total_a,assigned_occurrences,actual_occurrences,max_z,max_a,mean_z,mean_a,mean_a_times_log_mid,heuristic_sum,heuristic_times_log_mid,heuristic_times_log_mid_squared,norm_gap_q25,norm_gap_median,norm_gap_q75,scaled_gap_median\n";
    for(const Block&b:v)if(b.primes){long double mid=(b.lo+b.hi)/2.0L,L=std::log(mid),mz=(long double)b.totalz/b.primes,ma=(long double)b.totala/b.primes;
      o<<b.lo<<','<<b.hi<<','<<b.primes<<','<<b.qz<<','<<b.qa<<','<<b.qboth<<','<<b.totalz<<','<<b.totala<<','<<b.occ<<','<<b.actual<<','<<b.maxz<<','<<b.maxa<<','
       <<double(mz)<<','<<double(ma)<<','<<double(ma*L)<<','<<double(b.heuristic)<<','<<double(b.heuristic*L)<<','<<double(b.heuristic*L*L)<<','
       <<quantile(b.norm,.25)<<','<<quantile(b.norm,.5)<<','<<quantile(b.norm,.75)<<','<<quantile(b.scaled,.5)<<'\n';}
}

struct StateRow {
    u32 q=0,t=0;u64 rp=0,rm=0,ap=0,am=0;u32 repr_p=0;char repr_sign='?';Nearest near;u32 zcount=0;bool selected=false;
};

} // namespace

int main(int argc,char**argv){
  try{
    const Options opt=parse_options(argc,argv);fs::create_directories(opt.out_dir);auto start=std::chrono::steady_clock::now();
    std::cerr<<"Q4279_VALIDATION_BEGIN\n";const Validation val=validate();
    std::cerr<<"Q4279_VALIDATION_OK inverse_hits="<<val.inverse_hits<<" direct_hits="<<val.direct_hits<<'\n';

    const u32 qlo=opt.qmin_exclusive+1;
    const u32 pmin=static_cast<u32>((6ULL*qlo)/7ULL+1ULL);
    auto ts=std::chrono::steady_clock::now();const auto primes=primes_upto(opt.qmax);
    const double sieve_sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-ts).count();
    std::vector<i32> index(static_cast<std::size_t>(opt.qmax)+1,-1);
    for(std::size_t i=0;i<primes.size();++i)index[primes[i]]=static_cast<i32>(i);
    const std::size_t first=static_cast<std::size_t>(std::lower_bound(primes.begin(),primes.end(),pmin)-primes.begin());
    const std::size_t jobs=primes.size()-first;std::vector<std::vector<u32>> zeros(primes.size());
    std::atomic<std::size_t>cursor{0},done{0};std::atomic<u64>steps{0},zrecords{0};std::mutex pm;
    ts=std::chrono::steady_clock::now();
    auto worker=[&](){u64 ls=0,lz=0;for(;;){std::size_t task=cursor.fetch_add(1);if(task>=jobs)break;
        std::size_t i=primes.size()-1-task;u64 st=0;zeros[i]=zeros_mont(primes[i],&st);ls+=st;lz+=zeros[i].size();
        std::size_t d=done.fetch_add(1)+1;if(d%5000==0||d==jobs){steps.fetch_add(ls);zrecords.fetch_add(lz);ls=lz=0;
          std::lock_guard<std::mutex>g(pm);std::cerr<<"ZERO_PROGRESS "<<d<<'/'<<jobs<<" steps="<<steps.load()<<" elapsed="
          <<std::chrono::duration<double>(std::chrono::steady_clock::now()-ts).count()<<'\n';}}
        steps.fetch_add(ls);zrecords.fetch_add(lz);};
    std::vector<std::thread> th;for(u32 k=0;k<opt.threads;++k)th.emplace_back(worker);for(auto&x:th)x.join();
    const double zero_sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-ts).count();

    u64 refl=0,rfail=0,cons=0,cfail=0,elig3=0,elig4=0,bad3=0;
    for(std::size_t i=first;i<primes.size();++i){u32 p=primes[i];const auto&z=zeros[i];if(z.size()>=3)++elig3;if(z.size()>=4)++elig4;
      if(z.size()==3&&!std::binary_search(z.begin(),z.end(),(p-1)/2))++bad3;
      for(u32 r:z){++refl;if(!std::binary_search(z.begin(),z.end(),p-1-r))++rfail;}
      for(std::size_t j=1;j<z.size();++j){++cons;if(z[j]==z[j-1]+1)++cfail;}}
    if(rfail||cfail||bad3)throw std::runtime_error("zero-set structure validation failed");

    ts=std::chrono::steady_clock::now();ProjectionStats proj;auto raw=project_hits(primes,index,zeros,pmin,qlo,opt.qmax,&proj);auto occ=merge_hits(raw);
    const double projection_sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-ts).count();
    u64 rawp=0,rawm=0,overlap=0,assignedp=0,assignedm=0;
    for(const auto&o:occ){rawp+=o.rawp;rawm+=o.rawm;overlap+=o.rawp&&o.rawm;if(o.assigned()=='+')++assignedp;else++assignedm;}

    std::ofstream zo(fs::path(opt.out_dir)/"zero_records.csv"),zc(fs::path(opt.out_dir)/"prime_zero_counts.csv");
    zo<<"p,r\n";zc<<"p,zero_count\n";for(std::size_t i=first;i<primes.size();++i){zc<<primes[i]<<','<<zeros[i].size()<<'\n';for(u32 r:zeros[i])zo<<primes[i]<<','<<r<<'\n';}
    std::ofstream oo(fs::path(opt.out_dir)/"candidate_occurrences.csv");
    oo<<"assigned_sign,q,t,p,r,a_plus,a_minus,raw_plus,raw_minus,zq_count,selected,nearest_zero,linear_gap,cyclic_gap\n"<<std::setprecision(17);
    std::ofstream qs(fs::path(opt.out_dir)/"q_stats.csv");
    qs<<"q,z_count,a_union_count,a_plus_count,a_minus_count,raw_plus_occurrences,raw_minus_occurrences,raw_overlap_occurrences,assigned_occurrences,max_premark_fiber,max_selected_fiber,actual_plus_occurrences,actual_minus_occurrences,min_linear_gap,min_cyclic_gap,scaled_cyclic_gap,heuristic_a_times_z_over_q\n"<<std::setprecision(17);
    std::vector<Block>blocks((opt.qmax-opt.qmin_exclusive+opt.block_width-1)/opt.block_width);
    for(std::size_t b=0;b<blocks.size();++b){blocks[b].lo=u64(opt.qmin_exclusive)+b*opt.block_width;blocks[b].hi=std::min<u64>(opt.qmax,blocks[b].lo+opt.block_width);}
    std::vector<StateRow>states;std::vector<Occurrence>witnesses;std::size_t cur=0;u64 upper=0,qwith=0,totalA=0,maxpre=0,maxsel=0;long double heuristic_sum=0;
    auto qbegin=std::lower_bound(primes.begin(),primes.end(),qlo);
    for(auto it=qbegin;it!=primes.end()&&*it<=opt.qmax;++it){u32 q=*it;++upper;std::size_t qi=static_cast<std::size_t>(it-primes.begin());const auto&zq=zeros[qi];
      while(cur<occ.size()&&occ[cur].q<q)++cur;std::size_t b=cur;while(cur<occ.size()&&occ[cur].q==q)++cur;std::size_t e=cur;
      std::map<u32,StateRow> sm;u64 rqp=0,rqm=0,rqo=0,actp=0,actm=0;
      std::map<std::pair<u32,char>,u64> fibers;
      for(std::size_t k=b;k<e;++k){const auto&o=occ[k];bool selected=has_zero(zq,o.t);Nearest n=nearest(q,o.t,zq);
        rqp+=o.rawp;rqm+=o.rawm;rqo+=o.rawp&&o.rawm;++fibers[{o.t,o.assigned()}];if(selected){witnesses.push_back(o);if(o.assigned()=='+')++actp;else++actm;}
        auto&sr=sm[o.t];if(sr.q==0){sr.q=q;sr.t=o.t;sr.repr_p=o.p;sr.repr_sign=o.assigned();sr.near=n;sr.zcount=zq.size();sr.selected=selected;}
        sr.rp+=o.rawp;sr.rm+=o.rawm;if(o.assigned()=='+')++sr.ap;else++sr.am;
        oo<<o.assigned()<<','<<q<<','<<o.t<<','<<o.p<<','<<o.r<<','<<o.ap<<','<<o.am<<','<<o.rawp<<','<<o.rawm<<','<<zq.size()<<','<<selected<<','<<n.zero<<','<<n.linear<<','<<n.cyclic<<'\n';}
      u64 lp=0,lm=0;for(auto&[t,sr]:sm){states.push_back(sr);if(sr.ap)++lp;if(sr.am)++lm;}
      u64 localpre=0,localsel=0;for(auto&[key,c]:fibers){localpre=std::max(localpre,c);if(has_zero(zq,key.first))localsel=std::max(localsel,c);}
      maxpre=std::max(maxpre,localpre);maxsel=std::max(maxsel,localsel);if(!sm.empty())++qwith;totalA+=sm.size();
      i64 minlin=-1,mincyc=-1;for(auto&[t,sr]:sm){if(sr.near.linear>=0&&(minlin<0||sr.near.linear<minlin))minlin=sr.near.linear;if(sr.near.cyclic>=0&&(mincyc<0||sr.near.cyclic<mincyc))mincyc=sr.near.cyclic;}
      double scaled=(mincyc>=0&&!zq.empty())?2.0*zq.size()*mincyc/q:-1.0,heur=double(sm.size())*zq.size()/q;heuristic_sum+=heur;
      qs<<q<<','<<zq.size()<<','<<sm.size()<<','<<lp<<','<<lm<<','<<rqp<<','<<rqm<<','<<rqo<<','<<(e-b)<<','<<localpre<<','<<localsel<<','<<actp<<','<<actm<<','<<minlin<<','<<mincyc<<','<<scaled<<','<<heur<<'\n';
      std::size_t bi=std::min<std::size_t>(blocks.size()-1,(q-qlo)/opt.block_width);add(blocks[bi],q,zq.size(),sm.size(),e-b,actp+actm,mincyc);
      if(zq.size()==1&&zq[0]==(q-1)/2&&actm)throw std::runtime_error("singleton-minus validation failed");}
    write_blocks(fs::path(opt.out_dir)/"block_stats.csv",blocks);

    std::sort(states.begin(),states.end(),[](const StateRow&x,const StateRow&y){i64 a=x.near.cyclic<0?std::numeric_limits<i64>::max():x.near.cyclic,b=y.near.cyclic<0?std::numeric_limits<i64>::max():y.near.cyclic;
      return std::tie(a,x.q,x.t)<std::tie(b,y.q,y.t);});
    std::ofstream cs(fs::path(opt.out_dir)/"candidate_states_by_margin.csv");
    cs<<"rank,cyclic_gap,linear_gap,nearest_zero,q,t,raw_plus_count,raw_minus_count,assigned_plus_count,assigned_minus_count,representative_p,representative_sign,z_count,selected\n";
    for(std::size_t i=0;i<states.size();++i){const auto&s=states[i];cs<<i+1<<','<<s.near.cyclic<<','<<s.near.linear<<','<<s.near.zero<<','<<s.q<<','<<s.t<<','<<s.rp<<','<<s.rm<<','<<s.ap<<','<<s.am<<','<<s.repr_p<<','<<s.repr_sign<<','<<s.zcount<<','<<s.selected<<'\n';}
    std::ofstream wit(fs::path(opt.out_dir)/"witnesses.txt");if(witnesses.empty())wit<<"NO RAW MINUS-FIRST SELECTED q6 LEAF FOUND FOR "<<opt.qmin_exclusive<<" < q <= "<<opt.qmax<<"; forced p >= "<<pmin<<".\nFINITE EXHAUSTED-RANGE STATEMENT ONLY.\n";
    else{wit<<"ACTUAL RAW MINUS-FIRST SELECTED q6 LEAVES:\n";for(const auto&o:witnesses)wit<<o.assigned()<<" q="<<o.q<<" t="<<o.t<<" p="<<o.p<<" r="<<o.r<<" a_plus="<<o.ap<<" a_minus="<<o.am<<'\n';}

    double total_sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::ofstream mf(fs::path(opt.out_dir)/"manifest.json");mf<<std::setprecision(17)<<"{\n"
      <<"  \"qmin_exclusive\": "<<opt.qmin_exclusive<<",\n  \"qmax\": "<<opt.qmax<<",\n  \"forced_pmin\": "<<pmin<<",\n  \"threads\": "<<opt.threads<<",\n"
      <<"  \"prime_count_total\": "<<primes.size()<<",\n  \"lower_primes_scanned\": "<<jobs<<",\n  \"upper_primes_examined\": "<<upper<<",\n"
      <<"  \"recurrence_steps\": "<<steps.load()<<",\n  \"zero_records\": "<<zrecords.load()<<",\n  \"eligible_lower_primes_z_ge_3\": "<<elig3<<",\n  \"eligible_lower_primes_z_ge_4\": "<<elig4<<",\n"
      <<"  \"raw_plus_occurrences\": "<<rawp<<",\n  \"raw_minus_occurrences\": "<<rawm<<",\n  \"raw_sign_overlap\": "<<overlap<<",\n  \"assigned_plus_occurrences\": "<<assignedp<<",\n  \"assigned_minus_occurrences\": "<<assignedm<<",\n"
      <<"  \"candidate_upper_states\": "<<totalA<<",\n  \"q_with_candidates\": "<<qwith<<",\n  \"actual_leaf_occurrences\": "<<witnesses.size()<<",\n  \"max_premark_assigned_fiber\": "<<maxpre<<",\n  \"max_selected_assigned_fiber\": "<<maxsel<<",\n"
      <<"  \"naive_heuristic_sum\": "<<double(heuristic_sum)<<",\n  \"timing_seconds\": {\"sieve\": "<<sieve_sec<<", \"zero_phase\": "<<zero_sec<<", \"projection\": "<<projection_sec<<", \"total\": "<<total_sec<<"},\n"
      <<"  \"projection\": {\"p_seen\": "<<proj.p_seen<<", \"p_eligible\": "<<proj.p_eligible<<", \"ordered_pairs\": "<<proj.ordered_pairs<<", \"parity_pruned\": "<<proj.parity_pruned<<", \"gap_pruned\": "<<proj.gap_pruned<<", \"plus_integral\": "<<proj.plus_integral<<", \"plus_prime_hits\": "<<proj.plus_prime<<", \"minus_integral\": "<<proj.minus_integral<<", \"minus_prime_hits\": "<<proj.minus_prime<<"},\n"
      <<"  \"checks\": {\"montgomery_products\": "<<val.mont_products<<", \"montgomery_failures\": "<<val.mont_fail<<", \"full_recurrence_crosschecks\": "<<val.rec_checks<<", \"recurrence_failures\": "<<val.rec_fail<<", \"small_inverse_hits\": "<<val.inverse_hits<<", \"small_direct_hits\": "<<val.direct_hits<<", \"candidate_mismatches\": "<<val.candidate_fail<<", \"selected_states_le_5000\": "<<val.states5000<<", \"selected_leaves_le_5000\": "<<val.leaves5000<<", \"reflection_checks\": "<<refl<<", \"reflection_failures\": "<<rfail<<", \"consecutive_checks\": "<<cons<<", \"consecutive_failures\": "<<cfail<<", \"three_zero_primes_without_midpoint\": "<<bad3<<"},\n"
      <<"  \"complexity\": \"O(sum_{p in forced range} p/2 + sum_p |Z_p|^2)\"\n}\n";
    std::ofstream sm(fs::path(opt.out_dir)/"summary.md");sm<<"# Q4279 finite q6 leaf hunt\n\n- Exhausted: `"<<opt.qmin_exclusive<<" < q <= "<<opt.qmax<<"`.\n- Forced lower range: `p >= "<<pmin<<"`.\n- Cleared-recurrence steps: `"<<steps.load()<<"`.\n- Wall clock: `"<<total_sec<<"` seconds on `"<<opt.threads<<"` threads.\n- Lower-complete assigned occurrences: `"<<occ.size()<<"`.\n- Candidate upper states: `"<<totalA<<"`.\n- Actual selected leaves: `"<<witnesses.size()<<"`.\n- Maximum pre-mark assigned fibre: `"<<maxpre<<"`.\n- Naive sum `sum_q |A_q||Z_q|/q`: `"<<double(heuristic_sum)<<"`.\n\nFinite computation only; not a theorem.\n\n## Closest candidate states\n\n|rank|cyclic gap|q|t|nearest zero|representative p|sign|zcount|selected|\n|---:|---:|---:|---:|---:|---:|:---:|---:|:---:|\n";
    for(std::size_t i=0;i<std::min<std::size_t>(states.size(),20);++i){const auto&s=states[i];sm<<'|'<<i+1<<'|'<<s.near.cyclic<<'|'<<s.q<<'|'<<s.t<<'|'<<s.near.zero<<'|'<<s.repr_p<<'|'<<s.repr_sign<<'|'<<s.zcount<<'|'<<s.selected<<"|\n";}
    std::cout<<"Q4279_SCAN_COMPLETE qmin="<<opt.qmin_exclusive<<" qmax="<<opt.qmax<<" pmin="<<pmin<<" steps="<<steps.load()<<" candidates="<<occ.size()<<" leaves="<<witnesses.size()<<" max_premark_fiber="<<maxpre<<" seconds="<<total_sec<<'\n';
    return 0;
  }catch(const std::exception&e){std::cerr<<"Q4279_FATAL "<<e.what()<<'\n';return 2;}
}
