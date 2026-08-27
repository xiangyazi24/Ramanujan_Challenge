// Q4279 trend analyzer for the exact q6 scanner output.
// Standard-library only.  Input CSV files are unquoted and produced by
// q4279_q6_leaf_hunt.cpp.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using u64=std::uint64_t;

std::vector<std::string> split(const std::string&s){std::vector<std::string>v;std::stringstream ss(s);std::string x;while(std::getline(ss,x,','))v.push_back(x);return v;}
double qtile(std::vector<double>v,double q){if(v.empty())return std::numeric_limits<double>::quiet_NaN();std::sort(v.begin(),v.end());double p=q*(v.size()-1);std::size_t i=std::size_t(std::floor(p)),j=std::size_t(std::ceil(p));return i==j?v[i]:v[i]+(p-i)*(v[j]-v[i]);}
std::string hist_json(const std::map<u64,u64>&h){std::ostringstream o;o<<'{';bool first=true;for(auto[k,v]:h){if(!first)o<<',';first=false;o<<'"'<<k<<"\":"<<v;}o<<'}';return o.str();}

struct B {u64 lo=0,hi=0,n=0,qboth=0,totalz=0,totala=0;double meanz=0,meana=0,alog=0,hsum=0,hlog2=0,nmed=0,smed=0;};

int main(int argc,char**argv){try{
 if(argc!=4)throw std::runtime_error("usage: q4279_analyze Q_STATS BLOCK_STATS OUT_PREFIX");
 const std::string qpath=argv[1],bpath=argv[2],prefix=argv[3];
 std::ifstream qi(qpath);if(!qi)throw std::runtime_error("cannot open q_stats");std::string line;std::getline(qi,line);
 u64 N=0,sumz=0,suma=0,qz=0,qa=0,qboth=0,actual=0,maxz=0,maxa=0,maxpre=0;double hsum=0;std::map<u64,u64>hz,ha;
 std::vector<double>linear,cyclic,norm,scaled;u64 bestq=0;double bestgap=std::numeric_limits<double>::infinity();
 while(std::getline(qi,line)){auto c=split(line);if(c.size()<17)throw std::runtime_error("short q_stats row");
  u64 q=std::stoull(c[0]),z=std::stoull(c[1]),a=std::stoull(c[2]),pre=std::stoull(c[9]),ap=std::stoull(c[11]),am=std::stoull(c[12]);
  long long lg=std::stoll(c[13]),cg=std::stoll(c[14]);double sc=std::stod(c[15]),h=std::stod(c[16]);
  ++N;sumz+=z;suma+=a;qz+=z>0;qa+=a>0;qboth+=(z>0&&a>0);actual+=ap+am;maxz=std::max(maxz,z);maxa=std::max(maxa,a);maxpre=std::max(maxpre,pre);hsum+=h;++hz[z];++ha[a];
  if(z&&a&&cg>=0){linear.push_back(double(lg));cyclic.push_back(double(cg));norm.push_back(double(cg)/q);scaled.push_back(sc);if(cg<bestgap){bestgap=cg;bestq=q;}}
 }
 std::ifstream bi(bpath);if(!bi)throw std::runtime_error("cannot open block_stats");std::getline(bi,line);std::vector<B>bs;
 while(std::getline(bi,line)){auto c=split(line);if(c.size()<22)throw std::runtime_error("short block_stats row");B b;
  b.lo=std::stoull(c[0]);b.hi=std::stoull(c[1]);b.n=std::stoull(c[2]);b.qboth=std::stoull(c[5]);b.totalz=std::stoull(c[6]);b.totala=std::stoull(c[7]);b.meanz=std::stod(c[12]);b.meana=std::stod(c[13]);b.alog=std::stod(c[14]);b.hsum=std::stod(c[15]);b.hlog2=std::stod(c[17]);b.nmed=std::stod(c[19]);b.smed=std::stod(c[21]);bs.push_back(b);
 }
 auto finite_range=[&](auto f){double lo=std::numeric_limits<double>::infinity(),hi=-lo;for(const auto&b:bs){double x=f(b);if(std::isfinite(x)){lo=std::min(lo,x);hi=std::max(hi,x);}}return std::pair<double,double>(lo,hi);};
 auto ar=finite_range([](const B&b){return b.alog;});auto hr=finite_range([](const B&b){return b.hlog2;});auto zr=finite_range([](const B&b){return b.meanz;});auto sr=finite_range([](const B&b){return b.smed;});
 const std::size_t cut=bs.size()/2;u64 n1=0,n2=0,a1=0,a2=0,z1=0,z2=0;double h1=0,h2=0;
 for(std::size_t i=0;i<bs.size();++i){if(i<cut){n1+=bs[i].n;a1+=bs[i].totala;z1+=bs[i].totalz;h1+=bs[i].hsum;}else{n2+=bs[i].n;a2+=bs[i].totala;z2+=bs[i].totalz;h2+=bs[i].hsum;}}
 auto midpoint=[](const std::vector<B>&v,std::size_t i,std::size_t j){long double num=0,den=0;for(std::size_t k=i;k<j;++k){long double m=(v[k].lo+v[k].hi)/2.0L;num+=v[k].n*std::log(m);den+=v[k].n;}return std::exp(double(num/den));};
 double mid1=cut?midpoint(bs,0,cut):1,mid2=cut<bs.size()?midpoint(bs,cut,bs.size()):1;
 std::ofstream jo(prefix+".json");jo<<std::setprecision(17)<<"{\n"
  <<"  \"upper_primes\": "<<N<<",\n  \"total_z\": "<<sumz<<",\n  \"total_a\": "<<suma<<",\n  \"q_with_z\": "<<qz<<",\n  \"q_with_a\": "<<qa<<",\n  \"q_with_both\": "<<qboth<<",\n  \"actual_occurrences\": "<<actual<<",\n  \"max_z\": "<<maxz<<",\n  \"max_a\": "<<maxa<<",\n  \"max_premark_fiber\": "<<maxpre<<",\n  \"mean_z\": "<<double(sumz)/N<<",\n  \"mean_a\": "<<double(suma)/N<<",\n  \"heuristic_sum\": "<<hsum<<",\n  \"hist_z\": "<<hist_json(hz)<<",\n  \"hist_a\": "<<hist_json(ha)<<",\n"
  <<"  \"margin_count\": "<<norm.size()<<",\n  \"best_cyclic_gap\": "<<bestgap<<",\n  \"best_gap_q\": "<<bestq<<",\n  \"linear_gap_quantiles\": ["<<qtile(linear,.1)<<','<<qtile(linear,.25)<<','<<qtile(linear,.5)<<','<<qtile(linear,.75)<<','<<qtile(linear,.9)<<"],\n"
  <<"  \"cyclic_gap_quantiles\": ["<<qtile(cyclic,.1)<<','<<qtile(cyclic,.25)<<','<<qtile(cyclic,.5)<<','<<qtile(cyclic,.75)<<','<<qtile(cyclic,.9)<<"],\n"
  <<"  \"normalized_gap_quantiles\": ["<<qtile(norm,.1)<<','<<qtile(norm,.25)<<','<<qtile(norm,.5)<<','<<qtile(norm,.75)<<','<<qtile(norm,.9)<<"],\n"
  <<"  \"spacing_scaled_gap_quantiles\": ["<<qtile(scaled,.1)<<','<<qtile(scaled,.25)<<','<<qtile(scaled,.5)<<','<<qtile(scaled,.75)<<','<<qtile(scaled,.9)<<"],\n"
  <<"  \"block_ranges\": {\"mean_z\": ["<<zr.first<<','<<zr.second<<"], \"mean_a_times_log_mid\": ["<<ar.first<<','<<ar.second<<"], \"heuristic_times_log_mid_squared\": ["<<hr.first<<','<<hr.second<<"], \"scaled_margin_median\": ["<<sr.first<<','<<sr.second<<"]},\n"
  <<"  \"first_half\": {\"prime_count\": "<<n1<<", \"mean_z\": "<<double(z1)/n1<<", \"mean_a\": "<<double(a1)/n1<<", \"mean_a_times_log_mid\": "<<double(a1)/n1*std::log(mid1)<<", \"heuristic_sum\": "<<h1<<", \"heuristic_times_log_mid_squared\": "<<h1*std::log(mid1)*std::log(mid1)<<"},\n"
  <<"  \"second_half\": {\"prime_count\": "<<n2<<", \"mean_z\": "<<double(z2)/n2<<", \"mean_a\": "<<double(a2)/n2<<", \"mean_a_times_log_mid\": "<<double(a2)/n2*std::log(mid2)<<", \"heuristic_sum\": "<<h2<<", \"heuristic_times_log_mid_squared\": "<<h2*std::log(mid2)*std::log(mid2)<<"}\n}\n";
 std::ofstream md(prefix+".md");md<<std::setprecision(8)<<"# Q4279 trend analysis\n\n"
  <<"- Upper primes: `"<<N<<"`; raw candidate states: `"<<suma<<"`; actual hits: `"<<actual<<"`.\n"
  <<"- Mean `|Z_q|`: `"<<double(sumz)/N<<"`; mean raw `|A_q|`: `"<<double(suma)/N<<"`.\n"
  <<"- Max raw `|A_q|`: `"<<maxa<<"`; max pre-mark state/sign fibre: `"<<maxpre<<"`.\n"
  <<"- Naive sum `sum |A_q||Z_q|/q`: `"<<hsum<<"`.\n"
  <<"- Closest cyclic margin: `"<<bestgap<<"` at `q="<<bestq<<"`.\n"
  <<"- Normalized cyclic-margin quantiles `(10%,25%,50%,75%,90%)`: `"<<qtile(norm,.1)<<", "<<qtile(norm,.25)<<", "<<qtile(norm,.5)<<", "<<qtile(norm,.75)<<", "<<qtile(norm,.9)<<"`.\n"
  <<"- Local-spacing-scaled margin quantiles: `"<<qtile(scaled,.1)<<", "<<qtile(scaled,.25)<<", "<<qtile(scaled,.5)<<", "<<qtile(scaled,.75)<<", "<<qtile(scaled,.9)<<"`.\n\n"
  <<"## Block stability\n\n"
  <<"- Mean `|Z_q|` block range: `"<<zr.first<<"` to `"<<zr.second<<"`.\n"
  <<"- `mean |A_q| * log(mid)` block range: `"<<ar.first<<"` to `"<<ar.second<<"`.\n"
  <<"- `(sum |A_q||Z_q|/q) * log(mid)^2` block range: `"<<hr.first<<"` to `"<<hr.second<<"`.\n"
  <<"- Spacing-scaled median-margin block range: `"<<sr.first<<"` to `"<<sr.second<<"`.\n\n"
  <<"These are finite descriptive statistics.  The scaling comparisons are heuristic diagnostics, not asymptotic theorems.\n";
 return 0;
}catch(const std::exception&e){std::cerr<<"Q4279_ANALYZE_FATAL "<<e.what()<<'\n';return 2;}}
