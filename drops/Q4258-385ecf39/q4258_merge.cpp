#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;
namespace fs = std::filesystem;

struct Row {
    u32 q{}, t{}, p{}, x{}, y{}, h{}, zp_size{}, zq_size{};
    char sign{};
    std::string type;
    bool selected{};
    std::int64_t nearest_zero{-1}, linear_distance{-1}, cyclic_distance{-1};
};

static std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> v;
    std::size_t start = 0;
    while (true) {
        const std::size_t pos = line.find(',', start);
        if (pos == std::string::npos) { v.push_back(line.substr(start)); break; }
        v.push_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return v;
}

static Row parse_row(const std::string& line) {
    const auto v = split_csv(line);
    if (v.size() != 14) throw std::runtime_error("bad candidate CSV field count: " + std::to_string(v.size()));
    Row r;
    r.q = static_cast<u32>(std::stoul(v[0]));
    r.t = static_cast<u32>(std::stoul(v[1]));
    r.p = static_cast<u32>(std::stoul(v[2]));
    if (v[3].size() != 1 || (v[3][0] != '+' && v[3][0] != '-')) throw std::runtime_error("bad sign");
    r.sign = v[3][0];
    r.type = v[4];
    r.x = static_cast<u32>(std::stoul(v[5]));
    r.y = static_cast<u32>(std::stoul(v[6]));
    r.h = static_cast<u32>(std::stoul(v[7]));
    r.zp_size = static_cast<u32>(std::stoul(v[8]));
    r.zq_size = static_cast<u32>(std::stoul(v[9]));
    r.selected = std::stoul(v[10]) != 0;
    r.nearest_zero = std::stoll(v[11]);
    r.linear_distance = std::stoll(v[12]);
    r.cyclic_distance = std::stoll(v[13]);
    return r;
}

static bool row_less(const Row& a, const Row& b) {
    return std::tie(a.q, a.t, a.sign, a.p, a.type, a.x, a.y)
         < std::tie(b.q, b.t, b.sign, b.p, b.type, b.x, b.y);
}
static bool row_same(const Row& a, const Row& b) {
    return std::tie(a.q, a.t, a.sign, a.p)
         == std::tie(b.q, b.t, b.sign, b.p);
}

struct Options {
    fs::path input_root = ".";
    fs::path output_dir = ".";
    u32 q_min_exclusive = 1000000;
    u32 q_max = 10000000;
    u32 block_width = 1000000;
};
static Options parse_options(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&]() -> std::string {
            if (++i >= argc) throw std::runtime_error("missing value after " + a);
            return argv[i];
        };
        if (a == "--input-root") o.input_root = need();
        else if (a == "--output-dir") o.output_dir = need();
        else if (a == "--q-min-exclusive") o.q_min_exclusive = static_cast<u32>(std::stoul(need()));
        else if (a == "--q-max") o.q_max = static_cast<u32>(std::stoul(need()));
        else if (a == "--block-width") o.block_width = static_cast<u32>(std::stoul(need()));
        else throw std::runtime_error("unknown argument: " + a);
    }
    if (o.q_max <= o.q_min_exclusive || o.block_width == 0) throw std::runtime_error("bad range");
    return o;
}

static double extract_json_number(const std::string& text, const std::string& key) {
    const std::string needle = "\"" + key + "\"";
    std::size_t pos = text.find(needle);
    if (pos == std::string::npos) return 0.0;
    pos = text.find(':', pos + needle.size());
    if (pos == std::string::npos) return 0.0;
    ++pos;
    while (pos < text.size() && (text[pos] == ' ' || text[pos] == '\t')) ++pos;
    std::size_t end = pos;
    while (end < text.size() && (std::isdigit(static_cast<unsigned char>(text[end])) || text[end] == '.' || text[end] == '-' || text[end] == '+' || text[end] == 'e' || text[end] == 'E')) ++end;
    return std::stod(text.substr(pos, end - pos));
}

struct QSummary {
    u32 q{};
    u32 zq_size{};
    u64 occurrences{}, plus_occurrences{}, minus_occurrences{};
    u64 selected_occurrences{};
    std::set<u32> digits, selected_digits;
    std::int64_t min_linear{-1}, min_cyclic{-1};
    double min_normalized = std::numeric_limits<double>::infinity();
    u32 max_k_premark{}, max_k_selected{};
};

struct Block {
    u64 candidate_q{}, occurrences{}, plus_occurrences{}, minus_occurrences{};
    u64 unique_a{}, total_zq{}, selected_occurrences{}, selected_digits{};
    u64 q_zq_zero{}, q_zq_nonzero{};
    double heuristic_sum{};
    std::int64_t min_linear{-1}, min_cyclic{-1};
    double min_normalized = std::numeric_limits<double>::infinity();
    u32 max_k_premark{}, max_k_selected{};
};

static double quantile(std::vector<double> v, double q) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const double x = q * (v.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(x));
    const std::size_t hi = static_cast<std::size_t>(std::ceil(x));
    if (lo == hi) return v[lo];
    return v[lo] + (x - lo) * (v[hi] - v[lo]);
}

int main(int argc, char** argv) {
    try {
        const Options opt = parse_options(argc, argv);
        fs::create_directories(opt.output_dir);
        std::vector<fs::path> candidate_files, summary_files;
        for (const auto& entry : fs::recursive_directory_iterator(opt.input_root)) {
            if (!entry.is_regular_file()) continue;
            const std::string name = entry.path().filename().string();
            if (name.rfind("candidates-", 0) == 0 && entry.path().extension() == ".csv") candidate_files.push_back(entry.path());
            if (name.rfind("summary-", 0) == 0 && entry.path().extension() == ".json") summary_files.push_back(entry.path());
        }
        std::sort(candidate_files.begin(), candidate_files.end());
        std::sort(summary_files.begin(), summary_files.end());
        if (candidate_files.empty()) throw std::runtime_error("no candidate CSV files found");

        std::vector<Row> rows;
        for (const fs::path& path : candidate_files) {
            std::ifstream in(path);
            if (!in) throw std::runtime_error("cannot read " + path.string());
            std::string line;
            if (!std::getline(in, line) || line.rfind("q,t,p,", 0) != 0) throw std::runtime_error("bad CSV header in " + path.string());
            while (std::getline(in, line)) if (!line.empty()) rows.push_back(parse_row(line));
        }
        std::sort(rows.begin(), rows.end(), row_less);
        rows.erase(std::unique(rows.begin(), rows.end(), row_same), rows.end());

        double max_shard_wall = 0.0, sum_shard_wall = 0.0;
        u64 lower_steps = 0, upper_steps = 0, lower_primes = 0, admissible_pairs = 0;
        for (const fs::path& path : summary_files) {
            std::ifstream in(path);
            std::ostringstream ss; ss << in.rdbuf();
            const std::string text = ss.str();
            const double wall = extract_json_number(text, "wall_seconds");
            max_shard_wall = std::max(max_shard_wall, wall);
            sum_shard_wall += wall;
            lower_steps += static_cast<u64>(extract_json_number(text, "lower_recurrence_steps"));
            upper_steps += static_cast<u64>(extract_json_number(text, "upper_recurrence_steps"));
            lower_primes += static_cast<u64>(extract_json_number(text, "lower_primes_scanned"));
            admissible_pairs += static_cast<u64>(extract_json_number(text, "admissible_opposite_parity_zero_pairs"));
        }

        std::map<u32, QSummary> by_q;
        u64 total_plus = 0, total_minus = 0, selected_plus = 0, selected_minus = 0;
        for (const Row& r : rows) {
            if (r.q <= opt.q_min_exclusive || r.q > opt.q_max) throw std::runtime_error("row outside requested q range");
            QSummary& q = by_q[r.q];
            if (q.q == 0) { q.q = r.q; q.zq_size = r.zq_size; }
            else if (q.zq_size != r.zq_size) throw std::runtime_error("inconsistent Z_q size across shards");
            ++q.occurrences;
            if (r.sign == '+') { ++q.plus_occurrences; ++total_plus; }
            else { ++q.minus_occurrences; ++total_minus; }
            q.digits.insert(r.t);
            if (r.selected) {
                ++q.selected_occurrences;
                q.selected_digits.insert(r.t);
                if (r.sign == '+') ++selected_plus; else ++selected_minus;
            }
            if (r.linear_distance >= 0) {
                if (q.min_linear < 0 || r.linear_distance < q.min_linear) q.min_linear = r.linear_distance;
                if (q.min_cyclic < 0 || r.cyclic_distance < q.min_cyclic) q.min_cyclic = r.cyclic_distance;
                q.min_normalized = std::min(q.min_normalized, static_cast<double>(r.linear_distance) / r.q);
            }
        }

        for (std::size_t i = 0; i < rows.size();) {
            std::size_t j = i + 1;
            while (j < rows.size() && std::tie(rows[j].q, rows[j].t, rows[j].sign) == std::tie(rows[i].q, rows[i].t, rows[i].sign)) ++j;
            std::set<u32> ps, selected_ps;
            for (std::size_t k = i; k < j; ++k) {
                ps.insert(rows[k].p);
                if (rows[k].selected) selected_ps.insert(rows[k].p);
            }
            QSummary& q = by_q[rows[i].q];
            q.max_k_premark = std::max<u32>(q.max_k_premark, static_cast<u32>(ps.size()));
            q.max_k_selected = std::max<u32>(q.max_k_selected, static_cast<u32>(selected_ps.size()));
            i = j;
        }

        std::ofstream all(opt.output_dir / "all_candidates.csv");
        all << "q,t,p,sign,type,x,y,h,zp_size,zq_size,selected,nearest_zero,linear_distance,cyclic_distance\n";
        for (const Row& r : rows) {
            all << r.q << ',' << r.t << ',' << r.p << ',' << r.sign << ',' << r.type << ',' << r.x << ',' << r.y << ',' << r.h
                << ',' << r.zp_size << ',' << r.zq_size << ',' << (r.selected ? 1 : 0) << ',' << r.nearest_zero
                << ',' << r.linear_distance << ',' << r.cyclic_distance << '\n';
        }

        std::ofstream qcsv(opt.output_dir / "by_q.csv");
        qcsv << "q,occurrences,plus_occurrences,minus_occurrences,unique_Aq,zq_size,selected_occurrences,selected_digits,min_linear_distance,min_cyclic_distance,min_normalized_distance,heuristic_AZ_over_q,max_K_premark,max_K_selected\n";
        std::map<u64, Block> blocks;
        std::vector<double> normalized_distances, heuristic_values;
        double heuristic_total = 0.0;
        u32 global_max_k_premark = 0, global_max_k_selected = 0;
        std::int64_t global_min_linear = -1, global_min_cyclic = -1;
        double global_min_normalized = std::numeric_limits<double>::infinity();
        for (auto& [qv, q] : by_q) {
            const double heuristic = static_cast<double>(q.digits.size()) * q.zq_size / qv;
            heuristic_total += heuristic;
            heuristic_values.push_back(heuristic);
            if (std::isfinite(q.min_normalized)) normalized_distances.push_back(q.min_normalized);
            global_max_k_premark = std::max(global_max_k_premark, q.max_k_premark);
            global_max_k_selected = std::max(global_max_k_selected, q.max_k_selected);
            if (q.min_linear >= 0 && (global_min_linear < 0 || q.min_linear < global_min_linear)) global_min_linear = q.min_linear;
            if (q.min_cyclic >= 0 && (global_min_cyclic < 0 || q.min_cyclic < global_min_cyclic)) global_min_cyclic = q.min_cyclic;
            global_min_normalized = std::min(global_min_normalized, q.min_normalized);
            qcsv << qv << ',' << q.occurrences << ',' << q.plus_occurrences << ',' << q.minus_occurrences << ',' << q.digits.size()
                 << ',' << q.zq_size << ',' << q.selected_occurrences << ',' << q.selected_digits.size() << ',' << q.min_linear << ',' << q.min_cyclic << ',';
            if (std::isfinite(q.min_normalized)) qcsv << std::setprecision(17) << q.min_normalized;
            qcsv << ',' << std::setprecision(17) << heuristic << ',' << q.max_k_premark << ',' << q.max_k_selected << '\n';

            const u64 bi = (static_cast<u64>(qv) - (static_cast<u64>(opt.q_min_exclusive) + 1)) / opt.block_width;
            Block& b = blocks[bi];
            ++b.candidate_q;
            b.occurrences += q.occurrences;
            b.plus_occurrences += q.plus_occurrences;
            b.minus_occurrences += q.minus_occurrences;
            b.unique_a += q.digits.size();
            b.total_zq += q.zq_size;
            b.selected_occurrences += q.selected_occurrences;
            b.selected_digits += q.selected_digits.size();
            if (q.zq_size == 0) ++b.q_zq_zero; else ++b.q_zq_nonzero;
            b.heuristic_sum += heuristic;
            if (q.min_linear >= 0 && (b.min_linear < 0 || q.min_linear < b.min_linear)) b.min_linear = q.min_linear;
            if (q.min_cyclic >= 0 && (b.min_cyclic < 0 || q.min_cyclic < b.min_cyclic)) b.min_cyclic = q.min_cyclic;
            b.min_normalized = std::min(b.min_normalized, q.min_normalized);
            b.max_k_premark = std::max(b.max_k_premark, q.max_k_premark);
            b.max_k_selected = std::max(b.max_k_selected, q.max_k_selected);
        }

        std::ofstream bcsv(opt.output_dir / "blocks.csv");
        bcsv << "q_lo,q_hi,candidate_q,occurrences,plus_occurrences,minus_occurrences,unique_A,total_Zq,selected_occurrences,selected_digits,q_with_Zq_zero,q_with_Zq_nonzero,heuristic_sum,heuristic_per_candidate_q,min_linear_distance,min_cyclic_distance,min_normalized_distance,max_K_premark,max_K_selected\n";
        for (const auto& [bi, b] : blocks) {
            const u64 lo = static_cast<u64>(opt.q_min_exclusive) + 1 + bi * opt.block_width;
            const u64 hi = std::min<u64>(opt.q_max, lo + opt.block_width - 1);
            bcsv << lo << ',' << hi << ',' << b.candidate_q << ',' << b.occurrences << ',' << b.plus_occurrences << ',' << b.minus_occurrences
                 << ',' << b.unique_a << ',' << b.total_zq << ',' << b.selected_occurrences << ',' << b.selected_digits << ',' << b.q_zq_zero
                 << ',' << b.q_zq_nonzero << ',' << std::setprecision(17) << b.heuristic_sum << ','
                 << (b.candidate_q ? b.heuristic_sum / b.candidate_q : 0.0) << ',' << b.min_linear << ',' << b.min_cyclic << ',';
            if (std::isfinite(b.min_normalized)) bcsv << b.min_normalized;
            bcsv << ',' << b.max_k_premark << ',' << b.max_k_selected << '\n';
        }

        std::vector<Row> near;
        for (const Row& r : rows) if (!r.selected && r.zq_size > 0 && r.linear_distance >= 0) near.push_back(r);
        std::sort(near.begin(), near.end(), [](const Row& a, const Row& b) {
            const long double da = static_cast<long double>(a.linear_distance) / a.q;
            const long double db = static_cast<long double>(b.linear_distance) / b.q;
            if (da != db) return da < db;
            if (a.linear_distance != b.linear_distance) return a.linear_distance < b.linear_distance;
            return row_less(a, b);
        });
        std::ofstream ncsv(opt.output_dir / "near_misses.csv");
        ncsv << "rank,q,t,p,sign,type,zq_size,nearest_zero,linear_distance,cyclic_distance,normalized_distance,zp_size,x,y,h\n";
        for (std::size_t i = 0; i < near.size(); ++i) {
            const Row& r = near[i];
            ncsv << i + 1 << ',' << r.q << ',' << r.t << ',' << r.p << ',' << r.sign << ',' << r.type << ',' << r.zq_size
                 << ',' << r.nearest_zero << ',' << r.linear_distance << ',' << r.cyclic_distance << ','
                 << std::setprecision(17) << static_cast<double>(r.linear_distance) / r.q << ',' << r.zp_size << ',' << r.x << ',' << r.y << ',' << r.h << '\n';
        }

        std::ofstream wcsv(opt.output_dir / "witnesses.csv");
        wcsv << "q,t,p,sign,type,x,y,h,zp_size,zq_size\n";
        for (const Row& r : rows) if (r.selected) wcsv << r.q << ',' << r.t << ',' << r.p << ',' << r.sign << ',' << r.type << ',' << r.x << ',' << r.y << ',' << r.h << ',' << r.zp_size << ',' << r.zq_size << '\n';

        std::ofstream ccsv(opt.output_dir / "collisions.csv");
        ccsv << "q,t,sign,K_premark,K_selected,primes\n";
        for (std::size_t i = 0; i < rows.size();) {
            std::size_t j = i + 1;
            while (j < rows.size() && std::tie(rows[j].q, rows[j].t, rows[j].sign) == std::tie(rows[i].q, rows[i].t, rows[i].sign)) ++j;
            std::set<u32> ps, selected_ps;
            for (std::size_t k = i; k < j; ++k) { ps.insert(rows[k].p); if (rows[k].selected) selected_ps.insert(rows[k].p); }
            if (ps.size() >= 2) {
                ccsv << rows[i].q << ',' << rows[i].t << ',' << rows[i].sign << ',' << ps.size() << ',' << selected_ps.size() << ',';
                bool first = true;
                for (u32 p : ps) { if (!first) ccsv << ';'; first = false; ccsv << p; }
                ccsv << '\n';
            }
            i = j;
        }

        std::ofstream js(opt.output_dir / "summary.json");
        js << "{\n"
           << "  \"q_min_exclusive\": " << opt.q_min_exclusive << ",\n"
           << "  \"q_max\": " << opt.q_max << ",\n"
           << "  \"shard_files\": " << candidate_files.size() << ",\n"
           << "  \"summary_files\": " << summary_files.size() << ",\n"
           << "  \"candidate_rows\": " << rows.size() << ",\n"
           << "  \"candidate_upper_primes\": " << by_q.size() << ",\n"
           << "  \"assigned_plus\": " << total_plus << ",\n"
           << "  \"assigned_minus\": " << total_minus << ",\n"
           << "  \"selected_plus\": " << selected_plus << ",\n"
           << "  \"selected_minus\": " << selected_minus << ",\n"
           << "  \"first_witness_q\": ";
        u32 first_witness = 0;
        for (const Row& r : rows) if (r.selected) { first_witness = r.q; break; }
        if (first_witness) js << first_witness; else js << "null";
        js << ",\n"
           << "  \"max_K_premark\": " << global_max_k_premark << ",\n"
           << "  \"max_K_selected\": " << global_max_k_selected << ",\n"
           << "  \"heuristic_sum_AZ_over_q\": " << std::setprecision(17) << heuristic_total << ",\n"
           << "  \"heuristic_median_per_candidate_q\": " << quantile(heuristic_values, 0.5) << ",\n"
           << "  \"global_min_linear_distance\": " << global_min_linear << ",\n"
           << "  \"global_min_cyclic_distance\": " << global_min_cyclic << ",\n"
           << "  \"global_min_normalized_distance\": ";
        if (std::isfinite(global_min_normalized)) js << global_min_normalized; else js << "null";
        js << ",\n"
           << "  \"normalized_distance_q10\": " << quantile(normalized_distances, 0.1) << ",\n"
           << "  \"normalized_distance_median\": " << quantile(normalized_distances, 0.5) << ",\n"
           << "  \"normalized_distance_q90\": " << quantile(normalized_distances, 0.9) << ",\n"
           << "  \"lower_primes_scanned\": " << lower_primes << ",\n"
           << "  \"lower_recurrence_steps\": " << lower_steps << ",\n"
           << "  \"upper_recurrence_steps\": " << upper_steps << ",\n"
           << "  \"admissible_opposite_parity_zero_pairs\": " << admissible_pairs << ",\n"
           << "  \"max_shard_wall_seconds\": " << max_shard_wall << ",\n"
           << "  \"sum_shard_wall_seconds\": " << sum_shard_wall << "\n"
           << "}\n";

        std::cout << (opt.output_dir / "summary.json") << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "fatal: " << e.what() << '\n';
        return 2;
    }
}
