import java.io.IOException;
import java.math.BigInteger;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;
import java.util.function.Predicate;

/**
 * Post-processes the exact TSV emitted by Q4313IndependentCensus.
 * Uses only Java standard-library arithmetic and never reads prior census data.
 */
public final class Q4316CensusAnalysis {
    private static final BigInteger ONE = BigInteger.ONE;

    private Q4316CensusAnalysis() {}

    private static final class Config {
        Path records;
        Path out;
        int loExclusive = 5_000;
        int hiInclusive = 15_000;
    }

    private static final class Row {
        final int q;
        final int t;
        final char sign;
        final BigInteger g;
        final NavigableMap<Long, Integer> factors;
        final int omega;
        final boolean rawWindow;
        final boolean m5Capable;
        final int rawLeaves;
        final int assignedLeaves;
        final boolean sizeCertified;

        Row(int q, int t, char sign, BigInteger g,
            NavigableMap<Long, Integer> factors, int omega,
            boolean rawWindow, boolean m5Capable,
            int rawLeaves, int assignedLeaves, boolean sizeCertified) {
            this.q = q;
            this.t = t;
            this.sign = sign;
            this.g = g;
            this.factors = factors;
            this.omega = omega;
            this.rawWindow = rawWindow;
            this.m5Capable = m5Capable;
            this.rawLeaves = rawLeaves;
            this.assignedLeaves = assignedLeaves;
            this.sizeCertified = sizeCertified;
        }

        long largestFactor() {
            return factors.isEmpty() ? 0L : factors.lastKey();
        }
    }

    private static final class Regression {
        final double slope;
        final double r2;
        final int n;

        Regression(double slope, double r2, int n) {
            this.slope = slope;
            this.r2 = r2;
            this.n = n;
        }
    }

    private static final class Stats {
        final String name;
        final List<Row> rows;
        final int ones;
        final BigInteger median;
        final BigInteger medianPositive;
        final Row maxG;
        final int maxOmega;
        final Row maxFactor;
        final int sizeCertified;
        final int rawLeaves;
        final int assignedLeaves;
        final Regression power;

        Stats(String name, List<Row> rows, int ones, BigInteger median,
              BigInteger medianPositive, Row maxG, int maxOmega,
              Row maxFactor, int sizeCertified, int rawLeaves,
              int assignedLeaves, Regression power) {
            this.name = name;
            this.rows = rows;
            this.ones = ones;
            this.median = median;
            this.medianPositive = medianPositive;
            this.maxG = maxG;
            this.maxOmega = maxOmega;
            this.maxFactor = maxFactor;
            this.sizeCertified = sizeCertified;
            this.rawLeaves = rawLeaves;
            this.assignedLeaves = assignedLeaves;
            this.power = power;
        }
    }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.ROOT);
        Config cfg = parseArgs(args);
        List<Row> rows = readRows(cfg.records);
        String report = buildReport(rows, cfg);
        Files.createDirectories(cfg.out.getParent());
        Files.writeString(cfg.out, report, StandardCharsets.UTF_8);
        System.out.println(report);
    }

    private static Config parseArgs(String[] args) {
        Config cfg = new Config();
        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--records" -> cfg.records = Path.of(args[++i]);
                case "--out" -> cfg.out = Path.of(args[++i]);
                case "--lo-exclusive" -> cfg.loExclusive = Integer.parseInt(args[++i]);
                case "--hi-inclusive" -> cfg.hiInclusive = Integer.parseInt(args[++i]);
                default -> throw new IllegalArgumentException("unknown argument: " + args[i]);
            }
        }
        if (cfg.records == null || cfg.out == null) {
            throw new IllegalArgumentException("--records and --out are required");
        }
        if (cfg.loExclusive >= cfg.hiInclusive) {
            throw new IllegalArgumentException("invalid range");
        }
        return cfg;
    }

    private static List<Row> readRows(Path path) throws IOException {
        List<String> lines = Files.readAllLines(path, StandardCharsets.UTF_8);
        if (lines.isEmpty()) throw new IllegalStateException("empty records TSV");
        String[] header = lines.get(0).split("\\t", -1);
        Map<String, Integer> col = new HashMap<>();
        for (int i = 0; i < header.length; i++) col.put(header[i], i);
        String[] required = {"q", "t", "sign", "primitiveGcd", "factors", "omega",
            "rawWindow", "m5Capable", "rawLeaves", "assignedLeaves", "sizeCertified"};
        for (String key : required) {
            if (!col.containsKey(key)) throw new IllegalStateException("missing column " + key);
        }
        List<Row> out = new ArrayList<>();
        for (int lineNo = 1; lineNo < lines.size(); lineNo++) {
            String line = lines.get(lineNo);
            if (line.isBlank()) continue;
            String[] f = line.split("\\t", -1);
            out.add(new Row(
                Integer.parseInt(f[col.get("q")]),
                Integer.parseInt(f[col.get("t")]),
                f[col.get("sign")].charAt(0),
                new BigInteger(f[col.get("primitiveGcd")]),
                parseFactors(f[col.get("factors")]),
                Integer.parseInt(f[col.get("omega")]),
                Boolean.parseBoolean(f[col.get("rawWindow")]),
                Boolean.parseBoolean(f[col.get("m5Capable")]),
                Integer.parseInt(f[col.get("rawLeaves")]),
                Integer.parseInt(f[col.get("assignedLeaves")]),
                Boolean.parseBoolean(f[col.get("sizeCertified")])));
        }
        return out;
    }

    private static NavigableMap<Long, Integer> parseFactors(String text) {
        NavigableMap<Long, Integer> out = new TreeMap<>();
        if (text.equals("1") || text.isBlank()) return out;
        for (String token : text.split("\\*")) {
            int caret = token.indexOf('^');
            long p = Long.parseLong(caret < 0 ? token : token.substring(0, caret));
            int e = caret < 0 ? 1 : Integer.parseInt(token.substring(caret + 1));
            out.put(p, e);
        }
        return out;
    }

    private static String buildReport(List<Row> all, Config cfg) {
        StringBuilder out = new StringBuilder();
        out.append("# Q4316 extension analysis\n\n");
        out.append("Exact input: the independently generated `records.tsv`.\n\n");

        List<Row> cumulative = select(all, r -> r.q <= cfg.hiInclusive);
        List<Row> extension = select(all,
            r -> r.q > cfg.loExclusive && r.q <= cfg.hiInclusive);

        out.append("## Cumulative through q <= ").append(cfg.hiInclusive).append("\n\n");
        appendFourCategories(out, cumulative, cfg.hiInclusive / 4);

        out.append("## New band ").append(cfg.loExclusive)
            .append(" < q <= ").append(cfg.hiInclusive).append("\n\n");
        appendFourCategories(out, extension, cfg.loExclusive + 1);

        out.append("## Required point checks\n\n```text\n");
        Row q10000 = find(all, 5647, 4553, '-');
        if (q10000 != null) {
            out.append("minus state (5647,4553)   G=").append(q10000.g)
                .append(" = ").append(formatFactors(q10000.factors)).append('\n');
        } else {
            out.append("minus state (5647,4553)   MISSING\n");
        }
        Stats minus15000 = stats("minus", cumulative,
            r -> r.sign == '-', cfg.hiInclusive / 4);
        out.append("minus maximum through 15000 G=").append(minus15000.maxG.g)
            .append(" = ").append(formatFactors(minus15000.maxG.factors))
            .append(" at (").append(minus15000.maxG.q).append(',')
            .append(minus15000.maxG.t).append(")\n");
        appendLargestFactorCheck(out, "plus", stats("plus", cumulative,
            r -> r.sign == '+', cfg.hiInclusive / 4));
        appendLargestFactorCheck(out, "minus all", minus15000);
        appendLargestFactorCheck(out, "minus raw", stats("minus raw", cumulative,
            r -> r.sign == '-' && r.rawWindow, cfg.hiInclusive / 4));
        appendLargestFactorCheck(out, "minus M5", stats("minus M5", cumulative,
            r -> r.sign == '-' && r.m5Capable, cfg.hiInclusive / 4));
        out.append("```\n\n");

        out.append("## Largest retained prime factor by dyadic q block\n\n");
        appendDyadicTable(out, cumulative, cfg.hiInclusive);

        out.append("## Dyadic-maximum trend diagnostics\n\n");
        appendBlockTrend(out, "plus", cumulative, r -> r.sign == '+', cfg.hiInclusive);
        appendBlockTrend(out, "minus all", cumulative, r -> r.sign == '-', cfg.hiInclusive);
        appendBlockTrend(out, "minus raw-window", cumulative,
            r -> r.sign == '-' && r.rawWindow, cfg.hiInclusive);
        appendBlockTrend(out, "minus M5-capable", cumulative,
            r -> r.sign == '-' && r.m5Capable, cfg.hiInclusive);

        return out.toString();
    }

    private static void appendFourCategories(StringBuilder out, List<Row> source,
                                             int fitStart) {
        appendStats(out, stats("plus", source, r -> r.sign == '+', fitStart));
        appendStats(out, stats("minus all", source, r -> r.sign == '-', fitStart));
        appendStats(out, stats("minus raw-window", source,
            r -> r.sign == '-' && r.rawWindow, fitStart));
        appendStats(out, stats("minus M5-capable", source,
            r -> r.sign == '-' && r.m5Capable, fitStart));
    }

    private static Stats stats(String name, List<Row> source,
                               Predicate<Row> predicate, int fitStart) {
        List<Row> rows = select(source, predicate);
        if (rows.isEmpty()) throw new IllegalStateException("empty category " + name);
        List<BigInteger> values = new ArrayList<>();
        List<BigInteger> positive = new ArrayList<>();
        int ones = 0;
        Row maxG = rows.get(0);
        Row maxFactor = rows.get(0);
        int maxOmega = 0;
        int certified = 0;
        int rawLeaves = 0;
        int assignedLeaves = 0;
        List<Row> fit = new ArrayList<>();
        for (Row r : rows) {
            values.add(r.g);
            if (r.g.equals(ONE)) ones++; else positive.add(r.g);
            if (r.g.compareTo(maxG.g) > 0) maxG = r;
            if (r.largestFactor() > maxFactor.largestFactor()) maxFactor = r;
            maxOmega = Math.max(maxOmega, r.omega);
            if (r.sizeCertified) certified++;
            rawLeaves += r.rawLeaves;
            assignedLeaves += r.assignedLeaves;
            if (r.q >= fitStart && r.g.compareTo(ONE) > 0) fit.add(r);
        }
        Collections.sort(values);
        Collections.sort(positive);
        return new Stats(name, rows, ones, lowerMedian(values),
            positive.isEmpty() ? BigInteger.ZERO : lowerMedian(positive),
            maxG, maxOmega, maxFactor, certified, rawLeaves, assignedLeaves,
            regressionRows(fit));
    }

    private static void appendStats(StringBuilder out, Stats s) {
        out.append("### ").append(s.name).append("\n\n```text\n");
        out.append("records                 ").append(s.rows.size()).append('\n');
        out.append("G=1                     ").append(s.ones).append('\n');
        out.append("median G                ").append(s.median).append('\n');
        out.append("median nontrivial G     ").append(s.medianPositive).append('\n');
        out.append("maximum G               ").append(s.maxG.g).append(" = ")
            .append(formatFactors(s.maxG.factors)).append(" at (")
            .append(s.maxG.q).append(',').append(s.maxG.t).append(")\n");
        out.append("maximum omega           ").append(s.maxOmega).append('\n');
        out.append("largest prime factor    ").append(s.maxFactor.largestFactor())
            .append(" at (").append(s.maxFactor.q).append(',')
            .append(s.maxFactor.t).append(")\n");
        out.append(String.format(Locale.ROOT,
            "rho(log log G~log q)  %.9g  R2=%.9g  n=%d%n",
            s.power.slope, s.power.r2, s.power.n));
        out.append("size certified          ").append(s.sizeCertified)
            .append(" / ").append(s.rows.size()).append('\n');
        out.append("raw leaves              ").append(s.rawLeaves).append('\n');
        out.append("assigned leaves         ").append(s.assignedLeaves).append('\n');
        out.append("```\n\n");
    }

    private static void appendLargestFactorCheck(StringBuilder out, String label, Stats s) {
        out.append(String.format(Locale.ROOT, "%-24s %d at (%d,%d)%n",
            label + " largest factor", s.maxFactor.largestFactor(),
            s.maxFactor.q, s.maxFactor.t));
    }

    private static void appendDyadicTable(StringBuilder out, List<Row> source, int qMax) {
        out.append("| q block | plus | minus all | minus raw | minus M5 |\n");
        out.append("|---|---:|---:|---:|---:|\n");
        for (int lo = 16; lo <= qMax; lo *= 2) {
            int hi = lo * 2;
            long plus = blockMax(source, lo, hi, r -> r.sign == '+');
            long minus = blockMax(source, lo, hi, r -> r.sign == '-');
            long raw = blockMax(source, lo, hi, r -> r.sign == '-' && r.rawWindow);
            long m5 = blockMax(source, lo, hi, r -> r.sign == '-' && r.m5Capable);
            out.append('|').append('[').append(lo).append(',').append(hi).append(")|")
                .append(plus).append('|').append(minus).append('|')
                .append(raw).append('|').append(m5).append("|\n");
            if (lo > Integer.MAX_VALUE / 2) break;
        }
        out.append('\n');
    }

    private static long blockMax(List<Row> source, int lo, int hi,
                                 Predicate<Row> predicate) {
        long max = 0L;
        for (Row r : source) {
            if (r.q >= lo && r.q < hi && predicate.test(r)) {
                max = Math.max(max, r.largestFactor());
            }
        }
        return max;
    }

    private static void appendBlockTrend(StringBuilder out, String name,
                                         List<Row> source, Predicate<Row> predicate,
                                         int qMax) {
        List<Double> x = new ArrayList<>();
        List<Double> y = new ArrayList<>();
        for (int lo = 16; lo <= qMax; lo *= 2) {
            int hi = lo * 2;
            long max = blockMax(source, lo, hi, predicate);
            if (max > 1L) {
                x.add(Math.log(Math.sqrt((double) lo * hi)));
                y.add(Math.log(max));
            }
            if (lo > Integer.MAX_VALUE / 2) break;
        }
        Regression r = regressionXY(x, y);
        out.append(String.format(Locale.ROOT,
            "- %s: log(block max factor) slope %.6g, R^2 %.6g, blocks %d.%n",
            name, r.slope, r.r2, r.n));
    }

    private static Regression regressionRows(List<Row> rows) {
        List<Double> x = new ArrayList<>();
        List<Double> y = new ArrayList<>();
        for (Row r : rows) {
            double lg = logBig(r.g);
            x.add(Math.log(r.q));
            y.add(Math.log(lg));
        }
        return regressionXY(x, y);
    }

    private static Regression regressionXY(List<Double> x, List<Double> y) {
        int n = x.size();
        if (n < 2) return new Regression(Double.NaN, Double.NaN, n);
        double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
        for (int i = 0; i < n; i++) {
            double a = x.get(i), b = y.get(i);
            sx += a; sy += b; sxx += a * a; syy += b * b; sxy += a * b;
        }
        double vx = sxx - sx * sx / n;
        double vy = syy - sy * sy / n;
        double cov = sxy - sx * sy / n;
        double slope = cov / vx;
        double r2 = vx > 0.0 && vy > 0.0 ? cov * cov / (vx * vy) : 0.0;
        return new Regression(slope, r2, n);
    }

    private static Row find(List<Row> rows, int q, int t, char sign) {
        for (Row r : rows) if (r.q == q && r.t == t && r.sign == sign) return r;
        return null;
    }

    private static List<Row> select(List<Row> rows, Predicate<Row> predicate) {
        List<Row> out = new ArrayList<>();
        for (Row r : rows) if (predicate.test(r)) out.add(r);
        return out;
    }

    private static BigInteger lowerMedian(List<BigInteger> sorted) {
        return sorted.get((sorted.size() - 1) / 2);
    }

    private static String formatFactors(NavigableMap<Long, Integer> factors) {
        if (factors.isEmpty()) return "1";
        StringBuilder out = new StringBuilder();
        boolean first = true;
        for (Map.Entry<Long, Integer> e : factors.entrySet()) {
            if (!first) out.append('*');
            first = false;
            out.append(e.getKey());
            if (e.getValue() > 1) out.append('^').append(e.getValue());
        }
        return out.toString();
    }

    private static double logBig(BigInteger x) {
        int bits = x.bitLength();
        if (bits <= 62) return Math.log(x.longValueExact());
        int shift = bits - 53;
        long top = x.shiftRight(shift).longValueExact();
        return Math.log(top) + shift * Math.log(2.0);
    }
}
