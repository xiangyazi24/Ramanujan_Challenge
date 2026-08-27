import java.io.BufferedWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

/**
 * Independent exact census for Q4313.
 *
 * <p>This implementation was written from the defining Apéry recurrence and
 * the fixed-gcd conventions stated in the task/report. It does not read any
 * prior census records or source code. It uses only the Java standard
 * library; arbitrary precision is java.math.BigInteger.</p>
 */
public final class Q4313IndependentCensus {
    private static final BigInteger ZERO = BigInteger.ZERO;
    private static final BigInteger ONE = BigInteger.ONE;

    private static final long[] PLUS_LEDGER_LOWER = {
        2L, 3L, 5L, 7L, 17L, 315493L, 858433L, 2720027329L
    };
    private static final long[] PLUS_LEDGER_FULL = {
        2L, 3L, 5L, 7L, 17L, 315493L, 858433L, 2720027329L,
        12073365010564729L
    };
    private static final long[] MINUS_LEDGER = {
        2L, 3L, 5L, 7L, 858433L
    };

    private static final long[] MR_BASES = {
        2L, 325L, 9375L, 28178L, 450775L, 9780504L, 1795265022L
    };

    private Q4313IndependentCensus() {}

    private static final class Config {
        int minQ = 17;
        int qMax = 20_000;
        boolean geometryOnly = false;
        Path outDir = Path.of("q4313_results");
        final List<Integer> prefixes = new ArrayList<>();
    }

    private static final class State {
        final int q;
        final int t;
        final int n6;
        final int nPlus;
        final int nMinus;
        final Interval plusWindow;
        final Interval minusWindow;
        final Interval m5Window;
        BigInteger b6;
        BigInteger bPlus;
        BigInteger bMinus;

        State(int q, int t, List<Integer> primes) {
            this.q = q;
            this.t = t;
            this.n6 = 6 * q + t;
            this.nPlus = 12 * q + t;
            this.nMinus = q - 1 - t;
            this.plusWindow = plusWindow(q, t, primes);
            this.minusWindow = minusWindow(q, t, primes);
            this.m5Window = m5Window(q, t, primes);
        }
    }

    private static final class Interval {
        final int lo;
        final int hi;
        final int firstPrime;

        Interval(int lo, int hi, int firstPrime) {
            this.lo = lo;
            this.hi = hi;
            this.firstPrime = firstPrime;
        }

        boolean nonemptyPrime() {
            return firstPrime >= 0;
        }

        boolean contains(int p) {
            return lo <= p && p <= hi;
        }
    }

    private static final class Slot {
        static final int SIX = 0;
        static final int PLUS = 1;
        static final int MINUS = 2;
        final int stateIndex;
        final int kind;

        Slot(int stateIndex, int kind) {
            this.stateIndex = stateIndex;
            this.kind = kind;
        }
    }

    private static final class StripResult {
        final BigInteger value;
        final int valuation;

        StripResult(BigInteger value, int valuation) {
            this.value = value;
            this.valuation = valuation;
        }
    }

    private static final class Record {
        final State state;
        final char sign;
        final BigInteger rawGcd;
        final BigInteger qFreeGcd;
        final BigInteger primitiveGcd;
        final BigInteger alternativePrimitive;
        final int qValuation;
        final NavigableMap<Long, Integer> factors;
        final int rawLeafCount;
        final int m5LeafCount;
        final int assignedLeafCount;
        final boolean sizeCertified;
        final double thresholdRatio;

        Record(State state, char sign, BigInteger rawGcd,
               BigInteger qFreeGcd, BigInteger primitiveGcd,
               BigInteger alternativePrimitive, int qValuation,
               NavigableMap<Long, Integer> factors, int rawLeafCount,
               int m5LeafCount, int assignedLeafCount,
               boolean sizeCertified, double thresholdRatio) {
            this.state = state;
            this.sign = sign;
            this.rawGcd = rawGcd;
            this.qFreeGcd = qFreeGcd;
            this.primitiveGcd = primitiveGcd;
            this.alternativePrimitive = alternativePrimitive;
            this.qValuation = qValuation;
            this.factors = factors;
            this.rawLeafCount = rawLeafCount;
            this.m5LeafCount = m5LeafCount;
            this.assignedLeafCount = assignedLeafCount;
            this.sizeCertified = sizeCertified;
            this.thresholdRatio = thresholdRatio;
        }
    }

    private static final class CategoryStats {
        final String name;
        final int records;
        final int ones;
        final BigInteger median;
        final BigInteger medianNontrivial;
        final BigInteger maximum;
        final State maxState;
        final NavigableMap<Long, Integer> maxFactors;
        final int maxOmega;
        final long largestPrimeFactor;
        final State largestPrimeState;
        final double rho;
        final double rhoR2;
        final double linearSlope;
        final double linearR2;
        final double maxThresholdRatio;
        final State maxThresholdState;
        final int sizeCertified;

        CategoryStats(String name, int records, int ones, BigInteger median,
                      BigInteger medianNontrivial, BigInteger maximum,
                      State maxState, NavigableMap<Long, Integer> maxFactors,
                      int maxOmega, long largestPrimeFactor,
                      State largestPrimeState, double rho, double rhoR2,
                      double linearSlope, double linearR2,
                      double maxThresholdRatio, State maxThresholdState,
                      int sizeCertified) {
            this.name = name;
            this.records = records;
            this.ones = ones;
            this.median = median;
            this.medianNontrivial = medianNontrivial;
            this.maximum = maximum;
            this.maxState = maxState;
            this.maxFactors = maxFactors;
            this.maxOmega = maxOmega;
            this.largestPrimeFactor = largestPrimeFactor;
            this.largestPrimeState = largestPrimeState;
            this.rho = rho;
            this.rhoR2 = rhoR2;
            this.linearSlope = linearSlope;
            this.linearR2 = linearR2;
            this.maxThresholdRatio = maxThresholdRatio;
            this.maxThresholdState = maxThresholdState;
            this.sizeCertified = sizeCertified;
        }
    }

    private static final class Regression {
        final double slope;
        final double r2;

        Regression(double slope, double r2) {
            this.slope = slope;
            this.r2 = r2;
        }
    }

    private static final class LeafCounts {
        final int rawPlus;
        final int rawMinus;
        final int m5Minus;
        final int assignedPlus;
        final int assignedMinus;

        LeafCounts(int rawPlus, int rawMinus, int m5Minus,
                   int assignedPlus, int assignedMinus) {
            this.rawPlus = rawPlus;
            this.rawMinus = rawMinus;
            this.m5Minus = m5Minus;
            this.assignedPlus = assignedPlus;
            this.assignedMinus = assignedMinus;
        }
    }

    private static final class SizeResult {
        final boolean certified;
        final double ratio;

        SizeResult(boolean certified, double ratio) {
            this.certified = certified;
            this.ratio = ratio;
        }
    }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.ROOT);
        Config cfg = parseArgs(args);
        if (cfg.prefixes.isEmpty()) {
            if (cfg.qMax >= 10_000) cfg.prefixes.add(10_000);
            if (cfg.qMax >= 15_000) cfg.prefixes.add(15_000);
            if (!cfg.prefixes.contains(cfg.qMax)) cfg.prefixes.add(cfg.qMax);
        }
        Collections.sort(cfg.prefixes);

        Instant started = Instant.now();
        List<Integer> primes = sievePrimes(cfg.qMax);
        List<State> states = enumerateSelectedStates(cfg, primes);
        System.out.printf(Locale.ROOT,
            "GEOMETRY minQ=%d qMax=%d primes=%d selected=%d minusRaw=%d m5=%d%n",
            cfg.minQ, cfg.qMax, primes.size(), states.size(),
            countMinusRaw(states, cfg.qMax), countM5(states, cfg.qMax));
        printPrefixGeometry(states, cfg.prefixes);

        if (cfg.geometryOnly) return;

        Files.createDirectories(cfg.outDir);
        int maxIndex = attachExactAperyValues(states);
        int maxBits = states.stream()
            .mapToInt(s -> Math.max(s.b6.bitLength(),
                Math.max(s.bPlus.bitLength(), s.bMinus.bitLength())))
            .max().orElse(0);

        Map<Integer, BitSet> zeroCache = new HashMap<>();
        List<Record> plusRecords = new ArrayList<>(states.size());
        List<Record> minusRecords = new ArrayList<>(states.size());
        int plusLedgerDifferences = 0;
        int plusOneCopyDifferences = 0;
        int minusOneCopyDifferences = 0;

        for (int i = 0; i < states.size(); i++) {
            State s = states.get(i);
            BigInteger qBig = BigInteger.valueOf(s.q);

            BigInteger rawPlus = s.b6.gcd(s.bPlus);
            StripResult plusQ = stripAll(rawPlus, qBig);
            if (plusQ.valuation < 1) {
                throw new IllegalStateException("selected q missing from plus gcd");
            }
            BigInteger plusLower = stripLedger(plusQ.value, PLUS_LEDGER_LOWER);
            BigInteger plusFull = stripLedger(plusQ.value, PLUS_LEDGER_FULL);
            if (!plusLower.equals(plusFull)) plusLedgerDifferences++;
            if (!rawPlus.divide(qBig).equals(plusQ.value)) plusOneCopyDifferences++;
            NavigableMap<Long, Integer> plusFactors = factorExactLong(plusFull);

            BigInteger rawMinus = s.b6.gcd(s.bMinus);
            StripResult minusQ = stripAll(rawMinus, qBig);
            if (minusQ.valuation < 1) {
                throw new IllegalStateException("selected q missing from minus gcd");
            }
            BigInteger minusPrimitive = stripLedger(minusQ.value, MINUS_LEDGER);
            if (!rawMinus.divide(qBig).equals(minusQ.value)) minusOneCopyDifferences++;
            NavigableMap<Long, Integer> minusFactors = factorExactLong(minusPrimitive);

            LeafCounts leaves = auditLeaves(s, plusFactors, minusFactors, zeroCache);
            SizeResult plusSize = sizeResult(plusFull, s.plusWindow);
            SizeResult minusSize = sizeResult(minusPrimitive, s.minusWindow);

            plusRecords.add(new Record(s, '+', rawPlus, plusQ.value, plusFull,
                plusLower, plusQ.valuation, plusFactors, leaves.rawPlus, 0,
                leaves.assignedPlus, plusSize.certified, plusSize.ratio));
            minusRecords.add(new Record(s, '-', rawMinus, minusQ.value,
                minusPrimitive, minusPrimitive, minusQ.valuation, minusFactors,
                leaves.rawMinus, leaves.m5Minus, leaves.assignedMinus,
                minusSize.certified, minusSize.ratio));

            if ((i + 1) % 250 == 0 || i + 1 == states.size()) {
                System.out.printf(Locale.ROOT,
                    "GCD progress %d/%d elapsed=%s%n", i + 1, states.size(),
                    formatDuration(Duration.between(started, Instant.now())));
            }
        }

        writeRecords(cfg.outDir.resolve("records.tsv"), plusRecords, minusRecords);
        writeFactorRecords(cfg.outDir.resolve("factors.tsv"), plusRecords, minusRecords);
        Duration elapsed = Duration.between(started, Instant.now());
        String report = buildReport(cfg, states, plusRecords, minusRecords,
            primes.size(), maxIndex, maxBits, plusLedgerDifferences,
            plusOneCopyDifferences, minusOneCopyDifferences, elapsed);
        Files.writeString(cfg.outDir.resolve("report.md"), report,
            StandardCharsets.UTF_8);
        Files.writeString(cfg.outDir.resolve("summary.json"),
            buildJsonSummary(cfg, states, plusRecords, minusRecords,
                primes.size(), maxIndex, maxBits, elapsed), StandardCharsets.UTF_8);
        System.out.println(report);
    }

    private static Config parseArgs(String[] args) {
        Config cfg = new Config();
        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--min-q" -> cfg.minQ = Integer.parseInt(args[++i]);
                case "--qmax" -> cfg.qMax = Integer.parseInt(args[++i]);
                case "--out" -> cfg.outDir = Path.of(args[++i]);
                case "--geometry-only" -> cfg.geometryOnly = true;
                case "--prefixes" -> {
                    for (String x : args[++i].split(",")) {
                        if (!x.isBlank()) cfg.prefixes.add(Integer.parseInt(x));
                    }
                }
                default -> throw new IllegalArgumentException("unknown argument: " + args[i]);
            }
        }
        if (cfg.minQ < 2 || cfg.qMax < cfg.minQ) {
            throw new IllegalArgumentException("invalid q range");
        }
        return cfg;
    }

    private static List<Integer> sievePrimes(int limit) {
        boolean[] composite = new boolean[limit + 1];
        for (int p = 2; (long) p * p <= limit; p++) {
            if (!composite[p]) {
                for (int m = p * p; m <= limit; m += p) composite[m] = true;
            }
        }
        List<Integer> primes = new ArrayList<>();
        for (int n = 2; n <= limit; n++) if (!composite[n]) primes.add(n);
        return primes;
    }

    private static List<State> enumerateSelectedStates(Config cfg, List<Integer> primes) {
        List<State> states = new ArrayList<>();
        for (int q : primes) {
            if (q < cfg.minQ) continue;
            int[] residues = aperyResiduesModPrime(q);
            for (int t = 0; t < q; t++) {
                if (residues[t] == 0) states.add(new State(q, t, primes));
                if (residues[t] != residues[q - 1 - t]) {
                    throw new IllegalStateException("reflection failed at q=" + q + " t=" + t);
                }
            }
        }
        return states;
    }

    private static int[] aperyResiduesModPrime(int q) {
        int[] b = new int[q];
        b[0] = 1 % q;
        if (q > 1) b[1] = 5 % q;
        long[] inv = new long[q];
        inv[1] = 1L;
        for (int i = 2; i < q; i++) {
            inv[i] = (q - ((q / i) * inv[q % i]) % q) % q;
        }
        for (int n = 1; n <= q - 2; n++) {
            long coeff = floorMod(pCoeff(n), q);
            long n3 = floorMod(cube(n), q);
            long numerator = floorMod(coeff * b[n] - n3 * b[n - 1], q);
            long invDen = mulModSmall(mulModSmall(inv[n + 1], inv[n + 1], q),
                inv[n + 1], q);
            b[n + 1] = (int) mulModSmall(numerator, invDen, q);
        }
        return b;
    }

    private static long mulModSmall(long a, long b, long mod) {
        return (a * b) % mod;
    }

    private static long floorMod(long a, long mod) {
        long r = a % mod;
        return r < 0 ? r + mod : r;
    }

    private static long cube(long n) {
        return n * n * n;
    }

    private static long pCoeff(long n) {
        return (2L * n + 1L) * (17L * n * (n + 1L) + 5L);
    }

    private static Interval plusWindow(int q, int t, List<Integer> primes) {
        int lo = max(17,
            ceilDiv(6L * q + t + 1L, 7L),
            ceilDiv(12L * q + t + 1L, 14L));
        int hi = min(q - 1, (int) ((12L * q + t) / 13L));
        return interval(lo, hi, primes);
    }

    private static Interval minusWindow(int q, int t, List<Integer> primes) {
        int lo = max(17,
            ceilDiv(6L * q + t + 1L, 7L),
            ceilDiv(q - (long) t, 2L));
        int hi = min(q - 1, q - 1 - t);
        return interval(lo, hi, primes);
    }

    private static Interval m5Window(int q, int t, List<Integer> primes) {
        Interval raw = minusWindow(q, t, primes);
        int lo = max(raw.lo,
            ceilDiv(12L * q + 2L * t + 1L, 13L),
            ceilDiv(2L * q - 2L * t - 1L, 3L));
        return interval(lo, raw.hi, primes);
    }

    private static Interval interval(int lo, int hi, List<Integer> primes) {
        return new Interval(lo, hi, firstPrime(primes, lo, hi));
    }

    private static int firstPrime(List<Integer> primes, int lo, int hi) {
        if (lo > hi) return -1;
        int pos = Collections.binarySearch(primes, lo);
        if (pos < 0) pos = -pos - 1;
        if (pos >= primes.size()) return -1;
        int p = primes.get(pos);
        return p <= hi ? p : -1;
    }

    private static int ceilDiv(long a, long b) {
        if (a < 0 || b <= 0) throw new IllegalArgumentException("ceilDiv domain");
        return Math.toIntExact((a + b - 1L) / b);
    }

    private static int max(int... x) {
        int m = Integer.MIN_VALUE;
        for (int v : x) if (v > m) m = v;
        return m;
    }

    private static int min(int... x) {
        int m = Integer.MAX_VALUE;
        for (int v : x) if (v < m) m = v;
        return m;
    }

    private static int countMinusRaw(List<State> states, int cutoff) {
        int n = 0;
        for (State s : states) if (s.q <= cutoff && s.minusWindow.nonemptyPrime()) n++;
        return n;
    }

    private static int countM5(List<State> states, int cutoff) {
        int n = 0;
        for (State s : states) if (s.q <= cutoff && s.m5Window.nonemptyPrime()) n++;
        return n;
    }

    private static void printPrefixGeometry(List<State> states, List<Integer> prefixes) {
        for (int cutoff : prefixes) {
            int selected = 0;
            for (State s : states) if (s.q <= cutoff) selected++;
            System.out.printf(Locale.ROOT,
                "PREFIX q<=%d selected=%d signs=%d minusRaw=%d m5=%d%n",
                cutoff, selected, 2 * selected,
                countMinusRaw(states, cutoff), countM5(states, cutoff));
        }
    }

    private static int attachExactAperyValues(List<State> states) {
        Map<Integer, List<Slot>> targets = new HashMap<>();
        int maxIndex = 0;
        for (int i = 0; i < states.size(); i++) {
            State s = states.get(i);
            addSlot(targets, s.n6, new Slot(i, Slot.SIX));
            addSlot(targets, s.nPlus, new Slot(i, Slot.PLUS));
            addSlot(targets, s.nMinus, new Slot(i, Slot.MINUS));
            maxIndex = Math.max(maxIndex, Math.max(s.nPlus, Math.max(s.n6, s.nMinus)));
        }
        BigInteger bPrev = ONE;
        assignTarget(0, bPrev, targets, states);
        if (maxIndex == 0) return 0;
        BigInteger bCur = BigInteger.valueOf(5L);
        assignTarget(1, bCur, targets, states);
        Instant start = Instant.now();
        for (int n = 1; n < maxIndex; n++) {
            BigInteger numerator = bCur.multiply(BigInteger.valueOf(pCoeff(n)))
                .subtract(bPrev.multiply(BigInteger.valueOf(cube(n))));
            BigInteger[] qr = numerator.divideAndRemainder(BigInteger.valueOf(cube(n + 1L)));
            if (qr[1].signum() != 0) {
                throw new IllegalStateException("nonexact recurrence at n=" + n);
            }
            BigInteger bNext = qr[0];
            assignTarget(n + 1, bNext, targets, states);
            bPrev = bCur;
            bCur = bNext;
            if ((n + 1) % 10_000 == 0 || n + 1 == maxIndex) {
                System.out.printf(Locale.ROOT,
                    "APERY index=%d/%d bits=%d elapsed=%s%n",
                    n + 1, maxIndex, bNext.bitLength(),
                    formatDuration(Duration.between(start, Instant.now())));
            }
        }
        for (State s : states) {
            if (s.b6 == null || s.bPlus == null || s.bMinus == null) {
                throw new IllegalStateException("unfilled target q=" + s.q + " t=" + s.t);
            }
        }
        return maxIndex;
    }

    private static void addSlot(Map<Integer, List<Slot>> targets, int index, Slot slot) {
        targets.computeIfAbsent(index, ignored -> new ArrayList<>()).add(slot);
    }

    private static void assignTarget(int index, BigInteger value,
                                     Map<Integer, List<Slot>> targets,
                                     List<State> states) {
        List<Slot> slots = targets.get(index);
        if (slots == null) return;
        for (Slot slot : slots) {
            State s = states.get(slot.stateIndex);
            switch (slot.kind) {
                case Slot.SIX -> s.b6 = value;
                case Slot.PLUS -> s.bPlus = value;
                case Slot.MINUS -> s.bMinus = value;
                default -> throw new IllegalStateException("bad slot kind");
            }
        }
    }

    private static StripResult stripAll(BigInteger n, BigInteger p) {
        int v = 0;
        BigInteger x = n;
        while (x.mod(p).signum() == 0) {
            x = x.divide(p);
            v++;
        }
        return new StripResult(x, v);
    }

    private static BigInteger stripLedger(BigInteger n, long[] ledger) {
        BigInteger x = n;
        for (long p : ledger) {
            BigInteger bp = BigInteger.valueOf(p);
            while (x.mod(bp).signum() == 0) x = x.divide(bp);
        }
        return x;
    }

    private static NavigableMap<Long, Integer> factorExactLong(BigInteger n) {
        if (n.signum() <= 0) throw new IllegalArgumentException("factor nonpositive");
        if (n.bitLength() > 63) {
            throw new IllegalStateException(
                "primitive gcd exceeds signed-long exact factor backend: bits=" + n.bitLength());
        }
        long x = n.longValueExact();
        NavigableMap<Long, Integer> factors = new TreeMap<>();
        factorRec(x, factors);
        long check = 1L;
        for (Map.Entry<Long, Integer> e : factors.entrySet()) {
            for (int i = 0; i < e.getValue(); i++) {
                check = Math.multiplyExact(check, e.getKey());
            }
        }
        if (check != x) throw new IllegalStateException("factor reconstruction failed");
        return factors;
    }

    private static void factorRec(long n, NavigableMap<Long, Integer> out) {
        if (n == 1L) return;
        if (isPrime64(n)) {
            out.merge(n, 1, Integer::sum);
            return;
        }
        long d = pollardRho(n);
        factorRec(d, out);
        factorRec(n / d, out);
    }

    private static boolean isPrime64(long n) {
        if (n < 2L) return false;
        int[] small = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
        for (int p : small) {
            if (n == p) return true;
            if (n % p == 0L) return false;
        }
        long d = n - 1L;
        int s = Long.numberOfTrailingZeros(d);
        d >>= s;
        BigInteger mod = BigInteger.valueOf(n);
        BigInteger nMinusOne = mod.subtract(ONE);
        for (long base : MR_BASES) {
            long a0 = Math.floorMod(base, n);
            if (a0 == 0L) continue;
            BigInteger x = BigInteger.valueOf(a0).modPow(BigInteger.valueOf(d), mod);
            if (x.equals(ONE) || x.equals(nMinusOne)) continue;
            boolean composite = true;
            for (int r = 1; r < s; r++) {
                x = x.multiply(x).mod(mod);
                if (x.equals(nMinusOne)) {
                    composite = false;
                    break;
                }
            }
            if (composite) return false;
        }
        return true;
    }

    private static long pollardRho(long n) {
        if ((n & 1L) == 0L) return 2L;
        if (n % 3L == 0L) return 3L;
        BigInteger mod = BigInteger.valueOf(n);
        for (long c = 1L; ; c++) {
            long x = 2L;
            long y = 2L;
            long d = 1L;
            int steps = 0;
            while (d == 1L && steps < 2_000_000) {
                x = rhoStep(x, c, mod);
                y = rhoStep(rhoStep(y, c, mod), c, mod);
                d = gcdLong(Math.abs(x - y), n);
                steps++;
            }
            if (d > 1L && d < n) return d;
        }
    }

    private static long rhoStep(long x, long c, BigInteger mod) {
        return BigInteger.valueOf(x).multiply(BigInteger.valueOf(x))
            .add(BigInteger.valueOf(c)).mod(mod).longValueExact();
    }

    private static long gcdLong(long a, long b) {
        while (b != 0L) {
            long t = a % b;
            a = b;
            b = t;
        }
        return a;
    }

    private static LeafCounts auditLeaves(State s,
                                          NavigableMap<Long, Integer> plusFactors,
                                          NavigableMap<Long, Integer> minusFactors,
                                          Map<Integer, BitSet> zeroCache) {
        Map<Integer, Boolean> plusLeaf = new HashMap<>();
        Map<Integer, Boolean> minusLeaf = new HashMap<>();
        for (long lp : plusFactors.keySet()) {
            if (lp <= Integer.MAX_VALUE) {
                int p = (int) lp;
                plusLeaf.put(p, rawLeaf(s, p, '+', zeroCache));
            }
        }
        for (long lp : minusFactors.keySet()) {
            if (lp <= Integer.MAX_VALUE) {
                int p = (int) lp;
                minusLeaf.put(p, rawLeaf(s, p, '-', zeroCache));
            }
        }
        int rawPlus = 0;
        int rawMinus = 0;
        int m5Minus = 0;
        int assignedPlus = 0;
        int assignedMinus = 0;
        for (Map.Entry<Integer, Boolean> e : plusLeaf.entrySet()) {
            int p = e.getKey();
            if (Boolean.TRUE.equals(e.getValue())) {
                rawPlus++;
                if (!Boolean.TRUE.equals(minusLeaf.get(p))) assignedPlus++;
            }
        }
        for (Map.Entry<Integer, Boolean> e : minusLeaf.entrySet()) {
            int p = e.getKey();
            if (Boolean.TRUE.equals(e.getValue())) {
                rawMinus++;
                assignedMinus++;
                if (s.m5Window.contains(p)) m5Minus++;
            }
        }
        return new LeafCounts(rawPlus, rawMinus, m5Minus, assignedPlus, assignedMinus);
    }

    private static boolean rawLeaf(State s, int p, char sign,
                                   Map<Integer, BitSet> zeroCache) {
        if (p < 17 || p >= s.q || !isPrime64(p)) return false;
        Interval window = sign == '+' ? s.plusWindow : s.minusWindow;
        if (!window.contains(p)) return false;
        int r6 = s.n6 - 6 * p;
        int rOther = sign == '+' ? s.nPlus - 13 * p : s.nMinus - p;
        if (r6 < 0 || r6 >= p || rOther < 0 || rOther >= p) return false;
        BitSet zeros = zeroCache.computeIfAbsent(p, Q4313IndependentCensus::zeroSet);
        return zeros.get(r6) && zeros.get(rOther);
    }

    private static BitSet zeroSet(int p) {
        int[] residues = aperyResiduesModPrime(p);
        BitSet z = new BitSet(p);
        for (int i = 0; i < p; i++) if (residues[i] == 0) z.set(i);
        return z;
    }

    private static SizeResult sizeResult(BigInteger g, Interval window) {
        if (!window.nonemptyPrime()) return new SizeResult(true, Double.NaN);
        BigInteger p = BigInteger.valueOf(window.firstPrime);
        boolean certified = g.compareTo(p.multiply(p)) < 0;
        double ratio = logBig(g) / (2.0 * Math.log(window.firstPrime));
        return new SizeResult(certified, ratio);
    }

    private static CategoryStats stats(String name, List<Record> source,
                                       int cutoff,
                                       java.util.function.Predicate<Record> filter) {
        List<Record> records = new ArrayList<>();
        for (Record r : source) {
            if (r.state.q <= cutoff && filter.test(r)) records.add(r);
        }
        if (records.isEmpty()) {
            return new CategoryStats(name, 0, 0, ZERO, ZERO, ZERO, null,
                new TreeMap<>(), 0, 0L, null, Double.NaN, Double.NaN,
                Double.NaN, Double.NaN, Double.NaN, null, 0);
        }
        List<BigInteger> values = new ArrayList<>(records.size());
        List<BigInteger> positive = new ArrayList<>();
        int ones = 0;
        Record maxRecord = records.get(0);
        int maxOmega = 0;
        long largestPrime = 0L;
        State largestPrimeState = null;
        int certified = 0;
        double maxRatio = Double.NEGATIVE_INFINITY;
        State maxRatioState = null;
        for (Record r : records) {
            values.add(r.primitiveGcd);
            if (r.primitiveGcd.equals(ONE)) ones++;
            else positive.add(r.primitiveGcd);
            if (r.primitiveGcd.compareTo(maxRecord.primitiveGcd) > 0) maxRecord = r;
            maxOmega = Math.max(maxOmega, r.factors.size());
            if (!r.factors.isEmpty() && r.factors.lastKey() > largestPrime) {
                largestPrime = r.factors.lastKey();
                largestPrimeState = r.state;
            }
            if (r.sizeCertified) certified++;
            if (Double.isFinite(r.thresholdRatio) && r.thresholdRatio > maxRatio) {
                maxRatio = r.thresholdRatio;
                maxRatioState = r.state;
            }
        }
        Collections.sort(values);
        Collections.sort(positive);
        int fitStart = Math.max(17, cutoff / 4);
        List<Record> fit = new ArrayList<>();
        for (Record r : records) {
            if (r.state.q >= fitStart && r.primitiveGcd.compareTo(ONE) > 0) fit.add(r);
        }
        Regression power = regression(fit, true);
        Regression linear = regression(fit, false);
        return new CategoryStats(name, records.size(), ones, medianValue(values),
            positive.isEmpty() ? ZERO : medianValue(positive),
            maxRecord.primitiveGcd, maxRecord.state, maxRecord.factors,
            maxOmega, largestPrime, largestPrimeState, power.slope, power.r2,
            linear.slope, linear.r2,
            maxRatio == Double.NEGATIVE_INFINITY ? Double.NaN : maxRatio,
            maxRatioState, certified);
    }

    private static BigInteger medianValue(List<BigInteger> sorted) {
        if (sorted.isEmpty()) return ZERO;
        int n = sorted.size();
        return sorted.get((n - 1) / 2);
    }

    private static Regression regression(List<Record> records, boolean power) {
        int n = records.size();
        if (n < 2) return new Regression(Double.NaN, Double.NaN);
        double sx = 0.0;
        double sy = 0.0;
        double sxx = 0.0;
        double syy = 0.0;
        double sxy = 0.0;
        for (Record r : records) {
            double logG = logBig(r.primitiveGcd);
            double x = power ? Math.log(r.state.q) : r.state.q;
            double y = power ? Math.log(logG) : logG;
            sx += x;
            sy += y;
            sxx += x * x;
            syy += y * y;
            sxy += x * y;
        }
        double vx = sxx - sx * sx / n;
        double vy = syy - sy * sy / n;
        double cov = sxy - sx * sy / n;
        double slope = cov / vx;
        double r2 = vx > 0.0 && vy > 0.0 ? cov * cov / (vx * vy) : 0.0;
        return new Regression(slope, r2);
    }

    private static double logBig(BigInteger x) {
        if (x.signum() <= 0) throw new IllegalArgumentException("log nonpositive");
        int bits = x.bitLength();
        if (bits <= 62) return Math.log(x.longValueExact());
        int shift = bits - 53;
        long top = x.shiftRight(shift).longValueExact();
        return Math.log(top) + shift * Math.log(2.0);
    }

    private static String buildReport(Config cfg, List<State> states,
                                      List<Record> plus, List<Record> minus,
                                      int primeCount, int maxIndex, int maxBits,
                                      int plusLedgerDifferences,
                                      int plusOneCopyDifferences,
                                      int minusOneCopyDifferences,
                                      Duration duration) {
        StringBuilder out = new StringBuilder();
        out.append("# Q4313 independent fixed-gcd census\n\n");
        out.append("Implementation: independent Java standard-library exact arithmetic.\n\n");
        out.append("```text\n");
        out.append("min q                   ").append(cfg.minQ).append('\n');
        out.append("max q                   ").append(cfg.qMax).append('\n');
        out.append("primes through qmax     ").append(primeCount).append('\n');
        out.append("selected states         ").append(states.size()).append('\n');
        out.append("sign records            ").append(2 * states.size()).append('\n');
        out.append("largest Apery index     ").append(maxIndex).append('\n');
        out.append("largest bit length      ").append(maxBits).append('\n');
        out.append("elapsed                  ").append(formatDuration(duration)).append('\n');
        out.append("plus ledger differences ").append(plusLedgerDifferences).append('\n');
        out.append("plus extra q powers     ").append(plusOneCopyDifferences).append('\n');
        out.append("minus extra q powers    ").append(minusOneCopyDifferences).append('\n');
        out.append("```\n\n");
        for (int cutoff : cfg.prefixes) {
            int selected = 0;
            for (State s : states) if (s.q <= cutoff) selected++;
            out.append("## Prefix q <= ").append(cutoff).append("\n\n");
            out.append("Selected states: ").append(selected)
                .append("; sign records: ").append(2 * selected)
                .append("; minus raw-window states: ")
                .append(countMinusRaw(states, cutoff))
                .append("; M5-capable states: ").append(countM5(states, cutoff))
                .append(".\n\n");
            appendStats(out, stats("plus", plus, cutoff, ignored -> true));
            appendStats(out, stats("minus", minus, cutoff, ignored -> true));
            appendStats(out, stats("minus raw-window", minus, cutoff,
                r -> r.state.minusWindow.nonemptyPrime()));
            appendStats(out, stats("minus M5-capable", minus, cutoff,
                r -> r.state.m5Window.nonemptyPrime()));
        }
        int rawPlus = plus.stream().mapToInt(r -> r.rawLeafCount).sum();
        int rawMinus = minus.stream().mapToInt(r -> r.rawLeafCount).sum();
        int assignedPlus = plus.stream().mapToInt(r -> r.assignedLeafCount).sum();
        int assignedMinus = minus.stream().mapToInt(r -> r.assignedLeafCount).sum();
        int m5Leaves = minus.stream().mapToInt(r -> r.m5LeafCount).sum();
        out.append("## Direct raw-leaf audit\n\n```text\n");
        out.append("raw plus leaves        ").append(rawPlus).append('\n');
        out.append("raw minus leaves       ").append(rawMinus).append('\n');
        out.append("assigned plus leaves   ").append(assignedPlus).append('\n');
        out.append("assigned minus leaves  ").append(assignedMinus).append('\n');
        out.append("raw minus M5 leaves    ").append(m5Leaves).append('\n');
        out.append("```\n\n");
        return out.toString();
    }

    private static void appendStats(StringBuilder out, CategoryStats s) {
        out.append("### ").append(s.name).append("\n\n```text\n");
        out.append("records                 ").append(s.records).append('\n');
        out.append("G=1                     ").append(s.ones).append('\n');
        out.append("median G                ").append(s.median).append('\n');
        out.append("median nontrivial G     ").append(s.medianNontrivial).append('\n');
        out.append("maximum G               ").append(s.maximum)
            .append(" = ").append(formatFactors(s.maxFactors)).append('\n');
        if (s.maxState != null) {
            out.append("maximum state           (").append(s.maxState.q)
                .append(',').append(s.maxState.t).append(")\n");
        }
        out.append("maximum omega           ").append(s.maxOmega).append('\n');
        out.append("largest prime factor    ").append(s.largestPrimeFactor);
        if (s.largestPrimeState != null) {
            out.append(" at (").append(s.largestPrimeState.q).append(',')
                .append(s.largestPrimeState.t).append(')');
        }
        out.append('\n');
        out.append(String.format(Locale.ROOT,
            "rho(log log G~log q)  %.9g  R2=%.9g%n", s.rho, s.rhoR2));
        out.append(String.format(Locale.ROOT,
            "slope(log G~q)        %.9g  R2=%.9g%n",
            s.linearSlope, s.linearR2));
        out.append("size certified          ").append(s.sizeCertified)
            .append(" / ").append(s.records).append('\n');
        if (Double.isFinite(s.maxThresholdRatio)) {
            out.append(String.format(Locale.ROOT,
                "max logG/(2log pmin)   %.12f", s.maxThresholdRatio));
            if (s.maxThresholdState != null) {
                out.append(" at (").append(s.maxThresholdState.q).append(',')
                    .append(s.maxThresholdState.t).append(')');
            }
            out.append('\n');
        }
        out.append("```\n\n");
    }

    private static String buildJsonSummary(Config cfg, List<State> states,
                                           List<Record> plus, List<Record> minus,
                                           int primeCount, int maxIndex,
                                           int maxBits, Duration duration) {
        return "{\n"
            + "  \"min_q\": " + cfg.minQ + ",\n"
            + "  \"qmax\": " + cfg.qMax + ",\n"
            + "  \"prime_count\": " + primeCount + ",\n"
            + "  \"selected_states\": " + states.size() + ",\n"
            + "  \"sign_records\": " + 2 * states.size() + ",\n"
            + "  \"max_index\": " + maxIndex + ",\n"
            + "  \"max_bits\": " + maxBits + ",\n"
            + "  \"elapsed_seconds\": "
            + String.format(Locale.ROOT, "%.3f", duration.toMillis() / 1000.0) + ",\n"
            + "  \"raw_plus_leaves\": "
            + plus.stream().mapToInt(r -> r.rawLeafCount).sum() + ",\n"
            + "  \"raw_minus_leaves\": "
            + minus.stream().mapToInt(r -> r.rawLeafCount).sum() + "\n"
            + "}\n";
    }

    private static void writeRecords(Path path, List<Record> plus,
                                     List<Record> minus) throws IOException {
        try (BufferedWriter w = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            w.write("q\tt\tsign\tn6\tnOther\trawGcdBits\trawGcdSha256\tqVal\tqFreeGcd\tprimitiveGcd\tfactors\tomega\trawWindow\tm5Capable\trawLeaves\tassignedLeaves\tsizeCertified\tthresholdRatio\n");
            for (Record r : concat(plus, minus)) {
                w.write(Integer.toString(r.state.q));
                w.write('\t'); w.write(Integer.toString(r.state.t));
                w.write('\t'); w.write(Character.toString(r.sign));
                w.write('\t'); w.write(Integer.toString(r.state.n6));
                w.write('\t'); w.write(Integer.toString(r.sign == '+' ? r.state.nPlus : r.state.nMinus));
                w.write('\t'); w.write(Integer.toString(r.rawGcd.bitLength()));
                w.write('\t'); w.write(sha256(r.rawGcd.toByteArray()));
                w.write('\t'); w.write(Integer.toString(r.qValuation));
                w.write('\t'); w.write(r.qFreeGcd.toString());
                w.write('\t'); w.write(r.primitiveGcd.toString());
                w.write('\t'); w.write(formatFactors(r.factors));
                w.write('\t'); w.write(Integer.toString(r.factors.size()));
                w.write('\t'); w.write(Boolean.toString(r.sign == '+'
                    ? r.state.plusWindow.nonemptyPrime()
                    : r.state.minusWindow.nonemptyPrime()));
                w.write('\t'); w.write(Boolean.toString(r.sign == '-' && r.state.m5Window.nonemptyPrime()));
                w.write('\t'); w.write(Integer.toString(r.rawLeafCount));
                w.write('\t'); w.write(Integer.toString(r.assignedLeafCount));
                w.write('\t'); w.write(Boolean.toString(r.sizeCertified));
                w.write('\t'); w.write(Double.isFinite(r.thresholdRatio)
                    ? String.format(Locale.ROOT, "%.15g", r.thresholdRatio) : "NA");
                w.newLine();
            }
        }
    }

    private static void writeFactorRecords(Path path, List<Record> plus,
                                           List<Record> minus) throws IOException {
        try (BufferedWriter w = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            w.write("q\tt\tsign\tprime\texponent\n");
            for (Record r : concat(plus, minus)) {
                for (Map.Entry<Long, Integer> e : r.factors.entrySet()) {
                    w.write(r.state.q + "\t" + r.state.t + "\t" + r.sign
                        + "\t" + e.getKey() + "\t" + e.getValue());
                    w.newLine();
                }
            }
        }
    }

    private static List<Record> concat(List<Record> a, List<Record> b) {
        List<Record> all = new ArrayList<>(a.size() + b.size());
        all.addAll(a);
        all.addAll(b);
        all.sort(Comparator.comparingInt((Record r) -> r.state.q)
            .thenComparingInt(r -> r.state.t).thenComparingInt(r -> r.sign));
        return all;
    }

    private static String formatFactors(NavigableMap<Long, Integer> factors) {
        if (factors.isEmpty()) return "1";
        StringBuilder s = new StringBuilder();
        boolean first = true;
        for (Map.Entry<Long, Integer> e : factors.entrySet()) {
            if (!first) s.append('*');
            first = false;
            s.append(e.getKey());
            if (e.getValue() > 1) s.append('^').append(e.getValue());
        }
        return s.toString();
    }

    private static String sha256(byte[] data) {
        try {
            MessageDigest md = MessageDigest.getInstance("SHA-256");
            byte[] digest = md.digest(data);
            StringBuilder s = new StringBuilder(64);
            for (byte b : digest) {
                s.append(String.format(Locale.ROOT, "%02x", b & 0xff));
            }
            return s.toString();
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException(e);
        }
    }

    private static String formatDuration(Duration d) {
        long seconds = d.toSeconds();
        long h = seconds / 3600;
        long m = (seconds % 3600) / 60;
        long s = seconds % 60;
        return String.format(Locale.ROOT, "%02d:%02d:%02d", h, m, s);
    }
}