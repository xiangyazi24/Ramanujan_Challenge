import java.math.BigInteger;
import java.util.Locale;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

/** Exact from-definition checks for the load-bearing Q4316 states. */
public final class Q4316PointCheck {
    private static final BigInteger ONE = BigInteger.ONE;

    private static final long[] PLUS_LEDGER = {
        2L, 3L, 5L, 7L, 17L, 315493L, 858433L, 2720027329L,
        12073365010564729L
    };
    private static final long[] MINUS_LEDGER = {2L, 3L, 5L, 7L, 858433L};

    private static final class Check {
        final int q;
        final int t;
        final char sign;
        final long expectedG;
        final long expectedLargestFactor;
        BigInteger bt;
        BigInteger first;
        BigInteger second;

        Check(int q, int t, char sign, long expectedG, long expectedLargestFactor) {
            this.q = q;
            this.t = t;
            this.sign = sign;
            this.expectedG = expectedG;
            this.expectedLargestFactor = expectedLargestFactor;
        }

        int firstIndex() { return 6 * q + t; }
        int secondIndex() { return sign == '+' ? 12 * q + t : q - 1 - t; }
    }

    private Q4316PointCheck() {}

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        Check[] checks = {
            new Check(5647, 4553, '-', 10013L, 31L),
            new Check(13217, 385, '-', 33847L, 181L),
            new Check(12539, 3701, '+', -1L, 4283L),
            new Check(12979, 11781, '-', -1L, 709L),
            new Check(11437, 139, '-', -1L, 271L)
        };

        int maxIndex = 1;
        for (Check c : checks) {
            maxIndex = Math.max(maxIndex,
                Math.max(c.t, Math.max(c.firstIndex(), c.secondIndex())));
        }

        BigInteger previous = ONE;
        assign(0, previous, checks);
        BigInteger current = BigInteger.valueOf(5L);
        assign(1, current, checks);
        for (int n = 1; n < maxIndex; n++) {
            BigInteger numerator = current.multiply(BigInteger.valueOf(middle(n)))
                .subtract(previous.multiply(BigInteger.valueOf(cube(n))));
            BigInteger[] qr = numerator.divideAndRemainder(
                BigInteger.valueOf(cube(n + 1L)));
            if (qr[1].signum() != 0) {
                throw new IllegalStateException("nonexact recurrence at n=" + n);
            }
            BigInteger next = qr[0];
            assign(n + 1, next, checks);
            previous = current;
            current = next;
        }

        System.out.println("Q4316 EXACT POINT CHECKS");
        System.out.println("largest_index=" + maxIndex + " largest_bits=" + current.bitLength());
        for (Check c : checks) verifyAndPrint(c);
    }

    private static void assign(int index, BigInteger value, Check[] checks) {
        for (Check c : checks) {
            if (index == c.t) c.bt = value;
            if (index == c.firstIndex()) c.first = value;
            if (index == c.secondIndex()) c.second = value;
        }
    }

    private static void verifyAndPrint(Check c) {
        if (c.bt == null || c.first == null || c.second == null) {
            throw new IllegalStateException("missing target for q=" + c.q + " t=" + c.t);
        }
        BigInteger q = BigInteger.valueOf(c.q);
        if (c.bt.mod(q).signum() != 0) {
            throw new IllegalStateException("upper mark failed at q=" + c.q + " t=" + c.t);
        }
        BigInteger raw = c.first.gcd(c.second);
        int qVal = 0;
        BigInteger primitive = raw;
        while (primitive.mod(q).signum() == 0) {
            primitive = primitive.divide(q);
            qVal++;
        }
        long[] ledger = c.sign == '+' ? PLUS_LEDGER : MINUS_LEDGER;
        for (long p : ledger) {
            BigInteger bp = BigInteger.valueOf(p);
            while (primitive.mod(bp).signum() == 0) primitive = primitive.divide(bp);
        }
        if (primitive.bitLength() > 63) {
            throw new IllegalStateException("point primitive exceeds long backend");
        }
        long g = primitive.longValueExact();
        NavigableMap<Long, Integer> factors = factor(g);
        long largest = factors.isEmpty() ? 0L : factors.lastKey();
        if (c.expectedG >= 0L && g != c.expectedG) {
            throw new IllegalStateException("G mismatch at q=" + c.q + " t=" + c.t
                + ": expected " + c.expectedG + " got " + g);
        }
        if (largest != c.expectedLargestFactor) {
            throw new IllegalStateException("largest-factor mismatch at q=" + c.q
                + " t=" + c.t + ": expected " + c.expectedLargestFactor
                + " got " + largest);
        }
        System.out.printf(Locale.ROOT,
            "q=%d t=%d sign=%c qVal=%d G=%d factors=%s largest=%d selected=true%n",
            c.q, c.t, c.sign, qVal, g, format(factors), largest);
    }

    private static NavigableMap<Long, Integer> factor(long value) {
        if (value < 1L) throw new IllegalArgumentException("factor domain");
        NavigableMap<Long, Integer> out = new TreeMap<>();
        long x = value;
        for (long p = 2L; p <= x / p; p += (p == 2L ? 1L : 2L)) {
            while (x % p == 0L) {
                out.merge(p, 1, Integer::sum);
                x /= p;
            }
        }
        if (x > 1L) out.merge(x, 1, Integer::sum);
        long check = 1L;
        for (Map.Entry<Long, Integer> e : out.entrySet()) {
            for (int i = 0; i < e.getValue(); i++) check = Math.multiplyExact(check, e.getKey());
        }
        if (check != value) throw new IllegalStateException("factor reconstruction failed");
        return out;
    }

    private static String format(NavigableMap<Long, Integer> factors) {
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

    private static long cube(long n) { return n * n * n; }

    private static long middle(long n) {
        return (2L * n + 1L) * (17L * n * (n + 1L) + 5L);
    }
}
