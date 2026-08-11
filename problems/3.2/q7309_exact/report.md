# Q7309 exact p-adic high-load scan

Mechanical result only; no proof claim.

## Frozen scope and result

Preregistration SHA-256: `9e3b050e0faf6e7a8acc2ee6706f45772a9226efdf3902fc8091676fc494521a`.

```text
Q7309 relation screen v1
scales: X in {128,256,512,1024}
training: all high-load hits with X<=512
holdout: all high-load hits with X=1024
variables: q,rho,tau,eta,jet_u where jet_u=(b_n/p mod p)/u mod p
linear dictionary: 1 plus the five variables
quadratic dictionary: all monomials of total degree at most two
candidate coefficients: primitive, first nonzero positive, each nonzero in {-2,-1,1,2}
support bound: at most three monomials
false positive: exact zero on every training hit and nonzero on at least one holdout hit
matched control: replace q by (q+1) mod X and recompute the Gessel jet, keeping the local zero labels fixed
```

Independent two-digit formula self-test: PASS on 195 values for p=5,7,11.

| X | primes | active | sum |Z_p| | K>=3 rows | hits | K histogram | max K | C(K,3) sum |
|---:|---:|---:|---:|---:|---:|---|---:|---:|
| 128 | 23 | 13 | 30 | 10 | 30 | `{3: 10}` | 3 | 10 |
| 256 | 43 | 17 | 46 | 10 | 30 | `{3: 10}` | 3 | 10 |
| 512 | 75 | 30 | 70 | 25 | 75 | `{3: 25}` | 3 | 25 |
| 1024 | 137 | 47 | 120 | 78 | 235 | `{3: 77, 4: 1}` | 4 | 81 |

All 370 high-load hit records passed the recurrence, cleared-recurrence, reflection, integrality, nonzero-slope, oriented-lift, independent binomial, and full Gessel congruence assertions.

## Gessel first jet checked

For n=qp+r and p|b_r, the checked prediction is

$$
b_{qp+r}\equiv b_r b_q+p q b'_r b_q\pmod {p^2},
\qquad
\frac{b_{qp+r}}p\equiv b_q\left(\frac{b_r}p+q b'_r\right)\pmod p.
$$

The left side was computed independently from the defining binomial sum using the two-base-p-digit no-carry expansion, not from this congruence.

## Relation screen

| dictionary | columns | candidates | training survivors | held-out false positives | all-data survivors |
|---|---:|---:|---:|---:|---:|
| linear | 6 | 656 | 0 | 0 | 0 |
| quadratic | 21 | 38521 | 0 | 0 | 0 |

Ranks are reported per characteristic in `ranks.csv`; pooling rows over different prime fields into one numerical rank would be meaningless.


## Exact row certificates

### X=128

| n | K | primes |
|---:|---:|---|
| 321 | 3 | `179;193;211` |
| 743 | 3 | `181;197;211` |
| 1609 | 3 | `163;181;191` |
| 3901 | 3 | `179;191;251` |
| 4092 | 3 | `139;191;197` |
| 5048 | 3 | `179;181;193` |
| 9006 | 3 | `193;227;241` |
| 10452 | 3 | `163;193;241` |
| 12849 | 3 | `139;151;197` |
| 16270 | 3 | `181;227;233` |

### X=256

| n | K | primes |
|---:|---:|---|
| 1667 | 3 | `271;431;499` |
| 12678 | 3 | `379;443;499` |
| 24882 | 3 | `419;439;499` |
| 25874 | 3 | `307;367;463` |
| 30870 | 3 | `379;439;499` |
| 35271 | 3 | `271;443;503` |
| 40481 | 3 | `293;419;443` |
| 44131 | 3 | `257;271;443` |
| 51931 | 3 | `283;347;367` |
| 55717 | 3 | `293;331;499` |

### X=512

| n | K | primes |
|---:|---:|---|
| 18808 | 3 | `617;751;757` |
| 25554 | 3 | `571;761;821` |
| 54104 | 3 | `541;571;709` |
| 81786 | 3 | `769;857;907` |
| 92670 | 3 | `593;907;977` |
| 98926 | 3 | `659;857;953` |
| 101703 | 3 | `541;631;991` |
| 102334 | 3 | `631;887;907` |
| 108089 | 3 | `709;821;907` |
| 118762 | 3 | `593;617;821` |
| 133700 | 3 | `659;691;787` |
| 135777 | 3 | `631;761;857` |
| 140703 | 3 | `593;709;887` |
| 146606 | 3 | `541;571;769` |
| 151220 | 3 | `599;769;937` |
| 156106 | 3 | `659;733;953` |
| 189737 | 3 | `599;617;911` |
| 195680 | 3 | `701;821;947` |
| 219301 | 3 | `709;881;911` |
| 220182 | 3 | `541;659;881` |
| 229328 | 3 | `593;701;709` |
| 230427 | 3 | `631;769;977` |
| 240445 | 3 | `857;881;937` |
| 247731 | 3 | `691;733;937` |
| 259185 | 3 | `577;857;937` |

### X=1024

| n | K | primes |
|---:|---:|---|
| 1773 | 3 | `1229;1427;1823` |
| 4325 | 3 | `1193;1907;1951` |
| 5518 | 3 | `1193;1489;1823` |
| 15468 | 4 | `1069;1381;1847;2011` |
| 100215 | 3 | `1823;1949;1973` |
| 127256 | 3 | `1087;1259;1831` |
| 149393 | 3 | `1117;1229;1709` |
| 152257 | 3 | `1069;1087;1193` |
| 153013 | 3 | `1381;1499;1619` |
| 177538 | 3 | `1129;1453;1741` |
| 195824 | 3 | `1381;1531;1709` |
| 211520 | 3 | `1097;1637;1873` |
| 214186 | 3 | `1741;1847;1951` |
| 242400 | 3 | `1559;1723;1831` |
| 248190 | 3 | `1453;1847;1973` |
| 252786 | 3 | `1069;1201;1447` |
| 253313 | 3 | `1049;1559;1621` |
| 258633 | 3 | `1297;1427;1487` |
| 260710 | 3 | `1259;1499;1831` |
| 262667 | 3 | `1153;1381;1559` |
| 285882 | 3 | `1049;1069;1669` |
| 286480 | 3 | `1129;1427;1931` |
| 320387 | 3 | `1657;1741;1951` |
| 337237 | 3 | `1049;1069;1759` |
| 361360 | 3 | `1049;1873;2011` |
| 365697 | 3 | `1129;1229;1487` |
| 389719 | 3 | `1049;1381;2011` |
| 390732 | 3 | `1049;1097;1129` |
| 397093 | 3 | `1471;1657;1931` |
| 402977 | 3 | `1259;1471;1559` |
| 423121 | 3 | `1259;1657;1709` |
| 424409 | 3 | `1447;1669;1759` |
| 431679 | 3 | `1049;1069;1973` |
| 436191 | 3 | `1193;1531;1907` |
| 446679 | 3 | `1087;1097;1433` |
| 447833 | 3 | `1217;1621;2011` |
| 448477 | 3 | `1049;1069;1669` |
| 451314 | 3 | `1069;1129;1499` |
| 463149 | 3 | `1049;1741;2011` |
| 474739 | 3 | `1433;1559;1949` |
| 502306 | 3 | `1129;1523;1621` |
| 505919 | 3 | `1523;1559;1619` |
| 513596 | 3 | `1097;1129;1217` |
| 514272 | 3 | `1499;1531;1873` |
| 522628 | 3 | `1129;1427;1433` |
| 533933 | 3 | `1069;1559;1873` |
| 538631 | 3 | `1049;1129;1847` |
| 544021 | 3 | `1217;1499;1759` |
| 553594 | 3 | `1229;1741;1759` |
| 578361 | 3 | `1087;1381;1447` |
| 594278 | 3 | `1049;1433;1447` |
| 597526 | 3 | `1129;1489;1949` |
| 634606 | 3 | `1709;1723;1931` |
| 676110 | 3 | `1049;1069;2011` |
| 683888 | 3 | `1117;1129;1489` |
| 713864 | 3 | `1049;1487;1723` |
| 716629 | 3 | `1129;1487;1831` |
| 732895 | 3 | `1499;1621;1823` |
| 767845 | 3 | `1193;1867;1931` |
| 770605 | 3 | `1087;1259;1873` |
| 775221 | 3 | `1069;1201;1831` |
| 776660 | 3 | `1069;1907;2011` |
| 797969 | 3 | `1453;1621;1847` |
| 834884 | 3 | `1217;1823;1931` |
| 842502 | 3 | `1087;1117;1523` |
| 867698 | 3 | `1201;1217;1823` |
| 869166 | 3 | `1453;1559;2011` |
| 870816 | 3 | `1229;1427;1873` |
| 900586 | 3 | `1049;1117;1723` |
| 930394 | 3 | `1087;1129;1933` |
| 943882 | 3 | `1447;1831;1847` |
| 951202 | 3 | `1087;1709;1831` |
| 953361 | 3 | `1433;1657;1723` |
| 954011 | 3 | `1447;1489;1907` |
| 961456 | 3 | `1489;1621;1709` |
| 990684 | 3 | `1427;1847;1951` |
| 1032151 | 3 | `1069;1153;1873` |
| 1039094 | 3 | `1087;1117;1831` |

## Files

- `rows.csv`: every exact K>=3 row.
- `hits.csv`: every p-adic hit certificate and every checked first jet.
- `primes.csv`: exact zero sets used by the scatter.
- `ranks.csv`: per-characteristic ranks.
- `relations.csv`: every preregistered sparse candidate and all failure counts.
- `summary.json`: machine-readable summary.
- `MANIFEST.sha256`: hashes of code and exact outputs.

