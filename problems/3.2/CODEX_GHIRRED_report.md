# Absolute irreducibility certificates for $G_h$, $2\le h\le40$

## Verdict

The contiguous chain is **complete for every $h=2,\ldots,40$**. For each height, two independent exact methods certify that $G_h$ is absolutely irreducible over $\overline{\mathbf Q}$:

1. direct characteristic-zero absolute factorization by Singular `absfact.lib`;
2. the critical-value/monodromy criterion, certified modulo two independent good primes.

Therefore the primitive integral model of $G_h$ remains absolutely irreducible modulo $p$ for all but finitely many primes $p$, for every fixed $2\le h\le40$.

## Convention and $h=2,\ldots,6$ reproduction

The certified convention is exactly the one in the specification and ledger appendix AS.4:

```text
N_h(r) = K_(h-1)(r+1)
R_h(x,y) = N_h(x) prod_(j=1)^h (y+j)^3
           - N_h(y) prod_(j=1)^h (x+j)^3
G_h(x,y) = R_h(x,y)/(x-y).
```

For even $h$, the forced rational factor $(2X+h+1)$ of $N_h$ is **retained**. No cofactor normalization is made. Exact SymPy construction from the $K_m$ recurrence was compared coefficient-for-coefficient with the Sage construction.

| h | deg $N_h$ | deg $G_h$ | exact remainder | `factor_list` over $\mathbf Q$ | convention match | wall (s) |
|---:|---:|---:|:---:|:---:|:---:|---:|
| 2 | 3 | 8 | 0 | one factor | yes | 1.189 |
| 3 | 6 | 14 | 0 | one factor | yes | 2.307 |
| 4 | 9 | 20 | 0 | one factor | yes | 4.779 |
| 5 | 12 | 26 | 0 | one factor | yes | 7.591 |
| 6 | 15 | 32 | 0 | one factor | yes | 9.662 |

## Exact criteria

### Method AF: direct absolute factorization

Singular 4.4.1's `absFactorizeBCG` works over characteristic zero and returns all absolute factors, grouped with their conjugates. For every $h=2,\ldots,40$, `absolute_factors[4]` is 1, its multiplicity is 1. This directly proves absolute irreducibility and hence $\mathbf Q$-irreducibility. For $h\le6$, the returned factor was additionally parsed back over $\mathbf Q[x,y]$ and verified coefficientwise to be an associate of the exact $G_h$.

Two controls guard the interface/parser: $x^2+y^2$ returns two absolute factors, while $x^2+y^2+1$ returns one.

### Method BG: critical values and geometric monodromy

Put $q_h(X)=\prod_{j=1}^h(X+j)$, $\delta_h=N_h/q_h^3$, and

$$A_h=N_h'q_h-3N_hq_h',\qquad C_h(T)=\operatorname{Res}_X(N_h-Tq_h^3,A_h).$$

For a listed prime $q$, the certificate stores every coefficient of $C_h\bmod q$ and checks

$$\deg C_h=4h-4,\quad C_h(0)\ne0,\quad \gcd(C_h,C_h')=1,$$

together with all construction degrees, pole coprimality, three out-of-sample resultant evaluations, and (for $h\le4$) an independent exact-$\mathbf Z$ resultant reduction. Expected degree modulo one prime prevents leading-term loss, so these checks prove the same nonvanishing statements in characteristic zero.

Here is the exact implication used. The checks give the branch profile of the degree-$3h$ map $\delta_h$: an index-3 zero at infinity, $h$ index-3 poles, and $4h-4$ simple ramification points with pairwise distinct nonzero critical values. If $\delta_h=g\circ u$ with $\deg u=a>1$, every one of those $4h-4$ simple points must ramify in $u$, hence $a\ge2h-1$ by Riemann--Hurwitz; but the nontrivial outer factor gives $a\le3h/2$. This contradicts $h\ge3$, and the $h=2$ endpoint is ruled out by its index-3 zero. Thus the geometric monodromy is primitive. A simple branch point supplies a transposition, so the monodromy group is $S_{3h}$. Its action on ordered unequal pairs is transitive; consequently the off-diagonal factor $G_h$ is geometrically irreducible.

This is an exact criterion (method (c) in the specification), not a point-count heuristic.

## Per-height certificate chain

`AF` is the direct Singular certificate; `BG(q)` is the critical-value certificate at prime $q$; `SQ` is the independent SymPy $\mathbf Q$-factorization recorded for $h\le6$. Hashes are SHA-256 of each certificate's canonical compact JSON payload.

| h | deg $G_h$ | methods | primes | SHA-256 certificates | wall (s) |
|---:|---:|:---|:---|:---|---:|
| 2 | 8 | AF + BGx2 + SQ | 65537, 1000003 | AF=`60d2ade2c69c8c0c14a4003d64a96e7aedf657d6094802482d114607b2dc6dfa`<br>BG(65537)=`896acbda15ef41e201a0732a9ff05ee28201aa98d9e17013f65e061fcd06ce71`<br>BG(1000003)=`70acc1480b35a078558903d23eaeb07cb500a5687b6caad71505e5564bd88b4f`<br>SQ=`ad6356a001746e72542ce25d735ba3bde3faf08a4b44f9cfd6d3d5e51a65490e` | 1.422 |
| 3 | 14 | AF + BGx2 + SQ | 65537, 1000003 | AF=`cfda5cdd4858a822b5b5289003850748d94f9d7eb61cf1955b969a7c2e8b42c9`<br>BG(65537)=`2d14fcbbe1e94282e9e2198a744f5f9203990ac5b2b4532a0ce21e74a842a839`<br>BG(1000003)=`53a3f0d0181424db67a60b820ff254666a0533f2af4634da98c3897e1f47d0d4`<br>SQ=`4b3e36378e30225f4d2708c0e23318a725b7b410ca8960fc86ba5182bf191f79` | 0.037 |
| 4 | 20 | AF + BGx2 + SQ | 65537, 1000003 | AF=`e0f12798256cbf877ec0bc2aafc98916898109b3abd7a053d4d4e00576f5605d`<br>BG(65537)=`1fd7e7dcaa876106c73a05d479030e3bd794f2de5df4d2b18cb7dd4cf3daa0a9`<br>BG(1000003)=`5df795d95e9d48351de620e3c810ae0e7f13a2129a3fd1145b79def5e3ed49f5`<br>SQ=`c0ffdfa50adf3ae1ce88da1235741acd64d4f45556289f3c4734ed46b6a3fcba` | 0.151 |
| 5 | 26 | AF + BGx2 + SQ | 65537, 1000003 | AF=`58ecd13afa1814642d01261c1010645475ac19ddb0dbdc21aaf99b703a66ce23`<br>BG(65537)=`443230151e6d850bf68d941148b71a64a7b6087d56258695f03e7e5178250eef`<br>BG(1000003)=`41f904de4091876dbce737e3283e41f8a88594bb7f02946db96571b83120b4ef`<br>SQ=`342ed383faed41b72b0e3ca2fdcacb360d275918dd3d392ed6e59c1e8722bf80` | 0.231 |
| 6 | 32 | AF + BGx2 + SQ | 65537, 1000003 | AF=`d789f91daea0063baba6c7b352f880e97e344e663f958bf6c196766ff2f88088`<br>BG(65537)=`1cc2a25b9f384d16a374198d45ad76b94df0c3f4d1e747b0fd14c6a485a7a58c`<br>BG(1000003)=`f9d8ac629639f7c8877dade3433ea93a78fd70790705385c848efb28d1bbe595`<br>SQ=`a726c42e1ba29d82232cccbda7d1ed25cb9e9abdb4c2fa2318eb3fcc65606dd3` | 0.356 |
| 7 | 38 | AF + BGx2 | 65537, 1000003 | AF=`744e4c23e797ad9203c38a893b49766304b9a37c970ddc9894faec81f281b166`<br>BG(65537)=`eced7160017da8be13c8e21f75807c17915fbac2e27b3cf45bd736cc911c1421`<br>BG(1000003)=`8d0306ff0fbf7a7c864732840ff954b4a9d9859947f9190429c813d90c0a01ae` | 0.237 |
| 8 | 44 | AF + BGx2 | 65537, 1000003 | AF=`8d150be4cfeb3285d1e58e97e5e7dccb0133329d218d33c543c202bfb85f662e`<br>BG(65537)=`7d30a7cac5e842c18d1000646b3c9ddfc780229909c5c2e9c05d1de284031034`<br>BG(1000003)=`96715e0f8754ac0aebeb080faaa6988ffd84263ae5b5f2543b256ac37b104b7d` | 0.487 |
| 9 | 50 | AF + BGx2 | 65537, 1000003 | AF=`8dafbfa561abda796b64a80fb9cd9cd29d22fd26a88dfdb202a19a137f5945f4`<br>BG(65537)=`07734fef30ef9f2b324b4db2947ec24b8b5cb5c1b4eebec791c3327dc4092e83`<br>BG(1000003)=`5a58c826d81b3dfa8bd32015a24e52699f0bf409ae2e39861f236eb93d45451b` | 0.614 |
| 10 | 56 | AF + BGx2 | 65537, 1000003 | AF=`10b5b023404be5e02fb4a9e7fc69bab7921c5b497bc9fbc948c20279c14cd88d`<br>BG(65537)=`46699730db3e8f8f835ff89b43156130ea049a64534d01dc2d7f2de9df5010b7`<br>BG(1000003)=`cad1942e0761feda1afe4b7f1d335f63ba9c4ef618d9482a9ad33bedd275a3c1` | 0.619 |
| 11 | 62 | AF + BGx2 | 65537, 1000003 | AF=`fc593d8083e36a07970d5a6875b0706a3406e40212224383312aeeba09446e69`<br>BG(65537)=`421caedee4c56fc6a0a6027bca7606743e2ca45da936bfd93f567f654ffc6f58`<br>BG(1000003)=`193cd9146337108b67a38b5fa22e381e8d55ac258f339c3b34d89ff5a9c0ec4a` | 0.695 |
| 12 | 68 | AF + BGx2 | 65537, 1000003 | AF=`78a9c490c6a0e7e3ffb5c0293713daba350994c96a4f41b6fd7db42995a9741d`<br>BG(65537)=`c37373efb066b6ca63ff8a37e21023ca3899f03dcb37e18768e895d70cfffff0`<br>BG(1000003)=`39576593f010bda24984c1130d2381c9c0b023b5dc587a2bc7f81627fc718913` | 1.864 |
| 13 | 74 | AF + BGx2 | 65537, 1000003 | AF=`d8ebd68600f9c126eb028addc99ef20111146823108b9cddee37ac764bfffe17`<br>BG(65537)=`19310d6b9d69daec6c5d293ab15b76c27f10cd76bbe7b9413196a4ecc3876141`<br>BG(1000003)=`0692d6b2d25b2649635165d6eacd9c98c0a6d398ed93c4db406f1de2c009f662` | 1.009 |
| 14 | 80 | AF + BGx2 | 65537, 1000003 | AF=`0517c0be65ff87bd8da1407dde7b56db1dae4b0bd1fbba16c884a248c1de750b`<br>BG(65537)=`4fdd1ce3be4502d944c77bf44cce9ccdbe19fa4f3e5dab51d8ce136643726404`<br>BG(1000003)=`2e9f7bc66751876224cc301c32a08ee9a61635ee1a3c9432c189d5a15e7d10cf` | 1.869 |
| 15 | 86 | AF + BGx2 | 65537, 1000003 | AF=`10ab76670b5d4b811e18fad6bd5f7c8b6903a741cd3ae96e933d993e724dc869`<br>BG(65537)=`f036e7ef48108782012f0e5b6043fb8aa3ec88a86c84d68c92b6874a1684be31`<br>BG(1000003)=`1509176a8ea10674d5ea62686e39536a3285aa44e5bc67d822956a75dee2b808` | 1.674 |
| 16 | 92 | AF + BGx2 | 65537, 1000003 | AF=`b3bce2c89ee977504aa0fafaf09a16b719f6c2c71fc01ad0cae0ea6571252056`<br>BG(65537)=`3058ca64f11299c959315b0fb6cc4bc19890f7a2a4a91e95f555116b56912475`<br>BG(1000003)=`c3037d47a9c26983c74c36acca1e90cd80b688cc71c498b6656e33a866b2eee1` | 2.305 |
| 17 | 98 | AF + BGx2 | 65537, 1000003 | AF=`24e8fe180054ff71d826efb0e714d5aeb3c2e5c4a567cf2c4b2798501a9023e3`<br>BG(65537)=`fdd7b4544a9a6df97aafee70bd052a3ae4c83b084266e8546544eed8435fabf6`<br>BG(1000003)=`3fbf534dce907eabb479e6ae157779a202f952725d3d2bb6b63255833ec0f455` | 3.885 |
| 18 | 104 | AF + BGx2 | 65537, 1000003 | AF=`67cfb0f050279727ce37c346e82892515c31394a2e15ff2b75617021bbe6d4b3`<br>BG(65537)=`dbbbc5cb32c15e8da273ff232087109ca8f2ebb5111cda0fd4f8135b6c828b8a`<br>BG(1000003)=`7fa69824506e9832a290a5db2a0d11c3d1f67ccc76ccbcd5be35ac2b237d2a4f` | 3.488 |
| 19 | 110 | AF + BGx2 | 65537, 1000003 | AF=`0b23692a2163a5ab3853597a4472eab479f5aefd1942acb7d230d8c8635da1cf`<br>BG(65537)=`dfb6188425af2730b48690559d58fb4ff28709cbd0f37723a52e169ea6940834`<br>BG(1000003)=`c63d96810283cf97f151826a45d909f8007c2e70a447d7e22ca02acdf520773a` | 8.663 |
| 20 | 116 | AF + BGx2 | 65537, 1000003 | AF=`f1652878d82ff299692ebb6748957b3d6124a39017d0c3e44f4a5206868e5bda`<br>BG(65537)=`8a7ea66af93e74fe8f7006a1b51f58a89be102d52465677e73f5026be3e24143`<br>BG(1000003)=`ac67650a43862a64c1e96c23b8cc5cbf93886704bce58da0f344332363eb9cda` | 3.889 |
| 21 | 122 | AF + BGx2 | 65537, 1000003 | AF=`25f18dbbfc17c1a18cae25a8b171f52ffe29e52f529e97bd349dde32c02f2fdf`<br>BG(65537)=`8ccefc81d80a44a9375f1ffeccc00c2330002ea519740f39618fc2d0b5cc2a7f`<br>BG(1000003)=`e3255eead8e4458fe181e4cbe8207b68a3bcef39492f6ea2a327be2198b4d6de` | 4.850 |
| 22 | 128 | AF + BGx2 | 65537, 1000003 | AF=`3f9fcd25fc553e3a228d68bc4dcf46b515f1303f0631e5faedc834aa61c18fbb`<br>BG(65537)=`5886e9a5a23388cee94f2e4116a668344f03a45036b736c473887b88e3c4caec`<br>BG(1000003)=`942741aa1ee792deb96e655f43a05ff118acd6bc8ac8c1f4450996960cd1dd73` | 5.521 |
| 23 | 134 | AF + BGx2 | 65537, 1000003 | AF=`313fff3f272992e32d780753f940daf76315cc7f291194adb4a1f9d1d7b5c13a`<br>BG(65537)=`bda277a99befd1ca1a6c82ff43169eb9aa924cd9fe8b63e4787e200a17c1e873`<br>BG(1000003)=`9c571485fd55b362d34bdc749b2ef64ae153a8f6a236f24eb6106796d941f50c` | 20.337 |
| 24 | 140 | AF + BGx2 | 65537, 1000003 | AF=`502deb1b6b0443df7ca3fefd457e22fe2727dd118daca031821eed4472940610`<br>BG(65537)=`65345c040496dbce647700d55edbac704c28c7a1b13af844d639f9f154ad5176`<br>BG(1000003)=`b77f301ae02de4a7cd156b7a5d182e51e9c664359a2ff36abdc30de457f3f3c7` | 6.730 |
| 25 | 146 | AF + BGx2 | 65537, 1000003 | AF=`d22188608191e7c845efe5226be49277090caa92228dd4fb949dcfc5baff286d`<br>BG(65537)=`cfc2d833d8d308d1b9bc6f89d0815a679f031d5b3e60ade44000bf2464c5fa46`<br>BG(1000003)=`dc347672bb66bc4a3c1f6a7331a6cde441ec544db81e9e6938283a371b036c3e` | 4.453 |
| 26 | 152 | AF + BGx2 | 65537, 1000003 | AF=`d4fdad047329b8053ffb5ad3256997fac2e8686b0dc7076454e36ba069c5faaf`<br>BG(65537)=`1990d163fb6940ca2f5910754c0309cd86ebb53b83d6facefd636bd3c0612df8`<br>BG(1000003)=`919d6103273a05209c1a83f89d382e5a9cfceccdfe8f53a8a84bf9f77f7db128` | 8.704 |
| 27 | 158 | AF + BGx2 | 65537, 1000003 | AF=`19503cc1f41eb7d24ef1067731983bd4c007e157839d5e09550ba6812db8704d`<br>BG(65537)=`199674f3d0fa9549916f0df61309e0216ae0058ec1388b6ef0d6ebb8be399458`<br>BG(1000003)=`2c45dd5c8096de31de8c010a22a1110582584045eb5d442e6881877aae41de03` | 11.214 |
| 28 | 164 | AF + BGx2 | 65537, 1000003 | AF=`0239909d4fc1433451cc82c73570e6a9d5e5f20299537d22dcb3aac14d992dcc`<br>BG(65537)=`f72988dc2cec60c022338adfbe18f4081e566f36cf86c60cffc4f7167c6ffea7`<br>BG(1000003)=`70e152aa645594478246894a8c728ebc0bdb0ed41ae6b2f7b7f98092e1702fad` | 35.586 |
| 29 | 170 | AF + BGx2 | 65537, 1000003 | AF=`bd34340174575f8f61b825a35420541eaabe62bd1f4ecd32e27704b3d8712b27`<br>BG(65537)=`b871a31d8666aa57dba482cef65fc00e4a72338fd253da4e7f07f2e272460142`<br>BG(1000003)=`bf3b0795baff95efb44c10f2ee72f16f632217221cafa966f393b90f72ea3791` | 17.346 |
| 30 | 176 | AF + BGx2 | 65537, 1000003 | AF=`c186006e21c8b952c8cb4aeef2000481d527541de9312b1bf659595d63d7e788`<br>BG(65537)=`f2182413e9dce32860106ed60379ce4e3ef12fcd76b9dd66c2c1d62c2d085693`<br>BG(1000003)=`71ee7988b71adc76685f7531200ff2df298dc719978d52d1a0a0d1728885f2ff` | 32.790 |
| 31 | 182 | AF + BGx2 | 65537, 1000003 | AF=`cf345242a88a60d0c9fe2eedc02e46c8cd6b6605e826882152c998c92bc78366`<br>BG(65537)=`0314c6da9ea22b43c538c93e3d79b5fcfb30306a68fbcea813b9815163754273`<br>BG(1000003)=`dff749757d9620b672fbe8cead7086250a1497869baff9c2ee4649ec9d331332` | 19.342 |
| 32 | 188 | AF + BGx2 | 65537, 1000003 | AF=`e775179ebb9df9633269165584bf76f19801d37b74ce5262c9ac6f0ea257f3b2`<br>BG(65537)=`17c4ac5ff6efa39f44b882343e9c1e55c94fe92536a50407f54485f062a4607f`<br>BG(1000003)=`369d43b7a082787d67f4a96036d2b6e50649ee32490b6c8f19cd470a8994c511` | 34.942 |
| 33 | 194 | AF + BGx2 | 65537, 1000003 | AF=`414a095a811f9c56699497ab16bc1974048bddd2247d35f347273a894a146038`<br>BG(65537)=`46b4270e6db3168e9cfeebf6a53d8911692177e6e8987f71c8d8ac32f2bf5e54`<br>BG(1000003)=`65d513d0c50c4cd01c8edf818c354abbe133501d06425c3305f2f931777bd57d` | 30.286 |
| 34 | 200 | AF + BGx2 | 65537, 1000003 | AF=`2e1de1be7ea6a895af2b92da79e7472113a8280e98760b562b7127a5e4278169`<br>BG(65537)=`e7f855b45e8572e8b39a1b42c9d05c88032145332dd8c53fcf4091a2bd45a7dd`<br>BG(1000003)=`40c47520f3a6ad7b82cd94565775de7fc2a46813e9a96f177ed0c8b0d826bfea` | 26.700 |
| 35 | 206 | AF + BGx2 | 65537, 1000003 | AF=`be89e5a32f3df75e64e56efc076797e9cd2fdd22ec817604be44f546eeb700ae`<br>BG(65537)=`372ea9f518f4400189c9cf10d1dd95bb074969fe6dd36c8c7a7f4ac636ccbf4d`<br>BG(1000003)=`e4dd1a164e6d6d317f9703fa97cd227354f656b232c807e1856b4ad36c65e86c` | 25.866 |
| 36 | 212 | AF + BGx2 | 65537, 1000003 | AF=`20528b6167441ec39f93cd56ca85e0982888d6a10002166843beffa670e50448`<br>BG(65537)=`541f925bed107d1bfc6bb58821032c8a518a347b033a848046596a268fb23d37`<br>BG(1000003)=`da89964309951aedadf3c16e00214308a46537aedee4e1b0e62c7f505c6b5265` | 63.884 |
| 37 | 218 | AF + BGx2 | 65537, 1000003 | AF=`f4644cc2da15b3ee54529a6985f5fa10070cf1e72aba8f163b78f75e77efd5af`<br>BG(65537)=`58a1cbc932edfeb5a81c3be02657b6db0487330e4d6fd1a4a3e10e030778bf4c`<br>BG(1000003)=`0cfa85eb3a83261740fc12737dce4bf76e96c211d07df6cd7c895ace4cb91162` | 49.423 |
| 38 | 224 | AF + BGx2 | 65537, 1000003 | AF=`979c01529a0c212e79f270284b31279def424e56afc7d1d76a3af66e6c48691f`<br>BG(65537)=`d78d946105ae95048e7d74606b817d8365c455eff75c67fd77fe8f0d1d9b31f9`<br>BG(1000003)=`96ed38f442a55b4d49fc7d60791e960ba76618cc17db52d2c8c31f2d4d746c8a` | 57.128 |
| 39 | 230 | AF + BGx2 | 65537, 1000003 | AF=`e5de06fd5a16f2a2834ef096670dc17c5d24511c583239947b92bf7de92214d1`<br>BG(65537)=`0cbacfdb41c4fc0f3c38ccfa3cb5ccb607d28a188581aa12623ec671f82dc966`<br>BG(1000003)=`c3039577121ad4c08f5fc03de79f72421f990f3a0a4bee68e5e33d6be24dc5c7` | 61.283 |
| 40 | 236 | AF + BGx2 | 65537, 1000003 | AF=`006dcd864090f85c3013fe8c685eebc505224dce9acdb7dd7c9e3a4f61c3c26c`<br>BG(65537)=`936af25176f630691571f3836b5ed64c0a74b2ba377dd18dc3daceb372685216`<br>BG(1000003)=`c30aaaa391c4cdb0083e26a9581cea96058889be6eb0ed58986a5084de92e1b8` | 63.856 |

## Anomalies

None. No factorization, convention mismatch, or failed certificate gate occurred.

## Environment and reproducibility

The invoking system Python did not have `python-flint`; the run therefore used Sage's FLINT-backed prime-field univariate polynomials for modular resultants/interpolation and Singular's exact `absfact.lib` for bivariate absolute factorization. No floating-point arithmetic enters a certificate.

- Sage: `10.9`
- SymPy: `1.14.0`
- modular polynomial backend: `sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint`
- invoking Python `flint` importable: `false`
- Singular:

```text
Singular for arm64-Darwin version 4.4.1 (44100, 64 bit) Jan 2025
with
	GMP(6.3.0),NTL(11.6.0),FLINT(3.4.0),
	omalloc,static readline(8),Plural,DBM,
```

Canonical certificate encoding: `UTF-8 compact JSON; sort_keys=true; separators=(',', ':')`.

- chain SHA-256: `bddedd0b0068089cf41a284cd8dd61b9ea749c1be98a282a1db7c97d05698e62`
- generated JSON file SHA-256: `e657433996ce3bfc32d91357d4c6c8fa19994e2da97ca53484702e7b2ef4b3de`

Reproduce from this directory with:

```console
$ python3 CRON_ghirred.py
```

The script prints progress for every height, rewrites both deliverables, rereads the JSON, and verifies every stored certificate and bundle hash before returning success.
