import RamanujanChallenge.Problem25DualCertificateSemantics

noncomputable section

namespace RamanujanChallenge.P25

def dualCertPp0Terms : List DualCertTerm := [
    { nCoeffs := [-1326456, -4494852, -6601788, -5492772, -2832956, -927848, -188512, -21728, -1088], pExp := 3, qExp := 5, vExp := 4 },
    { nCoeffs := [-1326456, -4494852, -6601788, -5492772, -2832956, -927848, -188512, -21728, -1088], pExp := 3, qExp := 5, vExp := 2 },
    { nCoeffs := [1326456, 4494852, 6601788, 5492772, 2832956, 927848, 188512, 21728, 1088], pExp := 3, qExp := 3, vExp := 4 },
    { nCoeffs := [1326456, 4494852, 6601788, 5492772, 2832956, 927848, 188512, 21728, 1088], pExp := 3, qExp := 3, vExp := 2 },
    { nCoeffs := [-4384800, -14608800, -21076152, -17210200, -8705008, -2794080, -555968, -62720, -3072], pExp := 2, qExp := 4, vExp := 3 },
    { nCoeffs := [-1196640, -5611536, -11458032, -13390432, -9882672, -4782032, -1518528, -305408, -35328, -1792], pExp := 2, qExp := 4, vExp := 1 },
    { nCoeffs := [4384800, 14608800, 21076152, 17210200, 8705008, 2794080, 555968, 62720, 3072], pExp := 2, qExp := 2, vExp := 3 },
    { nCoeffs := [239328, 1026576, 1880976, 1925696, 1206256, 473904, 114144, 15424, 896], pExp := 2, qExp := 2, vExp := 1 },
    { nCoeffs := [884304, 3291336, 5203536, 4594016, 2487808, 848664, 178464, 21184, 1088], pExp := 1, qExp := 5, vExp := 4 },
    { nCoeffs := [884304, 3291336, 5203536, 4594016, 2487808, 848664, 178464, 21184, 1088], pExp := 1, qExp := 5, vExp := 2 },
    { nCoeffs := [14112, 52752, 83528, 73564, 39556, 13332, 2756, 320, 16], pExp := 1, qExp := 3, vExp := 6 },
    { nCoeffs := [-841968, -3133080, -4952952, -4373324, -2369140, -808668, -170196, -20224, -1040], pExp := 1, qExp := 3, vExp := 4 },
    { nCoeffs := [274896, 1232136, 2283560, 2296884, 1380844, 511236, 114412, 14208, 752], pExp := 1, qExp := 3, vExp := 2 },
    { nCoeffs := [6396192, 28917264, 56914120, 64156892, 45732740, 21409172, 6589764, 1287232, 144912, 7168], pExp := 1, qExp := 3, vExp := 0 },
    { nCoeffs := [159552, 1216224, 3340256, 4757280, 4021824, 2135200, 720544, 150240, 17664, 896], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [638208, 2950272, 5946176, 6871392, 5024928, 2414272, 762784, 152928, 17664, 896], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [-105840, -437976, -722976, -643848, -344224, -114256, -23136, -2624, -128], pExp := 0, qExp := 4, vExp := 5 },
    { nCoeffs := [3791232, 14199552, 22485984, 19830944, 10710496, 3640544, 762432, 90112, 4608], pExp := 0, qExp := 4, vExp := 3 },
    { nCoeffs := [-105840, -437976, -722976, -643848, -344224, -114256, -23136, -2624, -128], pExp := 0, qExp := 4, vExp := 1 },
    { nCoeffs := [49392, 212856, 357280, 320496, 171872, 57112, 11568, 1312, 64], pExp := 0, qExp := 2, vExp := 5 },
    { nCoeffs := [-12181248, -50762688, -91918336, -95465264, -62907264, -27344368, -7854368, -1439360, -152832, -7168], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [31959792, 141344376, 271673504, 298740208, 207594592, 94711320, 28411184, 5410080, 593984, 28672], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [5743872, 24318720, 44731968, 47105920, 31384768, 13747264, 3964800, 726848, 76928, 3584], pExp := 0, qExp := 0, vExp := 3 },
    { nCoeffs := [-19784448, -86033664, -162285632, -174809152, -118783232, -52903360, -15467968, -2866688, -305920, -14336], pExp := 0, qExp := 0, vExp := 1 }
  ]

def dualCertPp0Poly : DualCertPoly :=
  dualCertOfTerms dualCertPp0Terms

def dualCertPq0Terms : List DualCertTerm := [
    { nCoeffs := [957312, 4584960, 9577056, 11464736, 8676416, 4308128, 1404384, 289984, 34432, 1792], pExp := 4, qExp := 2, vExp := 1 },
    { nCoeffs := [-4786560, -22446144, -45832128, -53561728, -39530688, -19128128, -6074112, -1221632, -141312, -7168], pExp := 3, qExp := 1, vExp := 0 },
    { nCoeffs := [42336, 200592, 384144, 393072, 238000, 88192, 19680, 2432, 128], pExp := 2, qExp := 2, vExp := 5 },
    { nCoeffs := [-2238624, -8626608, -14152560, -12977280, -7300656, -2587024, -564992, -69632, -3712], pExp := 2, qExp := 2, vExp := 3 },
    { nCoeffs := [-436320, -1852560, -3377808, -3458320, -2174512, -859616, -208608, -28416, -1664], pExp := 2, qExp := 2, vExp := 1 },
    { nCoeffs := [8615808, 37435392, 70831968, 76869536, 52857920, 23924768, 7137120, 1354432, 148480, 7168], pExp := 2, qExp := 0, vExp := 3 },
    { nCoeffs := [-28719360, -128294784, -249319296, -277376512, -195097728, -90118784, -27373824, -5278208, -586752, -28672], pExp := 2, qExp := 0, vExp := 1 },
    { nCoeffs := [147384, 499428, 700780, 532076, 237276, 62352, 8960, 544], pExp := 1, qExp := 3, vExp := 4 },
    { nCoeffs := [147384, 499428, 700780, 532076, 237276, 62352, 8960, 544], pExp := 1, qExp := 3, vExp := 2 },
    { nCoeffs := [3528, 11424, 15170, 10806, 4486, 1090, 144, 8], pExp := 1, qExp := 1, vExp := 6 },
    { nCoeffs := [85464, 318060, 376418, 201046, 48390, 2262, -1040, -136], pExp := 1, qExp := 1, vExp := 4 },
    { nCoeffs := [-729144, -3263004, -6190790, -6906802, -5017090, -2432042, -777120, -156584, -17984, -896], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [1315224, 4896408, 7408090, 5809326, 2358366, 292802, -151296, -75032, -13376, -896], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [-52920, -218988, -361488, -321924, -172112, -57128, -11568, -1312, -64], pExp := 0, qExp := 2, vExp := 5 },
    { nCoeffs := [2211552, 8177760, 12757496, 11063216, 5865208, 1953664, 400288, 46208, 2304], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [-52920, -218988, -361488, -321924, -172112, -57128, -11568, -1312, -64], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [-3528, -9660, -10340, -5636, -1668, -256, -16], pExp := 0, qExp := 0, vExp := 5 },
    { nCoeffs := [-6757632, -28282560, -51340800, -53270960, -34914944, -15020560, -4247744, -762304, -78848, -3584], pExp := 0, qExp := 0, vExp := 3 },
    { nCoeffs := [22617576, 97621044, 182551436, 194558796, 130534204, 57286432, 16471344, 2995840, 313088, 14336], pExp := 0, qExp := 0, vExp := 1 }
  ]

def dualCertPq0Poly : DualCertPoly :=
  dualCertOfTerms dualCertPq0Terms

def dualCertPv0Terms : List DualCertTerm := [
    { nCoeffs := [478656, 2292480, 4788528, 5732368, 4338208, 2154064, 702192, 144992, 17216, 896], pExp := 3, qExp := 1, vExp := 0 },
    { nCoeffs := [2871936, 13116672, 26100000, 29763616, 21473152, 10177504, 3172320, 627584, 71552, 3584], pExp := 2, qExp := 0, vExp := 1 },
    { nCoeffs := [-283824, -1533384, -3605544, -4723080, -3808616, -1976000, -664336, -140320, -16960, -896], pExp := 1, qExp := 3, vExp := 0 },
    { nCoeffs := [-134064, -559944, -1040168, -1070696, -660056, -249760, -56896, -7168, -384], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [-159552, -657792, -1157648, -1139024, -686720, -260208, -60592, -7936, -448], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [-52920, -218988, -361488, -321924, -172112, -57128, -11568, -1312, -64], pExp := 0, qExp := 4, vExp := 3 },
    { nCoeffs := [-52920, -218988, -361488, -321924, -172112, -57128, -11568, -1312, -64], pExp := 0, qExp := 4, vExp := 1 },
    { nCoeffs := [-21168, -79128, -130584, -121308, -68712, -24244, -5208, -624, -32], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [-1440288, -7604136, -17462568, -22281804, -17477704, -8816980, -2882520, -592240, -69664, -3584], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [-10584, -32508, -40680, -27248, -10640, -2436, -304, -16], pExp := 0, qExp := 0, vExp := 3 },
    { nCoeffs := [-1392840, -5421564, -8758024, -7889840, -4389056, -1551492, -340784, -42512, -2304], pExp := 0, qExp := 0, vExp := 1 }
  ]

def dualCertPv0Poly : DualCertPoly :=
  dualCertOfTerms dualCertPv0Terms


def dualCertM00Poly : DualCertPoly :=
  (2 * dualCertN + 5) * (dualCertN + 3) ^ 2 *
    (136 * dualCertN ^ 4 + 1424 * dualCertN ^ 3 + 5548 * dualCertN ^ 2 +
      9551 * dualCertN + 6141)

def dualCertM01Poly : DualCertPoly :=
  384 * dualCertN ^ 6 + 6384 * dualCertN ^ 5 + 44168 * dualCertN ^ 4 +
    162698 * dualCertN ^ 3 + 336377 * dualCertN ^ 2 + 369933 * dualCertN + 169011

def dualCertM02Poly : DualCertPoly :=
  480 * dualCertN ^ 4 + 4980 * dualCertN ^ 3 + 19210 * dualCertN ^ 2 +
    32690 * dualCertN + 20730

set_option maxHeartbeats 0 in
set_option maxRecDepth 100000 in
theorem dualCert_identity_zero_poly :
    dualCertNormalize (dualCertDeltaPoly *
        (dualCertM00Poly * dualCertNextNumPoly 0 +
          dualCertM01Poly * dualCertNextNumPoly 1 +
          dualCertM02Poly * dualCertNextNumPoly 2 -
          dualCertLambdaPoly * dualCertCurNumPoly 0)) =
      dualCertNormalize (
      dualCertOpPPoly dualCertPp0Poly +
        dualCertOpQPoly dualCertPq0Poly +
        dualCertOpVPoly dualCertPv0Poly) := by
  native_decide

def dualCertPp0 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPp0Poly n p q v

def dualCertPq0 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPq0Poly n p q v

def dualCertPv0 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPv0Poly n p q v

def dualCertDPp0 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertPDeriv dualCertPp0Poly) n p q v

def dualCertDPq0 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertQDeriv dualCertPq0Poly) n p q v

def dualCertDPv0 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertVDeriv dualCertPv0Poly) n p q v

theorem dualCertPp0_hasDerivAt (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertPp0 n x q v) (dualCertDPp0 n p q v) p := by
  simpa [dualCertPp0, dualCertDPp0] using
    dualCertEval_hasDerivAt_p dualCertPp0Poly n p q v

theorem dualCertPq0_hasDerivAt (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertPq0 n p x v) (dualCertDPq0 n p q v) q := by
  simpa [dualCertPq0, dualCertDPq0] using
    dualCertEval_hasDerivAt_q dualCertPq0Poly n p q v

theorem dualCertPv0_hasDerivAt (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertPv0 n p q x) (dualCertDPv0 n p q v) v := by
  simpa [dualCertPv0, dualCertDPv0] using
    dualCertEval_hasDerivAt_v dualCertPv0Poly n p q v

macro "row0_eval_num" name:ident value:num : command =>
  `(@[simp] private theorem $name (n p q v : ℝ) :
      dualCertEval ($value : DualCertPoly) n p q v = ($value : ℝ) := by
    change dualCertEval (dualCertConst ($value : ℤ)) n p q v = _
    simp)

row0_eval_num dualCertEval_row0_136 136
row0_eval_num dualCertEval_row0_1424 1424
row0_eval_num dualCertEval_row0_5548 5548
row0_eval_num dualCertEval_row0_9551 9551
row0_eval_num dualCertEval_row0_6141 6141
row0_eval_num dualCertEval_row0_384 384
row0_eval_num dualCertEval_row0_6384 6384
row0_eval_num dualCertEval_row0_44168 44168
row0_eval_num dualCertEval_row0_162698 162698
row0_eval_num dualCertEval_row0_336377 336377
row0_eval_num dualCertEval_row0_369933 369933
row0_eval_num dualCertEval_row0_169011 169011
row0_eval_num dualCertEval_row0_480 480
row0_eval_num dualCertEval_row0_4980 4980
row0_eval_num dualCertEval_row0_19210 19210
row0_eval_num dualCertEval_row0_32690 32690
row0_eval_num dualCertEval_row0_20730 20730

theorem dualCert_identity_zero (n : ℕ) (p q v : ℝ) :
    (4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)) *
        ((positiveMatrix (n : ℤ) 0 0 : ℝ) * dualCertNextNum (n : ℝ) 0 p q v +
          (positiveMatrix (n : ℤ) 0 1 : ℝ) * dualCertNextNum (n : ℝ) 1 p q v +
          (positiveMatrix (n : ℤ) 0 2 : ℝ) * dualCertNextNum (n : ℝ) 2 p q v -
          dualCertLambda (n : ℝ) * dualCertCurNum (n : ℝ) 0 p q v) =
      dualCertOpP (n : ℝ) (dualCertPp0 (n : ℝ) p q v)
          (dualCertDPp0 (n : ℝ) p q v) p q v +
        dualCertOpQ (n : ℝ) (dualCertPq0 (n : ℝ) p q v)
          (dualCertDPq0 (n : ℝ) p q v) p q v +
        dualCertOpV (n : ℝ) (dualCertPv0 (n : ℝ) p q v)
          (dualCertDPv0 (n : ℝ) p q v) p q v := by
  have h := congrArg (fun P => dualCertEval P (n : ℝ) p q v)
    dualCert_identity_zero_poly
  dsimp only at h
  rw [dualCertEval_normalize, dualCertEval_normalize] at h
  simp only [dualCertDeltaPoly, dualCertM00Poly, dualCertM01Poly, dualCertM02Poly,
    dualCertNextNumPoly, dualCertCurNumPoly, dualCertLambdaPoly,
    dualCertOpPPoly, dualCertOpQPoly, dualCertOpVPoly,
    dualCertDPoly, dualCertSnumPoly] at h
  simp only [dualCertEval_add, dualCertEval_sub,
    dualCertEval_mul, dualCertEval_pow,
    dualCertEval_const,
    dualCertEval_N, dualCertEval_P, dualCertEval_Q, dualCertEval_V] at h
  simp at h
  simpa [dualCertD, dualCertSnum,
    dualCertNextNum, dualCertCurNum, dualCertLambda,
    dualCertOpP, dualCertOpQ, dualCertOpV,
    dualCertPp0, dualCertPq0, dualCertPv0,
    dualCertDPp0, dualCertDPq0, dualCertDPv0,
    positiveMatrix, Matrix.cons_val_two] using h

end RamanujanChallenge.P25

end
