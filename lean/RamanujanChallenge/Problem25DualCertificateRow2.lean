import RamanujanChallenge.Problem25DualCertificateSemantics

noncomputable section

namespace RamanujanChallenge.P25

def dualCertPp2Terms : List DualCertTerm := [
    { nCoeffs := [-1404864, -6624288, -13926672, -17193624, -13806192, -7535248, -2831296, -723256, -120224, -11744, -512], pExp := 3, qExp := 5, vExp := 4 },
    { nCoeffs := [-1404864, -6624288, -13926672, -17193624, -13806192, -7535248, -2831296, -723256, -120224, -11744, -512], pExp := 3, qExp := 5, vExp := 2 },
    { nCoeffs := [1404864, 6624288, 13926672, 17193624, 13806192, 7535248, 2831296, 723256, 120224, 11744, 512], pExp := 3, qExp := 3, vExp := 4 },
    { nCoeffs := [1404864, 6624288, 13926672, 17193624, 13806192, 7535248, 2831296, 723256, 120224, 11744, 512], pExp := 3, qExp := 3, vExp := 2 },
    { nCoeffs := [-4838400, -22570560, -46914720, -57227536, -45373112, -24435080, -9052944, -2278624, -372928, -35840, -1536], pExp := 2, qExp := 4, vExp := 3 },
    { nCoeffs := [-1157760, -6963264, -18785184, -30017968, -31580688, -22975168, -11797328, -4276624, -1072832, -177408, -17408, -768], pExp := 2, qExp := 4, vExp := 1 },
    { nCoeffs := [4838400, 22570560, 46914720, 57227536, 45373112, 24435080, 9052944, 2278624, 372928, 35840, 1536], pExp := 2, qExp := 2, vExp := 3 },
    { nCoeffs := [231552, 1300032, 3237024, 4708784, 4432624, 2821984, 1230672, 363056, 69344, 7744, 384], pExp := 2, qExp := 2, vExp := 1 },
    { nCoeffs := [936576, 4728384, 10548384, 13714608, 11523472, 6545312, 2547488, 671376, 114736, 11488, 512], pExp := 1, qExp := 5, vExp := 4 },
    { nCoeffs := [936576, 4728384, 10548384, 13714608, 11523472, 6545312, 2547488, 671376, 114736, 11488, 512], pExp := 1, qExp := 5, vExp := 2 },
    { nCoeffs := [-936576, -4728384, -10548384, -13714608, -11523472, -6545312, -2547488, -671376, -114736, -11488, -512], pExp := 1, qExp := 3, vExp := 4 },
    { nCoeffs := [144000, 926784, 2520992, 3862832, 3722192, 2370336, 1014528, 289104, 52624, 5536, 256], pExp := 1, qExp := 3, vExp := 2 },
    { nCoeffs := [6174720, 36108288, 94684160, 147066880, 150433664, 106460288, 53212672, 18793088, 4597376, 742144, 71168, 3072], pExp := 1, qExp := 3, vExp := 0 },
    { nCoeffs := [154368, 1381248, 4858304, 9399072, 11403808, 9238304, 5141952, 1977184, 516704, 87648, 8704, 384], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [617472, 3672576, 9804800, 15518976, 16189120, 11692064, 5967456, 2152960, 538208, 88800, 8704, 384], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [105840, 402696, 674004, 651454, 400784, 162782, 43656, 7456, 736, 32], pExp := 0, qExp := 4, vExp := 5 },
    { nCoeffs := [3492288, 17964000, 40633584, 53379720, 45202328, 25826568, 10096824, 2669936, 457440, 45888, 2048], pExp := 0, qExp := 4, vExp := 3 },
    { nCoeffs := [105840, 402696, 674004, 651454, 400784, 162782, 43656, 7456, 736, 32], pExp := 0, qExp := 4, vExp := 1 },
    { nCoeffs := [-63504, -237384, -389988, -369678, -222890, -88670, -23282, -3892, -376, -16], pExp := 0, qExp := 2, vExp := 5 },
    { nCoeffs := [-11912832, -65558880, -161363664, -234956776, -225248644, -149494248, -70171108, -23318056, -5380368, -821536, -74752, -3072], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [30810096, 177216696, 456520348, 695877874, 697955414, 483990946, 236920846, 81914828, 19613320, 3098608, 290816, 12288], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [5557248, 30892032, 76881408, 113208832, 109703872, 73520256, 34795456, 11636800, 2696448, 412480, 37504, 1536], pExp := 0, qExp := 0, vExp := 3 },
    { nCoeffs := [-19141632, -108601344, -275664896, -413571328, -407771328, -277620928, -133255296, -45117248, -10564800, -1630208, -149248, -6144], pExp := 0, qExp := 0, vExp := 1 }
  ]

def dualCertPp2Poly : DualCertPoly :=
  dualCertOfTerms dualCertPp2Terms

def dualCertPq2Terms : List DualCertTerm := [
    { nCoeffs := [926208, 5663232, 15548160, 25309184, 27148064, 20153184, 10566656, 3913568, 1003488, 169664, 17024, 768], pExp := 4, qExp := 2, vExp := 1 },
    { nCoeffs := [-4631040, -27853056, -75140736, -120071872, -126322752, -91900672, -47189312, -17106496, -4291328, -709632, -69632, -3072], pExp := 3, qExp := 1, vExp := 0 },
    { nCoeffs := [-42336, -186480, -356088, -387748, -265796, -119108, -34932, -6472, -688, -32], pExp := 2, qExp := 2, vExp := 5 },
    { nCoeffs := [-2083968, -10774080, -24653472, -32926640, -28459312, -16650208, -6683280, -1818608, -321248, -33280, -1536], pExp := 2, qExp := 2, vExp := 3 },
    { nCoeffs := [-505440, -2786544, -6830136, -9805316, -9131044, -5763076, -2496276, -732584, -139376, -15520, -768], pExp := 2, qExp := 2, vExp := 1 },
    { nCoeffs := [8335872, 47264256, 120162048, 181013504, 179690528, 123519456, 60033536, 20640992, 4922208, 775616, 72704, 3072], pExp := 2, qExp := 0, vExp := 3 },
    { nCoeffs := [-27786240, -160943616, -417823488, -642257536, -649956480, -454960384, -224897408, -78545536, -19001600, -3033600, -287744, -12288], pExp := 2, qExp := 0, vExp := 1 },
    { nCoeffs := [156096, 736032, 1512720, 1781528, 1326736, 648640, 208368, 42440, 4976, 256], pExp := 1, qExp := 3, vExp := 4 },
    { nCoeffs := [156096, 736032, 1512720, 1781528, 1326736, 648640, 208368, 42440, 4976, 256], pExp := 1, qExp := 3, vExp := 2 },
    { nCoeffs := [-378360, -1445160, -2487078, -2535236, -1686046, -756832, -228480, -44552, -5072, -256], pExp := 1, qExp := 1, vExp := 4 },
    { nCoeffs := [97848, -1038000, -5487154, -11445514, -13835758, -10908094, -5875888, -2188240, -555328, -91744, -8896, -384], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [1016208, 5374008, 12130124, 15297426, 11682552, 5304202, 1121264, -153832, -173040, -49760, -6848, -384], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [52920, 201348, 337002, 325727, 200392, 81391, 21828, 3728, 368, 16], pExp := 0, qExp := 2, vExp := 5 },
    { nCoeffs := [2037168, 10381992, 23236260, 30165014, 25209640, 14196006, 5462240, 1419520, 238656, 23456, 1024], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [52920, 201348, 337002, 325727, 200392, 81391, 21828, 3728, 368, 16], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [-6709248, -36877056, -90650112, -131703168, -125783424, -82983600, -38616192, -12682320, -2882112, -431808, -38400, -1536], pExp := 0, qExp := 0, vExp := 3 },
    { nCoeffs := [22138992, 124365456, 312733420, 464749112, 453638508, 305480832, 144866320, 48395744, 11165120, 1694592, 152320, 6144], pExp := 0, qExp := 0, vExp := 1 }
  ]

def dualCertPq2Poly : DualCertPoly :=
  dualCertOfTerms dualCertPq2Terms

def dualCertPv2Terms : List DualCertTerm := [
    { nCoeffs := [463104, 2831616, 7774080, 12654592, 13574032, 10076592, 5283328, 1956784, 501744, 84832, 8512, 384], pExp := 3, qExp := 1, vExp := 0 },
    { nCoeffs := [2778624, 16372224, 43280640, 67804672, 69986656, 49999200, 25238080, 9004192, 2225760, 363136, 35200, 1536], pExp := 2, qExp := 0, vExp := 1 },
    { nCoeffs := [-527472, -2881224, -7413300, -11697582, -12428360, -9255414, -4902096, -1840424, -479088, -82272, -8384, -384], pExp := 1, qExp := 3, vExp := 0 },
    { nCoeffs := [-134064, -559944, -1101908, -1343822, -1113628, -646070, -262148, -72688, -13088, -1376, -64], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [-154368, -840960, -2030720, -2864384, -2615088, -1615472, -684128, -196176, -36464, -3968, -192], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [52920, 201348, 337002, 325727, 200392, 81391, 21828, 3728, 368, 16], pExp := 0, qExp := 4, vExp := 3 },
    { nCoeffs := [52920, 201348, 337002, 325727, 200392, 81391, 21828, 3728, 368, 16], pExp := 0, qExp := 4, vExp := 1 },
    { nCoeffs := [-10584, -36036, -52986, -43951, -22498, -7279, -1454, -164, -8], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [-2647944, -14178420, -35810742, -55479583, -57841646, -42231495, -21910654, -8051876, -2050440, -344320, -34304, -1536], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [-578736, -4010976, -10862764, -16183788, -15106340, -9363572, -3934672, -1112112, -202944, -21632, -1024], pExp := 0, qExp := 0, vExp := 1 }
  ]

def dualCertPv2Poly : DualCertPoly :=
  dualCertOfTerms dualCertPv2Terms



def dualCertM20Poly : DualCertPoly :=
  (4 * dualCertN + 10) * (dualCertN + 2) ^ 2 * (dualCertN + 3) ^ 2 *
    (32 * dualCertN ^ 4 + 302 * dualCertN ^ 3 + 1037 * dualCertN ^ 2 +
      1530 * dualCertN + 813)

def dualCertM21Poly : DualCertPoly :=
  (dualCertN + 2) ^ 2 *
    (192 * dualCertN ^ 6 + 2984 * dualCertN ^ 5 + 19116 * dualCertN ^ 4 +
      64452 * dualCertN ^ 3 + 120256 * dualCertN ^ 2 + 117279 * dualCertN + 46476)

def dualCertM22Poly : DualCertPoly :=
  (dualCertN + 2) ^ 2 *
    (16 * dualCertN ^ 5 + 408 * dualCertN ^ 4 + 2912 * dualCertN ^ 3 +
      8884 * dualCertN ^ 2 + 12254 * dualCertN + 6240)

set_option maxHeartbeats 0 in
set_option maxRecDepth 100000 in
theorem dualCert_identity_two_poly :
    dualCertNormalize (dualCertDeltaPoly *
        (dualCertM20Poly * dualCertNextNumPoly 0 +
          dualCertM21Poly * dualCertNextNumPoly 1 +
          dualCertM22Poly * dualCertNextNumPoly 2 -
          dualCertLambdaPoly * dualCertCurNumPoly 2)) =
      dualCertNormalize (
      dualCertOpPPoly dualCertPp2Poly +
        dualCertOpQPoly dualCertPq2Poly +
        dualCertOpVPoly dualCertPv2Poly) := by
  native_decide

def dualCertPp2 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPp2Poly n p q v

def dualCertPq2 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPq2Poly n p q v

def dualCertPv2 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPv2Poly n p q v

def dualCertDPp2 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertPDeriv dualCertPp2Poly) n p q v

def dualCertDPq2 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertQDeriv dualCertPq2Poly) n p q v

def dualCertDPv2 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertVDeriv dualCertPv2Poly) n p q v

theorem dualCertPp2_hasDerivAt (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertPp2 n x q v) (dualCertDPp2 n p q v) p := by
  simpa [dualCertPp2, dualCertDPp2] using
    dualCertEval_hasDerivAt_p dualCertPp2Poly n p q v

theorem dualCertPq2_hasDerivAt (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertPq2 n p x v) (dualCertDPq2 n p q v) q := by
  simpa [dualCertPq2, dualCertDPq2] using
    dualCertEval_hasDerivAt_q dualCertPq2Poly n p q v

theorem dualCertPv2_hasDerivAt (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertPv2 n p q x) (dualCertDPv2 n p q v) v := by
  simpa [dualCertPv2, dualCertDPv2] using
    dualCertEval_hasDerivAt_v dualCertPv2Poly n p q v

macro "row2_eval_num" name:ident value:num : command =>
  `(@[simp] private theorem $name (n p q v : ℝ) :
      dualCertEval ($value : DualCertPoly) n p q v = ($value : ℝ) := by
    change dualCertEval (dualCertConst ($value : ℤ)) n p q v = _
    simp)

row2_eval_num dualCertEval_row2_32 32
row2_eval_num dualCertEval_row2_302 302
row2_eval_num dualCertEval_row2_1037 1037
row2_eval_num dualCertEval_row2_1530 1530
row2_eval_num dualCertEval_row2_813 813
row2_eval_num dualCertEval_row2_192 192
row2_eval_num dualCertEval_row2_2984 2984
row2_eval_num dualCertEval_row2_19116 19116
row2_eval_num dualCertEval_row2_64452 64452
row2_eval_num dualCertEval_row2_120256 120256
row2_eval_num dualCertEval_row2_117279 117279
row2_eval_num dualCertEval_row2_46476 46476
row2_eval_num dualCertEval_row2_16 16
row2_eval_num dualCertEval_row2_408 408
row2_eval_num dualCertEval_row2_2912 2912
row2_eval_num dualCertEval_row2_8884 8884
row2_eval_num dualCertEval_row2_12254 12254
row2_eval_num dualCertEval_row2_6240 6240

theorem dualCert_identity_two (n : ℕ) (p q v : ℝ) :
    (4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)) *
        ((positiveMatrix (n : ℤ) 2 0 : ℝ) * dualCertNextNum (n : ℝ) 0 p q v +
          (positiveMatrix (n : ℤ) 2 1 : ℝ) * dualCertNextNum (n : ℝ) 1 p q v +
          (positiveMatrix (n : ℤ) 2 2 : ℝ) * dualCertNextNum (n : ℝ) 2 p q v -
          dualCertLambda (n : ℝ) * dualCertCurNum (n : ℝ) 2 p q v) =
      dualCertOpP (n : ℝ) (dualCertPp2 (n : ℝ) p q v)
          (dualCertDPp2 (n : ℝ) p q v) p q v +
        dualCertOpQ (n : ℝ) (dualCertPq2 (n : ℝ) p q v)
          (dualCertDPq2 (n : ℝ) p q v) p q v +
        dualCertOpV (n : ℝ) (dualCertPv2 (n : ℝ) p q v)
          (dualCertDPv2 (n : ℝ) p q v) p q v := by
  have h := congrArg (fun P => dualCertEval P (n : ℝ) p q v)
    dualCert_identity_two_poly
  dsimp only at h
  rw [dualCertEval_normalize, dualCertEval_normalize] at h
  simp only [dualCertDeltaPoly, dualCertM20Poly, dualCertM21Poly, dualCertM22Poly,
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
    dualCertPp2, dualCertPq2, dualCertPv2,
    dualCertDPp2, dualCertDPq2, dualCertDPv2,
    positiveMatrix, Matrix.cons_val_two] using h

end RamanujanChallenge.P25

end
