import RamanujanChallenge.Problem25DualCertificateSemantics

noncomputable section

namespace RamanujanChallenge.P25

def dualCertPp1Terms : List DualCertTerm := [
    { nCoeffs := [-1518912, -6060960, -10667088, -10868952, -7066344, -3040072, -865512, -157248, -16544, -768], pExp := 3, qExp := 5, vExp := 4 },
    { nCoeffs := [-1518912, -6060960, -10667088, -10868952, -7066344, -3040072, -865512, -157248, -16544, -768], pExp := 3, qExp := 5, vExp := 2 },
    { nCoeffs := [1518912, 6060960, 10667088, 10868952, 7066344, 3040072, 865512, 157248, 16544, 768], pExp := 3, qExp := 3, vExp := 4 },
    { nCoeffs := [1518912, 6060960, 10667088, 10868952, 7066344, 3040072, 865512, 157248, 16544, 768], pExp := 3, qExp := 3, vExp := 2 },
    { nCoeffs := [-4959360, -19528704, -33896064, -34040672, -21799576, -9232360, -2585856, -461888, -47744, -2176], pExp := 2, qExp := 4, vExp := 3 },
    { nCoeffs := [-1486080, -7805952, -18162048, -24674816, -21695376, -12908416, -5266352, -1455424, -260864, -27392, -1280], pExp := 2, qExp := 4, vExp := 1 },
    { nCoeffs := [4959360, 19528704, 33896064, 34040672, 21799576, 9232360, 2585856, 461888, 47744, 2176], pExp := 2, qExp := 2, vExp := 3 },
    { nCoeffs := [297216, 1442304, 3055488, 3712768, 2853968, 1440096, 477232, 100192, 12096, 640], pExp := 2, qExp := 2, vExp := 1 },
    { nCoeffs := [1012608, 4378176, 8233248, 8868528, 6044512, 2707936, 798432, 149552, 16160, 768], pExp := 1, qExp := 5, vExp := 4 },
    { nCoeffs := [1012608, 4378176, 8233248, 8868528, 6044512, 2707936, 798432, 149552, 16160, 768], pExp := 1, qExp := 5, vExp := 2 },
    { nCoeffs := [-1012608, -4378176, -8233248, -8868528, -6044512, -2707936, -798432, -149552, -16160, -768], pExp := 1, qExp := 3, vExp := 4 },
    { nCoeffs := [374400, 1824192, 3813856, 4500304, 3311808, 1580256, 489984, 95376, 10592, 512], pExp := 1, qExp := 3, vExp := 2 },
    { nCoeffs := [7925760, 40310784, 90806272, 119493632, 101845760, 58802048, 23308288, 6267008, 1094400, 112128, 5120], pExp := 1, qExp := 3, vExp := 0 },
    { nCoeffs := [198144, 1622016, 4999936, 8251392, 8339424, 5486208, 2403584, 697088, 128672, 13696, 640], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [792576, 4110336, 9452032, 12708864, 11075712, 6542208, 2654048, 730496, 130592, 13696, 640], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [-211680, -805392, -1348008, -1302908, -801568, -325564, -87312, -14912, -1472, -64], pExp := 0, qExp := 4, vExp := 5 },
    { nCoeffs := [4558464, 19411776, 36052896, 38431536, 25961168, 11541104, 3379984, 629344, 67648, 3200], pExp := 0, qExp := 4, vExp := 3 },
    { nCoeffs := [-211680, -805392, -1348008, -1302908, -801568, -325564, -87312, -14912, -1472, -64], pExp := 0, qExp := 4, vExp := 1 },
    { nCoeffs := [127008, 474768, 779976, 739356, 445780, 177340, 46564, 7784, 752, 32], pExp := 0, qExp := 2, vExp := 5 },
    { nCoeffs := [-15088896, -71267904, -149119968, -182516592, -145001592, -78248592, -29081272, -7356688, -1213152, -117824, -5120], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [39755808, 198065808, 436637384, 561891356, 468085844, 264059068, 102251236, 26857064, 4582128, 458784, 20480], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [7133184, 34219008, 72597504, 89921024, 72139136, 39219008, 14647616, 3713856, 612160, 59264, 2560], pExp := 0, qExp := 0, vExp := 3 },
    { nCoeffs := [-24569856, -120683520, -261839872, -331076096, -270568576, -149509824, -56626048, -14526784, -2417408, -235776, -10240], pExp := 0, qExp := 0, vExp := 1 }
  ]

def dualCertPp1Poly : DualCertPoly :=
  dualCertOfTerms dualCertPp1Terms

def dualCertPq1Terms : List DualCertTerm := [
    { nCoeffs := [1188864, 6363648, 15106560, 20962048, 18841408, 11468320, 4789120, 1355232, 248768, 26752, 1280], pExp := 4, qExp := 2, vExp := 1 },
    { nCoeffs := [-5944320, -31223808, -72648192, -98699264, -86781504, -51633664, -21065408, -5821696, -1043456, -109568, -5120], pExp := 3, qExp := 1, vExp := 0 },
    { nCoeffs := [84672, 372960, 712176, 775496, 531592, 238216, 69864, 12944, 1376, 64], pExp := 2, qExp := 2, vExp := 5 },
    { nCoeffs := [-2674944, -11791872, -22654848, -24975616, -17450192, -8025888, -2432752, -469024, -52224, -2560], pExp := 2, qExp := 2, vExp := 3 },
    { nCoeffs := [-509760, -2511648, -5398800, -6650040, -5176344, -2641976, -884600, -187440, -22816, -1216], pExp := 2, qExp := 2, vExp := 1 },
    { nCoeffs := [10699776, 52517376, 114203136, 145212160, 119752000, 67003936, 25782784, 6741600, 1146944, 114688, 5120], pExp := 2, qExp := 0, vExp := 3 },
    { nCoeffs := [-35665920, -179417088, -399541248, -519563264, -437511552, -249544960, -97718912, -25958656, -4479488, -453632, -20480], pExp := 2, qExp := 0, vExp := 1 },
    { nCoeffs := [168768, 673440, 1147728, 1095512, 642248, 237240, 53992, 6928, 384], pExp := 1, qExp := 3, vExp := 4 },
    { nCoeffs := [168768, 673440, 1147728, 1095512, 642248, 237240, 53992, 6928, 384], pExp := 1, qExp := 3, vExp := 2 },
    { nCoeffs := [275760, 744816, 800988, 411904, 76372, -20856, -13768, -2704, -192], pExp := 1, qExp := 1, vExp := 4 },
    { nCoeffs := [-1292976, -5298432, -10471820, -13031700, -11067028, -6549436, -2688160, -747520, -134016, -13952, -640], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [1740384, 7211664, 12761896, 12412956, 6978472, 2011340, 18184, -199536, -69696, -10624, -640], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [-105840, -402696, -674004, -651454, -400784, -162782, -43656, -7456, -736, -32], pExp := 0, qExp := 2, vExp := 5 },
    { nCoeffs := [2659104, 11196912, 20533848, 21582596, 14355072, 6274148, 1803792, 329184, 34624, 1600], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [-105840, -402696, -674004, -651454, -400784, -162782, -43656, -7456, -736, -32], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [-8434944, -40041984, -83954304, -102631040, -81151440, -43423424, -15939376, -3966080, -640576, -60672, -2560], pExp := 0, qExp := 0, vExp := 3 },
    { nCoeffs := [27981216, 136898208, 295079816, 369875488, 299101832, 163259056, 60976224, 15399552, 2518272, 240896, 10240], pExp := 0, qExp := 0, vExp := 1 }
  ]

def dualCertPq1Poly : DualCertPoly :=
  dualCertOfTerms dualCertPq1Terms

def dualCertPv1Terms : List DualCertTerm := [
    { nCoeffs := [594432, 3181824, 7553280, 10481024, 9420704, 5734160, 2394560, 677616, 124384, 13376, 640], pExp := 3, qExp := 1, vExp := 0 },
    { nCoeffs := [3566592, 18298368, 41605632, 55291136, 47612864, 27784928, 11135168, 3027744, 534784, 55424, 2560], pExp := 2, qExp := 0, vExp := 1 },
    { nCoeffs := [-241056, -1875312, -5396184, -8381412, -8085784, -5157044, -2224696, -644848, -120640, -13184, -640], pExp := 1, qExp := 3, vExp := 0 },
    { nCoeffs := [-14112, -283248, -877720, -1255332, -1031664, -527380, -171120, -34368, -3904, -192], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [-198144, -928512, -1898752, -2227840, -1655008, -808048, -259488, -52880, -6208, -320], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [-105840, -402696, -674004, -651454, -400784, -162782, -43656, -7456, -736, -32], pExp := 0, qExp := 4, vExp := 3 },
    { nCoeffs := [-105840, -402696, -674004, -651454, -400784, -162782, -43656, -7456, -736, -32], pExp := 0, qExp := 4, vExp := 1 },
    { nCoeffs := [21168, 72072, 105972, 87902, 44996, 14558, 2908, 328, 16], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [-1184112, -9183960, -25997556, -39559762, -37322916, -23258274, -9798244, -2772728, -506352, -54016, -2560], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [-2050848, -8286336, -14985416, -15899928, -10885128, -4973960, -1513264, -294912, -33344, -1664], pExp := 0, qExp := 0, vExp := 1 }
  ]

def dualCertPv1Poly : DualCertPoly :=
  dualCertOfTerms dualCertPv1Terms



def dualCertM10Poly : DualCertPoly :=
  (dualCertN + 2) ^ 2 * (dualCertN + 3) ^ 2 * (4 * dualCertN + 10) *
    (48 * dualCertN ^ 3 + 386 * dualCertN ^ 2 + 1017 * dualCertN + 879)

def dualCertM11Poly : DualCertPoly :=
  (dualCertN + 2) ^ 2 *
    (272 * dualCertN ^ 5 + 3848 * dualCertN ^ 4 + 21732 * dualCertN ^ 3 +
      61184 * dualCertN ^ 2 + 85761 * dualCertN + 47808)

def dualCertM12Poly : DualCertPoly :=
  (dualCertN + 2) ^ 2 *
    (320 * dualCertN ^ 3 + 2540 * dualCertN ^ 2 + 6610 * dualCertN + 5640)

set_option maxHeartbeats 0 in
set_option maxRecDepth 100000 in
theorem dualCert_identity_one_poly :
    dualCertNormalize (dualCertDeltaPoly *
        (dualCertM10Poly * dualCertNextNumPoly 0 +
          dualCertM11Poly * dualCertNextNumPoly 1 +
          dualCertM12Poly * dualCertNextNumPoly 2 -
          dualCertLambdaPoly * dualCertCurNumPoly 1)) =
      dualCertNormalize (
      dualCertOpPPoly dualCertPp1Poly +
        dualCertOpQPoly dualCertPq1Poly +
        dualCertOpVPoly dualCertPv1Poly) := by
  rfl

def dualCertPp1 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPp1Poly n p q v

def dualCertPq1 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPq1Poly n p q v

def dualCertPv1 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPv1Poly n p q v

def dualCertDPp1 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertPDeriv dualCertPp1Poly) n p q v

def dualCertDPq1 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertQDeriv dualCertPq1Poly) n p q v

def dualCertDPv1 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertVDeriv dualCertPv1Poly) n p q v

theorem dualCertPp1_hasDerivAt (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertPp1 n x q v) (dualCertDPp1 n p q v) p := by
  simpa [dualCertPp1, dualCertDPp1] using
    dualCertEval_hasDerivAt_p dualCertPp1Poly n p q v

theorem dualCertPq1_hasDerivAt (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertPq1 n p x v) (dualCertDPq1 n p q v) q := by
  simpa [dualCertPq1, dualCertDPq1] using
    dualCertEval_hasDerivAt_q dualCertPq1Poly n p q v

theorem dualCertPv1_hasDerivAt (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertPv1 n p q x) (dualCertDPv1 n p q v) v := by
  simpa [dualCertPv1, dualCertDPv1] using
    dualCertEval_hasDerivAt_v dualCertPv1Poly n p q v

macro "row1_eval_num" name:ident value:num : command =>
  `(@[simp] private theorem $name (n p q v : ℝ) :
      dualCertEval ($value : DualCertPoly) n p q v = ($value : ℝ) := by
    change dualCertEval (dualCertConst ($value : ℤ)) n p q v = _
    simp)

row1_eval_num dualCertEval_row1_48 48
row1_eval_num dualCertEval_row1_386 386
row1_eval_num dualCertEval_row1_1017 1017
row1_eval_num dualCertEval_row1_879 879
row1_eval_num dualCertEval_row1_272 272
row1_eval_num dualCertEval_row1_3848 3848
row1_eval_num dualCertEval_row1_21732 21732
row1_eval_num dualCertEval_row1_61184 61184
row1_eval_num dualCertEval_row1_85761 85761
row1_eval_num dualCertEval_row1_47808 47808
row1_eval_num dualCertEval_row1_320 320
row1_eval_num dualCertEval_row1_2540 2540
row1_eval_num dualCertEval_row1_6610 6610
row1_eval_num dualCertEval_row1_5640 5640

theorem dualCert_identity_one (n : ℕ) (p q v : ℝ) :
    (4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)) *
        ((positiveMatrix (n : ℤ) 1 0 : ℝ) * dualCertNextNum (n : ℝ) 0 p q v +
          (positiveMatrix (n : ℤ) 1 1 : ℝ) * dualCertNextNum (n : ℝ) 1 p q v +
          (positiveMatrix (n : ℤ) 1 2 : ℝ) * dualCertNextNum (n : ℝ) 2 p q v -
          dualCertLambda (n : ℝ) * dualCertCurNum (n : ℝ) 1 p q v) =
      dualCertOpP (n : ℝ) (dualCertPp1 (n : ℝ) p q v)
          (dualCertDPp1 (n : ℝ) p q v) p q v +
        dualCertOpQ (n : ℝ) (dualCertPq1 (n : ℝ) p q v)
          (dualCertDPq1 (n : ℝ) p q v) p q v +
        dualCertOpV (n : ℝ) (dualCertPv1 (n : ℝ) p q v)
          (dualCertDPv1 (n : ℝ) p q v) p q v := by
  have h := congrArg (fun P => dualCertEval P (n : ℝ) p q v)
    dualCert_identity_one_poly
  dsimp only at h
  rw [dualCertEval_normalize, dualCertEval_normalize] at h
  simp only [dualCertDeltaPoly, dualCertM10Poly, dualCertM11Poly, dualCertM12Poly,
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
    dualCertPp1, dualCertPq1, dualCertPv1,
    dualCertDPp1, dualCertDPq1, dualCertDPv1,
    positiveMatrix, Matrix.cons_val_two] using h

end RamanujanChallenge.P25

end
