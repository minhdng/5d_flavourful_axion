(* ::Package:: *)

(* :Title: WarpFlavourAxion *)
(* :Context: WarpFlavourAxion` *)
(* :Author: Minh D. Nguyen, Department of Physics, University of Minnesota *)
(* :Summary: Supplementary Mathematica package *)
(* :Copyright: Copyright 2020, Minh D. Nguyen *)
(* :Mathematica Version: 12.0 *)
(* :History:
  V1.00 
  V1.01 Include leptons
  V1.01 Update definition for constraint 3 for PMNS
 *)
(* :Notes:

 *)

(****************************************************************
 * Begin package
 ****************************************************************)


BeginPackage["WarpFlavourAxion`"];


(****************************************************************
 * Messages
 ****************************************************************)


(* Usage messages *)


(* constants *)
$UpMass::usage = "$UpMass return a vector of {up, charm, top} quark masses of in GeV";
$UpMassError::usage = "$UpMassError returns uncertainty of UpMass in GeV";
$DownMass::usage = "$DownMass[1], $DownMass[2], $DownMass[3] return masses of down, strange and bottom quark in GeV";
$DownMassError::usage = "$DownMassError returns uncertainty of DownMass in GeV";
$LeptonMass::usage = "$LeptonMass[1], $LeptonMass[2], $LeptonMass[3] return masses of electron, muon and tau in GeV";
$LeptonMassError::usage = "$LeptonMassError returns uncertainty of LeptonMass in GeV";
$Yu::usage = "$Yu returns a vector of {up, charm, top} Yukawa couplings at 10^10 GeV for tan beta = 3";
$YuError::usage = "$YuError returns a vector of {up, charm, top} Yukawa coupling uncertainty at 10^10 GeV for tan beta = 3";
$Yd::usage = "$Yd returns a vector of {down, strange, bottom} Yukawa couplings at 10^10 GeV for tan beta = 3";
$YdError::usage = "$YdError returns a vector of {down, strange, bottom} Yukawa coupling uncertainty at 10^10 GeV for tan beta = 3";
$Ye::usage = "$Ye returns a vector of {electron, mu, tau} Yukawa couplings at 10^10 GeV for tan beta = 3";
$YeError::usage = "$YeError returns a vector of {electron, mu, tau} Yukawa coupling uncertainty at 10^10 GeV for tan beta = 3";
$CkmLambda::usage = "0.22453"; 
$CkmLambdaError::usage = "0.00044";
$CkmA::usage = "0.836"; 
$CkmAError::usage = "0.015";
$CkmRhoBar::usage = "0.122"; 
$CkmRhoBarError::usage = "0.018";
$CkmEtaBar::usage = "0.355"; 
$CkmEtaBarError::usage = "0.012";
$PmnsLambda::usage = "{0.555, 0.555}"; 
$PmnsLambdaError::usage = "{0.015, 0.015}";
$PmnsA::usage = "{2.15, 2.19}"; 
$PmnsAError::usage = "{0.18, 0.31}";
$PmnsRhoBar::usage = "{-0.16, -0.12}"; 
$PmnsRhoBarError::usage = "{0.14, 0.15}";
$PmnsEtaBar::usage = "{0.477, 0.488}"; 
$PmnsEtaBarError::usage = "{0.061, 0.053}";

(* Utitlities *)

RandomMatrix::usage = "randomMatrix[min_, max_, dim_]";
RandomUnitaryMatrix::usage = "randomUnitaryMatrix[min_, max_, dim_]";

(* Constraint 1 *)

RhoEtaBar::usage = "CkmRhoEtaBar[yU_,yD_,yUMinor_,yDMinor_]";
CkmQ::usage = "CkmQ[yU, yD, yUMinor, yDMinor]";
PmnsQ::usage = "PmnsQ[yU_,yD_,yUMinor_,yDMinor_]";

(* Constraint 2 *)

QuarkEffMass::usage = "QuarkEffMass[yU_,yD_,yUMinor_,yDMinor_,v_,\[Beta]_]";
LeptonEffMass::usage = "LeptonEffMass[yU_,yD_,yUMinor_,yDMinor_,v_,\[Beta]_]";

QuarkEffYukawa::usage = "QuarkEffYukawa[yU_, yD_, yUMinor_, yDMinor_]";
LeptonEffYukawa::usage = "LeptonEffYukawa[yN_, yE_, yNMinor_, yEMinor_]";

QuarkProfileBoundedQ::usage = "QuarkProfileBoundedQ[yU_,yD_,yUMinor_,yDMinor_,v_,\[Beta]_,threshold_]";
LeptonProfileBoundedQ::usage = "LeptonProfileBoundedQ[yU_,yD_,yUMinor_,yDMinor_,v_,\[Beta]_,threshold_]";

(* Constraint 3 *)

QuarkAMatrices::usage = "QuarkAMatrices[yU_, yD_, yUMinor_, yDMinor_, v_, \[Beta]_]";
LeptonAMatrices::usage = "LeptonAMatrices[yN_, yE_, yNMinor_, yEMinor_, v_, \[Beta]_]";
UnitaryQ::usage = "UnitaryQ[mat_,thres_]";
QuarkAMatricesUnitaryQ::usage = "QuarkAMatricesUnitaryQ[yU_,yD_,yUMinor_,yDMinor_,thres_, v_, \[Beta]_]";
LeptonAMatricesUnitaryQ::usage = "LeptonAMatricesUnitaryQ[yN_,yE_,yNMinor_,yEMinor_,thres_, v_, \[Beta]_]";
QuarkConstraintSummary::usage = "QuarkConstraintSummary[yU_, yD_, yUMinor_, yDMinor_]";

(* Analytical profile and profile overlap functions *)

FermionProfile::usage = "FermionProfile[c_?NumericQ, z_, zir_:10^8]";
FermionProfileUVOverlap::usage = "FermionProfileUVOverlap[cL_?NumericQ, cR_?NumericQ, zir_:10^8]";
FermionProfileBulkOverlap::usage = "FermionProfileBulkOverlap[cL_?NumericQ,cR_?NumericQ,zir_:10^8]";
TurnRight::usage = "TurnRight[listLR_]";
TurnLeft::usage = "TurnLeft[listPM_]";
FermionProfileUVOverlapCM::usage = "FermionProfileUVOverlapCM[cP_?NumericQ, cM_?NumericQ, zir_:10^8]";
FermionProfileBulkOverlapCM::usage = "FermionProfileBulkOverlapCM[cP_?NumericQ, cM_?NumericQ, zir_:10^8]";
FermionAxionOverlap::usage = "FermionAxionOverlap[c_,\[CapitalDelta]_,zir_:10^8]";

(* Plot Utilities *)

FixedRange::usage = "FixedRange[min_, max_, length_]";
ListMirror::usage = "ListMirror[list_]";


(****************************************************************
 * Begin private context
 ****************************************************************)


Begin["`Private`"];


(****************************************************************
 * Constants and experimental values
 ****************************************************************)


(* Quark and lepton masses in GeV *)


$UpMass = {0.0023, 1.275, 173.21};
$UpMassError = {0.0012, 0.025, 1.22};
$DownMass = {0.0048 , 0.095, 4.18}; 
$DownMassError = {0.0008, 0.005, 0.03};
$LeptonMass = {0.00051099895000, 0.1056583755, 1.77686};
$LeptonMassError = {0.00000000000015, 0.0000000023, 0.00012};


(* Quark and lepton Yukawa at 10^10 GeV, tan \[Beta] = 3 
	See running.nb for detail calculation
*)


$Yu = {3.29456*10^-6, 0.00165737, 0.497757};
$YuError = {1.41*10^-9, 7.4*10^-7, 2.*10^-6};
$Yd = {0.0000244146, 0.000486184, 0.0237974};
$YdError = {5.*10^-10, 2.21*10^-7, 4.*10^-7};
$Ye = {2.96535*10^-6, 0.000624451, 0.0106073};
$YeError = {3.*10^-11, 3.*10^-9, 1.*10^-7};


(* Wolfenstefan parameters for CKM *)


$CkmLambda = 0.22453; $CkmLambdaError = 0.00044;
$CkmA = 0.836; $CkmAError = 0.015;
$CkmRhoBar = 0.122; $CkmRhoBarError = 0.018;
$CkmEtaBar = 0.355; $CkmEtaBarError = 0.012;


(* Wolfenstefan parameters for PMNS in {normal, inverted} hierarchy, see ???.nb for details *)


$PmnsLambda = {0.555, 0.555}; $PmnsLambdaError = {0.015, 0.015};
$PmnsA = {2.15, 2.19}; $PmnsAError = {0.18, 0.31};
$PmnsRhoBar = {-0.16, -0.12}; $PmnsRhoBarError = {0.14, 0.15};
$PmnsEtaBar = {0.477, 0.488}; $PmnsEtaBarError = {0.061, 0.053};


(****************************************************************
 * Constraint functions
 ****************************************************************)


RandomMatrix[dim_, min_:0.001, max_:3.] := RandomReal[{min, max},{dim, dim}] Exp[I RandomReal[{0., 2\[Pi]},{dim, dim}]]


RandomUnitaryMatrix[dim_, min_:0.001, max_:3.] := QRDecomposition[RandomMatrix[dim, min, max]][[1]];


(* Constraints 1. *)


RhoEtaBar[yU_,yD_,yUMinor_,yDMinor_]:= 
Through[
	{Re,Im}
	[ 
		(yD[[3,3]] yUMinor[[3,1]]-yD[[2,3]] yUMinor[[2,1]]+yD[[1,3]] yUMinor[[1,1]])/
		(
			yD[[3,3]] yUMinor[[1,1]] 
			(yD[[2,3]]/yD[[3,3]]-yU[[2,3]]/yU[[3,3]]) 
			(yDMinor[[2,1]]/yDMinor[[1,1]]-yUMinor[[2,1]]/yUMinor[[1,1]])
		) 
	]
];


CkmQ[yU_,yD_,yUMinor_,yDMinor_]:=(
(Norm[ RhoEtaBar[yU, yD, yUMinor, yDMinor][[1]] - $CkmRhoBar] < $CkmRhoBarError) 
&& 
(Norm[ RhoEtaBar[yU, yD, yUMinor, yDMinor][[2]] + $CkmEtaBar] < $CkmEtaBarError)
);


PmnsQ[yN_, yE_, yNMinor_, yEMinor_, ordering_]:=(
(Norm[ RhoEtaBar[yN, yE, yNMinor, yEMinor][[1]] - $PmnsRhoBar[[ordering]]] < $PmnsRhoBarError[[ordering]])
&&
(Norm[RhoEtaBar[yN, yE, yNMinor, yEMinor][[2]]+ $PmnsEtaBar[[ordering]]] < $PmnsEtaBarError[[ordering]])
);


(* Constraint 2.  *)


QuarkEffMass[yU_, yD_, yUMinor_, yDMinor_, v_:246, \[Beta]_:ArcTan[3]]:= <|
"u"->(Sqrt[2]$UpMass[[1]])/(v Sin[\[Beta]]) Norm[yUMinor[[1,1]]]/Norm[Det[yU]],
"c"->(Sqrt[2]$UpMass[[2]])/(v Sin[\[Beta]]) Norm[yU[[3,3]]]/Norm[yUMinor[[1,1]]],
"t"->(Sqrt[2]$UpMass[[3]])/(v Sin[\[Beta]]) 1/Norm[yU[[3,3]]],
"d"->(Sqrt[2]$DownMass[[1]])/(v Cos[\[Beta]]) Norm[yDMinor[[1,1]]]/Norm[Det[yD]],
"s"->(Sqrt[2]$DownMass[[2]])/(v Cos[\[Beta]]) Norm[yD[[3,3]]]/Norm[yDMinor[[1,1]]],
"b"->(Sqrt[2]$DownMass[[3]])/(v Cos[\[Beta]]) 1/Norm[yD[[3,3]]]
|>;


LeptonEffMass[yN_, yE_, yNMinor_, yEMinor_, v_:246, \[Beta]_:ArcTan[3]]:= <|
"e"->(Sqrt[2]$LeptonMass[[1]])/(v Cos[\[Beta]]) Norm[yEMinor[[1,1]]]/Norm[Det[yE]],
"mu"->(Sqrt[2]$LeptonMass[[2]])/(v Cos[\[Beta]]) Norm[yE[[3,3]]]/Norm[yEMinor[[1,1]]],
"tau"->(Sqrt[2]$LeptonMass[[3]])/(v Cos[\[Beta]]) 1/Norm[yE[[3,3]]]
|>;


QuarkEffYukawa[yU_, yD_, yUMinor_, yDMinor_]:= <|
"u"-> $Yu[[1]] Norm[yUMinor[[1,1]]]/Norm[Det[yU]],
"c"-> $Yu[[2]] Norm[yU[[3,3]]]/Norm[yUMinor[[1,1]]],
"t"-> $Yu[[3]]/Norm[yU[[3,3]]],
"d"-> $Yd[[1]] Norm[yDMinor[[1,1]]]/Norm[Det[yD]],
"s"-> $Yd[[2]] Norm[yD[[3,3]]]/Norm[yDMinor[[1,1]]],
"b"-> $Yd[[3]]/Norm[yD[[3,3]]]
|>;


LeptonEffYukawa[yN_, yE_, yNMinor_, yEMinor_]:= <|
"e"-> $Ye[[1]]Norm[yEMinor[[1,1]]]/Norm[Det[yE]],
"mu"-> $Ye[[2]] Norm[yE[[3,3]]]/Norm[yEMinor[[1,1]]],
"tau"-> $Ye[[3]]/Norm[yE[[3,3]]]
|>;


(* Old Boolean constraint 2 using QuarkEffMass[] *)
(* 
QuarkProfileBoundedQ[yU_, yD_, yUMinor_, yDMinor_, v_:246, \[Beta]_:ArcTan[3], threshold_:0.95]:= 
	AllTrue[ Values[QuarkEffMass[yU, yD, yUMinor, yDMinor, v, \[Beta]]], (#<threshold)& ];
	*)


(* Old Boolean constraint 2 using LeptonEffMass[] *)
(*
LeptonProfileBoundedQ[yN_, yE_, yNMinor_, yEMinor_, v_:246, \[Beta]_:ArcTan[3], threshold_:0.95]:= 
	AllTrue[ Values[LeptonEffMass[yN, yE, yNMinor, yEMinor, v, \[Beta]]], (#<threshold)& ];
	*)


QuarkProfileBoundedQ[yU_, yD_, yUMinor_, yDMinor_, threshold_:0.95]:= 
	AllTrue[ Values[QuarkEffYukawa[yU, yD, yUMinor, yDMinor]], (#<threshold)& ];
LeptonProfileBoundedQ[yN_, yE_, yNMinor_, yEMinor_, threshold_:0.95]:= 
	AllTrue[ Values[LeptonEffYukawa[yN, yE, yNMinor, yEMinor]], (#<threshold)& ];


(* Constraint 3 *)


QuarkAMatrices[yU_, yD_, yUMinor_, yDMinor_, v_:246, \[Beta]_:ArcTan[3]]:= Module[{fQ, fu, fd, AuL, AuR, AdL, AdR, phaseU, phaseD},
fQ = {
		$CkmLambda/Norm[yDMinor[[2,1]]/yDMinor[[1,1]]-yUMinor[[2,1]]/yUMinor[[1,1]]], 
		1, 
		Norm[yD[[2,3]]/yD[[3,3]]-yU[[2,3]]/yU[[3,3]]]/($CkmA $CkmLambda^2)
	};
fu = (QuarkEffYukawa[yU, yD, yUMinor, yDMinor][#]&/@{"u","c","t"})/fQ;
fd = (QuarkEffYukawa[yU, yD, yUMinor, yDMinor][#]&/@{"d","s","b"})/fQ;
phaseU = ({
 {E^(-I(Arg[Det[yU]]-Arg[yUMinor[[1,1]]])), 0, 0},
 {0, E^(-I(Arg[yUMinor[[1,1]]]- Arg[yU[[3,3]]])), 0},
 {0, 0, E^(-I Arg[yU[[3,3]]])}
});
phaseD = ({
 {E^(-I(Arg[Det[yD]]-Arg[yDMinor[[1,1]]])), 0, 0},
 {0, E^(-I(Arg[yDMinor[[1,1]]]- Arg[yD[[3,3]]])), 0},
 {0, 0, E^(-I Arg[yD[[3,3]]])}
});
AuL = ({
 {1, 
 yUMinor[[2,1]]/yUMinor[[1,1]] fQ[[1]]/fQ[[2]], 
 yU[[1,3]]/yU[[3,3]] fQ[[1]]/fQ[[3]]
 },
 {-((yUMinor[[2,1]]//Conjugate)/(yUMinor[[1,1]]//Conjugate)) fQ[[1]]/fQ[[2]], 
 1, 
 yU[[2,3]]/yU[[3,3]] fQ[[2]]/fQ[[3]]
 },
 {(yUMinor[[3,1]]//Conjugate)/(yUMinor[[1,1]]//Conjugate) fQ[[1]]/fQ[[3]], 
 -((yU[[2,3]]//Conjugate)/(yU[[3,3]]//Conjugate)) fQ[[2]]/fQ[[3]], 
 1}
});
AuR = phaseU.({
 {1, 
 (yUMinor[[1,2]]//Conjugate)/(yUMinor[[1,1]]//Conjugate) fu[[1]]/fu[[2]], 
 (yU[[3,1]]//Conjugate)/(yU[[3,3]]//Conjugate) fu[[1]]/fu[[3]]
 },
 {-(yUMinor[[1,2]]/yUMinor[[1,1]]) fu[[1]]/fu[[2]], 
 1, 
 (yU[[3,2]]//Conjugate)/(yU[[3,3]]//Conjugate) fu[[2]]/fu[[3]]
 },
 {yUMinor[[1,3]]/yUMinor[[1,1]] fu[[1]]/fu[[3]], 
 -(yU[[3,2]]/yU[[3,3]]) fu[[2]]/fu[[3]], 
 1}
});
AdL = ({
 {1, 
 yDMinor[[2,1]]/yDMinor[[1,1]] fQ[[1]]/fQ[[2]], 
 yD[[1,3]]/yD[[3,3]] fQ[[1]]/fQ[[3]]
 },
 {-((yDMinor[[2,1]]//Conjugate)/(yDMinor[[1,1]]//Conjugate)) fQ[[1]]/fQ[[2]], 
 1, 
 yD[[2,3]]/yD[[3,3]] fQ[[2]]/fQ[[3]]
 },
 {(yDMinor[[3,1]]//Conjugate)/(yDMinor[[1,1]]//Conjugate) fQ[[1]]/fQ[[3]], 
 -((yD[[2,3]]//Conjugate)/(yD[[3,3]]//Conjugate)) fQ[[2]]/fQ[[3]], 
 1}
});
AdR = phaseD.({
 {1, 
 (yDMinor[[1,2]]//Conjugate)/(yDMinor[[1,1]]//Conjugate) fd[[1]]/fd[[2]], 
 (yD[[3,1]]//Conjugate)/(yD[[3,3]]//Conjugate) fd[[1]]/fd[[3]]
 },
 {-(yDMinor[[1,2]]/yDMinor[[1,1]]) fd[[1]]/fd[[2]], 
 1, 
 (yD[[3,2]]//Conjugate)/(yD[[3,3]]//Conjugate) fd[[2]]/fd[[3]]
 },
 {yDMinor[[1,3]]/yDMinor[[1,1]] fd[[1]]/fd[[3]], 
 -(yD[[3,2]]/yD[[3,3]]) fd[[2]]/fd[[3]], 
 1}
});
{AuL, AuR, AdL, AdR}
];


LeptonAMatrices[yN_, yE_, yNMinor_, yEMinor_, ordering_, v_:246, \[Beta]_:ArcTan[3]]:= Module[{fl, fe, AeL, AeR, phaseE},
fl = {
		$PmnsLambda[[ordering]]/Norm[yEMinor[[2,1]]/yEMinor[[1,1]]-yNMinor[[2,1]]/yNMinor[[1,1]]], 
		1, 
		Norm[yE[[2,3]]/yE[[3,3]]-yN[[2,3]]/yN[[3,3]]]/($PmnsA[[ordering]] $PmnsLambda[[ordering]]^2)
	};
fe = (LeptonEffYukawa[yN, yE, yNMinor, yEMinor][#]&/@{"e","mu","tau"}) / fl;
phaseE = ({
 {E^(-I(Arg[Det[yE]]-Arg[yEMinor[[1,1]]])), 0, 0},
 {0, E^(-I(Arg[yEMinor[[1,1]]]- Arg[yE[[3,3]]])), 0},
 {0, 0, E^(-I Arg[yE[[3,3]]])}
});
AeL = ({
 {1, 
 yEMinor[[2,1]]/yEMinor[[1,1]] fl[[1]]/fl[[2]], 
 yE[[1,3]]/yE[[3,3]] fl[[1]]/fl[[3]]
 },
 {-((yEMinor[[2,1]]//Conjugate)/(yEMinor[[1,1]]//Conjugate)) fl[[1]]/fl[[2]], 
 1, 
 yE[[2,3]]/yE[[3,3]] fl[[2]]/fl[[3]]
 },
 {(yEMinor[[3,1]]//Conjugate)/(yEMinor[[1,1]]//Conjugate) fl[[1]]/fl[[3]], 
 -((yE[[2,3]]//Conjugate)/(yE[[3,3]]//Conjugate)) fl[[2]]/fl[[3]], 
 1}
});
AeR = phaseE.({
 {1, 
 (yEMinor[[1,2]]//Conjugate)/(yEMinor[[1,1]]//Conjugate) fe[[1]]/fe[[2]], 
 (yE[[3,1]]//Conjugate)/(yE[[3,3]]//Conjugate) fe[[1]]/fe[[3]]
 },
 {-(yEMinor[[1,2]]/yEMinor[[1,1]]) fe[[1]]/fe[[2]], 
 1, 
 (yE[[3,2]]//Conjugate)/(yE[[3,3]]//Conjugate) fe[[2]]/fe[[3]]
 },
 {yEMinor[[1,3]]/yEMinor[[1,1]] fe[[1]]/fe[[3]], 
 -(yE[[3,2]]/yE[[3,3]]) fe[[2]]/fe[[3]], 
 1}
});
{AeL, AeR}
];


(* Different from UnitaryMatrixQ of Mathematica *)
UnitaryQ[mat_,thres_]:= AllTrue[ Abs[Flatten[mat.ConjugateTranspose[mat]-IdentityMatrix[3]]], (#<thres)& ] 


QuarkAMatricesUnitaryQ[yU_, yD_, yUMinor_, yDMinor_, v_:246, \[Beta]_:ArcTan[3], thres_:0.2]:= 
	AllTrue[ {1,2,3,4}, UnitaryQ[QuarkAMatrices[yU,yD,yUMinor,yDMinor, v, \[Beta]][[#]], thres]& ]


LeptonAMatricesUnitaryQ[yN_, yE_, yNMinor_, yEMinor_, ordering_, v_:246, \[Beta]_:ArcTan[3], thres_:0.2]:= 
	AllTrue[ {1,2}, UnitaryQ[LeptonAMatrices[yN, yE, yNMinor, yEMinor, ordering, v, \[Beta]][[#]], thres]& ]


(* Script running through all three constraint *)


QuarkConstraintSummary[yU_, yD_, yUMinor_, yDMinor_]:= Module[{AdR},
	Print["Constraint 1 result for quark sector"];
	Print["\!\(\*SubscriptBox[\(Y\), \(u\)]\) = ",MatrixForm[yU]];
	Print["\!\(\*SubscriptBox[\(Y\), \(d\)]\) = ",MatrixForm[yD]];
	Print["The resulted \[Rho]bar = ",RhoEtaBar[yU,yD,yUMinor,yDMinor][[1]], ", \[Eta]bar = ", - RhoEtaBar[yU,yD,yUMinor,yDMinor][[2]]];
	Print["Experimental \[Rho]bar = ",$CkmRhoBar, " \[PlusMinus] ", $CkmRhoBarError,", \[Eta]bar = ",$CkmEtaBar, " \[PlusMinus] ", $CkmEtaBarError];
	Print["CkmQ[] returns ",CkmQ[yU,yD,yUMinor,yDMinor]];
	Print["Constraint 2 result for quark sector"];
	Print["\!\(\*SubscriptBox[\(Y\), \(u\)]\) = ",MatrixForm[yU]];
	Print["\!\(\*SubscriptBox[\(Y\), \(d\)]\) = ",MatrixForm[yD]];
	Print["The resulted effective quark masses ", QuarkEffYukawa[yU,yD,yUMinor,yDMinor]];
	Print["QuarkProfileBoundedQ[] returns ",QuarkProfileBoundedQ[yU,yD,yUMinor,yDMinor]];
	Print["Constraint 3 result for quark sector"];
	Print["\!\(\*SubscriptBox[\(Y\), \(u\)]\) = ",MatrixForm[yU]];
	Print["\!\(\*SubscriptBox[\(Y\), \(d\)]\) = ",MatrixForm[yD]];
	AdR = QuarkAMatrices[yU, yD, yUMinor, yDMinor][[4]];
	Print["The resulted \!\(\*SubscriptBox[SuperscriptBox[\(A\), \(d\)], \(R\)]\) = ", MatrixForm[AdR]];
	Print[ " with norm ", MatrixForm[AdR.ConjugateTranspose[AdR]]];
	Print["QuarkAMatricesUnitaryQ[] returns ", QuarkAMatricesUnitaryQ[yU, yD, yUMinor, yDMinor]];
]


(****************************************************************
 * Fermion profiles and their overlap
 ****************************************************************)


FermionProfile[c_?NumericQ, z_, zir_:10^8]:=
	Piecewise[{{Sqrt[(1-2c)/(zir^(1-2c) - 1)]z^(2-c),c>1/2||c<1/2},{1/Sqrt[Log[zir]],c==1/2}}]


FermionProfileUVOverlap[cL_?NumericQ, cR_?NumericQ, zir_:10^8]:=Module[{zuv=1.1}, 
	FermionProfile[cL, zuv, zir] FermionProfile[cR, zuv, zir]
]


FermionProfileBulkOverlap[cL_?NumericQ,cR_?NumericQ,zir_:10^8]:=
	Piecewise[{{Sqrt[(1-2cL)/(zir^(1-2cL) - 1)],cL!=1/2},{Sqrt[1/Log[zir]],cL==1/2}}]*
	Piecewise[{{Sqrt[(1-2cR)/(zir^(1-2cR) - 1)],cR!=1/2},{Sqrt[1/Log[zir]],cR==1/2}}]* 
	Piecewise[{{(1-zir^-(cL+cR))/(cL+cR),cL+cR!=0},{Log[zir],cL+cR==0}}
];


TurnRight[listLR_]:=Module[{cL=Transpose[listLR][[1]], cR=Transpose[listLR][[2]]}, 
	Transpose[{cL-cR, cL+cR}]
];


TurnLeft[listPM_]:=Module[{cM=Transpose[listPM][[1]], cP=Transpose[listPM][[2]]},
	Transpose[{cP+cM, cP-cM}]/2
];


(* Axion-fermion-fermion overlap integral (non-flat part) 
	c either cL or -cR
*)


FermionAxionOverlap[c_,\[CapitalDelta]_,zir_:10^8]:=(1-2c)/(zir^(1-2c) - 1) 1/(4\[CapitalDelta](\[CapitalDelta]-1)) ((zir^(1-2 c)-zir^(-2 \[CapitalDelta]))/(1-2 c+2 \[CapitalDelta])+\[CapitalDelta] (1/zir^2-zir^(1-2 c))/(3-2 c));


(* fermion overlap functions in Subscript[c, P], Subscript[c, M] *)


FermionProfileUVOverlapCM[cP_?NumericQ, cM_?NumericQ, zir_:10^8]:= Module[{zuv=1.1}, 
	FermionProfile[(cP-cM)/2, zuv, zir] FermionProfile[(cP+cM)/2, zuv, zir]
];


FermionProfileBulkOverlapCM[cP_?NumericQ, cM_?NumericQ, zir_:10^8]:=
	Piecewise[{{Sqrt[(1-(cP+cM))/(zir^(1-(cP+cM)) - 1)],cP+cM!=1},{Sqrt[1/Log[zir]],cP+cM==1}}]*
	Piecewise[{{Sqrt[(1-(cP-cM))/(zir^(1-(cP-cM)) - 1)],cP-cM!=1},{Sqrt[1/Log[zir]],cP-cM==1}}]* 
	Piecewise[{{(1-zir^-cP)/cP,cP!=0},{Log[zir],cP==0}}
];


FixedRange[min_, max_, length_]:=Module[{step=(max-min)/(length-1)},Range[min, max+1/2 step, step]];


ListMirror[list_]:= Reverse[Reverse/@list]~Join~list;


(****************************************************************
 * Graphing functions
 ****************************************************************)


(****************************************************************
 * End private context
 ****************************************************************)


End[];


(****************************************************************
 * End package
 ****************************************************************)


EndPackage[];
