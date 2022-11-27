(* ::Package:: *)

(* Wolfram Language Package *)

(* Created by the Wolfram Workbench Jun 18, 2021 *)

BeginPackage["QTAIM`"]
(* Exported symbols added here with SymbolName::usage *) 

LocateCriticalPoints::"Locate Critical Points in a 3D Field.";
ReadWavefunctionFromWFX::"Read a Wavefunction from a WFX file.";
(*ReadWavefunctionFromMOLDEN::"";*)
$QTAIMContours::"Default Contours in QTAIM.";
$QTAIMContoursPositive::"Default Contours in QTAIM (Positive Valued).";
ElectronDensityDerivativeInternal::"Compute Spatial Derivative of Electron Density. (Native Routines)";
ElectronDensityDerivativeCompiled::"Compute Spatial Derivative of Electron Density. (Compiled)";
ElectronDensityCompiled::"Return Electron Density (rho).";
ElectronDensityGradientCompiled::"Return Electron Density Gradient.";
KineticEnergyDensityGCompiled::"Return Kinetic Energy Density G.";
ElectronDensityGradientXCompiled::"Return Electron Density Gradient (x component).";
ElectronDensityGradientYCompiled::"Return Electron Density Gradient (y component).";
ElectronDensityGradientZCompiled::"Return Electron Density Gradient (z component).";
ReadWavefunctionFromPySCF::"Read Wavefunvtion from the result of running a PySCF input.";
LocateNuclearCriticalPoints::"Locate nuclear critical points.";
AssociatedNuclearAttractor::"Index of nuclear attractor in supplied nuclear critical point list that is associated with point via path of gradient ascent.";;
SteepestAscentPath::"Return the path of steepest ascent through the electron density."
SteepestDescentPath::"Return the path of steepest descent through the electron density."
BondPath::"Return the paths of steepest ascent in the electron density forward and backward the lowest Hessian eigenvector."
ElectronDensity::"Return (potentially symbolic) electron density"
ElectronDensityDerivative::"Electron Density Derivative Internal Function."
InteratomicSurface::"Locate an interatomic surface radially from a nuclear critical point."
CriticalPointGradientPlus3::""
CriticalPointGradientPlus1::""
CriticalPointGradientMinus1::""
CriticalPointGradientMinus3::""
ParallelNInt::""
rho::""
g::""
gx::""
gy::""
gz::""
H::""
KineticEnergyDensityG::""
KineticEnergyDensityK::""
StressTensorTrace:""
BS32::""
BS54::""
DOPRI54amat::""
DOPRI54bvec::""
DOPRI54cvec::""
DOPRI54evec::""
DOPRI54Coefficients::""
CashKarp45amat::""
CashKarp45bvec::""
CashKarp45bhatvec::""
CashKarp45cvec::""
CashKarp45evec::""
CashKarp45Coefficients::""
CashKarp45::"Cash-Karp"
DOPRIamat::""
DOPRI54bvec::""
DOPRI54cvec::""
DOPRI54evec::""
DOPRI54Coefficients::""
DOPRI54::"Doermand-Prince"
Fehlberg45amat::""
Fehlberg45bvec::""
Fehlberg45cvec::""
Fehlberg45evec::""
Fehlberg45Coefficients::""
Fehlberg45::"Fehlberg"

Begin["`Private`"]
(* Implementation of the package *)

$QTAIMContours = Sort[Join[{0}, Flatten[Table[{-8,-4,-2,2,4,8}*10^n, {n, -3,3}]]]];
$QTAIMContoursPositive = Sort[Flatten[Table[{2,4,8}*10^n, {n, -3,8}]]];

(* ODE integrators and Butcher Tables *)
AdamsBMCoefficients[hlist_List] := 
 Module[{k, h, \[CapitalDelta]h, brat, \[Beta], \[Alpha], \[Sigma], c},
  k = Length[hlist];
  h = Last[hlist];
  \[CapitalDelta]h = Drop[FoldList[Plus, 0, Reverse[hlist]], 1];
  brat = Drop[\[CapitalDelta]h, -1]/(Drop[\[CapitalDelta]h, 1] - h);
  \[Beta] = FoldList[Times, 1, brat];
  \[Alpha] = h/\[CapitalDelta]h;
  \[Sigma] = FoldList[Times, 1, \[Alpha] Range[Length[\[Alpha]]]];
  c[0] = Table[1/q, {q, 1, k}];
  c[1] = Table[1/(q (q + 1)), {q, 1, k}]; 
  Do[c[j] = 
    Drop[c[j - 1], -1] - (
     Drop[c[j - 1], 1] h)/\[CapitalDelta]h[[j]], {j, 2, 
    k}]; {(First[c[#1]] &) /@ Range[0, k], \[Beta], \[Sigma]}];
Moulton[0] = 1;
Moulton[m_] := Moulton[m] = -Sum[Moulton[k]/(1 + m - k), {k, 0, m - 1}];
AdamsBM /: NDSolve`InitializeMethod[AdamsBM, {Automatic, __}, sd_, rhs_, 
   ndstate_, opts___] := 
  Module[{prec, norm, hlist, \[CapitalPhi], mord},
   mord = MaxDifferenceOrder /. Flatten[{opts, Options[AdamsBM]}];
   If[mord != \[Infinity] && ! (IntegerQ[mord] && mord > 0), 
    Return[$Failed]];
   prec = ndstate["WorkingPrecision"];
   norm = ndstate["Norm"];
   hlist = {};
   \[CapitalPhi] = {ndstate["SolutionDerivativeVector"["Active"]]}; 
   AdamsBM[{{hlist, \[CapitalPhi], 
      N[0, prec] \[CapitalPhi][[1]]}, {norm, prec, mord, 0, True}}]];
Options[AdamsBM] = {MaxDifferenceOrder -> \[Infinity]};
AdamsBM[___]["StepMode"] = Automatic;
AdamsBM[data_]["DifferenceOrder"] := Length[data[[1, 2]]];      
AdamsBM[data_]["Step"[rhs_, t_, h_, y_, yp_]] := 
  Module[{prec, norm, hlist, \[CapitalPhi], \[CapitalPhi]1, ns, 
    starting, k, zero, g, \[Beta], \[Sigma], p, f, \[CapitalDelta]y, 
    normh, ev, err, PE, knew, hnew, 
    temp}, {{hlist, \[CapitalPhi], \[CapitalPhi]1}, {norm, prec, mord,
       ns, starting}} = data;
   (* Norm scaling will be based on current solution y. *)
   normh = (Abs[h] temp[#1, y] &) /. {temp -> norm};
   k = Length[\[CapitalPhi]];
   zero = N[0, prec];
   (* Keep track of number of steps at this stepsize h. *)
   
   If[Length[hlist] > 0 && Last[hlist] == h, ns++, ns = 1]; 
   hlist = Join[hlist, {h}];
   {g, \[Beta], \[Sigma]} = AdamsBMCoefficients[hlist];
   (* Convert \[CapitalPhi] to \[CapitalPhi]^* *)
   \[CapitalPhi] = \[CapitalPhi] Reverse[\[Beta]];
   (* PE: Predict and evaluate *)
   
   p = Reverse[Drop[g, -1]] . \[CapitalPhi];
   f = rhs[h + t, h p + y];
   (* Update divided differences *)
   \[CapitalPhi] = 
    FoldList[Plus, zero \[CapitalPhi]1, \[CapitalPhi]];
   (* Compute scaled error estimate *)
   ev = f - Last[\[CapitalPhi]];
   err = (g[[-2]] - g[[-1]]) normh[ev];
   (* First order check: determines if order should be lowered
   even in the case of a rejected step *)
   
   knew = OrderCheck[PE, k, \[CapitalPhi], ev, normh, \[Sigma]];
   If[err > 1,
    (* Rejected step: reduce h by half, 
    make sure starting mode flag is unset and reset \[CapitalPhi] to previous values *)
    hnew = h/2; \[CapitalDelta]y = $Failed; 
    f = None; starting = False; \[CapitalPhi] = data[[1, 2]],
    (* Sucessful step:
    CE: Correct and evaluate *)
    \[CapitalDelta]y = h (p + ev Last[g]);
    f = rhs[h + t, y + \[CapitalDelta]y]; 
    temp = f - Last[\[CapitalPhi]];
    (* Update the divided differences *)
    \[CapitalPhi] = (temp + #1 &) /@ \[CapitalPhi];
    (* Determine best order and stepsize for the next step *)
\[CapitalPhi]1 = temp - \[CapitalPhi]1; 
    knew = 
     ChooseNextOrder[starting, PE, k, knew, \[CapitalPhi]1, 
      normh, \[Sigma], mord, ns]; 
    hnew = ChooseNextStep[PE, knew, h]];
   (* Truncate hlist and \[CapitalPhi] to the appropriate length for 
the chosen order. *)
   hlist = Take[hlist, 1 - knew];
   If[Length[\[CapitalPhi]] > 
     knew, \[CapitalPhi]1 = \[CapitalPhi][[
      Length[\[CapitalPhi]] - knew]]; \[CapitalPhi] = 
     Take[\[CapitalPhi], -knew];];
   (* Return step data along with updated method data *)
   {hnew, \[CapitalDelta]y, f, 
    AdamsBM[{{hlist, \[CapitalPhi], \[CapitalPhi]1}, {norm, prec, 
       mord, ns, starting}}]}];
OrderCheck[PE_, k_, \[CapitalPhi]_, ev_, normh_, \[Sigma]_] := 
  Module[{knew = k},
   Subscript[PE, k] = Abs[\[Sigma][[k + 1]] Moulton[k] normh[ev]]; 
   If[k > 1,
    Subscript[PE, k - 1] = 
     Abs[\[Sigma][[k]] Moulton[k - 1] normh[
        ev + \[CapitalPhi][[2]]]];
    If[k > 2,
     Subscript[PE, k - 2] = 
      Abs[\[Sigma][[k - 1]] Moulton[k - 2] normh[
         ev + \[CapitalPhi][[3]]]]; 
     If[Max[Subscript[PE, k - 1], Subscript[PE, k - 2]] < Subscript[
       PE, k], knew = k - 1]],
    If[Subscript[PE, k - 1] < Subscript[PE, k]/2, knew = k - 1];
    ];
   knew
   ];
SetAttributes[ChooseNextOrder, HoldFirst];
ChooseNextOrder[starting_, PE_, k_, knw_, \[CapitalPhi]1_, 
   normh_, \[Sigma]_, mord_, ns_] := Module[{knew = knw},
   starting = starting && knew >= k && k < mord;
   If[starting,
    knew = k + 1; Subscript[PE, k + 1] = 0,
    If[knew >= k && ns >= k + 1,
      Subscript[PE, k + 1] = Abs[Moulton[k + 1] normh[\[CapitalPhi]1]];
      If[k > 1,
       If[
        Subscript[PE, k - 1] <= 
         Min[Subscript[PE, k], Subscript[PE, k + 1]],
        knew = k - 1,
        If[Subscript[PE, k + 1] < Subscript[PE, k] && k < mord, 
         knew = k + 1]
        ],
       If[Subscript[PE, k + 1] < Subscript[PE, k]/2, knew = k + 1]
       ];
      ];
    ];
   knew
   ];
ChooseNextStep[PE_, k_, h_] :=
  If[Subscript[PE, k] < 2^-(k + 2),
   2 h,
   If[Subscript[PE, k] < 1/2, h, 
    h Max[1/2, Min[9/10, (1/(2 Subscript[PE, k]))^(1/(k + 1))]]]
   ];

(* Bogacki-Shampine 5(4) FSAL *) 
BS54 = {"ExplicitRungeKutta", 
 "Coefficients" -> "EmbeddedExplicitRungeKuttaCoefficients", 
 "DifferenceOrder" -> 5, "StiffnessTest" -> False}

(* Fehlberg *)
Fehlberg45amat = {
       {1/4},
       {3/32, 9/32},
       {1932/2197, -7200/2197, 7296/2197},
       {439/216, -8, 3680/513, -845/4104},
       {-8/27, 2, -3544/2565, 1859/4104, -11/40}};
Fehlberg45bvec = {25/216, 0, 1408/2565, 2197/4104, -1/5, 0};
Fehlberg45cvec = {1/4, 3/8, 12/13, 1, 1/2};
Fehlberg45evec = {-1/360, 0, 128/4275, 2197/75240, -1/50, -2/55};
Fehlberg45Coefficients[4, p_] := N[{Fehlberg45amat, Fehlberg45bvec, Fehlberg45cvec, Fehlberg45evec}, p];
Fehlberg45 = {"ExplicitRungeKutta", "Coefficients" -> Fehlberg45Coefficients, "DifferenceOrder" -> 4, "EmbeddedDifferenceOrder" -> 5, "StiffnessTest" -> False};

(* Cash-Karp *)
CashKarp45amat = {
      {1/5},
      {3/40, 9/40},
      {3/10, -9/10, 6/5},
      {-11/54, 5/2, -70/27, 35/27},
      {1631/55296, 175/512, 575/13824, 44275/110592, 253/4096}};
CashKarp45bvec = {2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 1/4};
CashKarp45bhatvec = {37/378, 0, 250/621, 125/594, 0, 512/1771};
CashKarp45cvec = {1/5, 3/10, 3/5, 1, 7/8};
CashKarp45evec = CashKarp45bvec - CashKarp45bhatvec;
CashKarp45Coefficients[4, p_] := N[{CashKarp45amat, CashKarp45bvec, CashKarp45cvec, CashKarp45evec}, p];
CashKarp45 = {"ExplicitRungeKutta", "Coefficients" -> CashKarp45Coefficients, "DifferenceOrder" -> 4, "EmbeddedDifferenceOrder" -> 5, "StiffnessTest" -> False};

(* Doermand-Prince *)
DOPRI54amat = {
      {1/5},
      {3/40, 9/40},
      {44/45, -56/15, 32/9},
      {19372/6561, -25360/2187, 64448/6561, -212/729},
      {9017/3168, -355/33, 46732/5247, 49/176, -5103/18656},
      {35/384, 0, 500/1113, 125/192, -2187/6784, 11/84}};
DOPRI54bvec = {35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0};
DOPRI54cvec = {1/5, 3/10, 4/5, 8/9, 1, 1};
DOPRI54evec = {71/57600, 0, -71/16695, 71/1920, -17253/339200, 22/525, -1/40};
DOPRI54Coefficients[5, p_] := N[{DOPRI54amat, DOPRI54bvec, DOPRI54cvec, DOPRI54evec}, p];
DOPRI54 = {"ExplicitRungeKutta", "Coefficients" -> DOPRI54Coefficients, "DifferenceOrder" -> 5, "StiffnessTest" -> False};

(* From: https://mathematica.stackexchange.com/questions/11364/how-do-i-find-all-the-solutions-of-three-simultaneous-equations-within-a-given-b *)
Options[LocateCriticalPoints] = 
  Sort[Join[Options[ContourPlot3D], Options[FindRoot], {MaxRecursion -> Automatic, 
       PerformanceGoal :> $PerformanceGoal, PlotPoints -> Automatic}]];

LocateCriticalPoints[funcs_?VectorQ,
                   {x_, xmin_, xmax_}, {y_, ymin_, ymax_}, {z_, zmin_, zmax_}, opts___] := 
 Module[{contourData, seeds, fz = Quiet[Compile[{x,y,z}, Evaluate[funcs[[3]]]]]}, 
  contourData = Quiet[Cases[Normal[ContourPlot3D[
      Evaluate[Most[funcs]], {x, xmin, xmax}, {y, ymin, ymax}, {z, zmin, zmax}, 
      BoundaryStyle -> {1 -> None, 2 -> None, {1, 2} -> {}}, 
      ContourStyle -> None, Mesh -> None, Method -> Automatic, 
      Evaluate[Sequence @@ FilterRules[Join[{opts}, Options[LocateCriticalPoints]], 
         Options[ContourPlot3D]]]]], Line[l_] :> l, Infinity]];
  seeds = Flatten[Pick[Rest[#], 
                  Most[#] Rest[#] &@Sign[Apply[fz, #, 2]], -1] & /@ contourData, 1];
  If[seeds === {}, seeds, 
   Select[Union[Map[{x, y, z} /. FindRoot[funcs, Transpose[{{x, y, z}, #}], 
         Evaluate[Sequence @@ FilterRules[Join[{opts}, Options[LocateCriticalPoints]], 
            Options[FindRoot]]] ] &, seeds]],
          (xmin < #[[1]] < xmax && ymin < #[[2]] < ymax && zmin < #[[3]] < zmax) &]]];

Options[ReadWavefunctionFromWFX] =  {UnusedOption -> Automatic};

(* Reading Wavefunctions *)

ReadWavefunctionFromWFX[filename_?StringQ] := Module[
  {wfx, lx, ly, lz, l,
   numberOfNuclei, numberOfOccupiedMolecularOrbitals, 
   numberOfPrimitives,
   atomicNumbers, nuclearCartesianCoordinates, 
   primitiveTypes, primitiveCenters, primitiveExponents,
   molecularOrbitalOccupationNumbers, 
   molecularOrbitalPrimitiveCoefficients,
   energy, virialRatio},
  wfx = Import[filename, "Text"];
  wfx = StringReplace[wfx, 
    "<" ~~ Shortest[x___] ~~ ">" :> 
     "<" <> StringReplace[x, {" " -> ""}] <> ">"];
  wfx = StringReplace[wfx,
    {
     "<Energy=T+Vne+Vee+Vnn>" -> "<Energy>",
     "<VirialRatio(-V/T)>" -> "<VirialRatio>",
     "</Energy=T+Vne+Vee+Vnn>" -> "</Energy>",
     "</VirialRatio(-V/T)>" -> "</VirialRatio>"
     }
    ];
  wfx = StringReplace[wfx,
    {
     "Numberof" -> "NumberOf"}
    ];
  wfx = "<wfx>
" <>
    wfx <> "
</wfx>
";
  wfx = ImportString[
    wfx
    , "XML"
    ];
  
  numberOfNuclei = 
   Cases[wfx, XMLElement["NumberOfNuclei", _, _], Infinity];
  numberOfNuclei = ToExpression[
    numberOfNuclei // First // Last // First
    ];
  
  numberOfOccupiedMolecularOrbitals = 
   Cases[wfx, XMLElement["NumberOfOccupiedMolecularOrbitals", _, _], 
    Infinity];
  numberOfOccupiedMolecularOrbitals = ToExpression[
    numberOfOccupiedMolecularOrbitals // First // Last // First
    ];
  
  atomicNumbers = 
   Cases[wfx, XMLElement["AtomicNumbers", _, _], Infinity];
  atomicNumbers = Map[ToExpression, StringSplit[
     atomicNumbers // First // Last // First
     ]
    ];
  
  nuclearCartesianCoordinates = 
   Cases[wfx, XMLElement["NuclearCartesianCoordinates", _, _], 
    Infinity];
  nuclearCartesianCoordinates = Partition[
    Map[
     ToExpression,
     StringSplit[
      StringReplace[
       nuclearCartesianCoordinates // First // Last // First,
       {
        "e" -> "*^",
        "E" -> "*^"
        }
       ]
      ]
     ]
    , 3];
  
  numberOfPrimitives = 
   Cases[wfx, XMLElement["NumberOfPrimitives", _, _], Infinity];
  numberOfPrimitives = ToExpression[
    numberOfPrimitives // First // Last // First
    ];
  
  primitiveTypes = 
   Cases[wfx, XMLElement["PrimitiveTypes", _, _], Infinity];
  primitiveTypes = Map[ToExpression, StringSplit[
     primitiveTypes // First // Last // First
     ]
    ];
  
  lx = Table[0, {j, 1, numberOfPrimitives}];
  ly = Table[0, {j, 1, numberOfPrimitives}];
  lz = Table[0, {j, 1, numberOfPrimitives}];
  l = Table[
    Which[ 
     primitiveTypes[[j]] == 1, {0, 0, 0},
     primitiveTypes[[j]] == 2, {1, 0, 0},
     primitiveTypes[[j]] == 3, {0, 1, 0},
     primitiveTypes[[j]] == 4, {0, 0, 1},
     primitiveTypes[[j]] == 5, {2, 0, 0},
     primitiveTypes[[j]] == 6, {0, 2, 0},
     primitiveTypes[[j]] == 7, {0, 0, 2},
     primitiveTypes[[j]] == 8, {1, 1, 0},
     primitiveTypes[[j]] == 9, {1, 0, 1},
     primitiveTypes[[j]] == 10, {0, 1, 1},
     primitiveTypes[[j]] == 11, {3, 0, 0},
     primitiveTypes[[j]] == 12, {0, 3, 0},
     primitiveTypes[[j]] == 13, {0, 0, 3},
     primitiveTypes[[j]] == 14, {2, 1, 0},
     primitiveTypes[[j]] == 15, {2, 0, 1},
     primitiveTypes[[j]] == 16, {0, 2, 1},
     primitiveTypes[[j]] == 17, {1, 2, 0},
     primitiveTypes[[j]] == 18, {1, 0, 2},
     primitiveTypes[[j]] == 19, {0, 1, 2},
     primitiveTypes[[j]] == 20, {1, 1, 1},
     primitiveTypes[[j]] == 21, {4, 0, 0},
     primitiveTypes[[j]] == 22, {0, 4, 0},
     primitiveTypes[[j]] == 23, {0, 0, 4},
     primitiveTypes[[j]] == 24, {3, 1, 0},
     primitiveTypes[[j]] == 25, {3, 0, 1},
     primitiveTypes[[j]] == 26, {1, 3, 0},
     primitiveTypes[[j]] == 27, {0, 3, 1},
     primitiveTypes[[j]] == 28, {1, 0, 3},
     primitiveTypes[[j]] == 29, {0, 1, 3},
     primitiveTypes[[j]] == 30, {2, 2, 0},
     primitiveTypes[[j]] == 31, {2, 0, 2},
     primitiveTypes[[j]] == 32, {0, 2, 2},
     primitiveTypes[[j]] == 33, {2, 1, 1},
     primitiveTypes[[j]] == 34, {1, 2, 1},
     primitiveTypes[[j]] == 35, {1, 1, 2},
     (* H and higher:
            # Do IX=0,L
            # Do IY=0,(L-IX)
            # IZ=L-IX-IY
     *)
     primitiveTypes[[j]] == 36, {0, 0, 5},
     primitiveTypes[[j]] == 37, {0, 1, 4},
     primitiveTypes[[j]] == 38, {0, 2, 3},
     primitiveTypes[[j]] == 39 , {0, 3, 2},
     primitiveTypes[[j]] == 40, {0, 4, 1},
     primitiveTypes[[j]] == 41, {0, 5, 0},
     primitiveTypes[[j]] == 42, {1, 0, 4},
     primitiveTypes[[j]] == 43, {1, 1, 3},
     primitiveTypes[[j]] == 44, {1, 2, 2},
     primitiveTypes[[j]] == 45, {1, 3, 1},
     primitiveTypes[[j]] == 46, {1, 4, 0},
     primitiveTypes[[j]] == 47, {2, 0, 3},
     primitiveTypes[[j]] == 48, {2, 1, 2},
     primitiveTypes[[j]] == 49, {2, 2, 1},
     primitiveTypes[[j]] == 50, {2, 3, 0},
     primitiveTypes[[j]] == 51, {3, 0, 2},
     primitiveTypes[[j]] == 52, {3, 1, 1},
     primitiveTypes[[j]] == 53, {3, 2, 0},
     primitiveTypes[[j]] == 54, {4, 0, 1},
     primitiveTypes[[j]] == 55, {4, 1, 0},
     primitiveTypes[[j]] == 56, {5, 0, 0}
     ]
    , {j, 1, numberOfPrimitives}];
  lx = l[[All, 1]];
  ly = l[[All, 2]];
  lz = l[[All, 3]];
  
  primitiveCenters = 
   Cases[wfx, XMLElement["PrimitiveCenters", _, _], Infinity];
  primitiveCenters = Map[ToExpression, StringSplit[
     primitiveCenters // First // Last // First
     ]
    ];
  
  primitiveExponents = 
   Cases[wfx, XMLElement["PrimitiveExponents", _, _], Infinity];
  primitiveExponents = Map[ToExpression,
    StringSplit[
     StringReplace[
      primitiveExponents // First // Last // First,
      {
       "e" -> "*^",
       "E" -> "*^"
       }
      ]
     ]
    ];
  
  molecularOrbitalOccupationNumbers = 
   Cases[wfx, XMLElement["MolecularOrbitalOccupationNumbers", _, _], 
    Infinity];
  molecularOrbitalOccupationNumbers = Map[ToExpression,
    StringSplit[
     StringReplace[
      molecularOrbitalOccupationNumbers // First // Last // First,
      {
       "e" -> "*^",
       "E" -> "*^"
       }
      ]
     ]
    ];
  
  molecularOrbitalPrimitiveCoefficients = 
   Cases[wfx, 
    XMLElement["MolecularOrbitalPrimitiveCoefficients", _, _], 
    Infinity];
  molecularOrbitalPrimitiveCoefficients = 
   Partition[molecularOrbitalPrimitiveCoefficients[[1, 3]], 2][[All, 
     2]];
  molecularOrbitalPrimitiveCoefficients = Map[
    Map[
      ToExpression,
      StringSplit[
       StringReplace[
        #,
        {
         "e" -> "*^",
         "E" -> "*^"
         }
        ]
       ]
      ] &
    , molecularOrbitalPrimitiveCoefficients];
  
  energy = Cases[wfx, XMLElement["Energy", _, _], Infinity];
  energy = ToExpression[StringReplace[
     energy // First // Last // First,
     {
      "e" -> "*^",
      "E" -> "*^"
      }
     ]
    ];
  
  virialRatio = 
   Cases[wfx, XMLElement["VirialRatio", _, _], Infinity];
  virialRatio = ToExpression[StringReplace[
     virialRatio // First // Last // First,
     {
      "e" -> "*^",
      "E" -> "*^"
      }
     ]
    ];
  
  Association[
   "NumberOfNuclei" -> numberOfNuclei,
   "NumberOfOccupiedMolecularOrbitals" -> 
    numberOfOccupiedMolecularOrbitals,
   "AtomicNumbers" -> atomicNumbers,
   "NuclearCartesianCoordinates" -> nuclearCartesianCoordinates,
   "NumberOfPrimitives" -> numberOfPrimitives,
   "PrimitiveTypes" -> primitiveTypes,
   "PrimitiveCenters" -> primitiveCenters,
   "PrimitiveExponents" -> primitiveExponents,
   "xp" -> Table[
     (nuclearCartesianCoordinates[[primitiveCenters[[j]]]])[[1]], {j, 
      1, numberOfPrimitives}],
    "yp" -> Table[
     (nuclearCartesianCoordinates[[primitiveCenters[[j]]]])[[2]], {j, 
      1, numberOfPrimitives}],
   "zp" -> Table[
     (nuclearCartesianCoordinates[[primitiveCenters[[j]]]])[[3]], {j, 
      1, numberOfPrimitives}],
   "l" -> l,
   "lx" -> lx,
   "ly" -> ly,
   "lz" -> lz,
   "MolecularOrbitalOccupationNumbers" -> 
    molecularOrbitalOccupationNumbers,
   "MolecularOrbitalPrimitiveCoefficients" -> 
    molecularOrbitalPrimitiveCoefficients,
    "cFlatten" -> Flatten[molecularOrbitalPrimitiveCoefficients],
    "cTranspose" -> Transpose[molecularOrbitalPrimitiveCoefficients],
    "cFlattenTranspose" -> Flatten[Transpose[molecularOrbitalPrimitiveCoefficients]],
   "Energy" -> energy,
   "VirialRatio" -> virialRatio
   ]
  
  
  ];

ElectronDensity[wfn_,{x_,y_,z_}] := Module[
{
nmo=wfn["NumberOfOccupiedMolecularOrbitals"],
np=wfn["NumberOfPrimitives"],
xp=wfn["xp"],
yp=wfn["yp"],
zp=wfn["zp"],
lx=wfn["lx"],
ly=wfn["ly"],
lz=wfn["lz"],
a=wfn["PrimitiveExponents"],
c=wfn["MolecularOrbitalPrimitiveCoefficients"],
o=wfn["MolecularOrbitalOccupationNumbers"]
},

   xxp = x-xp;
   yyp = y-yp;
   zzp = z-zp;
   expa = -a*(xxp^2 + yyp^2 + zzp^2);

   p = Table[
   	If[lx[[j]]==0, 1, xxp[[j]]^lx[[j]]] * 
   	If[ly[[j]]==0, 1, yyp[[j]]^ly[[j]]] * 
   	If[lz[[j]]==0, 1, zzp[[j]]^lz[[j]]] * 
   	Exp[-a[[j]] * (xxp[[j]]^2 + yyp[[j]]^2 + zzp[[j]]^2)]
   	,{j, 1, np}];
   cp = c . p;
   Total[Table[o[[i]]*cp[[i]]^2, {i,1,nmo}]]
];

(*
ElectronDensityDerivativeInternal[wfn_, 
	                      {x_?NumericQ, y_?NumericQ, z_?NumericQ}, 
                          {nx_?IntegerQ, ny_?IntegerQ, nz_?IntegerQ}] := Module[
  {nmo=wfn["NumberOfOccupiedMolecularOrbitals"],
   nprim=wfn["NumberOfPrimitives"],
   l=wfn["l"],
   xp=wfn["xp"],
   yp=wfn["yp"],
   zp=wfn["zp"],   
   a=wfn["PrimitiveExponents"],
   o=wfn["MolecularOrbitalOccupationNumbers"],
   c=wfn["MolecularOrbitalPrimitiveCoefficients"],
   expa, pos, gamma, p,cp, cpos, apos, xpos, ypos, zpos, xpow, ypow, zpow, 
   sqrta, binomialx, binomialy, binomialz, exp, angx, expx, angy, expy, angz, 
   expz, xxp, yyp, zzp, pa, pb},

   xxp = x-xp;
   yyp = y-yp;
   zzp = z-zp;

   expa = -N@a * (xxp^2 + yyp^2 + zzp^2);
   
   pos = Flatten[Position[expa, _?(# > Log[1.*^-12] &)]];
     
   If[ Length[pos] == 0, Return[0];];
   
   xxp = N@xxp[[pos]];
   yyp = N@yyp[[pos]];
   zzp = N@zzp[[pos]];
   l = l[[pos]];
   expa = N@expa[[pos]];
   cpos = N@c[[All,pos]];
   apos = N@a[[pos]];
   

   (*
   If[ nx == 0 && ny == 0 && nz == 0,
  	
     p = Table[
      myPower[xxp[[j]], l[[j, 1]]]*
      myPower[yyp[[j]], l[[j, 2]]]*
      myPower[zzp[[j]], l[[j, 3]]]*
      Exp[-apos[[j]]*(xxp[[j]]^2 + yyp[[j]]^2 + zzp[[j]]^2)], {j, 1, Length[pos]}];
      cp = cpos . p;
      Return[ Sum[o[[i]]*cp[[i]]^2, {i, 1, nmo}] ];
    ];
    *)
    
sqrta = N@Sqrt[a[[pos]]];
binomialx = N@Table[Table[Binomial[n, k], {k, 0, nx}], {n, 0, nx}];
binomialy = N@Table[Table[Binomial[n, k], {k, 0, ny}], {n, 0, ny}];
binomialz = N@Table[Table[Binomial[n, k], {k, 0, nz}], {n, 0, nz}];
gamma = N@Table[Gamma[n], {n, 1, Max[nx,ny,nz]}];

angx = Table[
  Table[If[l[[j, 1]] < n, 0., If[(l[[j, 1]] - n)==0,1.,xxp[[j]]^(l[[j, 1]] - n)]*
     gamma[[l[[j, 1]] + 1]]/gamma[[l[[j, 1]] - n + 1]]], {j, 1, Length[pos]}], {n, 0, nx}];
angy = Table[
  Table[If[l[[j, 2]] < n, 0., If[(l[[j, 2]] - n)==0,1.,yyp[[j]]^(l[[j, 2]] - n)]*
     gamma[[l[[j, 2]] + 1]]/gamma[[l[[j, 2]] - n + 1]]], {j, 1, Length[pos]}], {n, 0, ny}];
angz = Table[
  Table[If[l[[j, 3]] < n, 0., If[(l[[j, 3]] - n)==0,1.,zzp[[j]]^(l[[j, 3]] - n)]*
     gamma[[l[[j, 3]] + 1]]/gamma[[l[[j, 3]] - n + 1]]], {j, 1, Length[pos]}], {n, 0, nz}];

expx = Table[Table[(-1.)^n*sqrta[[j]]^n*HermiteH[n, sqrta[[j]]*xxp[[j]]], {j, 1, Length[pos]}], {n, 0, nx}];
expy = Table[Table[(-1.)^n*sqrta[[j]]^n*HermiteH[n, sqrta[[j]]*yyp[[j]]], {j, 1, Length[pos]}], {n, 0, ny}];
expz = Table[Table[(-1.)^n*sqrta[[j]]^n*HermiteH[n, sqrta[[j]]*zzp[[j]]], {j, 1, Length[pos]}], {n, 0, nz}];
exp = Table[Exp[expa[[j]]], {j, 1, Length[pos]}];



pa = Table[
       Table[
         Table[

             Table[ 
            Total@Table[
             binomialx[[nx - kx + 1, KX + 1]]*
              angx[[(nx - kx) - KX + 1, j]]*
              expx[[KX + 1, j]]*
              Total@Table[
               binomialy[[ny - ky + 1, KY + 1]]*
                angy[[(ny - ky) - KY + 1, j]]*
                expy[[KY + 1, j]]*
                Total@Table[
                 binomialz[[nz - kz + 1, KZ + 1]]*
                  angz[[(nz - kz) - KZ + 1, j]]*
                  expz[[KZ + 1, j]]
                 , {KZ, 0, nz - kz}]
               , {KY, 0, ny - ky}]
             , {KX, 0, nx - kx}]*
            exp[[j]]
           , {j, 1, Length[pos]}]

        , {kz, 0, nz}]
      , {ky, 0, ny}]
    , {kx, 0, nx}];


pb = Table[
       Table[
         Table[
               	
          Table[
         	Total@Table[
             binomialx[[kx + 1, KX + 1]]*
              angx[[(kx - KX) + 1, j]]*
              expx[[KX + 1, j]]*
              Total@Table[
               binomialy[[ky + 1, KY + 1]]*
                angy[[(ky - KY) + 1, j]]*
                expy[[KY + 1, j]]*
                Total@Table[
                 binomialz[[kz + 1, KZ + 1]]*
                  angz[[(kz - KZ) + 1, j]]*
                  expz[[KZ + 1, j]]
                 , {KZ, 0, kz}]
               , {KY, 0, ky}]
             , {KX, 0, kx}]*
            exp[[j]]
           , {j, 1, Length[pos]}]

        , {kz, 0, nz}]
      , {ky, 0, ny}]
    , {kx, 0, nx}];
Return[0]; 

Sum[
 (
  o[[i]]*
   Sum[
    binomialx[[nx + 1, kx + 1]]*
     Sum[
      binomialy[[ny + 1, ky + 1]]*
       Sum[
        binomialz[[nz + 1, kz + 1]]*
         (

         (*Sum[ 
           cpos[[i, j]]*pa[[j]]
            ,{j, 1, Length[pos]}] *)
            cpos[[i]].pa[[kx+1,ky+1,kz+1]]
          )
         *
         (

          (*Sum[
           cpos[[i, j]]*pb[[j]]
           , {j, 1, Length[pos]}]*)
           cpos[[i]].pb[[kx+1,ky+1,kz+1]]
          )
        , {kz, 0, nz}]
      , {ky, 0, ny}]
    , {kx, 0, nx}]
  )
 , {i, 1, nmo}]
 
  ];
*)

(*
ElectronDensityDerivativeInternal[wfn_, 
	                      xyz_?MatrixQ, 
                          deriv_?MatrixQ] := Module[
  {nmo=wfn["NumberOfOccupiedMolecularOrbitals"],
   nprim=wfn["NumberOfPrimitives"],
   l=wfn["l"],
   xp=wfn["xp"],
   yp=wfn["yp"],
   zp=wfn["zp"],   
   a=wfn["PrimitiveExponents"],
   o=wfn["MolecularOrbitalOccupationNumbers"],
   c=wfn["MolecularOrbitalPrimitiveCoefficients"],
   myPower, expa, pos, p,cp, cpos, apos,
   sqrta, binomialx, binomialy, binomialz, exp, angx, expx, angy, expy, angz, 
   expz, xxp, yyp, zzp, pa, pb, nxmax, nymax, nzmax, nx, ny, nz, lx, ly, lz,
    gamma, constant, px, py, pz},

   nxmax=Max[deriv[[All,1]]];
   nymax=Max[deriv[[All,2]]];
   nzmax=Max[deriv[[All,3]]];

Table[

   xxp = xyz[[nn,1]]-xp;
   yyp = xyz[[nn,2]]-yp;
   zzp = xyz[[nn,3]]-zp;
   expa = -a*(xxp^2 + yyp^2 + zzp^2);

   pos = Flatten[Position[expa, _?(# > Log[1.*^-12] &)]];   
   
   (*If[ Length[pos] == 0, Return[ConstantArray[0, Length[deriv]]]; ];*)

   xxp = xxp[[pos]];
   yyp = yyp[[pos]];
   zzp = zzp[[pos]];
   lx = l[[pos,1]];
   ly = l[[pos,2]];
   lz = l[[pos,3]];
   expa = expa[[pos]];
   cpos = c[[All,pos]];
   apos = a[[pos]];
   sqrta = Sqrt[apos];

(*
   If[ nxmax == 0 && nymax == 0 && nzmax == 0,
  	
   p = Table[
      Quiet[Check[xxp[[j]]^lx[[j]], 1]]*
      Quiet[Check[yyp[[j]]^ly[[j]], 1]]*
      Quiet[Check[zzp[[j]]^lz[[j]], 1]]*
      Exp[-apos[[j]]*(xxp[[j]]^2 + yyp[[j]]^2 + zzp[[j]]^2)], {j, 1, Length[pos]}];
   cp = cpos . p;
   Return[ { Sum[o[[i]]*cp[[i]]^2, {i, 1, nmo}] } ];
   
  ];
*)



binomialx = Table[Binomial[n, k], {n, 0, nxmax}, {k, 0, nxmax}];
binomialy = Table[Binomial[n, k], {n, 0, nymax}, {k, 0, nymax}];
binomialz = Table[Binomial[n, k], {n, 0, nzmax}, {k, 0, nzmax}];
gamma = Table[Gamma[n], {n,1,Max[Join[lx,ly,lz]]+1}];

angx = Table[If[lx[[j]] < n, 0, 
    Quiet[Check[xxp[[j]]^(lx[[j]] - n), 1]] *
     gamma[[lx[[j]] + 1]]/gamma[[lx[[j]] - n + 1]]]
     , {n, 0, nxmax}, {j, 1, Length[pos]}];

angy = Table[If[ly[[j]] < n, 0, 
    Quiet[Check[yyp[[j]]^(ly[[j]] - n), 1]] *
     gamma[[ly[[j]] + 1]]/gamma[[ly[[j]] - n + 1]]]
     , {n, 0, nymax}, {j, 1, Length[pos]}];

angz = Table[If[lz[[j]] < n, 0, 
    Quiet[Check[zzp[[j]]^(lz[[j]] - n), 1]] *
     gamma[[lz[[j]] + 1]]/gamma[[lz[[j]] - n + 1]]]
     , {n, 0, nzmax}, {j, 1, Length[pos]}];

constant = Table[
	(-1)^n*sqrta[[j]]^n
     , {n, 0, Max[nxmax,nymax,nzmax]}, {j, 1, Length[pos]}];

expx = Table[
    constant[[n+1, j]] *
    HermiteH[n, sqrta[[j]]*xxp[[j]]], {n, 0, nxmax}, {j, 1, Length[pos]}];
expy = Table[
    constant[[n+1, j]] *
    HermiteH[n, sqrta[[j]]*yyp[[j]]], {n, 0, nymax}, {j, 1, Length[pos]}];
expz = Table[
    constant[[n+1, j]] *
    HermiteH[n, sqrta[[j]]*zzp[[j]]], {n, 0, nzmax}, {j, 1, Length[pos]}];
exp = Exp[expa];


Table[

nx=deriv[[d,1]];
ny=deriv[[d,2]];
nz=deriv[[d,3]];

(* These loop orders need reversed, so indexes aren't recomputed *)

pa = Table[
            Sum[
             binomialx[[nx - kx + 1, KX + 1]]*
              angx[[(nx - kx) - KX + 1, j]]*
              expx[[KX + 1, j]]*
              Sum[
               binomialy[[ny - ky + 1, KY + 1]]*
                angy[[(ny - ky) - KY + 1, j]]*
                expy[[KY + 1, j]]*
                Sum[
                 binomialz[[nz - kz + 1, KZ + 1]]*
                  angz[[(nz - kz) - KZ + 1, j]]*
                  expz[[KZ + 1, j]]
                 , {KZ, 0, nz - kz}]
               , {KY, 0, ny - ky}]
             , {KX, 0, nx - kx}]*
            exp[[j]]
            , {kx, 0, nx}
           , {ky, 0, ny}
        , {kz, 0, nz}
       , {j, 1, Length[pos]}
    ];


pb = Table[
         	Sum[
             binomialx[[kx + 1, KX + 1]]*
              angx[[(kx - KX) + 1, j]]*
              expx[[KX + 1, j]]*
              Sum[
               binomialy[[ky + 1, KY + 1]]*
                angy[[(ky - KY) + 1, j]]*
                expy[[KY + 1, j]]*
                Sum[
                 binomialz[[kz + 1, KZ + 1]]*
                  angz[[(kz - KZ) + 1, j]]*
                  expz[[KZ + 1, j]]
                 , {KZ, 0, kz}]
               , {KY, 0, ky}]
             , {KX, 0, kx}]*
            exp[[j]]
                , {kx, 0, nx}
             , {ky, 0, ny}
            , {kz, 0, nz}
         , {j, 1, Length[pos]}];

Sum[
 (
  o[[i]]*
   Sum[
    binomialx[[nx + 1, kx + 1]]*
     Sum[
      binomialy[[ny + 1, ky + 1]]*
       Sum[
        binomialz[[nz + 1, kz + 1]]*
         (

         (*Sum[ 
           cpos[[i, j]]*pa[[j]]
            ,{j, 1, Length[pos]}] *)
           cpos[[i]].pa[[kx+1,ky+1,kz+1]]
           (*cpa[[i, kx+1, ky+1, kz+1]]*)
          )
         *
         (

          (*Sum[
           cpos[[i, j]]*pb[[j]]
           , {j, 1, Length[pos]}]*)
           cpos[[i]].pb[[kx+1,ky+1,kz+1]]
           (*cpb[[i, kx+1, ky+1, kz+1]]*)
          )
        , {kz, 0, nz}]
      , {ky, 0, ny}]
    , {kx, 0, nx}]
  )
 , {i, 1, nmo}]

 , {d,1,Length[deriv]}]

, {nn, 1, Length[xyz]}]
];
*)

ElectronDensityCompiled=Compile[
	{
	{nmo, _Integer},
	{np, _Integer},
	{x, _Real},
	{y, _Real},
	{z, _Real},
	{xp, _Real, 1},
	{yp, _Real, 1},
	{zp, _Real, 1},
	{lx, _Integer, 1},
	{ly, _Integer, 1},
	{lz, _Integer, 1},
	{a, _Real, 1},
	{c, _Real, 2},
	{o, _Real, 1}		
	}, Module[
		{xxp, yyp, zzp, xxppos, yyppos, zzppos, expa, poscut, pos, poscount, lxpos, lypos, lzpos, expapos, cpos, apos, sqrtapos, p, cp, ctr},

   xxp = x-xp;
   yyp = y-yp;
   zzp = z-zp;
   expa = -a*(xxp^2 + yyp^2 + zzp^2);

(*
   pos = Flatten[Position[expa, _?(# > Log[1.*^-12] &)]];   
*)   

   poscut = Table[If[expa[[j]] > Log[1.*^-12],1,0], {j, 1, np}]; 
   poscount = Total[poscut];
   pos = Table[1, {j, 1, poscount}];
   ctr=1;
   Do[
     If[poscut[[j]]==1, pos[[ctr]] = j; ctr=ctr+1;];
     , {j, 1, np}] ;

   xxppos = xxp[[pos]];
   yyppos = yyp[[pos]];
   zzppos = zzp[[pos]];
   lxpos = lx[[pos]];
   lypos = ly[[pos]];
   lzpos = lz[[pos]];
   expapos = expa[[pos]];
   cpos = c[[All,pos]];
   apos = a[[pos]];

   p = Table[
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	Exp[-apos[[j]] * (xxppos[[j]]^2 + yyppos[[j]]^2 + zzppos[[j]]^2)]
   	,{j, 1, Length[pos]}]; 
   cp = cpos . p;
   Total[Table[o[[i]]*cp[[i]]^2, {i,1,nmo}]]

	]   
, CompilationTarget -> "C", CompilationOptions -> {"ExpressionOptimization" -> True}];

ElectronDensityGradientCompiled=Compile[
	{
	{nmo, _Integer},
	{np, _Integer},
	{x, _Real},
	{y, _Real},
	{z, _Real},
	{xp, _Real, 1},
	{yp, _Real, 1},
	{zp, _Real, 1},
	{lx, _Integer, 1},
	{ly, _Integer, 1},
	{lz, _Integer, 1},
	{a, _Real, 1},
	{c, _Real, 2},
	{o, _Real, 1}		
	}, Module[
		{xxp, yyp, zzp, xxppos, yyppos, zzppos, expa, poscut, pos, poscount, lxpos, lypos, lzpos, expapos, exppos, cpos, apos, sqrtapos, 
			p000, p100, p010, p001, cp000, cp100, cp010, cp001, ctr},

   xxp = x-xp;
   yyp = y-yp;
   zzp = z-zp;
   expa = -a*(xxp^2 + yyp^2 + zzp^2);

(*
   pos = Flatten[Position[expa, _?(# > Log[1.*^-12] &)]];   
*)   

   poscut = Table[If[expa[[j]] > Log[1.*^-12],1,0], {j, 1, np}]; 
   poscount = Total[poscut];
   pos = Table[1, {j, 1, poscount}];
   ctr=1;
   Do[
     If[poscut[[j]]==1, pos[[ctr]] = j; ctr=ctr+1;];
     , {j, 1, np}] ;
   
   xxppos = xxp[[pos]];
   yyppos = yyp[[pos]];
   zzppos = zzp[[pos]];
   lxpos = lx[[pos]];
   lypos = ly[[pos]];
   lzpos = lz[[pos]];
   expapos = expa[[pos]];
   cpos = c[[All,pos]];
   apos = a[[pos]];
   exppos = Exp[-apos * (xxppos^2 + yyppos^2 + zzppos^2)];
   
   p000 = Table[
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	exppos[[j]]
   	,{j, 1, Length[pos]}]; 
   p100 = Table[(
   	If[lxpos[[j]] - 1 < 0, 0, If[(lxpos[[j]] - 1) == 0, 1, lxpos[[j]] * xxppos[[j]]^(lxpos[[j]]-1)]]* 
   	                          If[(lypos[[j]]    ) == 0, 1,              yyppos[[j]]^(lypos[[j]]  )] * 
   	                          If[(lzpos[[j]]    ) == 0, 1,              zzppos[[j]]^(lzpos[[j]]  )] + 
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	(-2 apos[[j]] * xxppos[[j]])) * 
   	exppos[[j]] 	
   	,{j, 1, Length[pos]}]; 
   p010 = Table[(
   	                          If[(lxpos[[j]]    ) == 0, 1,              xxppos[[j]]^(lxpos[[j]]  )] * 
   	If[lypos[[j]] - 1 < 0, 0, If[(lypos[[j]] - 1) == 0, 1, lypos[[j]] * yyppos[[j]]^(lypos[[j]]-1)]]* 
   	                          If[(lzpos[[j]]    ) == 0, 1,              zzppos[[j]]^(lzpos[[j]]  )] + 
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	(-2 apos[[j]] * yyppos[[j]])) * 
    exppos[[j]]
   	,{j, 1, Length[pos]}];
   p001 = Table[(
                              If[(lxpos[[j]]    ) == 0, 1,              xxppos[[j]]^(lxpos[[j]]  )] * 
   	                          If[(lypos[[j]]    ) == 0, 1,              yyppos[[j]]^(lypos[[j]]  )] * 
    If[lzpos[[j]] - 1 < 0, 0, If[(lzpos[[j]] - 1) == 0, 1, lzpos[[j]] * zzppos[[j]]^(lzpos[[j]]-1)]]
   	+ 
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	(-2 apos[[j]] * zzppos[[j]])) * 
   	exppos[[j]]
   	,{j, 1, Length[pos]}]; 
   	
   {cp000,
    cp100,
    cp010,
    cp001} = Transpose[cpos . Transpose[{p000,
                                       p100,
                                       p010,
                                       p001}]];
   {
     Total[Table[o[[i]]*2*cp100[[i]]*cp000[[i]], {i,1,nmo}]],
     Total[Table[o[[i]]*2*cp010[[i]]*cp000[[i]], {i,1,nmo}]],
     Total[Table[o[[i]]*2*cp001[[i]]*cp000[[i]], {i,1,nmo}]] 
   }
	]   
, CompilationTarget -> "C", CompilationOptions -> {"ExpressionOptimization" -> True}];

KineticEnergyDensityGCompiled=Compile[
	{
	{nmo, _Integer},
	{np, _Integer},
	{x, _Real},
	{y, _Real},
	{z, _Real},
	{xp, _Real, 1},
	{yp, _Real, 1},
	{zp, _Real, 1},
	{lx, _Integer, 1},
	{ly, _Integer, 1},
	{lz, _Integer, 1},
	{a, _Real, 1},
	{c, _Real, 2},
	{o, _Real, 1}		
	}, Module[
		{xxp, yyp, zzp, xxppos, yyppos, zzppos, expa, poscut, pos, poscount, lxpos, lypos, lzpos, expapos, exppos, cpos, apos, sqrtapos, 
			p000, p100, p010, p001, cp000, cp100, cp010, cp001, ctr},

   xxp = x-xp;
   yyp = y-yp;
   zzp = z-zp;
   expa = -a*(xxp^2 + yyp^2 + zzp^2);

(*
   pos = Flatten[Position[expa, _?(# > Log[1.*^-12] &)]];   
*)   

   poscut = Table[If[expa[[j]] > Log[1.*^-12],1,0], {j, 1, np}]; 
   poscount = Total[poscut];
   pos = Table[1, {j, 1, poscount}];
   ctr=1;
   Do[
     If[poscut[[j]]==1, pos[[ctr]] = j; ctr=ctr+1;];
     , {j, 1, np}] ;
   
   xxppos = xxp[[pos]];
   yyppos = yyp[[pos]];
   zzppos = zzp[[pos]];
   lxpos = lx[[pos]];
   lypos = ly[[pos]];
   lzpos = lz[[pos]];
   expapos = expa[[pos]];
   cpos = c[[All,pos]];
   apos = a[[pos]];
   exppos = Exp[-apos * (xxppos^2 + yyppos^2 + zzppos^2)];
   
   p000 = Table[
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	exppos[[j]]
   	,{j, 1, Length[pos]}]; 
   p100 = Table[(
   	If[lxpos[[j]] - 1 < 0, 0, If[(lxpos[[j]] - 1) == 0, 1, lxpos[[j]] * xxppos[[j]]^(lxpos[[j]]-1)]]* 
   	                          If[(lypos[[j]]    ) == 0, 1,              yyppos[[j]]^(lypos[[j]]  )] * 
   	                          If[(lzpos[[j]]    ) == 0, 1,              zzppos[[j]]^(lzpos[[j]]  )] + 
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	(-2 apos[[j]] * xxppos[[j]])) * 
   	exppos[[j]] 	
   	,{j, 1, Length[pos]}]; 
   p010 = Table[(
   	                          If[(lxpos[[j]]    ) == 0, 1,              xxppos[[j]]^(lxpos[[j]]  )] * 
   	If[lypos[[j]] - 1 < 0, 0, If[(lypos[[j]] - 1) == 0, 1, lypos[[j]] * yyppos[[j]]^(lypos[[j]]-1)]]* 
   	                          If[(lzpos[[j]]    ) == 0, 1,              zzppos[[j]]^(lzpos[[j]]  )] + 
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	(-2 apos[[j]] * yyppos[[j]])) * 
    exppos[[j]]
   	,{j, 1, Length[pos]}];
   p001 = Table[(
                              If[(lxpos[[j]]    ) == 0, 1,              xxppos[[j]]^(lxpos[[j]]  )] * 
   	                          If[(lypos[[j]]    ) == 0, 1,              yyppos[[j]]^(lypos[[j]]  )] * 
    If[lzpos[[j]] - 1 < 0, 0, If[(lzpos[[j]] - 1) == 0, 1, lzpos[[j]] * zzppos[[j]]^(lzpos[[j]]-1)]]
   	+ 
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	(-2 apos[[j]] * zzppos[[j]])) * 
   	exppos[[j]]
   	,{j, 1, Length[pos]}]; 
   	
   {cp000,
    cp100,
    cp010,
    cp001} = Transpose[cpos . Transpose[{p000,
                                       p100,
                                       p010,
                                       p001}]];
   
   (1/2) Total[Table[ o[[m]] *
{
        cp100[[m]],
        cp010[[m]],
        cp001[[m]]
}
.
{
        cp100[[m]],
        cp010[[m]],
        cp001[[m]]
}
      ,{m,nmo}]]
	]
, CompilationTarget -> "C", CompilationOptions -> {"ExpressionOptimization" -> True}];

ElectronDensityGradientXCompiled=Compile[
	{
	{nmo, _Integer},
	{np, _Integer},
	{x, _Real},
	{y, _Real},
	{z, _Real},
	{xp, _Real, 1},
	{yp, _Real, 1},
	{zp, _Real, 1},
	{lx, _Integer, 1},
	{ly, _Integer, 1},
	{lz, _Integer, 1},
	{a, _Real, 1},
	{c, _Real, 2},
	{o, _Real, 1}		
	}, Module[
		{xxp, yyp, zzp, xxppos, yyppos, zzppos, expa, poscut, pos, poscount, lxpos, lypos, lzpos, expapos, exppos, cpos, apos, sqrtapos, 
			p000, p100, p010, p001, cp000, cp100, cp010, cp001, ctr},

   xxp = x-xp;
   yyp = y-yp;
   zzp = z-zp;
   expa = -a*(xxp^2 + yyp^2 + zzp^2);

(*
   pos = Flatten[Position[expa, _?(# > Log[1.*^-12] &)]];   
*)   

   poscut = Table[If[expa[[j]] > Log[1.*^-12],1,0], {j, 1, np}]; 
   poscount = Total[poscut];
   pos = Table[1, {j, 1, poscount}];
   ctr=1;
   Do[
     If[poscut[[j]]==1, pos[[ctr]] = j; ctr=ctr+1;];
     , {j, 1, np}] ;
   
   xxppos = xxp[[pos]];
   yyppos = yyp[[pos]];
   zzppos = zzp[[pos]];
   lxpos = lx[[pos]];
   lypos = ly[[pos]];
   lzpos = lz[[pos]];
   expapos = expa[[pos]];
   cpos = c[[All,pos]];
   apos = a[[pos]];
   exppos = Exp[-apos * (xxppos^2 + yyppos^2 + zzppos^2)];
   
   p000 = Table[
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	exppos[[j]]
   	,{j, 1, Length[pos]}]; 
   p100 = Table[(
   	If[lxpos[[j]] - 1 < 0, 0, If[(lxpos[[j]] - 1) == 0, 1, lxpos[[j]] * xxppos[[j]]^(lxpos[[j]]-1)]]* 
   	                          If[(lypos[[j]]    ) == 0, 1,              yyppos[[j]]^(lypos[[j]]  )] * 
   	                          If[(lzpos[[j]]    ) == 0, 1,              zzppos[[j]]^(lzpos[[j]]  )] + 
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	(-2 apos[[j]] * xxppos[[j]])) * 
   	exppos[[j]] 	
   	,{j, 1, Length[pos]}]; 
   	
   {cp000,
    cp100} = Transpose[cpos . Transpose[{p000,
                                       p100}]];


     Total[Table[o[[i]]*2*cp100[[i]]*cp000[[i]], {i,1,nmo}]]
     
	]   
, CompilationTarget -> "C", CompilationOptions -> {"ExpressionOptimization" -> True}];

ElectronDensityGradientYCompiled=Compile[
	{
	{nmo, _Integer},
	{np, _Integer},
	{x, _Real},
	{y, _Real},
	{z, _Real},
	{xp, _Real, 1},
	{yp, _Real, 1},
	{zp, _Real, 1},
	{lx, _Integer, 1},
	{ly, _Integer, 1},
	{lz, _Integer, 1},
	{a, _Real, 1},
	{c, _Real, 2},
	{o, _Real, 1}		
	}, Module[
		{xxp, yyp, zzp, xxppos, yyppos, zzppos, expa, poscut, pos, poscount, lxpos, lypos, lzpos, expapos, exppos, cpos, apos, sqrtapos, 
			p000, p100, p010, p001, cp000, cp100, cp010, cp001, ctr},

   xxp = x-xp;
   yyp = y-yp;
   zzp = z-zp;
   expa = -a*(xxp^2 + yyp^2 + zzp^2);

(*
   pos = Flatten[Position[expa, _?(# > Log[1.*^-12] &)]];   
*)   

   poscut = Table[If[expa[[j]] > Log[1.*^-12],1,0], {j, 1, np}]; 
   poscount = Total[poscut];
   pos = Table[1, {j, 1, poscount}];
   ctr=1;
   Do[
     If[poscut[[j]]==1, pos[[ctr]] = j; ctr=ctr+1;];
     , {j, 1, np}] ;
   
   xxppos = xxp[[pos]];
   yyppos = yyp[[pos]];
   zzppos = zzp[[pos]];
   lxpos = lx[[pos]];
   lypos = ly[[pos]];
   lzpos = lz[[pos]];
   expapos = expa[[pos]];
   cpos = c[[All,pos]];
   apos = a[[pos]];
   exppos = Exp[-apos * (xxppos^2 + yyppos^2 + zzppos^2)];
   
   p000 = Table[
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	exppos[[j]]
   	,{j, 1, Length[pos]}]; 
   p010 = Table[(
   	                          If[(lxpos[[j]]    ) == 0, 1,              xxppos[[j]]^(lxpos[[j]]  )] * 
   	If[lypos[[j]] - 1 < 0, 0, If[(lypos[[j]] - 1) == 0, 1, lypos[[j]] * yyppos[[j]]^(lypos[[j]]-1)]]* 
   	                          If[(lzpos[[j]]    ) == 0, 1,              zzppos[[j]]^(lzpos[[j]]  )] + 
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	(-2 apos[[j]] * yyppos[[j]])) * 
    exppos[[j]]
   	,{j, 1, Length[pos]}];
 
   {cp000,
    cp010} = Transpose[cpos . Transpose[{p000,
                                       p010}]];
   
     Total[Table[o[[i]]*2*cp010[[i]]*cp000[[i]], {i,1,nmo}]]

	]   
, CompilationTarget -> "C", CompilationOptions -> {"ExpressionOptimization" -> True}];


ElectronDensityGradientZCompiled=Compile[
	{
	{nmo, _Integer},
	{np, _Integer},
	{x, _Real},
	{y, _Real},
	{z, _Real},
	{xp, _Real, 1},
	{yp, _Real, 1},
	{zp, _Real, 1},
	{lx, _Integer, 1},
	{ly, _Integer, 1},
	{lz, _Integer, 1},
	{a, _Real, 1},
	{c, _Real, 2},
	{o, _Real, 1}		
	}, Module[
		{xxp, yyp, zzp, xxppos, yyppos, zzppos, expa, poscut, pos, poscount, lxpos, lypos, lzpos, expapos, exppos, cpos, apos, sqrtapos, 
			p000, p100, p010, p001, cp000, cp100, cp010, cp001, ctr},

   xxp = x-xp;
   yyp = y-yp;
   zzp = z-zp;
   expa = -a*(xxp^2 + yyp^2 + zzp^2);

(*
   pos = Flatten[Position[expa, _?(# > Log[1.*^-12] &)]];   
*)   

   poscut = Table[If[expa[[j]] > Log[1.*^-12],1,0], {j, 1, np}]; 
   poscount = Total[poscut];
   pos = Table[1, {j, 1, poscount}];
   ctr=1;
   Do[
     If[poscut[[j]]==1, pos[[ctr]] = j; ctr=ctr+1;];
     , {j, 1, np}] ;
   
   xxppos = xxp[[pos]];
   yyppos = yyp[[pos]];
   zzppos = zzp[[pos]];
   lxpos = lx[[pos]];
   lypos = ly[[pos]];
   lzpos = lz[[pos]];
   expapos = expa[[pos]];
   cpos = c[[All,pos]];
   apos = a[[pos]];
   exppos = Exp[-apos * (xxppos^2 + yyppos^2 + zzppos^2)];
   
   p000 = Table[
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	exppos[[j]]
   	,{j, 1, Length[pos]}]; 
   p001 = Table[(
                              If[(lxpos[[j]]    ) == 0, 1,              xxppos[[j]]^(lxpos[[j]]  )] * 
   	                          If[(lypos[[j]]    ) == 0, 1,              yyppos[[j]]^(lypos[[j]]  )] * 
    If[lzpos[[j]] - 1 < 0, 0, If[(lzpos[[j]] - 1) == 0, 1, lzpos[[j]] * zzppos[[j]]^(lzpos[[j]]-1)]]
   	+ 
   	If[lxpos[[j]]==0, 1, xxppos[[j]]^lxpos[[j]]] * 
   	If[lypos[[j]]==0, 1, yyppos[[j]]^lypos[[j]]] * 
   	If[lzpos[[j]]==0, 1, zzppos[[j]]^lzpos[[j]]] * 
   	(-2 apos[[j]] * zzppos[[j]])) * 
   	exppos[[j]]
   	,{j, 1, Length[pos]}]; 
   	
   {cp000,
    cp001} = Transpose[cpos . Transpose[{p000,
                                       p001}]];

     Total[Table[o[[i]]*2*cp001[[i]]*cp000[[i]], {i,1,nmo}]] 

	]   
, CompilationTarget -> "C", CompilationOptions -> {"ExpressionOptimization" -> True}];

ElectronDensityDerivativeCompiled=Compile[
	{
	{nn, _Integer},
	{nd, _Integer},
	{nmo, _Integer},
	{np, _Integer},
	{x, _Real, 1},
	{y, _Real, 1},
	{z, _Real, 1},
	{dx, _Integer, 1},
	{dy, _Integer, 1},
	{dz, _Integer, 1},
	{xp, _Real, 1},
	{yp, _Real, 1},
	{zp, _Real, 1},
	{lx, _Integer, 1},
	{ly, _Integer, 1},
	{lz, _Integer, 1},
	{a, _Real, 1},
	{c, _Real, 2},
	{o, _Real, 1},
	{cutoff, _Real}		
	},
Module[
{nxmax,nymax, nzmax, gamma, pasum, pbsum, pxa, pya, pza, pxb, pyb, pzb,
sqrta, mosum, ax, ay, az, bx, by, bz, expa, binx, biny, binz,
nxkx, nyky, nzkz,
xx0, yy0, zz0,
pa, pb, pat, pbt, px, py, pz, 
ce,
n, d, i, j, k, kx, ky, kz, kkx, kky, kkz, nx, ny, nz, psum, KX, KY, KZ, mynp, tmp,
res,
xxp, yyp, zzp, xxppos, yyppos, zzppos, poscut, pos, poscount, lxpos, lypos, lzpos, expapos, exppos, cpos, apos, sqrtapos, 
			p000, p100, p010, p001, cp000, cp100, cp010, cp001, ctr, npos, nidx, cpa, cpb, kcx, kcy, kcz},



nxmax = Max[dx];
nymax = Max[dy];
nzmax = Max[dz];

binx = Table[Binomial[n,k], {n,0,nxmax}, {k, 0, nxmax}];
biny = Table[Binomial[n,k], {n,0,nymax}, {k, 0, nymax}];
binz = Table[Binomial[n,k], {n,0,nzmax}, {k, 0, nzmax}];
gamma = Table[Gamma[n], {n, 1, Max[Max[lx], Max[ly], Max[lz]]+1}];

res=Table[

   xxp = x[[nidx]] - xp;
   yyp = y[[nidx]] - yp;
   zzp = z[[nidx]] - zp;
   expa = -a*(xxp^2 + yyp^2 + zzp^2);

   poscut = Table[If[expa[[j]] > cutoff,1,0], {j, 1, np}]; 
   poscount = Total[poscut];
   pos = Table[1, {j, 1, poscount}];
   ctr=1;
   Do[
     If[poscut[[j]]==1, pos[[ctr]] = j; ctr=ctr+1;];
     , {j, 1, np}] ;

   npos = Length[pos];

   xxppos = xxp[[pos]];
   yyppos = yyp[[pos]];
   zzppos = zzp[[pos]];
   lxpos = lx[[pos]];
   lypos = ly[[pos]];
   lzpos = lz[[pos]];
   expapos = expa[[pos]];
   cpos = c[[All,pos]];
   apos = a[[pos]];
   sqrtapos = Sqrt[apos];
   exppos = Exp[expapos];
 
ax = Table[
	If[lxpos[[j]] - n < 0, 0, If[(lxpos[[j]] - n)==0, 1, xxppos[[j]]^(lxpos[[j]] - n)] *
     gamma[[lxpos[[j]] + 1]]/gamma[[lxpos[[j]] - n + 1]]]
	, {n, 0, nxmax}, {j, 1, npos}];
ay = Table[
	If[lypos[[j]] - n < 0, 0, If[(lypos[[j]] - n)==0, 1, yyppos[[j]]^(lypos[[j]] - n)] *
     gamma[[lypos[[j]] + 1]]/gamma[[lypos[[j]] - n + 1]]]
	, {n, 0, nymax}, {j, 1, npos}];
az = Table[
	If[lzpos[[j]] - n < 0, 0, If[(lzpos[[j]] - n)==0, 1, zzppos[[j]]^(lzpos[[j]] - n)] *
     gamma[[lzpos[[j]] + 1]]/gamma[[lzpos[[j]] - n + 1]]]
	, {n, 0, nzmax}, {j, 1, npos}];

bx = Table[
(-1)^n*sqrtapos[[j]]^n*
    HermiteH[n, sqrtapos[[j]]*xxppos[[j]]]
	, {n, 0, nxmax}, {j, 1, npos}];
by = Table[
(-1)^n*sqrtapos[[j]]^n*
    HermiteH[n, sqrtapos[[j]]*yyppos[[j]]]
	, {n, 0, nymax}, {j, 1, npos}];
bz = Table[
(-1)^n*sqrtapos[[j]]^n*
    HermiteH[n, sqrtapos[[j]]*zzppos[[j]]]
	, {n, 0, nzmax}, {j, 1, npos}];

Table[

nx=dx[[d]];
ny=dy[[d]];
nz=dz[[d]];

pa=Table[Table[Table[
     Table[
	    Total@Table[
             binx[[(nx - kx) + 1, KX + 1]]*
              ax[[(nx - kx) - KX + 1, j]]*
              bx[[KX + 1, j]]
              , {KX, 0, nx - kx}] *
        Total@Table[
               biny[[(ny-ky) + 1, KY + 1]]*
                ay[[(ny-ky) - KY + 1, j]]*
                by[[KY + 1, j]]
              , {KY, 0, ny - ky}] *
        Total@Table[
               binz[[(nz-kz) + 1, KZ + 1]]*
                az[[(nz-kz) - KZ + 1, j]]*
                bz[[KZ + 1, j]]
              , {KZ, 0, nz - kz}] *
                exppos[[j]]
          , {j, 1, npos}]
        , {kz, 0, nz}]
      , {ky, 0, ny}]
    , {kx, 0, nx}];

pb=Table[Table[Table[
	Table[
		Total@Table[
             binx[[kx + 1, KX + 1]]*
              ax[[(kx - KX) + 1, j]]*
              bx[[KX + 1, j]]
              , {KX, 0, kx}] *
        Total@Table[
               biny[[ky + 1, KY + 1]]*
                ay[[(ky - KY) + 1, j]]*
                by[[KY + 1, j]]
              , {KY, 0, ky}] *
        Total@Table[
                 binz[[kz + 1, KZ + 1]]*
                  az[[(kz - KZ) + 1, j]]*
                  bz[[KZ + 1, j]]
             , {KZ, 0, kz}]*
              exppos[[j]]
          , {j, 1, npos}]
        , {kz, 0, nz}]
      , {ky, 0, ny}]
    , {kx, 0, nx}];

cpa = Table[Table[Table[
	cpos . pa[[kx+1,ky+1,kz+1]]
	    , {kz, 0, nz}]
      , {ky, 0, ny}]
    , {kx, 0, nx}];

cpb = Table[Table[Table[
	cpos . pb[[kx+1,ky+1,kz+1]]
	     , {kz, 0, nz}]
      , {ky, 0, ny}]
    , {kx, 0, nx}];
 
Total@Table[
  o[[i]]*
   Total@Table[
    binx[[nx + 1, kx + 1]]*
     Total@Table[
       biny[[ny + 1, ky + 1]]*
       Total@Table[
         binz[[nz + 1, kz + 1]]*
	  cpa[[kx+1, ky+1, kz+1, i]]
         *
	  cpb[[kx+1, ky+1, kz+1, i]]
        , {kz, 0, nz}]
      , {ky, 0, ny}]
    , {kx, 0, nx}]
 , {i, 1, nmo}]
, {d, 1, nd}]
, {nidx, 1, nn}];

res

]

, CompilationTarget -> "C", CompilationOptions -> {"ExpressionOptimization" -> True}];

Options[ElectronDensityDerivative] = {Cutoff -> Log[1.*^-8]};
ElectronDensityDerivative[
   wfn_,
   {x_?NumericQ, y_?NumericQ, z_?NumericQ}, {nx_?IntegerQ, 
    ny_?IntegerQ, nz_?IntegerQ}, opts : OptionsPattern[]] := Module[
   {nmo = wfn["NumberOfOccupiedMolecularOrbitals"], 
    np = wfn["NumberOfPrimitives"], lx = wfn["lx"], ly = wfn["ly"], 
    lz = wfn["lz"], xp = wfn["xp"], yp = wfn["yp"], zp = wfn["zp"], 
    a = wfn["PrimitiveExponents"], 
    o = wfn["MolecularOrbitalOccupationNumbers"], 
    c = wfn["MolecularOrbitalPrimitiveCoefficients"], link, res},
   First@
    First@
     ElectronDensityDerivativeCompiled[1, 1, nmo, 
      np, {x}, {y}, {z}, {nx}, {ny}, {nz}, xp, yp, zp, lx, ly, lz, a, 
      c, o, OptionValue[Cutoff]]
   ];
ElectronDensityDerivative[
   wfn_,
   xyz_?MatrixQ, deriv_?MatrixQ, opts : OptionsPattern[]] := Module[
   {nn = Length[xyz], nd = Length[deriv], 
    nmo = wfn["NumberOfOccupiedMolecularOrbitals"], 
    np = wfn["NumberOfPrimitives"], lx = wfn["lx"], ly = wfn["ly"], 
    lz = wfn["lz"], xp = wfn["xp"], yp = wfn["yp"], zp = wfn["zp"], 
    a = wfn["PrimitiveExponents"], 
    o = wfn["MolecularOrbitalOccupationNumbers"], 
    c = wfn["MolecularOrbitalPrimitiveCoefficients"], res},
   ElectronDensityDerivativeCompiled[nn, nd, nmo, np, xyz[[All, 1]], 
     xyz[[All, 2]], xyz[[All, 3]], deriv[[All, 1]], deriv[[All, 2]], 
     deriv[[All, 3]], xp, yp, zp, lx, ly, lz, a, c, o, OptionValue[Cutoff]]
   ];
(*
ReadWavefunctionFromMOLDEN[filename_] := Module[
      {fn = filename, fc, lenfc, atomsMatches, atomsStart, atomsEnd, 
        gtoMatches, gtoStart, gtoEnd, gtoSubstring, matchString, 
        matches,
        start, end,
        moStart, moEnd, moMatches, moSubstring,
        natom,
        atomNames, atomNumbers, atomCharges, atomXCoordinates, 
        atomYCoordinates, atomZCoordinates,
        atomName, atomNumber, atomCharge, atomXCoordinate, 
        atomYCoordinate, atomZCoordinate, atomCoordinates,
        substring, substringStream, 
        atomgtoEnd,
        dummy, amomLetter, nprimsInBlock, ctr, w, z, mos, number,
        nmo, nprim, orbe, occno, mocoef, coef, amom, 
        norm, amoms, aset, ws, zs, zl, wl, currentAtom, x0y0z0, m,
        energyMatch, energyLine, energy
        },
      
      fc = ReadList[fn, Record, RecordSeparators -> "\n", NullRecords -> True];
      lenfc = Length[fc];
      
      matchString = Table["_ENERGD=*", {lenfc}];
      energyMatch = MapThread[StringMatchQ, {fc, matchString}];
      energyLine = First[First[Position[energyMatch, True]]];
      substringStream = StringToStream[ fc[[energyLine]] ];
      {dummy, energy} = Read[substringStream, {Word, Number}];
      
      (* [Atoms] *)
      
      (* start *)
      matchString = Table["[Atoms]*", {lenfc}];
      atomsMatches = MapThread[StringMatchQ, {fc, matchString}];
      atomsStart = First[First[Position[atomsMatches, True]]] + 1;
      
      (* end *)
      matchString = Table["[GTO]", {lenfc}];
      atomsMatches = MapThread[StringMatchQ, {fc, matchString}];
      atomsEnd = First[First[Position[atomsMatches, True]]] - 1;
      
      natom = atomsEnd - atomsStart + 1;
      
      atomNames = {};
      atomNumbers = {};
      atomCharges = {};
      atomXCoordinates = {};
      atomYCoordinates = {};
      atomZCoordinates = {};
      
      substring = fc[[atomsStart ;; atomsEnd]];
      
      
      Do[
        substringStream = StringToStream[substring[[i]] ];
        {atomName, atomNumber, atomCharge, atomXCoordinate, 
            atomYCoordinate, atomZCoordinate} = 
          Read[substringStream, {Word,   Number, Number,       Number, 
              Number, Number}];
        
        AppendTo[atomNames, atomName];
        AppendTo[atomNumbers, atomNumber];
        AppendTo[atomCharges, atomCharge];
        AppendTo[atomXCoordinates, atomXCoordinate];
        AppendTo[atomYCoordinates, atomYCoordinate];
        AppendTo[atomZCoordinates, atomZCoordinate];
        , {i, natom}
        ];
      atomCoordinates = 
        1.889725989*
          
     Transpose[{atomXCoordinates, atomYCoordinates, atomZCoordinates}];
      
      (* [MO] *)
      matchString = Table["[MO]", {lenfc}];
      moMatches = MapThread[StringMatchQ, {fc, matchString}];
      moStart = First[First[Position[moMatches, True]]] + 1;
      moEnd = Length[fc];
      
      moSubstring = fc[[moStart ;; moEnd]];
      
      (* count the number of molecular orbitals *)
      
      matchString = Table["*Ene=*", {Length[moSubstring]}];
      matches = MapThread[StringMatchQ, {moSubstring, matchString}];
      
      nmo = Count[matches, True];
      
      orbe = {};
      occno = {};
      mocoef = {};
      
      ctr = 1;
      Do[
        
        substringStream = StringToStream[moSubstring[[ctr]]];
        {dummy, number} = Read[substringStream, {Word, Number}];
        AppendTo[orbe, number];
        ctr++;
        
        (*skip spin*)
        
        substringStream = StringToStream[moSubstring[[ctr]]];
        {dummy, dummy} = Read[substringStream, {Word, Word}];
        ctr++;
        
        substringStream = StringToStream[moSubstring[[ctr]]];
        {dummy, number} = Read[substringStream, {Word, Number}];
        AppendTo[occno, number];
        ctr++;
        
        mos = {};
        Do[
          substringStream = StringToStream[moSubstring[[ctr]]];
          {dummy, number} = Read[substringStream, {Number, Number}];
          AppendTo[mos, number];
          ctr++;
          , {i, nmo}];
        
        AppendTo[mocoef, mos];
        (*ctr++;*)
        
        , {i, nmo}];
      
      coef = Table[{}, {m, 1, nmo}];
      ws = {};
      zs = {};
      
      (* [GTO] *)
      
      (* start *)
      matchString = Table["[GTO]", {lenfc}];
      gtoMatches = MapThread[StringMatchQ, {fc, matchString}];
      gtoStart = First[First[Position[atomsMatches, True]]] + 1;
      
      (* end *)
      matchString = Table["[MO]", {lenfc}];
      gtoMatches = MapThread[StringMatchQ, {fc, matchString}];
      gtoEnd = First[First[Position[gtoMatches, True]]] - 2;
      
      gtoSubstring = fc[[gtoStart ;; gtoEnd]];
      
      matchString = Table[Whitespace, {Length[gtoSubstring]}];
      atomgtoEnd = 
        Flatten[Position[
            MapThread[StringMatchQ, {gtoSubstring, matchString}], 
      True]];
      
      norm[a_, {lx_, ly_, 
            lz_}] := (2 a / Pi)^(3/
                4) Sqrt[(  (4 a)^(lx + ly + 
                      lz)  )/(   (2 lx - 1)!! (2 ly - 1)!! (2 lz - 
                        1)!!  )];(* ((2/Pi)^(3/4) (2^(lx+ly+
      lz) a^((2 lx+2 ly+2 lz+3)/
      4))/((2 lx-1)!! (2 ly-1)!! (2 lz-1)!!));*)
      
      amom = {};
      x0y0z0 = {};
      
      ctr = 1;
      start = 1;
      m = 1; (* mo counter *)
      nprim = 0; (* 
   primitive counter *)
      
      
      Do[
        end = atomgtoEnd[[a]];
        substring = gtoSubstring[[start ;; end]];
        substringStream = StringToStream[substring[[ctr]]];
        {currentAtom, dummy} = Read[substringStream, {Number, Number}];
        ctr++;
        While[ctr < end - start + 1,
          substringStream = StringToStream[substring[[ctr]]];
          {amomLetter, nprimsInBlock, dummy} = 
            Read[substringStream, {Word, Number, Number}];
          ctr++;
          
          nprim = 
            nprim + nprimsInBlock*
                
        Switch[amomLetter, "s", 1, "p", 3, "d", 6, "f", 10, "g", 
                  15, _, Print["Angular Momentum Too Large."]; Null];
          
          amoms[letter_] := Switch[
              letter,
              "s", {{0, 0, 0}},
              "p", {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
              "d", {{2, 0, 0}, {0, 2, 0}, {0, 0, 2}, {1, 1, 0}, {1, 0, 
                  1}, {0, 1, 1}},
              (*{{2,0,0},{1,1,0},{0,2,0},{1,0,1},{0,1,1},{0,0,
       2}},*)
              
              "f", {{3, 0, 0}, {0, 3, 0}, {0, 0, 3}, {1, 2, 0}, {2, 1, 
                  0}, {2, 0, 1}, {1, 0, 2}, {0, 1, 2}, {0, 2, 1}, {1, 
         1, 1}},
              "g", {{4, 0, 0}, {0, 4, 0}, {0, 0, 4}, {3, 1, 0}, {3, 0, 
                  1}, {1, 3, 0}, {0, 3, 1}, {1, 0, 3}, {0, 1, 3}, {2, 
         2, 
                  0}, {2, 0, 2}, {0, 2, 2}, {2, 1, 1}, {1, 2, 1}, {1, 
         1, 2}},
              _, Print["Angular Momentum Too Large."]; Null];
          
          zl = {};
          wl = {};
          
          Do[
            substringStream = StringToStream[substring[[ctr]]];
            {z, w} = Read[substringStream, {Number, Number}];
            AppendTo[wl, w];
            AppendTo[zl, z];
            ctr++;
            , {p, nprimsInBlock}];
          
          aset = amoms[amomLetter];
          Do[
            Do[
              (* Print[aset[[s]]]; *)
              
              AppendTo[x0y0z0, atomCoordinates[[currentAtom]]];
              AppendTo[amom, aset[[s]]];
              
              AppendTo[ws, wl[[p]]];
              AppendTo[zs, zl[[p]]];
              
              (*
              Print[{mocoef[[1,m]],ws[[p]],norm[
       zs[[p]],
              aset[[s]]]} ];
              Print[mocoef[[1,m]]*ws[[p]]*norm[zs[[p]],aset[[s]]] ];
              *)
              
              Do[
                AppendTo[coef[[mo]], 
                  mocoef[[mo, m]]*wl[[p]]*norm[zl[[p]], aset[[s]]]  ]
                , {mo, 1, nmo}];
              
              , {p, nprimsInBlock}];
            m++;
            , {s, Length[aset]}
            ];
          
          ];
        
        start = end + 1;
        ctr = 1;
        , {a, natom}
        ];
      
        Association[
   "NumberOfNuclei" -> natom,
   "NumberOfOccupiedMolecularOrbitals" -> nmo,
   "AtomicNumbers" -> atomCharges,
   "NuclearCartesianCoordinates" -> atomCoordinates,
   "NumberOfPrimitives" -> nprim,
   "PrimitiveTypes" -> Missing[],
   "PrimitiveCenters" -> x0y0z0,
   "PrimitiveExponents" -> zs,
   "xp" -> x0y0z0[[All,1]],
   "yp" ->x0y0z0[[All,2]],
   "zp" -> x0y0z0[[All,3]],
   "l" -> amom,
   "lx" -> amom[[All,1]],
   "ly" -> amom[[All,2]],
   "lz" -> amom[[All,3]],
   "MolecularOrbitalOccupationNumbers" -> 
    occno,
   "MolecularOrbitalPrimitiveCoefficients" -> 
    coef,
   "cFlatten" -> Flatten[coef],
   "cTranspose" -> Transpose[coef],
   "cFlattenTranspose" -> Flatten[Transpose[coef]],
   "Energy" -> N[energy],
   "VirialRatio" -> Missing[]
   ]
     
];
*)

ReadWavefunctionFromPySCF[session_] := Module[
	{extractWavefunctionInput,
	nprim,
	nmo,
	nnuc,
	XYZ,
	Znuc,
	centers,
	center,
	t,
	lx, ly, lz, l,
	a,
	o,
	c
	},


extractWavefunctionInput = "
#!/usr/bin/env python
# Copyright 2014-2020 The PySCF Developers. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the \"License\");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an \"AS IS\" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or \
implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Qiming Sun <osirpt.sun@gmail.com>
#         Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>
#         Modifications by: Eric Brown <ecbrown@ericcbrown.com>
#

# types
# 1 S
# 2 PX
# 3 PY
# 4 PZ
# 5 DXX
# 6 DYY
# 7 DZZ
# 8 DXY
# 9 DXZ
# 10 DYZ
# 11 FXXX
# 12 FYYY
# 13 FZZZ
# 14 FXXY
# 15 FXXZ
# 16 FYYZ
# 17 FXYY
# 18 FXZZ
# 19 FYZZ
# 20 FXYZ
# 21 GXXXX
# 22 GYYYY
# 23 GZZZZ
# 24 GXXXY
# 25 GXXXZ
# 26 GXYYY
# 27 GYYYZ
# 28 GXZZZ
# 29 GYZZZ
# 30 GXXYY
# 31 GXXZZ
# 32 GYYZZ
# 33 GXXYZ
# 34 GXYYZ
# 35 GXYZZ
# 36 HZZZZZ
# 37 HYZZZZ
# 38 HYYZZZ
# 39 HYYYZZ
# 40 HYYYYZ
# 41 HYYYYY
# 42 HXZZZZ
# 43 HXYZZZ
# 44 HXYYZZ
# 45 HXYYYZ
# 46 HXYYYY
# 47 HXXZZZ
# 48 HXXYZZ
# 49 HXXYYZ
# 50 HXXYYY
# 51 HXXXZZ
# 52 HXXXYZ
# 53 HXXXYY
# 54 HXXXXZ
# 55 HXXXXY
# 56 HXXXXX
TYPE_MAP = [
      [1],  # S
        [2, 3, 4],  # P
        [5, 8, 9, 6, 10, 7],  # D
        [11, 14, 15, 17, 20, 18, 12, 16, 19, 13],  # F
        [21, 24, 25, 30, 33, 31, 26, 34, 35, 28, 22, 27, 32, 29, 23],  # G
        [56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36],  # H
  ]


if mol.cart:
    raise NotImplementedError('Cartesian basis not available')

#FIXME: Duplicated primitives may lead to problems.x2c._uncontract_mol
# is the workaround at the moment to remove duplicated primitives.
nmo = mc.ncore + mc.ncas
nelecas = mc.nelecas[0] + mc.nelecas[1]
casdm1, casdm2 = mc.fcisolver.make_rdm12(mc.ci, mc.ncas, mc.nelecas)
rdm1, rdm2 = mcscf.addons._make_rdm12_on_mo(casdm1, casdm2, mc.ncore, mc.ncas, nmo)
orbsym =  symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, mc.mo_coeff[:, : nmo])
natocc, natorb = symm.eigh(-rdm1, orbsym)
for i, k in enumerate(numpy.argmax(abs(natorb), axis = 0)) :
      if natorb[k, i] < 0 :
          natorb[:, i] *= -1
natorb = numpy.dot(mc.mo_coeff[:, :nmo], natorb)
natocc = -natocc

mol, ctr = x2c._uncontract_mol(mol, True, 0.)
mo_coeff=natorb
mo_coeff = numpy.dot(ctr, mo_coeff)

nmo = mo_coeff.shape[1]
mo_cart = []
centers = []
types = []
exps = []
p0 = 0
for ib in range(mol.nbas):
    ia = mol.bas_atom(ib)
    l = mol.bas_angular(ib)
    es = mol.bas_exp(ib)
    c = mol._libcint_ctr_coeff(ib)
    np, nc = c.shape
    nd = nc*(2*l+1)
    mosub = mo_coeff[p0:p0+nd].reshape(-1,nc,nmo)
    c2s = gto.cart2sph(l)
    mosub = numpy.einsum('yki,cy,pk->pci', mosub, c2s, c)
    mo_cart.append(mosub.transpose(1,0,2).reshape(-1,nmo))
    for t in TYPE_MAP[l]:
        types.append([t]*np)
    ncart = mol.bas_len_cart(ib)
    exps.extend([es]*ncart)
    centers.extend([ia+1]*(np*ncart))
    p0 += nd
mo_cart = numpy.vstack(mo_cart)
centers = numpy.hstack(centers)
types = numpy.hstack(types)
exps = numpy.hstack(exps)
nprim, nmo = mo_cart.shape
";

ExternalEvaluate[session, extractWavefunctionInput];

nprim = ExternalEvaluate[session, "
nprim
"
   ] // Normal;
   
nmo = ExternalEvaluate[session, "
nmo
"
   ] // Normal;

nnuc = ExternalEvaluate[session, "
mol.natm
"
   ] // Normal;
 
XYZ = Table[
  ExternalEvaluate[session, "
mol.atom_coord(" <> ToString[i - 1] <> ")
"
    ] // Normal
  , {i, 1, nnuc}];
  
Znuc = Table[
  ExternalEvaluate[session, "
mol.atom_charge(" <> ToString[i - 1] <> ")
"
    ] // Normal
  , {i, 1, nnuc}];
  
centers = ExternalEvaluate[session, "
centers
"
    ] // Normal;
    
center = Table[
   {XYZ[[centers[[i]], 1]],
    XYZ[[centers[[i]], 2]],
    XYZ[[centers[[i]], 3]]
    },
   {i, 1, nprim}];

t = ExternalEvaluate[session, "
types
"
    ] // Normal;
    
lx = Table[0, {j, 1, nprim}];
ly = Table[0, {j, 1, nprim}];
lz = Table[0, {j, 1, nprim}];
l = Table[
   Which[ 
    t[[j]] == 1, {0, 0, 0},
    t[[j]] == 2, {1, 0, 0},
    t[[j]] == 3, {0, 1, 0},
    t[[j]] == 4, {0, 0, 1},
    t[[j]] == 5, {2, 0, 0},
    t[[j]] == 6, {0, 2, 0},
    t[[j]] == 7, {0, 0, 2},
    t[[j]] == 8, {1, 1, 0},
    t[[j]] == 9, {1, 0, 1},
    t[[j]] == 10, {0, 1, 1},
    t[[j]] == 11, {3, 0, 0},
    t[[j]] == 12, {0, 3, 0},
    t[[j]] == 13, {0, 0, 3},
    t[[j]] == 14, {2, 1, 0},
    t[[j]] == 15, {2, 0, 1},
    t[[j]] == 16, {0, 2, 1},
    t[[j]] == 17, {1, 2, 0},
    t[[j]] == 18, {1, 0, 2},
    t[[j]] == 19, {0, 1, 2},
    t[[j]] == 20, {1, 1, 1}
    ]
   , {j, 1, nprim}];
lx = l[[All, 1]];
ly = l[[All, 2]];
lz = l[[All, 3]];

a = ExternalEvaluate[session, "
exps
"
    ] // Normal;
    
o = ExternalEvaluate[session, "
natocc
#occ
"
    ] // Normal;

c = ExternalEvaluate[session, "
mo_cart
"
    ] // Normal;
    
 
  Association[
   "NumberOfNuclei" -> nnuc,
   "NumberOfOccupiedMolecularOrbitals" -> nmo,
   "AtomicNumbers" -> Znuc,
   "NuclearCartesianCoordinates" -> XYZ,
   "NumberOfPrimitives" -> nprim,
   "PrimitiveTypes" -> t,
   "PrimitiveCenters" -> center,
   "PrimitiveExponents" -> a,
   "xp" -> Table[
     (XYZ[[centers[[j]]]])[[1]], {j, 
      1, nprim}],
    "yp" -> Table[
     (XYZ[[centers[[j]]]])[[2]], {j, 
      1, nprim}],
   "zp" -> Table[
     (XYZ[[centers[[j]]]])[[3]], {j, 
      1, nprim}],
   "l" -> l,
   "lx" -> lx,
   "ly" -> ly,
   "lz" -> lz,
   "MolecularOrbitalOccupationNumbers" -> 
    o,
   "MolecularOrbitalPrimitiveCoefficients" -> 
    Transpose[c],
   "cFlatten" -> Flatten[c],
   "cTranspose" -> Transpose[c],
   "cFlattenTranspose" -> Flatten[Transpose[c]],
   "Energy" -> Missing[],
   "VirialRatio" -> Missing[]
   ]
];

Options[ElectronDensityDerivative] = {Cutoff -> Log[1.*^-12]};
rho[w_, {x_?NumericQ, y_?NumericQ, z_?NumericQ}] := 
  ElectronDensityCompiled[w["NumberOfOccupiedMolecularOrbitals"], 
   w["NumberOfPrimitives"], x, y, z, w["xp"], w["yp"], w["zp"], 
   w["lx"], w["ly"], w["lz"], w["PrimitiveExponents"], 
   w["MolecularOrbitalPrimitiveCoefficients"], 
   w["MolecularOrbitalOccupationNumbers"]];
g[w_, {x_?NumericQ, y_?NumericQ, z_?NumericQ}] := 
  ElectronDensityGradientCompiled[
   w["NumberOfOccupiedMolecularOrbitals"], w["NumberOfPrimitives"], x,
    y, z, w["xp"], w["yp"], w["zp"], w["lx"], w["ly"], w["lz"], 
   w["PrimitiveExponents"], 
   w["MolecularOrbitalPrimitiveCoefficients"], 
   w["MolecularOrbitalOccupationNumbers"]];
gx[w_, {x_?NumericQ, y_?NumericQ, z_?NumericQ}] := 
  ElectronDensityGradientXCompiled[
   w["NumberOfOccupiedMolecularOrbitals"],
   w["NumberOfPrimitives"], x, y, z, w["xp"], w["yp"], w["zp"], 
   w["lx"], w["ly"], w["lz"], w["PrimitiveExponents"], 
   w["MolecularOrbitalPrimitiveCoefficients"], 
   w["MolecularOrbitalOccupationNumbers"]];
gy[w_, {x_?NumericQ, y_?NumericQ, z_?NumericQ}] := 
  ElectronDensityGradientYCompiled[
   w["NumberOfOccupiedMolecularOrbitals"], w["NumberOfPrimitives"], x,
    y, z, w["xp"], w["yp"], w["zp"], w["lx"], w["ly"], w["lz"], 
   w["PrimitiveExponents"], 
   w["MolecularOrbitalPrimitiveCoefficients"], 
   w["MolecularOrbitalOccupationNumbers"]];
gz[w_, {x_?NumericQ, y_?NumericQ, z_?NumericQ}] := 
  ElectronDensityGradientZCompiled[
   w["NumberOfOccupiedMolecularOrbitals"], w["NumberOfPrimitives"], x,
    y, z, w["xp"], w["yp"], w["zp"], w["lx"], w["ly"], w["lz"], 
   w["PrimitiveExponents"], 
   w["MolecularOrbitalPrimitiveCoefficients"], 
   w["MolecularOrbitalOccupationNumbers"]];
H[w_, {x_?NumericQ, y_?NumericQ, z_?NumericQ}] := Module[
   {h},
   h = First[
     ElectronDensityDerivative[
      w, {{x, y, z}}, {{2, 0, 0}, {1, 1, 0}, {1, 0, 1}, {0, 2, 0}, {0,
         1, 1}, {0, 0, 2}}]];
   {{h[[1]], h[[2]], h[[3]]}, {h[[2]], h[[4]], h[[5]]}, {h[[3]], 
     h[[5]], h[[6]]}}];

KineticEnergyDensityG[w_, {x_?NumericQ, y_?NumericQ, z_?NumericQ}] :=   
KineticEnergyDensityGCompiled[
   w["NumberOfOccupiedMolecularOrbitals"], w["NumberOfPrimitives"], x,
    y, z, w["xp"], w["yp"], w["zp"], w["lx"], w["ly"], w["lz"], 
   w["PrimitiveExponents"], 
   w["MolecularOrbitalPrimitiveCoefficients"], 
   w["MolecularOrbitalOccupationNumbers"]];

KineticEnergyDensityK[w_, {x_?NumericQ, y_?NumericQ, z_?NumericQ}] := Module[
	{G,L},
	G=KineticEnergyDensityG[w,{x,y,z}];
	L = Total[First[
     ElectronDensityDerivative[
      w, {{x, y, z}}, {{2, 0, 0}, {0, 2, 0}, {0, 0, 2}}]]];
	G - L
];

StressTensorTrace[w_, {x_?NumericQ, y_?NumericQ, z_?NumericQ}] := Module[
	{G,K},
	G=KineticEnergyDensityG[w,{x,y,z}];
	K=KineticEnergyDensityK[w,{x,y,z}];
	-K - G
];

LocateNuclearCriticalPoints[wfn_] := Module[
{soln},

  ParallelTable[
  	
(* If maximization fails, it is almost certainly
   because of starting at essentially the maximum, 
   e.g. heavy atoms *)

Check[Quiet[
    soln = FindMaximum[rho[wfn,{x,y,z}],{{x, xyz[[1]]}, {y, xyz[[2]]}, {z, xyz[[3]]}},
                Gradient :> g[wfn,{x,y,z}]
                ,Method -> {"Newton", "Hessian" :> H[wfn,{x,y,z}]}];
    {x,y,z} /. soln[[2]]
], xyz]

   ,{xyz, wfn["NuclearCartesianCoordinates"]}
   , Method -> "FinestGrained"]
];

Options[AssociatedNuclearAttractor] = 
  Sort[Join[Options[NDSolveValue], {Method -> "Adams", MaxPathLength -> 100000, BetaSphereRadii -> {} }]];
AssociatedNuclearAttractor[wfn_, ncps_, {x0_?NumericQ, y0_?NumericQ, z0_?NumericQ}, opts : OptionsPattern[]] := Module[
{soln, nearestTable, end, endpoint, cutoffs, distances, ctr, res},
 
 If[ rho[wfn, {x0,y0,z0}] < Min[$QTAIMContoursPositive], Return[0] ];
 
 cutoffs=If[ Length[OptionValue[BetaSphereRadii]] == 0, 
 	Table[ 
 	If[ Round[(wfn["AtomicNumbers"])[[n]]] == 1,
 		0.1,
 		0.3] 
 	,{n,1,Length[ncps]}],
    OptionValue[BetaSphereRadii]
 ];

 distances=Table[EuclideanDistance[ ncps[[i]], {x0, y0, z0}], {i,1,Length[ncps]}];

 If[ AnyTrue[  distances - cutoffs,  (# < 0) &],
   res = Ordering[distances-cutoffs, 1];
   If[Length[res]>1, 0, First@res]
 ]; 


ctr = 0;

Quiet[
 	Check[
    soln = NDSolveValue[
     {
       X'[s] == g[wfn, X[s]],
       X[0] == {x0, y0, z0},
 WhenEvent[
 	If[ ctr < 2, False,
    distances=Table[EuclideanDistance[ ncps[[i]], X[s]], {i,1,Length[ncps]}];
    AnyTrue[  distances - cutoffs,  (# < 0) &]
 	]
,
end=s;
"StopIntegration"
]
     }, 
     X[s],
     {s, 0, 100000 (*smax*) }
     , Method -> Switch[OptionValue[Method],
     "ABM", AdamsBM,
     "BS32", BS32,
     "BS54", BS54,
     "CashKarp45", CashKarp45,
     "DOPRI54", DOPRI54,
     "Fehlberg45", Fehlberg45,
     "RK4", RK4,   
     _, OptionValue[Method]]
     , MaxSteps -> 1000 
     (*, StartingStepSize -> 0.1*)
     ,StepMonitor :> (ctr++(*; Print[ {x[s], y[s], z[s]} ] *) )
     (*, AccuracyGoal->3, PrecisionGoal->Infinity*)
     ];

     endpoint=soln /. {s->end};
     
     distances=Table[EuclideanDistance[ ncps[[i]], endpoint], {i,1,Length[ncps]}];
     res = Ordering[distances-cutoffs, 1];
     If[Length[res]>1, 0, First@res]

, 0]
]

];

Options[SteepestAscentPath] = 
  Sort[Join[Options[NDSolveValue], {Method -> "Adams", AccuracyGoal -> Automatic, PrecisionGoal -> Automatic}]];
SteepestAscentPath[wfn_, ncps_, {x0_?NumericQ, y0_?NumericQ, z0_?NumericQ}, opts : OptionsPattern[]] := Module[
{soln, end, distances, cutoffs, ctr},
	
 ctr = 0;
 cutoffs=Table[ 
 	
 	If[ Round[(wfn["AtomicNumbers"])[[n]]] == 1,
 		0.01,
 		0.01] 
 	
 	,{n,1,Length[ncps]}];

 
    soln = Reap[
    	NDSolveValue[
     {
       X'[s]==g[wfn, X[s]],
       X[0]=={x0,y0,z0},
 WhenEvent[
 	If[ ctr < 2, False,
    distances=Table[EuclideanDistance[ ncps[[i]], X[s]], {i,1,Length[ncps]}];
    AnyTrue[  distances - cutoffs,  (# < 0) &]
 	]
,
end=s;
Sow[{s,X[s]}];
"StopIntegration"
]
     }, 
     X[s],
     {s, 0, 100000 (*smax*) }
     , Method -> Switch[OptionValue[Method],
     "ABM", AdamsBM,
     "BS32", BS32,
     "BS54", BS54,
     "CashKarp45", CashKarp45,
     "DOPRI54", DOPRI54,
     "Fehlberg45", Fehlberg45,
     "RK4", RK4,   
     _, OptionValue[Method]]
     (*, StartingStepSize -> 0.1*)
     , StepMonitor :> (ctr++; Sow[{s,X[s]}])
     , AccuracyGoal -> OptionValue[AccuracyGoal], 
     PrecisionGoal -> OptionValue[PrecisionGoal]
     ]
    ];

    soln

];

Options[SteepestDescentPath] = 
  Sort[Join[Options[NDSolveValue], {Method -> "Adams"}]];
SteepestDescentPath[wfn_, ncps_, {x0_?NumericQ, y0_?NumericQ, z0_?NumericQ}, opts : OptionsPattern[]] := Module[
{soln, end},

    soln = Reap[
    	NDSolveValue[
     {
       X'[s]==-g[wfn,X[s]],
       X[0]=={x0,y0,z0},
 WhenEvent[
 
     rho[wfn, X[s]]  < Min[$QTAIMContoursPositive]

,
end=s;
Sow[{s,X[s]}];
"StopIntegration"
]
     }, 
     X[s],
     {s, 0, 100000 (*smax*) }
     , Method -> Switch[OptionValue[Method],
     "ABM", AdamsBM,
     "BS32", BS32,
     "BS54", BS54,
     "CashKarp45", CashKarp45,
     "DOPRI54", DOPRI54,
     "Fehlberg45", Fehlberg45,
     "RK4", RK4,
     _, OptionValue[Method]]
     (*, StartingStepSize -> 0.1*)
     , StepMonitor :> Sow[{s,X[s]}]
     (*, AccuracyGoal->3, PrecisionGoal->Infinity*)
     ]
    ];

    soln

];

BondPath[wfn_, ncps_, {x0_?NumericQ, y0_?NumericQ, z0_?NumericQ}] := Module[
{},
(* Two steepest ascent paths, forward and backward along the lowest eigenvector of the Hessian *)
{
SteepestAscentPath[wfn, ncps, {x0,y0,z0} + First[-0.001*Eigenvectors[H[wfn,{x0,y0,z0}], -1]]],
SteepestAscentPath[wfn, ncps, {x0,y0,z0} + First[+0.001*Eigenvectors[H[wfn,{x0,y0,z0}], -1]]]
}
]

Options[InteratomicSurface] = 
  Sort[Join[
    Options[FindRoot], {MaxIterations -> 100, StoppingCriterion -> 10^-3, BetaSphereRadii -> {}}]];
InteratomicSurface[w_, ncps_, ncp_, {theta_?NumericQ, phi_?NumericQ}, 
  opts : OptionsPattern[]] := Module[
  {x0, y0, z0, xr, yr, zr, rmax, r0, r1, ctr, rprev},
  {x0, y0, z0} = ncps[[ncp]];

  If[Length[OptionValue[BetaSphereRadii]] > 0,
   r0 = (3/4) (OptionValue[BetaSphereRadii])[[ncp]];
   r1 = (OptionValue[BetaSphereRadii])[[ncp]];
   ,
   r0 = 0.1;
   r1 = 0.2;
   ];
  rprev = 0;
  ctr=0;
  rmax = Quiet[Catch[
    r /. FindRoot[
      If[AssociatedNuclearAttractor[w, ncps,
         {r Cos[phi] Sin[theta], r Sin[phi] Sin[theta], 
           r Cos[theta]} + {x0, y0, z0}, 
         BetaSphereRadii -> OptionValue[BetaSphereRadii]] == ncp,
       rho[
        w, {r Cos[phi] Sin[theta], r Sin[phi] Sin[theta], 
          r Cos[theta]} + {x0, y0, z0}],
       -1
       ]
      (* TODO Implement Beta Sphere *)
      , {r, r0, r1}
      , MaxIterations -> OptionValue[MaxIterations] + 1
      , StepMonitor :> (ctr=ctr+1; If[ctr==OptionValue[MaxIterations] || ((r-rprev) < OptionValue[StoppingCriterion]), Throw[r]]; rprev=r; ) ]
  ]];
  {xr, yr, zr} = {r Cos[phi] Sin[theta], r Sin[phi] Sin[theta],  r Cos[theta]} + {x0, y0, z0} /. {r -> rmax};
  {
   {xr, yr, zr},
   rmax,
   If[rmax < Infinity,
    {
     AssociatedNuclearAttractor[w, ncps,
      {xr, yr, zr} + 
       First[-0.01*Eigenvectors[H[w, {xr, yr, zr}], -1]]
      , BetaSphereRadii -> OptionValue[BetaSphereRadii]],
     AssociatedNuclearAttractor[w, ncps,
      {xr, yr, zr} + 
       First[+0.01*Eigenvectors[H[w, {xr, yr, zr}], -1]]
      , BetaSphereRadii -> OptionValue[BetaSphereRadii]]
     },
    {ncp, Infinity}
    ]
   }
  ];

CriticalPointGradientMinus3[{g_?(VectorQ[#,NumericQ]&),H_}]:=Module[
{i,j,b,U,F,A,eval,lambda,h,denom,tiny},
tiny=1.*^-12;
(* Method of Cerjan, Miller, and Baker *)
{b,U}=Transpose[Sort[Transpose[Eigensystem[H]]]];
F=U . g;
(* -3 *)
A={
{  b[[1]]   ,     0       ,  0        , F[[1]] },
{      0    ,     b[[2]]  ,  0        , F[[2]] },
{      0    ,     0       ,  b[[3]]   , F[[3]] },
{  F[[1]]   ,     F[[2]]  ,  F[[3]]   ,    0   }
};
eval = First@Eigenvalues[A,1];
lambda={
eval,
eval,
eval
};
denom=(b-lambda);
Do[
	If[ Abs[denom[[i]]]<tiny,
 	denom[[i]]=denom[[i]]+Sign[denom[[i]]]*tiny
	];
 ,{i,1,3} 
 ];
h=Table[Sum[-F[[i]]*U[[i,j]]/denom[[i]],{i,1,3}],{j,1,3}];
h
];

CriticalPointGradientMinus1[{g_?(VectorQ[#,NumericQ]&),H_}]:=Module[
{i,j,b,U,F,A,eval,lambda,h,denom,tiny},
tiny=1.*^-12;
(* Method of Cerjan, Miller, and Baker *)
{b,U}=Transpose[Sort[Transpose[Eigensystem[H]]]];
F=U . g;
(* -1 *)
A={
{    b[[1]] ,       0   ,  F[[1]] },
{      0    ,     b[[2]],  F[[2]] },
{    F[[1]] ,     F[[2]],  0      }
};
eval = First@Eigenvalues[A,1];
lambda={ 
eval,   
eval, 
(1/2)(b[[3]]  - Sqrt[b[[3]]^2+ 4 F[[3]]^2])
};
denom=(b-lambda);
Do[
	If[ Abs[denom[[i]]]<tiny,
 	denom[[i]]=denom[[i]]+Sign[denom[[i]]]*tiny
	];
 ,{i,1,3} 
 ];
h=Table[Sum[-F[[i]]*U[[i,j]]/denom[[i]],{i,1,3}],{j,1,3}];
h
];

CriticalPointGradientPlus1[{g_?(VectorQ[#,NumericQ]&),H_}] := Module[
{i, j, b, U, F, A, eval, lambda, h, denom, tiny},
tiny=1.*^-12;
(* Method of Cerjan, Miller, and Baker *)
{b, U} = Transpose[Sort[Transpose[Eigensystem[H]]]];
F = U . g;
(* +1 *)
A ={
  {  b[[2]] ,  0      ,  F[[2]]   },
  {    0    ,  b[[3]] ,  F[[3]]   },
  {  F[[2]] ,  F[[3]] ,  0        }
  };
eval = First@Eigenvalues[A,-1];
lambda ={  
  (1/2) (b[[1]]  + Sqrt[b[[1]]^2 + 4 F[[1]]^2]), 
  eval, 
  eval   
  };
 denom=(b-lambda);
 Do[
	If[ Abs[denom[[i]]]<tiny,
 	denom[[i]]=denom[[i]]+Sign[denom[[i]]]*tiny
	];
 ,{i,1,3} 
 ];
h=Table[Sum[-F[[i]]*U[[i,j]]/denom[[i]],{i,1,3}],{j,1,3}];
h
];

CriticalPointGradientPlus3[{g_?(VectorQ[#,NumericQ]&),H_}]:=Module[
{i,j,b,U,F,A,eval,lambda,h,denom,tiny},
tiny=1.*^-12;
(* Method of Cerjan, Miller and Baker *)
{b,U}=Transpose[Sort[Transpose[Eigensystem[H]]]];
F=U . g;
(* +3 *)
A={
{  b[[1]] ,     0     ,  0      , F[[1]] },
{      0  ,     b[[2]],  0      , F[[2]] },
{      0  ,     0     ,  b[[3]] , F[[3]] },
{  F[[1]] ,     F[[2]],  F[[3]] ,    0   }
};
eval=First@Eigenvalues[A,-1];
lambda={ 
eval,   
eval, 
eval
};
denom=(b-lambda);
Do[
	If[ Abs[denom[[i]]]<tiny,
 	denom[[i]]=denom[[i]]+Sign[denom[[i]]]*tiny
	];
 ,{i,1,3} 
 ];
h=Table[Sum[-F[[i]]*U[[i,j]]/denom[[i]],{i,1,3}],{j,1,3}];
h
];

ParallelNInt[f_, {xmin_?NumericQ, xmax_?NumericQ}, maxr_?NumericQ] := Module[
   {memoize, mtemp, res, fake, mem, xnew, known},
   
   (* keep a list of known x values *)
   known = {};
   
   (* ensure numeric evaluation inside NIntegrate *)
   
   mem[x_?NumericQ] := memoize[x];
   fake[x_?NumericQ] := 0;
   
   (* increment MaxRecursion from 0 to maxr *)
   Do[
    
    (* First Invocation: "trick" NIntegrate 
    so that it either looks up values from previous
    lower recursive levels, or quickly yields a fake value *)
    
    (*
    Print["Recursion r = " <> ToString[r]];
    Print["---------"];
    *)
    
    xnew = Quiet[
      Reap[
        NIntegrate[
         If[MemberQ[known, x], mem[x], fake[x]],
         {x, xmin, xmax},
         Method -> {"LocalAdaptive", "SymbolicProcessing" -> 0},
         MinRecursion -> r,
         MaxRecursion -> r,
         EvaluationMonitor :> (
           (*
           Print[{"1st: ",x,If[MemberQ[known,x],mem[x],
           fake[x]],MemberQ[known,x]}];
           *)
           Sow[x])
         ]
        ][[2, 1]]
      ];
    
    (* return the previous Second Invocation result if no new data *)

        If[Length[xnew] == 0, Return[res]];
    
    (* for values that we have not seen, evaluate in parallel *)
    
    xnew = Complement[xnew, known];
    mtemp = ParallelMap[
      f[#] &,
      xnew
      , Method -> "FinestGrained"];
    
    (* memoize those values that were computed in parallel *)
    Do[
     memoize[xnew[[i]]] = mtemp[[i]]
     , {i, Length[mtemp]}];
    
    (* append those x-values to the known list *)
    
    known = Join[known, xnew];
    
    (* Second Invocation: NIntegrate with memoized values *)
    
    res = Quiet[
      NIntegrate[
       mem[x],
       {x, xmin, xmax},
       Method -> {"LocalAdaptive", "SymbolicProcessing" -> 0},
       MinRecursion -> r,
       MaxRecursion -> r
       (*,EvaluationMonitor:>Print[{"2nd: ", x,mem[x]}]*)
       ]
      ];
    
    (* increment up to maxr *)
    , {r, 0, maxr}];
   res
   ];

ParallelNInt[f_, {xmin_?NumericQ, xmax_?NumericQ}, {ymin_?NumericQ, ymax_?NumericQ}, maxr_?NumericQ] := Module[
   {memoize, mtemp, res, fake, mem, xnew, known},
   
   (* keep a list of known x values *)
   known = {};
   
   (* ensure numeric evaluation inside NIntegrate *)
   
   mem[{x_?NumericQ, y_?NumericQ}] := memoize[{x,y}];
   fake[{x_?NumericQ, y_?NumericQ}] := 0;
   
   (* increment MaxRecursion from 0 to maxr *)
   Do[
    
    (* First Invocation: "trick" NIntegrate 
    so that it either looks up values from previous
    lower recursive levels, or quickly yields a fake value *)
    
    
    Print["Recursion r = " <> ToString[r]];
    Print["---------"];
    
    
    xnew = Quiet[
      Reap[
        NIntegrate[
         If[MemberQ[known, {x,y}], mem[{x,y}], fake[{x,y}]],
         {x, xmin, xmax},
         {y, ymin, ymax},
         Method -> {"LocalAdaptive", Method -> "LobattoKronrodRule", "SymbolicProcessing" -> 0},
         MinRecursion -> r,
         MaxRecursion -> r,
         EvaluationMonitor :> (
           (*
           Print[{"1st: ",{x,y},If[MemberQ[known,{x,y}],mem[{x,y}],
           fake[{x,y}]],MemberQ[known,{x,y}]}];
           *)
           Sow[{x,y}])
         ]
        ][[2, 1]]
      ];
    
    (* return the previous Second Invocation result if no new data *)

        If[Length[xnew] == 0, Return[res]];
    
    (* for values that we have not seen, evaluate in parallel *)
    
    xnew = Complement[xnew, known];
    mtemp = ParallelMap[
      f[#] &,
      xnew
      , Method -> "FinestGrained"];
    
    (* memoize those values that were computed in parallel *)
    Do[
     memoize[xnew[[i]]] = mtemp[[i]]
     , {i, Length[mtemp]}];
    
    (* append those x-values to the known list *)
    
    known = Join[known, xnew];
    
    (* Second Invocation: NIntegrate with memoized values *)
    
    res = Quiet[
      NIntegrate[
       mem[{x,y}],
       (*f[{x,y}],*)
       {x, xmin, xmax},
       {y, ymin, ymax},
       Method -> {"LocalAdaptive", Method -> "LobattoKronrodRule", "SymbolicProcessing" -> 0},
       MinRecursion -> r,
       MaxRecursion -> r
       (*,EvaluationMonitor:>Print[{"2nd: ", {x,y} ,mem[{x,y}]}]*)
       ]
      ];
    
    (* increment up to maxr *)
    , {r, 0, maxr}];
   res
   ];

ParallelNInt[f_, {xmin_?NumericQ, xmax_?NumericQ}, {ymin_?NumericQ, ymax_?NumericQ}, {zmin_?NumericQ, zmax_?NumericQ}, maxr_?NumericQ] := Module[
   {memoize, mtemp, res, fake, mem, xnew, known},
   
   (* keep a list of known x values *)
   known = {};
   
   (* ensure numeric evaluation inside NIntegrate *)
   
   mem[{x_?NumericQ, y_?NumericQ, z_?NumericQ}] := memoize[{x,y,z}];
   fake[{x_?NumericQ, y_?NumericQ, z_?NumericQ}] := 0;
   
   (* increment MaxRecursion from 0 to maxr *)
   Do[
    
    (* First Invocation: "trick" NIntegrate 
    so that it either looks up values from previous
    lower recursive levels, or quickly yields a fake value *)
    
    (*
    Print["Recursion r = " <> ToString[r]];
    Print["---------"];
    *)
    
    xnew = Quiet[
      Reap[
        NIntegrate[
         If[MemberQ[known, {x,y,z}], mem[{x,y,z}], fake[{x,y,z}]],
         {x, xmin, xmax},
         {y, ymin, ymax},
         {z, zmin, zmax},
         Method -> {"LocalAdaptive", Method -> "LobattoKronrodRule", "SymbolicProcessing" -> 0},
         MinRecursion -> r,
         MaxRecursion -> r,
         EvaluationMonitor :> (
           (*
           Print[{"1st: ",x,If[MemberQ[known,x],mem[x],
           fake[x]],MemberQ[known,x]}];
           *)
           Sow[{x,y,z}])
         ]
        ][[2, 1]]
      ];
    
    (* return the previous Second Invocation result if no new data *)

    If[Length[xnew] == 0, Return[res]];
    
    (* for values that we have not seen, evaluate in parallel *)
    
    xnew = Complement[xnew, known];
    mtemp = ParallelMap[
      f[#] &,
      xnew
      , Method -> "FinestGrained"];
    
    (* memoize those values that were computed in parallel *)
    Do[
     memoize[xnew[[i]]] = mtemp[[i]]
     , {i, Length[mtemp]}];
    
    (* append those x-values to the known list *)
    
    known = Join[known, xnew];
    
    (* Second Invocation: NIntegrate with memoized values *)
    
    res = Quiet[
      NIntegrate[
       mem[{x,y,z}],
       {x, xmin, xmax},
       {y, ymin, ymax},
       {z, zmin, zmax},
       Method -> {"LocalAdaptive", Method -> "LobattoKronrodRule", "SymbolicProcessing" -> 0},
       MinRecursion -> r,
       MaxRecursion -> r
       (*,EvaluationMonitor:>Print[{"2nd: ", x,mem[x]}]*)
       ]
      ];
    
    (* increment up to maxr *)
    , {r, 0, maxr}];
   res
   ];

End[]

EndPackage[]
