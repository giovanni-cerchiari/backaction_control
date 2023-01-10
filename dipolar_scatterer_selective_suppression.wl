(* ::Package:: *)

(*
Copyright: Giovanni Cerchiari, Lorenzo Dania, Yannick Weiser, Dmitry Bykov
e-mail: giovanni.cerchiari@uibk.ac.at
date : 08/2022
*)
(*This file is a free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This file is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the other repository files.
    If not, it can be found at <https://www.gnu.org/licenses/>.
*)
(*
-----------------------------------------------------------------------
This file contains the script-program that has been used to evaluate
the results presented in the article:

If you wish to cite this work, we prepared a citation file "selective_suppression.bib"
in bibtex format in the repository.

This file was modified from 
doi: 10.5281/zenodo.4545692
*)
(*
-----------------------------------------------------------------------
The script is written for Mathematica 11.3 .
It was executed sucessfully on a machine with the following specifications
- operating system : Windows 10
- processor : Intel(R) Core(TM) i5-7300U CPU @ 2.60GHz 2.71 GHz
- RAM : 16 GB
-----------------------------------------------------------------------
*)
(*---------------------------------------------------------------------------------------------*)
(*---------------------------------------------------------------------------------------------*)
(*---------------------------------------------------------------------------------------------*)
ClearAll["Global`*"]
Needs["PlotLegends`"]
SetCoordinates[Spherical]
(*Generic direction in space*)
nr[\[Theta]_,\[Phi]_]:={Cos[\[Phi]]*Sin[\[Theta]],Sin[\[Phi]]*Sin[\[Theta]],Cos[\[Theta]]};
Print["Generic direction in space (\[Theta],\[Phi]) = ", MatrixForm[nr[\[Theta],\[Phi]]]]
(*---------------------------------------------------------------------------------------------*)
(*---------------------------------------------------------------------------------------------*)
(*RADIATED POWER PER UNIT SOLID ANGLE*)
(*differential scattered power into the solid andgle d\[CapitalOmega] = Sin[\[Theta]]d\[Theta]d\[Phi] *)
(*eq. 6 of the main text*)
(*Power radiated by a dipole:
- Pdip integral power of the dipole
- \[Theta]0,\[Phi]0 direction of linear polarization
- \[Theta],\[Phi] direction scattered light*)
(*Please remember the Sin[\[Theta]] missing factor at integration*)
(*linear polarization general case*)
dpdipd\[CapitalOmega][\[Theta]_,\[Phi]_,\[Theta]0_,\[Phi]0_]:=3/(8*\[Pi])*Pdip*(1-(nr[\[Theta],\[Phi]].nr[\[Theta]0,\[Phi]0])^2);
Print["(dpdipd\[CapitalOmega]) differential power radiate by a dipole
 with linear polarization: ",dpdipd\[CapitalOmega][\[Theta],\[Phi],\[Theta]p,\[Phi]p]]


(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*BACK ACTION CENTER OF MASS WITH SELECTIVE SUPPRESSION*)
(*power spectral density of the radiated power along the direction defined by \[Theta] and \[Phi]*)
dspp[\[Theta]_,\[Phi]_,\[Theta]p_,\[Phi]p_]:=hbar*(c/\[Lambda])*dpdipd\[CapitalOmega][\[Theta],\[Phi],\[Theta]p,\[Phi]p];
(*power spectral density fluctuation of the radiation pressure force*)
dsba[\[Theta]_,\[Phi]_,\[Theta]p_,\[Phi]p_,\[Theta]0_,\[Phi]0_]:=(nr[\[Theta],\[Phi]].nr[\[Theta]0,\[Phi]0])^2*dspp[\[Theta],\[Phi],\[Theta]p,\[Phi]p]/c^2;
(*Back action. Integral of the force fluctuations.*)
(*We select the light polarization aligned to x-axis (\[Theta]p = \[Pi]/2, \[Phi]p = 0)*)
Sba[\[Theta]_,\[Theta]0_,\[Phi]0_]:=Evaluate[Simplify[Integrate[Simplify[Integrate[dsba[\[Theta]\[Theta],\[Phi],\[Pi]/2,0,\[Theta]0,\[Phi]0]*Sin[\[Theta]\[Theta]],{\[Phi],0,2*\[Pi]}]], {\[Theta]\[Theta],0,\[Theta]}]]];
Sbavet[\[Theta]_]:=Evaluate[Simplify[{Sba[\[Theta],\[Pi]/2,0], Sba[\[Theta],\[Pi]/2,\[Pi]/2], Sba[\[Theta],0,0]}]];
(*Ratio between controlled back action and maximum back action in free space.*)
(*The factor of 2 in front account for both sides of the uncontrolled region.*)
Sbar[\[Theta]_,\[Theta]0_,\[Phi]0_]:=Evaluate[Simplify[2*Sba[\[Theta], \[Theta]0, \[Phi]0]/Sba[\[Pi], 0, 0]]];
Sbarvet[\[Theta]_]:=Evaluate[Simplify[{Sbar[\[Theta],\[Pi]/2,0], Sbar[\[Theta],\[Pi]/2,\[Pi]/2], Sbar[\[Theta],0,0]}]];
(*Print*)
Print["BackAction(x, free space) = ", Sbavet[\[Theta]][[1]]]
Print["BackAction(y, free space) = ", Sbavet[\[Theta]][[2]]]
Print["BackAction(z, free space) = ", Sbavet[\[Theta]][[3]]]
Print["BackAction(x,controlled/free) = ", Sbarvet[\[Theta]][[1]]]
Print["BackAction(y,controlled/free) = ", Sbarvet[\[Theta]][[2]]]
Print["BackAction(z,controlled/free) = ", Sbarvet[\[Theta]][[3]]]
Print["BackAction(x,controlled/free) = ", Series[Sbarvet[\[Theta]][[1]],{\[Theta],0,4}]]
Print["BackAction(y,controlled/free) = ", Series[Sbarvet[\[Theta]][[2]],{\[Theta],0,4}]]
Print["BackAction(z,controlled/free) = ", Series[Sbarvet[\[Theta]][[3]],{\[Theta],0,4}]]


(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*MEASUREMENT IMPRECISION OF THE SELF-HOMODYNE METHOD. LINEAR POLARIZATION*)
(*----------------------------------------------------------------------------------------------*)
(*The label "mirror" in the names of the variables refers to the self-homodyne method*)
(*----------------------------------------------------------------------------------------------*)
(*The imprecision s and S are described in the article Eqs 10 and 11: 
s = differential power spectral density of the imprecision noise associated with a differential 
    detector under direction (\[Theta], \[Phi])
S = integral of s in the angular integral \[Phi]=[0, 2\[Pi]] and \[Theta]=[0,\[Theta]]. Note the free parameter \[Theta].
With SS and ss natural units are introduced*)
(*Assumptions: shot noise dominated by the constant term*)
(*-------------------*)
(*differential imprecision for the self-homodyne configuration*)
smirror[\[Theta]_,\[Phi]_,\[Theta]p_,\[Phi]p_,\[Theta]0_,\[Phi]0_]:= Evaluate[Simplify[
	(hbar*c/\[Lambda])*(1+R^2)/( (4*R*(2*\[Pi]/\[Lambda])*(nr[\[Theta]0,\[Phi]0].nr[\[Theta],\[Phi]]))^2 *dpdipd\[CapitalOmega][\[Theta],\[Phi],\[Theta]p,\[Phi]p]),
	Assumptions->{0<=\[Theta]<=\[Pi], 0<=\[Phi]<2*\[Pi]}]];
Print["smirror = ", smirror[\[Theta],\[Phi],\[Pi]/2, 0, \[Theta]0,\[Phi]0]]
(*Imprecision for the self-homodyne configuration*)
Smirror[\[Theta]_,\[Theta]0_,\[Phi]0_]:=Evaluate[Simplify[
	(Integrate[Sin[\[Theta]\[Theta]]/smirror[\[Theta]\[Theta],\[Phi],\[Pi]/2,0,\[Theta]0,\[Phi]0],{\[Theta]\[Theta],0,\[Theta]},{\[Phi],0,2*\[Pi]}])^-1,
	 Assumptions->{0<=\[Theta]<=\[Pi]/2}]];
Smirrorvet[\[Theta]_]:=Evaluate[Simplify[{Smirror[\[Theta],\[Pi]/2,0], Smirror[\[Theta],\[Pi]/2,\[Pi]/2], Smirror[\[Theta],0,0]},{Assumptions->{R==1}}]];
Print["Imprecision self-homodyne(x) = ", Smirrorvet[\[Theta]][[1]]]
Print["Imprecision self-homodyne(y) = ", Smirrorvet[\[Theta]][[2]]]
Print["Imprecision self-homodyne(z) = ", Smirrorvet[\[Theta]][[3]]]
(*Enforce reference system with light polarization on the x-axis (\[Theta]p = \[Pi]/2, \[Phi]p = 0) and introducing natural units.*)
SSmirror[\[Theta]_,\[Theta]0_,\[Phi]0_]:= Evaluate[Simplify[Smirror[\[Theta],\[Pi]/2,0,\[Theta]0,\[Phi]0],
	Assumptions->{0<=\[Theta]<=\[Pi]/2, 0<=\[Phi]<2*\[Pi], 0<=\[Theta]0<=\[Pi], 0<=\[Phi]0<2*\[Pi], hbar==1, \[Lambda]==1, c==1, R==1, Pdip==1}]];
ssmirror[\[Theta]_,\[Phi]_,\[Theta]0_,\[Phi]0_]:= Evaluate[FullSimplify[smirror[\[Theta],\[Phi],\[Pi]/2,0,\[Theta]0,\[Phi]0],
	Assumptions->{0<=\[Theta]<=\[Pi]/2, 0<=\[Phi]<2*\[Pi], 0<=\[Theta]0<=\[Pi], 0<=\[Phi]0<2*\[Pi], hbar==1, \[Lambda]==1, c==1, R==1, Pdip==1}]];


(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*PRODUCT OF IMPRECISION AND BACK ACTION FOR HYBRID SYSTEM COMPOSED BY THE HOLLOW MIRROR AND SELF-HOMODYNE*)
Print["detection limit x = ", Simplify[Smirrorvet[\[Theta]][[1]]*Sbavet[\[Theta]][[1]]]]
Print["detection limit y = ", Simplify[Smirrorvet[\[Theta]][[2]]*Sbavet[\[Theta]][[2]]]]
Print["detection limit z = ", Simplify[Smirrorvet[\[Theta]][[3]]*Sbavet[\[Theta]][[3]]]]


Subscript[r, 1][\[Alpha]_,d_]:={d Cos[\[Alpha]],0,d Sin[\[Alpha]]};
Subscript[r, 2][\[Alpha]_,d_]:=-Subscript[r, 1][\[Alpha],d];
Subscript[EF, sc,1][\[Theta]_,\[Phi]_,\[Alpha]_,d_]:=Exp[(I (2 \[Pi]) nr[\[Theta],\[Phi]].Subscript[r, 1][\[Alpha],d])/\[Lambda]];
Subscript[EF, sc,2][\[Theta]_,\[Phi]_,\[Alpha]_,d_]:=Exp[(I (2 \[Pi]) nr[\[Theta],\[Phi]].Subscript[r, 2][\[Alpha],d])/\[Lambda]];
Subscript[EF, i,1][\[Theta]_,\[Phi]_,\[Alpha]_,d_]:=-\[Rho] Conjugate[Subscript[EF, sc,1][\[Theta],\[Phi],\[Alpha],d]] Exp[-((I (2 \[Pi]) 2 R)/\[Lambda])];
Subscript[EF, i,2][\[Theta]_,\[Phi]_,\[Alpha]_,d_]:=-\[Rho] Conjugate[Subscript[EF, sc,1][\[Theta],\[Phi],\[Alpha],d]] Exp[-((I (2 \[Pi]) 2 R)/\[Lambda])];
IN[\[Theta]_,\[Phi]_,\[Alpha]_,d_]:=Abs[Subscript[EF, sc,1][\[Theta],\[Phi],\[Alpha],d]+Subscript[EF, sc,2][\[Theta],\[Phi],\[Alpha],d]+Subscript[EF, i,1][\[Theta],\[Phi],\[Alpha],d]+Subscript[EF, i,2][\[Theta],\[Phi],\[Alpha],d]]^2;
assumptions={d==1,d\[Element]Reals,0<=\[Alpha]<=2 \[Pi],\[Alpha]\[Element]Reals,Subscript[EF, o,1]\[Element]Vectors[3,Complexes],Subscript[EF, o,2]\[Element]Vectors[3,Complexes],\[Lambda]>0,\[Lambda]\[Element]Reals,0<=\[Theta]<=2 \[Pi],\[Theta]\[Element]Reals,0<=\[Phi]<=2 \[Pi],\[Phi]\[Element]Reals,R==20 d,R\[Element]Reals,\[Rho]>0,\[Rho]\[Element]Reals,Abs[Subscript[EF, 0,1]]==1,Abs[Subscript[EF, 0,2]]==1};
Print[Evaluate[FullSimplify[IN[\[Theta],\[Phi],\[Alpha],d],Assumptions->assumptions]]]


(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*GENERAL FORMULAS FOR SELF-HOMDYNE DETECTION*)
Print[""]
(*Maximum and minimum imprecision*)
minS = FullSimplify[SSmirror[ArcSin[1],\[Theta]pmax,\[Theta]0max,\[Phi]0max]];
maxS = FullSimplify[SSmirror[ArcSin[1],\[Theta]pmax,\[Pi]/2,0]];
Print["Minimum imprecision (natural units): ", minS];
Print["Minimum imprecision (natural units): ", maxS];
Print[""]
Print["differential imprecision"]
Print["smirror = ", smirror[\[Theta],\[Phi],\[Pi]/2,0,\[Theta]0,\[Phi]0]]
Print[""]
Print["Total imprecision \[Theta]0=0"]
Print["Smirror(\[Theta]0=0)= ", FullSimplify[Smirror[\[Theta]D,\[Pi]/2,0,0,0]]]
Print[""]
Print["Total imprecision \[Theta]0=\[Pi]/2"]
Print["Smirror(\[Theta]0=\[Pi]/2)= ", FullSimplify[Smirror[\[Theta]D,\[Pi]/2,0,\[Pi]/2,\[Phi]0]]]
Print[""]
Print["Total imprecision \[Theta]D=\[Pi]/2"]
Print["Smirror(\[Theta]=\[Pi]/2)= ", FullSimplify[Smirror[\[Pi]/2,\[Pi]/2,0,\[Theta]0,\[Phi]0]]]
Print[""]
Print["SmirrorQPD =", SmirrorQPD[\[Theta],\[Theta]0,\[Phi]0]]


Needs["PlotLegends`"]
(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*PLOTS*)
(*----------------------------------------------------------------------------------------------*)
(*standard font size for all plots*)
lgdfontsize = 16;
(*frame flag that can be used for all plots*)
frameflg = False;
(*----------------------------------------------------------------------------------------------*)
(*3D axes label*)
xyzlabl = {Style["x",Bold,Black,lgdfontsize],Style["y",Bold,Black,lgdfontsize],Style["z",Bold,Black,lgdfontsize]};
(*2D axes label*)
xylabl = {Style["NA",Bold,Black,lgdfontsize],Style["Sba",Bold,Black,lgdfontsize]};
(*plot legend spec*)
pltlgdxyz100={Style["suppressed x",Black,lgdfontsize],Style["suppressed y",Black,lgdfontsize],Style["suppressed z",Black,lgdfontsize], Style["free space y-z",Black,lgdfontsize], Style["free space x",Black,lgdfontsize]};
(*line style spec*)
xyzlinestyle = {Directive[Red, Thick, Dashing[None]], Directive[Green, Thick, Dashing[None]],Directive[Blue, Thick, Dashing[None]],
Directive[Black, Thick, Dashed], Directive[Black, Thick, Dashed],Directive[Blue, Thick, Dashed],Directive[Orange, Thick, Dashing[None]], Directive[Purple, Thick, Dashing[None]]};
(*plotting*)
Plot[{Sbarvet[ArcSin[x]][[1]],Sbarvet[ArcSin[x]][[2]],Sbarvet[ArcSin[x]][[3]], Sba[\[Pi]/2,0,0]/Sba[\[Pi]/2,0,0], Sba[\[Pi]/2,\[Pi]/2,0]/Sba[\[Pi]/2,0,0]},{x,0,1},
	PlotLegend->pltlgdxyz100,
	LegendPosition->{-0.4,-0.2},LegendShadow->None,LegendSize->1,LegendBorder->Black,
	ImageSize->Large, AxesLabel -> xylabl , Axes -> True, Frame->frameflg, PlotStyle->xyzlinestyle, BaseStyle->{FontSize->16},
	 PlotTheme->"Monochrome", AspectRatio -> 3/4, Background->White]


(* ::InheritFromParent:: *)
(**)


Subscript[n, theta][\[Theta]_,\[Phi]_]:={Cos[\[Theta]] Cos[\[Phi]],Cos[\[Theta]] Sin[\[Phi]],-Sin[\[Theta]]};
Subscript[n, phi][\[Phi]_]:={-Sin[\[Phi]],Cos[\[Phi]],0};
p[\[CapitalDelta]\[Delta]_,\[CapitalDelta]\[Epsilon]_]:=Subscript[EF, 0] Subscript[\[Alpha], 0] {1,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon]};
Subscript[u, i][\[Theta]_,\[Phi]_,i_]:=Sqrt[3/(8 \[Pi])] (Subscript[n, theta][\[Theta],\[Phi]] Subscript[n, theta][\[Theta],\[Phi]][[i]]+Subscript[n, phi][\[Phi]] Subscript[n, phi][\[Phi]][[i]]);
U[\[Theta]_,\[Phi]_]:={Subscript[u, i][\[Theta],\[Phi],1],Subscript[u, i][\[Theta],\[Phi],2],Subscript[u, i][\[Theta],\[Phi],3]};
u[i_,j_]:=\!\(
\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(2\ \[Pi]\)]\(\((
\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Pi]\)]\(Sin[\[Theta]]\ \(U[\[Theta], \[Phi]]\)[[i]] . \(U[\[Theta], \[Phi]]\)[[j]]\) \[DifferentialD]\[Theta])\) \[DifferentialD]\[Phi]\)\);
G[\[Theta]_,\[Phi]_]:=Sqrt[(8 \[Pi])/3] U[\[Theta],\[Phi]];
Subscript[EF, sc][r_,\[Theta]_,\[Phi]_]:=(\!\(
\*SubsuperscriptBox[\(\[Omega]\), \(0\), \(2\)]\ 
\*SuperscriptBox[\(E\), \(I\ k\ r\)]\ \(G[\[Theta], \[Phi]] . p[\[CapitalDelta]\[Delta], \[CapitalDelta]\[Epsilon]]\)\))/((4 \[Pi] r) (c^2 Subscript[\[Epsilon], 0]));
Subscript[EF, lo][r_,\[Theta]_,\[Phi]_]:=-((\[Rho] Subscript[EF, sc][r,\[Theta],\[Phi]])/E^(I 2 k Subscript[R, s]));
Subscript[IN, tot][r_,\[Theta]_,\[Phi]_]:=FullSimplify[1/2 c Subscript[\[Epsilon], 0] Norm[Subscript[EF, lo][r,\[Theta],\[Phi]]+Subscript[EF, sc][r,\[Theta],\[Phi]]]^2,Assumptions->{0<=\[Theta]<=\[Pi]/2,0<=\[Phi]<=\[Pi],Subscript[R, s]>0,Subscript[R, s]\[Element]Reals,k\[Element]Reals,r>0,r\[Element]Reals,r>0,\[CapitalDelta]\[Delta]\[Element]Reals,\[CapitalDelta]\[Epsilon]\[Element]Reals,\[Rho]>0,\[Rho]\[Element]Reals,c\[Element]Reals,c>0,Subscript[\[Epsilon], 0]>0,Subscript[\[Omega], 0]>0,Subscript[\[Omega], 0]\[Element]Reals,Subscript[EF, 0]\[Element]Reals,Subscript[EF, 0]>0,Subscript[\[Alpha], 0]\[Element]Reals,Subscript[\[Alpha], 0]>0}];
Subscript[P, det][\[CapitalTheta]_]:=FullSimplify[\!\(
\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(2\ \[Pi]\)]\(\((
\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[CapitalTheta]\)]
\(\*SubscriptBox[\(IN\), \(tot\)]\)[r, \[Theta], \[Phi]] \[DifferentialD]\[Theta])\) \[DifferentialD]\[Phi]\)\),Assumptions->{\[Rho]==1,k Subscript[R, s]==\[Pi]/2}];


Print["Intensity: ",Subscript[IN, tot][r,\[Theta],\[Phi]]];
Print["Power: ",Subscript[P, det][\[CapitalTheta]]];
preFactor=(Subscript[EF, 0] Subscript[\[Alpha], 0] \!\(\*SubsuperscriptBox[\(\[Omega]\), \(0\), \(2\)]\))/(4 c^(3/2) Sqrt[\[Pi] Subscript[\[Epsilon], 0]] r);
Print["Normalised Power: ",Simplify[Subscript[P, det][\[CapitalTheta]]/preFactor^2]];
Subscript[P, \[CapitalDelta]\[Epsilon]][\[CapitalTheta]_]:=4 \[CapitalDelta]\[Epsilon]^2 \[CapitalTheta]-2 \[CapitalDelta]\[Epsilon]^2 Sin[2 \[CapitalTheta]];
Subscript[P, \[CapitalDelta]\[Delta]][\[CapitalTheta]_]:=6 \[CapitalDelta]\[Delta]^2 \[CapitalTheta]+\[CapitalDelta]\[Delta]^2 Sin[2 \[CapitalTheta]];
Plot[{Subscript[P, \[CapitalDelta]\[Epsilon]][\[CapitalTheta]]/. {\[CapitalDelta]\[Epsilon]->1},Subscript[P, \[CapitalDelta]\[Delta]][\[CapitalTheta]]/. {\[CapitalDelta]\[Delta]->1}},{\[CapitalTheta],0,\[Pi]/2},AxesLabel->{\[CapitalTheta],Subscript[P, \[CapitalDelta]]},PlotLegends->{Subscript[P, \[CapitalDelta]\[Epsilon]][\[CapitalTheta]],Subscript[P, \[CapitalDelta]\[Delta]][\[CapitalTheta]]}]
