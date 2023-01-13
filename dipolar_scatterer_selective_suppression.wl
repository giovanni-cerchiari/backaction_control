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


nr[\[Theta]_,\[Phi]_]:={Cos[\[Phi]] Sin[\[Theta]],Sin[\[Phi]] Sin[\[Theta]],Cos[\[Theta]]};
f[\[Theta]_,\[Phi]_,r_,a_,b_,\[Theta]1_,\[Phi]1_,r1_,\[Theta]2_,\[Phi]2_,r2_]:=a DiracDelta[\[Theta]-\[Theta]1] DiracDelta[\[Phi]-\[Phi]1] DiracDelta[r-r1]+b DiracDelta[\[Theta]-\[Theta]2] DiracDelta[\[Phi]-\[Phi]2] DiracDelta[r-r2];
g[\[Theta]_,\[Phi]_,\[Theta]0_,\[Phi]0_,r0_,\[Theta]k_,\[Phi]k_]:=Exp[(I r0 (2 \[Pi]) nr[\[Theta]k,\[Phi]k].nr[\[Theta]0,\[Phi]0])/\[Lambda]] Exp[(I r0 (2 \[Pi]) nr[\[Theta],\[Phi]].nr[\[Theta]0,\[Phi]0])/\[Lambda]];
gi[\[Theta]_,\[Phi]_,\[Theta]0_,\[Phi]0_,r0_,\[Theta]k_,\[Phi]k_,\[Rho]_,R_]:=-\[Rho] Exp[(I r0 (2 \[Pi]) nr[\[Theta]k,\[Phi]k].nr[\[Theta]0,\[Phi]0])/\[Lambda]] Exp[-((I r0 (2 \[Pi]) nr[\[Theta],\[Phi]].nr[\[Theta]0,\[Phi]0])/\[Lambda])-(I 2 (2 \[Pi]) R)/\[Lambda]];
assump={k\[Element]Reals,k>0,0<\[Theta]k<\[Pi],0<\[Phi]k<2 \[Pi],R\[Element]Reals,R>0,\[Rho]\[Element]Reals,\[Rho]>0,\[Lambda]\[Element]Reals,\[Lambda]>0,x0\[Element]Reals,x0>0,y\[Element]Reals,z\[Element]Reals,\[CapitalDelta]\[Element]Reals,x0>\[CapitalDelta]>0,0<=\[Theta]<\[Pi],0<\[Phi]<=2 \[Pi],r0>0,n>0,n\[Element]Integers};
gb[\[Theta]_,\[Phi]_,x0_,\[CapitalDelta]_,R_]:=Simplify[g[\[Theta],\[Phi],\[Pi]/2,0,x0+\[CapitalDelta],\[Pi]/2,0] (x0+\[CapitalDelta])^2+g[\[Theta],\[Phi],\[Pi]/2,\[Pi],x0-\[CapitalDelta],\[Pi]/2,0] (x0-\[CapitalDelta])^2+gi[\[Theta],\[Phi],\[Pi]/2,\[Pi],x0+\[CapitalDelta],\[Pi]/2,0,1,R] (x0+\[CapitalDelta])^2+gi[\[Theta],\[Phi],\[Pi]/2,0,x0-\[CapitalDelta],\[Pi]/2,0,1,R] (x0-\[CapitalDelta])^2];
gm=Simplify[Normal[Series[gb[\[Theta],\[Phi],x0,\[CapitalDelta],R],{x0,0,2}]],Assumptions->assump];

Ef[\[Theta]_,\[Phi]_,x0_,\[CapitalDelta]_,R_]:=gb[\[Theta],\[Phi],x0,\[CapitalDelta],R];
(*Ef[\[Theta]_,\[Phi]_,x0_,\[CapitalDelta]_,R_]:=gm*)
Print["Electric field E(\[Theta],\[Phi],\!\(\*SubscriptBox[\(x\), \(0\)]\),\[CapitalDelta]) \[Proportional] ",Ef[\[Theta],\[Phi],x0,\[CapitalDelta],R]];
Int[\[Theta]_,\[Phi]_,x0_,\[CapitalDelta]_,R_]:=FullSimplify[ComplexExpand[(Conjugate[Ef[\[Theta],\[Phi],x0,\[CapitalDelta],R]])*(Ef[\[Theta],\[Phi],x0,\[CapitalDelta],R])],Assumptions->assump];
Print["Intensity I(\[Theta],\[Phi],\!\(\*SubscriptBox[\(x\), \(0\)]\),\[CapitalDelta]) \[Proportional] ",Int[\[Theta],\[Phi],\[Lambda]/2,\[CapitalDelta],21/20 \[Lambda]]];

(*
Int2nd[\[Theta]_,\[Phi]_,x0_,\[CapitalDelta]_,R_]:=FullSimplify[Normal[Series[Int[\[Theta],\[Phi],x0,\[CapitalDelta],R],{\[CapitalDelta],0,2}]],Assumptions\[Rule]assump];
Print["Intensity up to and including second order in \[CapitalDelta]: I(\[Theta],\[Phi],Subscript[x, 0],\[CapitalDelta],R) \[Proportional] ", Int2nd[\[Theta],\[Phi],x0,\[CapitalDelta],R]]

P[\[CapitalDelta]_]:=Integrate[Integrate[(Int[\[Theta],\[Phi],\[Lambda]/2,\[CapitalDelta],21/20 \[Lambda]]/.{\[Lambda]\[Rule]1})*Sin[\[Theta]],{\[Theta],0,\[Pi]}],{\[Phi],0,\[Pi]}];
Plot[P[\[CapitalDelta]],{\[CapitalDelta],1/8,1/2}]*)

Print["Intensity in first order of \[CapitalDelta]: ", Series[Int[\[Theta],\[Phi],\[Lambda]/2,\[CapitalDelta],21/20 \[Lambda]],{\[CapitalDelta],0,1}]/.{\[Lambda]->1}]

Pphi[\[CapitalDelta]_,\[Phi]_]:=Simplify[Evaluate[Integrate[Normal[Series[Int[\[Theta],\[Phi],\[Lambda]/2,\[CapitalDelta],21/20 \[Lambda]]/.{\[Lambda]->1},{\[CapitalDelta],0,1}]],{\[Theta],0,\[Pi]}]],Assumptions->assump];

Plot[NIntegrate[Pphi[\[CapitalDelta],\[Phi]],{\[Phi],0,\[Pi]}],{\[CapitalDelta],0,1/2}]




