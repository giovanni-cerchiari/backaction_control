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


(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*IDEAL MEASUREMENT ROTATION*)



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
(*/.{\[CapitalDelta]\[Delta]->0,\[CapitalDelta]\[Epsilon]->0}*)


(* --------------------------------------------*)
(* ORIENTATION OF A DIPOLE
The calculation is based on the technique described in
----
F. Tebbenjohanns, A. Militaru, A. Norrman, F. van der Laan, L. Novotny and M. Frimmer, "Optimal orientation detection of an anisotropic dipolar scatterer",
Phys. Rev. A 105, 053504 (2022), doi:10.1103/PhysRevA.100.043821
----
*)
dipole[\[Alpha]0_,E0_,\[CapitalDelta]\[Delta]_,\[CapitalDelta]\[Epsilon]_]:= \[Alpha]0*E0*{1,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon]};
ux[\[Theta]_,\[Phi]_]:=Sqrt[3/(8*\[Pi])]*{nr[\[Theta],\[Phi]][[1]],0,0};
uy[\[Theta]_,\[Phi]_]:=Sqrt[3/(8*\[Pi])]*{0, nr[\[Theta],\[Phi]][[2]], 0};
uz[\[Theta]_,\[Phi]_]:=Sqrt[3/(8*\[Pi])]*{0,0,nr[\[Theta],\[Phi]][[3]]};
(*checking orthogonality condition*)
uiuj[ui_,uj_,\[Theta]_,\[Phi]_]:= Integrate[Integrate[ui[\[Theta],\[Phi]].uj[\[Theta],\[Phi]]*Sin[\[Theta]],{\[Theta],0,\[Pi]}],{\[Phi],0,2*\[Pi]}];
uxuy = Integrate[Integrate[ux[\[Theta],\[Phi]].uy[\[Theta],\[Phi]]*Sin[\[Theta]],{\[Theta],0,\[Pi]}],{\[Phi],0,2*\[Pi]}];
uxuz = Integrate[Integrate[ux[\[Theta],\[Phi]].uz[\[Theta],\[Phi]]*Sin[\[Theta]],{\[Theta],0,\[Pi]}],{\[Phi],0,2*\[Pi]}];
uyuy = Integrate[Integrate[uy[\[Theta],\[Phi]].uy[\[Theta],\[Phi]]*Sin[\[Theta]],{\[Theta],0,\[Pi]}],{\[Phi],0,2*\[Pi]}];
uyuz = Integrate[Integrate[uy[\[Theta],\[Phi]].uz[\[Theta],\[Phi]]*Sin[\[Theta]],{\[Theta],0,\[Pi]}],{\[Phi],0,2*\[Pi]}];
uzuz = Integrate[Integrate[uz[\[Theta],\[Phi]].uz[\[Theta],\[Phi]]*Sin[\[Theta]],{\[Theta],0,\[Pi]}],{\[Phi],0,2*\[Pi]}];
Print["check ortoghonality conditions"]
Print["ux.ux = ", uiuj[ux,ux,\[Theta],\[Phi]], "   ux.uy = ", uiuj[ux,uy,\[Theta],\[Phi]],  "   ux.uz = ", uiuj[ux,uz,\[Theta],\[Phi]]]
Print["uy.uy = ", uiuj[uy,uy,\[Theta],\[Phi]], "   uy.uz = ", uiuj[uy,uz,\[Theta],\[Phi]]]
Print["uz.uz = ", uiuj[uz,uz,\[Theta],\[Phi]]]
(*Green's tensor*)
Gt[\[Theta]_,\[Phi]_]:=Sqrt[8*\[Pi]/3]*{ux[\[Theta],\[Phi]], uy[\[Theta],\[Phi]],uz[\[Theta],\[Phi]]};
Print["Green's tensor = ", MatrixForm[Gt[\[Theta],\[Phi]]]]
(*Scattered electric field*)
Esc[\[Alpha]0_,E0_,\[CapitalDelta]\[Delta]_,\[CapitalDelta]\[Epsilon]_,\[Theta]_,\[Phi]_]:=(\[Omega]0^2/(4*\[Pi]*\[Epsilon]0*c^2))*Exp[((I*2*\[Pi])/\[Lambda])*r]*Gt[\[Theta],\[Phi]].dipole[\[Alpha]0,E0,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon]];
(*Scattered intensity*)
Intensity[\[Theta]D_]:=Evaluate[Integrate[Integrate[(\[Epsilon]0*c/2)*Conjugate[Esc[\[Alpha]0,E0,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon],\[Theta],\[Phi]]].Esc[\[Alpha]0,E0,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon],\[Theta],\[Phi]]*Sin[\[Theta]],{\[Phi],0,2*\[Pi]}],{\[Theta],0,\[Theta]D}, Assumptions->{\[Theta]D>0}]];
(*Integrated scattered power. The integral is doubled to consider the regions 0<\[Theta]<\[Theta]D and \[Pi]-\[Theta]D<\[Theta]<\[Pi] simultanously. Allowed values of \[Theta]D are in the interval [0,\[Pi]/2]*)
P0[\[Theta]D_]:= Simplify[2*Intensity[\[Theta]D],Assumptions->{c>0, \[Epsilon]0>0, \[Alpha]0 > 0, \[Omega]0 > 0, \[CapitalDelta]\[Delta] > 0, \[CapitalDelta]\[Epsilon] > 0, E0 > 0, \[Lambda]>0, r>0}];
Print["Total power free space (\[CapitalDelta]\[Delta]\[Rule]0,\[CapitalDelta]\[Epsilon]\[Rule]0). P0[\[Pi]/2] = ",P0[\[Pi]/2]/.{\[CapitalDelta]\[Delta]->0,\[CapitalDelta]\[Epsilon]->0}]
Print["Total power in presence of the hemipherical mirror (\[CapitalDelta]\[Delta]\[Rule]0,\[CapitalDelta]\[Epsilon]\[Rule]0). P0[\[Theta]D] = ",P0[\[Theta]D]/.{\[CapitalDelta]\[Delta]->0,\[CapitalDelta]\[Epsilon]->0}]
Simplify[Esc[\[Alpha]0,E0,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon],\[Theta],\[Phi]]/Sqrt[P0[\[Theta]D]/.{\[CapitalDelta]\[Delta]->0,\[CapitalDelta]\[Epsilon]->0}],Assumptions->{c>0, \[Epsilon]0>0, \[Alpha]0 > 0, \[Omega]0 > 0, \[CapitalDelta]\[Delta] > 0, \[CapitalDelta]\[Epsilon] > 0, E0 > 0, \[Lambda]>0, r>0}]







