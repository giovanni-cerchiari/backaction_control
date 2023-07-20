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

(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*PLOTS*)
(*----------------------------------------------------------------------------------------------*)
(*standard font size for all plots*)
lgdfontsize = 16;
(*frame flag that can be used for all plots*)
frameflg = False;
(*line style spec*)
xyzlinestyle = {Directive[Red, Thick, Dashing[None]], Directive[Green, Thick, Dashing[None]],Directive[Blue, Thick, Dashing[None]],
Directive[Black, Thick, Dashed], Directive[Black, Thick, Dashed],Directive[Blue, Thick, Dashed],
Directive[Orange, Thick, Dashing[None]], Directive[Purple, Thick, Dashing[None]]};
(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*SPHERICAL COORDINATES*)
SetCoordinates[Spherical]
(*Generic direction in space*)
nr[\[Theta]_,\[Phi]_]:={Cos[\[Phi]]*Sin[\[Theta]],Sin[\[Phi]]*Sin[\[Theta]],Cos[\[Theta]]};
n\[Theta][\[Theta]_,\[Phi]_]:={Cos[\[Phi]]*Cos[\[Theta]],Sin[\[Phi]]*Cos[\[Theta]],-Sin[\[Theta]]};
n\[Phi][\[Theta]_,\[Phi]_]:={-Sin[\[Phi]],Cos[\[Phi]],0};
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


(*---------------------------------------------------------------------------------------------*)
(*---------------------------------------------------------------------------------------------*)
(* ORIENTATION OF A DIPOLE
The calculation is based on the technique described in
----
F. Tebbenjohanns, A. Militaru, A. Norrman, F. van der Laan, L. Novotny and M. Frimmer, "Optimal orientation detection of an anisotropic dipolar scatterer",
Phys. Rev. A 105, 053504 (2022), doi:10.1103/PhysRevA.100.043821
----
*)
dipole[\[Alpha]0_,E0_,\[CapitalDelta]\[Delta]_,\[CapitalDelta]\[Epsilon]_]:= \[Alpha]0*E0*{1,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon]};
ui[\[Theta]_,\[Phi]_,n_]:=Sqrt[3/(8*\[Pi])]*(n\[Theta][\[Theta],\[Phi]]*n\[Theta][\[Theta],\[Phi]][[n]]+n\[Phi][\[Theta],\[Phi]]*n\[Phi][\[Theta],\[Phi]][[n]]);
(*checking orthogonality condition*)
uiuj[i_,j_]:= Integrate[Integrate[ui[\[Theta],\[Phi],i].ui[\[Theta],\[Phi],j]*Sin[\[Theta]],{\[Theta],0,\[Pi]}],{\[Phi],0,2*\[Pi]}];
Print["check ortoghonality conditions"]
Print["ux.ux = ", uiuj[1,1], "   ux.uy = ", uiuj[1,2],  "   ux.uz = ", uiuj[1,3]]
Print["uy.uy = ", uiuj[2,2], "   uy.uz = ", uiuj[2,3]]
Print["uz.uz = ", uiuj[3,3]]
(*-----------------------------------------------*)
(*Green's tensor*)
Gt[\[Theta]_,\[Phi]_]:=Sqrt[8*\[Pi]/3]*{ui[\[Theta],\[Phi],1], ui[\[Theta],\[Phi],2],ui[\[Theta],\[Phi],3]};
Print["Green's tensor = ", MatrixForm[Gt[\[Theta],\[Phi]]]]
(*-----------------------------------------------*)
(*Scattered electric field*)
Esc[\[Alpha]0_,E0_,\[CapitalDelta]\[Delta]_,\[CapitalDelta]\[Epsilon]_,\[Theta]_,\[Phi]_]:=Sqrt[\[Epsilon]0*c/2]*(\[Omega]0^2/(4*\[Pi]*\[Epsilon]0*c^2))*Exp[((I*2*\[Pi])/\[Lambda])*r]*Gt[\[Theta],\[Phi]].dipole[\[Alpha]0,E0,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon]];
(*-------------------*)
(*Scattered intensity*)
simplifyconditions = {c>0, \[Epsilon]0>0, \[Alpha]0 > 0, \[Omega]0 > 0, \[CapitalDelta]\[Delta] > 0, \[CapitalDelta]\[Epsilon] > 0, E0 > 0, \[Lambda]>0, r>0, Plo>0, PP0>0, \[Eta]>0, Element[\[Theta], Reals], Element[\[Phi], Reals]};
inte\[Theta]\[Phi] = Simplify[Conjugate[Esc[\[Alpha]0,E0,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon],\[Theta],\[Phi]]].Esc[\[Alpha]0,E0,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon],\[Theta],\[Phi]]*Sin[\[Theta]], Assumptions->simplifyconditions];
inte\[Theta] = Integrate[inte\[Theta]\[Phi],{\[Phi],0,2*\[Pi]}, Assumptions->simplifyconditions];
Intensity[\[Theta]DD_,\[CapitalDelta]\[CapitalDelta]\[Delta]_,\[CapitalDelta]\[CapitalDelta]\[Epsilon]_,\[Alpha]\[Alpha]0_,EE0_]:= Evaluate[Integrate[inte\[Theta],{\[Theta],0,\[Theta]D}, Assumptions->simplifyconditions]/.{\[Theta]D->\[Theta]DD,\[CapitalDelta]\[Delta]->\[CapitalDelta]\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon]->\[CapitalDelta]\[CapitalDelta]\[Epsilon],\[Alpha]0->\[Alpha]\[Alpha]0, E0->EE0}];
(*Integrated scattered power. The integral is doubled to consider the regions 0<\[Theta]<\[Theta]D and \[Pi]-\[Theta]D<\[Theta]<\[Pi] simultanously. Allowed values of \[Theta]D are in the interval [0,\[Pi]/2]*)
P0full[\[Theta]D_,\[CapitalDelta]\[Delta]_,\[CapitalDelta]\[Epsilon]_,\[Alpha]\[Alpha]0_,EE0_]:=Simplify[2*Intensity[\[Theta]D,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon],\[Alpha]\[Alpha]0,EE0],Assumptions->simplifyconditions];
P0[\[Theta]D_,\[CapitalDelta]\[Delta]_,\[CapitalDelta]\[Epsilon]_,\[Alpha]\[Alpha]0_,EE0_]:=Simplify[Normal[Series[P0full[\[Theta]D,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon],\[Alpha]\[Alpha]0,EE0],{\[CapitalDelta]\[Delta],0,1},{\[CapitalDelta]\[Epsilon],0,1}]],Assumptions->simplifyconditions];
Print["Scattered power in free space. P0[\[Pi]/2] = ", Simplify[P0[\[Pi]/2,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon],\[Alpha]0,E0]]]
Print["Scattered power in the cone with condition \[Theta]<\[Theta]D to first order in \[CapitalDelta]\[Delta] and \[CapitalDelta]\[Epsilon]. I[\[Theta]D] = ", Simplify[P0[\[Theta]D,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon],\[Alpha]0,E0]]]
(*-------------------*)
(*Interference*)
(*Local oscillator field*)
Elo[Plo_,\[Theta]_,\[Phi]_,n_]:=Sqrt[Plo]*ui[\[Theta],\[Phi],2]*Exp[I*((2*\[Pi]*r)/\[Lambda])]
(*Converting the scattered power expression from \[Alpha]0,E0 to PP0*)
EscP0[PPP0_,\[CapitalDelta]\[Delta]_,\[CapitalDelta]\[Epsilon]_,\[Theta]_,\[Phi]_]:= Evaluate[Simplify[Esc[\[Alpha]0,E0,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon],\[Theta],\[Phi]]/.{E0->E0*Sqrt[PP0/P0[\[Pi]/2,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon],\[Alpha]0,E0]]},Assumptions->simplifyconditions]/.{PP0->PPP0}];
(*Intensity of the interference.*)
(*WARNING: In this part, powers Plo and PP0 are squared for easier Taylor expansion*)
intdet\[Theta]\[Phi][Plo_,PP0_,\[CapitalDelta]\[Delta]_,\[CapitalDelta]\[Epsilon]_,\[Theta]_,\[Phi]_]:=Evaluate[Simplify[Conjugate[(Elo[Plo^2,\[Theta],\[Phi],2]+EscP0[PP0^2,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon],\[Theta],\[Phi]])].(Elo[Plo^2,\[Theta],\[Phi],2]+EscP0[PP0^2,\[CapitalDelta]\[Delta],\[CapitalDelta]\[Epsilon],\[Theta],\[Phi]]), Assumptions->simplifyconditions]];
(*-------------------*)
(*Measurement imprecision*)
(*Taylor expansion to order 0 and 1 in the variable of the measurement. In this case \[CapitalDelta]\[Delta], which is now called x*)
ty0[Plo_,PP0_,x_,\[Theta]_,\[Phi]_]:=Evaluate[Simplify[intdet\[Theta]\[Phi][Plo,PP0,0,0,\[Theta],\[Phi]], Assumptions->simplifyconditions]];
ty1[Plo_,PP0_,x_,\[Theta]_,\[Phi]_]:=Evaluate[Simplify[Normal[Series[intdet\[Theta]\[Phi][Plo,PP0,x,0,\[Theta],\[Phi]],{x,0,1}]], Assumptions->simplifyconditions]];
(*Extracting the part of the detected signal that is actually the signal connect to the variable to measure*)
d\[Beta][Plo_,PP0_,x_,\[Theta]_,\[Phi]_]:=Evaluate[Simplify[Normal[Series[ty1[Plo,PP0,x,\[Theta],\[Phi]]-ty0[Plo,PP0,x,\[Theta],\[Phi]],{PP0,0,1}]], Assumptions->simplifyconditions]];
Print["Actual signal cahnging with x. d\[Beta] = ", Simplify[d\[Beta][(Plo)^(1/2),(PP0)^(1/2),x,\[Theta],\[Phi]], Assumptions->simplifyconditions]]
(*differential power spectral density of the noise*)
d\[Sigma][Plo_,PP0_,x_,\[Theta]_,\[Phi]_]:=(hbar*c/\[Lambda])*Evaluate[Simplify[Normal[Series[ty0[Plo,PP0,x,\[Theta],\[Phi]],{PP0,0,0}]], Assumptions->simplifyconditions]];
Print["Power spectral density of the noise. d\[Sigma] = ", Simplify[d\[Sigma][(Plo)^(1/2),sqPP0,x,\[Theta],\[Phi]], Assumptions->simplifyconditions]]
(*Inverse differential imprecision for an ideal spherical differential detector.*)
sm1[Plo_,PP0_,x_,\[Theta]_,\[Phi]_]:=Evaluate[Simplify[d\[Beta][Plo,PP0,x,\[Theta],\[Phi]]^2/d\[Sigma][Plo,PP0,x,\[Theta],\[Phi]], Assumptions->simplifyconditions]];
Print["Inverse differential imprecision. d\[Beta]^2/d\[Sigma] = ",Simplify[sm1[(Plo)^(1/2),(PP0)^(1/2),x,\[Theta],\[Phi]], Assumptions->simplifyconditions]]
(*Integrating the differential imprecision to obtain the total imprecision*)
S\[Delta]\[Delta]\[Theta]m1[PP0_,x_,\[Theta]_]:= Evaluate[Integrate[sm1[PP0,x,\[Theta],\[Phi]]*Sin[\[Theta]],{\[Phi],0,2\[Pi]},Assumptions->simplifyconditions]];
S\[Delta]\[Delta]tmp[PP0_,x_,\[Theta]D_]:= Evaluate[1/Integrate[S\[Delta]\[Delta]\[Theta]m1[PP0,x,\[Theta]],{\[Theta],0,\[Theta]D},Assumptions->simplifyconditions]];
(*WARNING: Restoring the correct powers of Plo and PP0*)
S\[Delta]\[Delta][PP0_,x_,\[Theta]D_]:=Evaluate[Simplify[S\[Delta]\[Delta]tmp[Sqrt[PP0],x,\[Theta]D],Assumptions->simplifyconditions]];
Print["Total imprecision in free space. S\[Delta]\[Delta][\[Pi]] = ", S\[Delta]\[Delta][PP0,x,\[Pi]]/.{\[Lambda]->(2*\[Pi]*c)/\[Omega]0}]
Print["Total imprecision in the cone with condition \[Theta]<\[Theta]D. S\[Delta]\[Delta][\[Theta]D] = ", S\[Delta]\[Delta][PP0,x,\[Theta]D]]
(*-------------------*)
(*Back action*)
(*Torque fluctuation restricted to the emission cone.*)
Syy\[Tau]\[Tau][\[Theta]D_]:=(hbar/(2*\[Pi]*\[Omega]0))*(x^2)*P0[\[Theta]D,x,y,\[Alpha]0,E0];
Print["Back action free space. Syy\[Tau]\[Tau][\[Pi]] = ", Simplify[Syy\[Tau]\[Tau][\[Pi]], Assumptions->simplifyconditions]]
Print["Back action attribute to the scattering in the cone with condition \[Theta]<\[Theta]D. Syy\[Tau]\[Tau][\[Theta]D] = ", Syy\[Tau]\[Tau][\[Theta]D]]
(*-------------------*)
(*Product of back action and imprecision to test Heisenberg limit*)
Print["Product between imprecision and back action = ", Simplify[Syy\[Tau]\[Tau][\[Theta]D]*S\[Delta]\[Delta][P0[\[Pi],x,y,\[Alpha]0,E0],x,\[Theta]D], Assumptions->simplifyconditions]/.{\[Lambda]->(2*\[Pi]*c)/\[Omega]0}]


(*---------------------------------------------------------------------------------------------*)
(*---------------------------------------------------------------------------------------------*)
(*RADIATED POWER BY A LEVITATED SPHERE + SPHERICAL EMITTER*)
(*The sphere has radius R0 and it is located in the center of the reference system*)
(*---------------------------------------------------------------------------------------------*)
(*-----------------------------------------------*)
(*Optical Green's function in the measurement region (free space)*)
gm[\[Theta]k_,\[Phi]k_,\[Theta]0_,\[Phi]0_,\[Theta]_,\[Phi]_,r0_]:=Exp[I*(2*\[Pi]*r0/\[Lambda])*nr[\[Theta]k,\[Phi]k].nr[\[Theta]0,\[Phi]0]]*Exp[I*(2*\[Pi]*r0/\[Lambda])*nr[\[Theta],\[Phi]].nr[\[Theta]0,\[Phi]0]];
(*Optical Green's function of the image generated by the hemipherical mirror. It cannot be used directly*)
gi[\[Theta]k_,\[Phi]k_,\[Theta]0_,\[Phi]0_,\[Theta]_,\[Phi]_,\[Rho]_,R_,r0_]:=-\[Rho]*Exp[I*(2*\[Pi]*r0/\[Lambda])*nr[\[Theta]k,\[Phi]k].nr[\[Theta]0,\[Phi]0]]*Exp[-I*(2*\[Pi]*r0/\[Lambda])*(nr[\[Theta],\[Phi]].nr[\[Theta]0,\[Phi]0]-2*R/r0)];
(*Optical Green's function in the control region: direct field + image field*)
gc[\[Theta]k_,\[Phi]k_,\[Theta]0_,\[Phi]0_,\[Theta]_,\[Phi]_,r0_]:=Evaluate[Simplify[gm[\[Theta]k,\[Phi]k,\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0]+gi[\[Theta]k,\[Phi]k,\[Theta]0,\[Phi]0,\[Theta],\[Phi],1,\[Lambda],r0]]];
Print["Green's function free space/measurement region gm = ", gm[\[Theta]k,\[Phi]k,\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0]]
Print["Green's function control region gc = ", gc[\[Theta]k,\[Phi]k,\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0]]
(*Simplify conditions: angles real, wavelenght positive, radius positive*)
simplfyGreenConditions={\[Lambda]>0, \[Theta]>0, \[Phi]>0, \[Phi]0>0, \[Theta]0>0, \[Phi]k>0, \[Theta]k>0, r0>0, R0>0};
(*Green's functions to n-th order*)
Print["Expanding Greens' functions ..."]
gmn1[\[Theta]0_,\[Phi]0_,\[Theta]_,\[Phi]_,r0_]:= Evaluate[Simplify[Normal[Series[gm[\[Pi]/2,\[Pi]/2,\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0],{r0,0,1}]], Assumptions->simplfyGreenConditions]];
gmno[\[Theta]0_,\[Phi]0_,\[Theta]_,\[Phi]_,r0_]:= Evaluate[Simplify[Normal[Series[gm[\[Pi]/2,\[Pi]/2,\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0],{r0,0,2}]], Assumptions->simplfyGreenConditions]];
gcn1[\[Theta]0_,\[Phi]0_,\[Theta]_,\[Phi]_,r0_]:= Evaluate[Simplify[Normal[Series[gc[\[Pi]/2,\[Pi]/2,\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0],{r0,0,1}]], Assumptions->simplfyGreenConditions]];
gcno[\[Theta]0_,\[Phi]0_,\[Theta]_,\[Phi]_,r0_]:= Evaluate[Simplify[Normal[Series[gc[\[Pi]/2,\[Pi]/2,\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0],{r0,0,4}]], Assumptions->simplfyGreenConditions]];
Print["gm 1st = ", gmn1[\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0]]
Print["gc 1st = ", gcn1[\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0]]
Print["... Greens' functions expanded."]


fr[c_,r_,\[Theta]_,\[Phi]_]:=(ss*SphericalHarmonicY[0,0, \[Theta], \[Phi]]+c[[1]]*SphericalHarmonicY[1,-1, \[Theta], \[Phi]]+c[[2]]*SphericalHarmonicY[1, 0, \[Theta], \[Phi]]+c[[3]]*SphericalHarmonicY[1, 1, \[Theta], \[Phi]]);
mr[r_,\[Theta]_,\[Phi]_]:=r*{Sin[\[Theta]]*Exp[I*\[Phi]]/Sqrt[2],Cos[\[Theta]],Sin[\[Theta]]*Exp[-I*\[Phi]]/Sqrt[2]};
simplfyGreenintConditions={\[Lambda]>0, \[Theta]>0, \[Phi]>0, \[Phi]0>0, \[Theta]0>0, \[Phi]k>0, \[Theta]k>0, r0>0, R0>0, \[Theta]1>0, \[Phi]1>0};
Integrate[Integrate[fr[mr[c,\[Theta]1,\[Phi]1],r,\[Theta],\[Phi]]*gmn1[\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0],{\[Phi]0,0,2*\[Pi]}, Assumptions->simplfyGreenintConditions],{\[Theta]0,0,\[Pi]}, Assumptions->simplfyGreenintConditions]


(*The integral of the Green's function over the levitated sphere gives the electric field on the unitary sphere where the detector is located.*)
(*Normalization factors are disregarded because we are interested only in relative variations between the measurement region and the control region.*)
Print["Calculating the electric field ..."]
Emndet[\[Theta]_,\[Phi]_,R0_]:=Evaluate[Integrate[Integrate[Integrate[gmno[\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0]*r0^2*Sin[\[Theta]0],{\[Theta]0,0,\[Pi]}],{\[Phi]0,0,2*\[Pi]}],{r0,0,R0}]];
Ecndet[\[Theta]_,\[Phi]_,R0_]:=Evaluate[Integrate[Integrate[Integrate[gcno[\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0]*r0^2*Sin[\[Theta]0],{\[Theta]0,0,\[Pi]}],{\[Phi]0,0,2*\[Pi]}],{r0,0,R0}]];
Print["... electric field calculated"]
(*The power is proportional to the square of the electric field*)
Print["Calculating the radiated power ..."]
pmndet[\[Theta]_,\[Phi]_,R0_]:= Evaluate[Simplify[Conjugate[Emndet[\[Theta],\[Phi],R0]]*Emndet[\[Theta],\[Phi],R0], Assumptions->simplfyGreenConditions]];
pcndet[\[Theta]_,\[Phi]_,R0_]:= Evaluate[Simplify[Conjugate[Ecndet[\[Theta],\[Phi],R0]]*Ecndet[\[Theta],\[Phi],R0], Assumptions->simplfyGreenConditions]];
(*Integrals of the power in the measurement region and in the control region. Warning: different domains of integration.*)
Pmndet\[Theta][\[Theta]_,R0_]:= Evaluate[Integrate[pmndet[\[Theta],\[Phi],R0]*Sin[\[Theta]],{\[Phi],0,2*\[Pi]}]];
Pmndet[\[Theta]D_,RR0_]:=Evaluate[Simplify[Integrate[Pmndet\[Theta][\[Theta],R0],{\[Theta],0,\[Theta]D}]+Integrate[Pmndet\[Theta][\[Theta],R0],{\[Theta],\[Pi]-\[Theta]D,\[Pi]}], Assumptions->simplfyGreenConditions]/.{R0->RR0}];
Pmndetout[\[Theta]D_,RR0_]:=Evaluate[Simplify[Integrate[Pmndet\[Theta][\[Theta],R0],{\[Theta],\[Theta]D,\[Pi]-\[Theta]D}], Assumptions->simplfyGreenConditions]/.{R0->RR0}];
Pcndet[\[Theta]D_,RR0_]:=Evaluate[Integrate[Integrate[pcndet[\[Theta],\[Phi],R0]*Sin[\[Theta]],{\[Phi],0,2*\[Pi]}],{\[Theta],\[Theta]D,\[Pi]/2}]/.{R0->RR0}];
Pc2det[\[Theta]D_,R0_]:=Evaluate[Simplify[Normal[Series[Pcndet[\[Theta]D,R0],{R0,0,10}]], Assumptions->simplfyGreenConditions]];
Print["... radiated power calculated"]
Print["Power in the control region 4th order : ", Pcndet[\[Theta]D,R0]]
Print["Power in the control region 2nd order : ", Pc2det[\[Theta]D,R0]]
Print["Power in the measurement region 2nd order : ", Pmndet[\[Theta]D,R0]]
(*Estimate of how good is the series approximation as a function of the radius R0 the levitated sphere. This is done by comparing the calcualtion
in series expansion to order 2 and order 4.*)
relativeuncertainty[R0_]:=Evaluate[Simplify[(Pc2det[\[Theta]D,R0]-Pcndet[\[Theta]D,R0])/Pcndet[\[Theta]D,R0], Assumptions->{\[Theta]D>0,R0>0,\[Lambda]>0}]];
(*Power ratio between the scattered light in the measurement region and in the control region in presence of the hemipsherical mirror*)
Pratio2[\[Theta]D_,R0_]:=Evaluate[Simplify[Pc2det[\[Theta]D,R0]/Pmndet[\[Theta]D,R0]]];
Pratio4[\[Theta]D_,R0_]:=Evaluate[Simplify[Pcndet[\[Theta]D,R0]/Pmndet[\[Theta]D,R0]]];
(*Power ratio between the scattered light in the measurement region and in the control region without the hemipsherical mirror*)
Pratio2out[\[Theta]D_,R0_]:=Evaluate[Simplify[Pmndetout[\[Theta]D,R0]/Pmndet[\[Theta]D,R0]]];
(*------------------------------------------------------------------*)
(*Plots and values*)
Print["relative deviation of the scattered power between 2nd order to 4th order = ", relativeuncertainty[R0], "   (it does not depend on the angles!)"]
Print["relative deviation of the scattered power between 2nd order to 4th order at the condition R0=\[Lambda]/5 -> ",N[Simplify[relativeuncertainty[\[Lambda]/5]]]]
(*2D axes label*)
xyzlinestyle = {Directive[Black, Thick, Dashing[None]],
Directive[Black, Thick, Dashed], Directive[Red, Thick, Dashing[None]], Directive[Red, Thick, Dashed]};
xylabl = {Style["NA",Bold,Black,lgdfontsize],Style["(P2-P4)/P2",Bold,Black,lgdfontsize]};
pltlgdxyz={Style["relative",Black,lgdfontsize]};
Plot[relativeuncertainty[R0]/.{\[Lambda]->1},{R0,0,0.3},PlotLegend->pltlgdxyz, LegendPosition->{-0.4,-0.2},LegendShadow->None,LegendSize->1,LegendBorder->Black,
	ImageSize->Large, AxesLabel -> xylabl , Axes -> True, Frame->frameflg, PlotStyle->xyzlinestyle, BaseStyle->{FontSize->16},
	 PlotTheme->"Monochrome", AspectRatio -> 3/4, Background->White]
xylabl = {Style["R0/\[Lambda]",Bold,Black,lgdfontsize],Style["Pc/Pm",Bold,Black,lgdfontsize]};
pltlgdxyz={Style["full NA - 2nd order",Black,lgdfontsize],Style["full NA - 4th order",Black,lgdfontsize],Style["NA=0.4 - 2nd order",Black,lgdfontsize], Style["NA=0.4 - 4th order",Black,lgdfontsize]};
LogPlot[{(Pc2det[ArcSin[0],R0]/Pm2det[ArcSin[1],R0])/.{\[Lambda]->1},(Pc4det[ArcSin[0],R0]/Pm2det[ArcSin[1],R0])/.{\[Lambda]->1},Pratio2[ArcSin[0.4],R0]/.{\[Lambda]->1},Pratio4[ArcSin[0.4],R0]/.{\[Lambda]->1}},{R0,0,0.3},PlotRange->{0.0001,1},PlotLegend->pltlgdxyz,
	LegendPosition->{-0.4,-0.2},LegendShadow->None,LegendSize->1,LegendBorder->Black,
	ImageSize->Large, AxesLabel -> xylabl , Axes -> True, Frame->frameflg, PlotStyle->xyzlinestyle, BaseStyle->{FontSize->16},
	 PlotTheme->"Monochrome", AspectRatio -> 3/4, Background->White]
xylabl = {Style["NA",Bold,Black,lgdfontsize],Style["Pc/Pm",Bold,Black,lgdfontsize]};
pltlgdxyz={Style["\[Lambda]/10 - 2nd order",Black,lgdfontsize],Style["\[Lambda]/10 - 4th order",Black,lgdfontsize],Style["\[Lambda]/5 - 2nd order",Black,lgdfontsize],Style["\[Lambda]/5 - 4th order",Black,lgdfontsize]};
LogPlot[{Pratio2[ArcSin[NA],\[Lambda]/10]/.{\[Lambda]->1},Pratio4[ArcSin[NA],\[Lambda]/10]/.{\[Lambda]->1},Pratio2[ArcSin[NA],\[Lambda]/5]/.{\[Lambda]->1},Pratio4[ArcSin[NA],\[Lambda]/5]/.{\[Lambda]->1}},{NA,0,1},
    PlotRange->{0.001,1},PlotLegend->pltlgdxyz,
	LegendPosition->{-0.4,-0.2},LegendShadow->None,LegendSize->1,LegendBorder->Black,
	ImageSize->Large, AxesLabel -> xylabl , Axes -> True, Frame->frameflg, PlotStyle->xyzlinestyle, BaseStyle->{FontSize->16},
	 PlotTheme->"Monochrome", AspectRatio -> 3/4, Background->White]
Print["scattered power ratio at conditions: with mirror, NA=0.4, R0=\[Lambda]/10 (2nd order) = ", N[Simplify[Pratio2[ArcSin[0.4],\[Lambda]/10]]]]
Print["scattered power ratio at conditions: with mirror, NA=0.4, R0=\[Lambda]/10 (4th order) = ", N[Simplify[Pratio4[ArcSin[0.4],\[Lambda]/10]]]]
Print["scattered power ratio at conditions: with mirror, NA=0.4, R0=\[Lambda]/5  (2nd order) = ", N[Simplify[Pratio2[ArcSin[0.4],\[Lambda]/5]]]]
Print["scattered power ratio at conditions: with mirror, NA=0.4, R0=\[Lambda]/5  (4th order) = ", N[Simplify[Pratio4[ArcSin[0.4],\[Lambda]/5]]]]
Print["scattered power ratio at conditions: without mirror, NA=0.4, R0=\[Lambda]/10 (2nd order) = ", N[Simplify[Pratio2out[ArcSin[0.4],\[Lambda]/10]]]]
Print["scattered power ratio at conditions: without mirror, NA=0.4, R0=\[Lambda]/5  (2nd order) = ", N[Simplify[Pratio2out[ArcSin[0.4],\[Lambda]/5]]]]


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
(*PRODUCT OF IMPRECISION AND BACK ACTION FOR HYBRID SYSTEM COMPOSED BY THE HOLLOW MIRROR AND SELF-HOMODYNE*)
Print["detection limit x = ", Simplify[Smirrorvet[\[Theta]][[1]]*Sbavet[\[Theta]][[1]]]]
Print["detection limit y = ", Simplify[Smirrorvet[\[Theta]][[2]]*Sbavet[\[Theta]][[2]]]]
Print["detection limit z = ", Simplify[Smirrorvet[\[Theta]][[3]]*Sbavet[\[Theta]][[3]]]]
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
(*----------------------------------------------------------------------------------------------*)
(*2D axes label*)
xylabl = {Style["NA",Bold,Black,lgdfontsize],Style["Sba",Bold,Black,lgdfontsize]};
(*curve specs*)
basupplinestyle = {Directive[Red, Thick, Dashing[None]],
Directive[Green, Thick, Dashing[None]], Directive[Blue, Thick, Dashing[None]], Directive[Black, Thick, Dashed],  Directive[Black, Thick, Dashed]};
(*plotting*)
Plot[{Sbarvet[ArcSin[x]][[1]],Sbarvet[ArcSin[x]][[2]],Sbarvet[ArcSin[x]][[3]], Sba[\[Pi]/2,0,0]/Sba[\[Pi]/2,0,0], Sba[\[Pi]/2,\[Pi]/2,0]/Sba[\[Pi]/2,0,0]},{x,0,1},
	PlotLegend->pltlgdxyz,
	LegendPosition->{-0.4,-0.2},LegendShadow->None,LegendSize->1,LegendBorder->Black,
	ImageSize->Large, AxesLabel -> xylabl , Axes -> True, Frame->frameflg, PlotStyle->basupplinestyle, BaseStyle->{FontSize->16},
	 PlotTheme->"Monochrome", AspectRatio -> 3/4, Background->White]


(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*DISPLACEMENT - VARIANCE operator*)
(*radially symmetric Gaussian density distribution*)
fgauss[x_,y_,z_,\[Sigma]_]:=\[Rho]*Exp[-(x^2+y^2+z^2)/\[Sigma]^2];
(*Gradient in cartesian coordinates*)
Delfgauss[x_,y_,z_,\[Sigma]_]:= Grad[fgauss[x,y,z,\[Sigma]],{x,y,z}];
(*Geometric factor*)
geofgauss[\[Theta]_,\[Phi]_,\[Sigma]_]:=Evaluate[Simplify[Integrate[Integrate[Integrate[Delfgauss[x,y,z,\[Sigma]]*(nr[\[Theta],\[Phi]].{x,y,z}),{x,-Infinity,Infinity}, Assumptions->{Re[\[Sigma]^2]>0}],{y,-Infinity,Infinity}, Assumptions->{Re[\[Sigma]^2]>0}],{z,-Infinity,Infinity}, Assumptions->{Re[\[Sigma]^2]>0}],Assumptions->{\[Sigma]>0}]];
(*---*)
Print["fgauss = ", fgauss[x,y,z,\[Sigma]]]
Print["normalization fgauss = ", Simplify[Integrate[Integrate[Integrate[fgauss[x,y,z,\[Sigma]],{x,-Infinity,Infinity}, Assumptions->{Re[\[Sigma]^2]>0}],{y,-Infinity,Infinity}, Assumptions->{Re[\[Sigma]^2]>0}],{z,-Infinity,Infinity}, Assumptions->{Re[\[Sigma]^2]>0}],Assumptions->{\[Sigma]>0}]];
Print["Delfgauss = ", Delfgauss[x,y,z,\[Sigma]]]
Print["Geometric factors fgauss = ", geofgauss[\[Theta],\[Phi],\[Sigma]]]
(*----------------------------------------------*)
(*radially symmetric homogenous density*)
fsphere[r_,r0_]:=\[Rho]*HeavisideTheta[r0-r];
(*Gradient in spherical coordinates*)
Delfsphere[r_,\[Theta]_,\[Phi]_,r0_]:= Evaluate[D[fsphere[r,r0],{r,1}]*nr[\[Theta],\[Phi]]];
(*Geometric factor*)
fgeofsphere[r0_,\[Theta]n_,\[Phi]n_]:=Evaluate[Integrate[Integrate[Integrate[r^3*Delsphere[r,\[Theta],\[Phi],r0]*Sin[\[Theta]]*(nr[\[Theta],\[Phi]].nr[\[Theta]n,\[Phi]n]),{r,0,Infinity}, Assumptions->{r0>0}],{\[Phi],0,2*\[Pi]}],{\[Theta],0,\[Pi]/2}]];
(*---*)
Print["Normalization sphere = ", Evaluate[Integrate[Integrate[Integrate[r^2*Sin[\[Theta]]*fsphere[r,r0],{r,0,Infinity}, Assumptions->{r0>0}],{\[Phi],0,2*\[Pi]}],{\[Theta],0,\[Pi]}]]]
Print["fsphere = ", fsphere[r,r0]]
Print["Delsphere = ", Delsphere[r,\[Theta],\[Phi],r0]]
Print["Geometric factor fsphere = ", fgeofsphere[r0,\[Theta],\[Phi]]]


(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*Asymmetries - ROTATION*)
(*Approximating asymmetries with spherical harmonics*)
fr[c_,r_,\[Theta]_,\[Phi]_]:=(c[[1]]*SphericalHarmonicY[1,-1, \[Theta], \[Phi]]+c[[2]]*SphericalHarmonicY[1, 0, \[Theta], \[Phi]]+c[[3]]*SphericalHarmonicY[1, 1, \[Theta], \[Phi]]);
(*Radial vector in the spherical harmonics basis*)
mr[r_,\[Theta]_,\[Phi]_]:=r*{Sin[\[Theta]]*Exp[I*\[Phi]]/Sqrt[2],Cos[\[Theta]],Sin[\[Theta]]*Exp[-I*\[Phi]]/Sqrt[2]};
(*assumptions for integrals*)
asyassumptions = {c0>0, \[Theta]1>0, \[Phi]>0, \[Phi]1>0, r>0, \[Theta]d>0, \[Theta]a>0, \[Theta]>0}
(*Show that the zero order of the expansion of the object's shape does not contribute*)
Print["contribution from Y00 = ", Integrate[Integrate[SphericalHarmonicY[0,0, \[Theta], \[Phi]]*Sin[\[Theta]]*(nr[\[Theta]n,\[Phi]n].nr[\[Theta],\[Phi]]),{\[Theta],0,\[Pi]}],{\[Phi],0,2*\[Pi]}]]
(*Scattered electric field*)
Er[rr_,\[Theta]\[Theta]_,\[Phi]\[Phi]_,\[Theta]1_,\[Phi]1_]:=Evaluate[Integrate[Integrate[Integrate[r*fr[mr[c0,\[Theta]1,\[Phi]1],r,\[Theta],\[Phi]]*Sin[\[Theta]]*(nr[\[Theta]n,\[Phi]n].nr[\[Theta],\[Phi]]),{\[Theta],0,\[Pi]}],{\[Phi],0,2*\[Pi]}],{r,0,R0}]/.{R0->rr, \[Theta]n->\[Theta]\[Theta], \[Phi]n->\[Phi]\[Phi]}];
(*Detected intensity by a differential pixel detector occupying the emtpy half-sphere*)
Irplot[\[Theta]_,\[Phi]_,\[Theta]1_,\[Phi]1_]:=Evaluate[Simplify[(Conjugate[Er[1,\[Theta],\[Phi],\[Theta]1,\[Phi]1]]*Er[1,\[Theta],\[Phi],\[Theta]1,\[Phi]1])/.{c0->1}, Assumptions->asyassumptions]];
(*How the intensity in the direction of the detector looks like depending on rotation*)
(*----------------*)
(*3D axes label*)
xyzlabl = {Style["x",Bold,Black,lgdfontsize],Style["y",Bold,Black,lgdfontsize],Style["I",Bold,Black,lgdfontsize]};
(*The object is rotating in a plane including the optical axis*)
Animate[RevolutionPlot3D[Irplot[\[Theta],\[Phi],\[Theta]a,0],{\[Theta],0,\[Pi]/2},{\[Phi],0,2*\[Pi]}, AxesLabel ->xyzlabl],{\[Theta]a,0,\[Pi]},AnimationRunning->False]
(*The object is rotating in the plane orthoghonal to the optical axis*)
Animate[RevolutionPlot3D[Irplot[\[Theta],\[Phi],\[Pi]/2,\[Phi]a],{\[Theta],0,\[Pi]/2},{\[Phi],0,2*\[Pi]}, AxesLabel ->xyzlabl],{\[Phi]a,0,\[Pi]},AnimationRunning->False]



(*------------------------------------------------------*)
(*Examples for angle detection*)
(*Detection of rotations orthogonal to the optical axis. The object is fixed in the x-y plane*)
Ir\[Phi][\[Theta]_,\[Phi]_,\[Phi]1_]:=Evaluate[Simplify[Irplot[\[Theta],\[Phi],\[Pi]/2,\[Phi]1],Assumptions->asyassumptions]];
(*Scalar product. The phase of the reference rotation was selected to be \[Pi]/4 to simplify the calculation*)
Ir\[Phi]prj[\[Phi]1_]:=Evaluate[Integrate[Integrate[Ir\[Phi][\[Theta],\[Phi],\[Phi]1]*Ir\[Phi][\[Theta],\[Phi],\[Pi]/4]*Sin[\[Theta]],{\[Phi],0,2*\[Pi]}],{\[Theta],0,\[Pi]/2}]];
Print[Ir\[Phi]prj[\[Phi]1]]


Plot[Ir\[Phi]prj[\[Phi]1],{\[Phi]1,-\[Pi],\[Pi]}]




