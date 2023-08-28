(* ::Package:: *)

(*
Copyright: Giovanni Cerchiari, Yannick Weiser
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
(* GREEN'S TENSORS*)
(*-----------------------------------------------*)
(*The electric field is calculated by
E = -\[Omega]^2 (Gmat.pversor) * Integrate[ p(x)*g(n,x) d3x]
The full calculation to arrive to tensorial Green's function is presented in file "Green_function.wl". Here we just use the results
*)
(*Matrix part of the Green's tensor (shortcut to calculate it)*)
Gmat[\[Theta]_,\[Phi]_]:=Evaluate[FullSimplify[IdentityMatrix[3]+FullSimplify[Grad[Div[IdentityMatrix[3]*Exp[I*(x*x0+y*y0+z*z0)],{x,y,z},"Cartesian"],{x,y,z},"Cartesian"]/Exp[I*(x*x0+y*y0+z*z0)]]/.{x0->Sin[\[Theta]]*Cos[\[Phi]], y0->Sin[\[Theta]]*Sin[\[Phi]], z0->Cos[\[Theta]]}]];
Print["Matrix part of the Green's function = ", MatrixForm[Gmat[\[Theta],\[Phi]]]]
(*Scalar part of the Green's function in the measurement region (free space)*)
gm[\[Theta]k_,\[Phi]k_,\[Theta]0_,\[Phi]0_,\[Theta]_,\[Phi]_,r0_]:=I*\[Omega]*\[Mu]0*(Exp[I*(2*\[Pi]*Rdet/\[Lambda])]/(4*\[Pi]*Rdet))*Exp[I*(2*\[Pi]*r0/\[Lambda])*nr[\[Theta]k,\[Phi]k].nr[\[Theta]0,\[Phi]0]]*Exp[I*(2*\[Pi]*r0/\[Lambda])*nr[\[Theta],\[Phi]].nr[\[Theta]0,\[Phi]0]];
(*Scalar part of the Green's function of the image generated by the hemipherical mirror. It cannot be used directly*)
gi[\[Theta]k_,\[Phi]k_,\[Theta]0_,\[Phi]0_,\[Theta]_,\[Phi]_,\[Rho]_,R_,r0_]:=-I*\[Omega]*\[Mu]0*(Exp[I*(2*\[Pi]*Rdet/\[Lambda])]/(4*\[Pi]*Rdet))*\[Rho]*Exp[I*(2*\[Pi]*r0/\[Lambda])*nr[\[Theta]k,\[Phi]k].nr[\[Theta]0,\[Phi]0]]*Exp[-I*(2*\[Pi]*r0/\[Lambda])*(nr[\[Theta],\[Phi]].nr[\[Theta]0,\[Phi]0]-2*R/r0)];
(*Scalar part of the Green's function in the control region: direct field + image field*)
gc[\[Theta]k_,\[Phi]k_,\[Theta]0_,\[Phi]0_,\[Theta]_,\[Phi]_,r0_]:=Evaluate[Simplify[gm[\[Theta]k,\[Phi]k,\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0]+gi[\[Theta]k,\[Phi]k,\[Theta]0,\[Phi]0,\[Theta],\[Phi],1,\[Lambda],r0]]];
Print["Green's function free space/measurement region gm = ", gm[\[Theta]k,\[Phi]k,\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0]]
Print["Green's function control region gc = ", gc[\[Theta]k,\[Phi]k,\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0]]
(*Simplify conditions: angles real, wavelenght positive, radius positive*)
simplfyGreenConditions={\[Lambda]>0, \[Theta]>0, \[Phi]>0, \[Phi]0>0, \[Theta]0>0, \[Phi]k>0, \[Theta]k>0, r0>0, R0>0, Rdet>0, \[Omega]>0, \[Mu]0>0};
(*Green's functions to n-th order*)
Print["Expanding the scalar part  of the Greens' functions in series..."]
gmn1[\[Theta]0_,\[Phi]0_,\[Theta]_,\[Phi]_,r0_]:= Evaluate[Simplify[Normal[Series[gm[\[Pi]/2,\[Pi]/2,\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0],{r0,0,1}]], Assumptions->simplfyGreenConditions]];
gmno[\[Theta]0_,\[Phi]0_,\[Theta]_,\[Phi]_,r0_]:= Evaluate[Simplify[Normal[Series[gm[\[Pi]/2,\[Pi]/2,\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0],{r0,0,2}]], Assumptions->simplfyGreenConditions]];
gcn1[\[Theta]0_,\[Phi]0_,\[Theta]_,\[Phi]_,r0_]:= Evaluate[Simplify[Normal[Series[gc[\[Pi]/2,\[Pi]/2,\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0],{r0,0,1}]], Assumptions->simplfyGreenConditions]];
gcno[\[Theta]0_,\[Phi]0_,\[Theta]_,\[Phi]_,r0_]:= Evaluate[Simplify[Normal[Series[gc[\[Pi]/2,\[Pi]/2,\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0],{r0,0,4}]], Assumptions->simplfyGreenConditions]];
Print["gm 1st = ", gmn1[\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0]]
Print["gc 1st = ", gcn1[\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0]]
Print["... Greens' functions expanded."]



(*---------------------------------------------------------------------------------------------*)
(*---------------------------------------------------------------------------------------------*)
(*RADIATED POWER BY A LEVITATED SPHERE + SPHERICAL EMITTER*)
(*The sphere has radius R0 and it is located in the center of the reference system*)
(*---------------------------------------------------------------------------------------------*)
(*The integral of the Green's function over the levitated sphere gives the electric field on the unitary sphere where the detector is located.*)
(*Normalization factors are disregarded because we are interested only in relative variations between the measurement region and the control region.*)
Print["Calculating the electric field (scalar part)..."]
Emndet[\[Theta]_,\[Phi]_,R0_]:=Evaluate[Integrate[Integrate[Integrate[gmno[\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0]*r0^2*Sin[\[Theta]0],{\[Theta]0,0,\[Pi]}],{\[Phi]0,0,2*\[Pi]}],{r0,0,R0}]];
Ecndet[\[Theta]_,\[Phi]_,R0_]:=Evaluate[Integrate[Integrate[Integrate[gcno[\[Theta]0,\[Phi]0,\[Theta],\[Phi],r0]*r0^2*Sin[\[Theta]0],{\[Theta]0,0,\[Pi]}],{\[Phi]0,0,2*\[Pi]}],{r0,0,R0}]];
Print["... electric field calculated"]
(*The power is proportional to the square of the electric field*)
Print["Calculating the radiated power ..."]
(*Align polarization to the x-axis and calculating the power term deriving from Gmat.*)
Pmat[\[Theta]_,\[Phi]_]:=(Gmat[\[Theta],\[Phi]].{p,0,0}).(Gmat[\[Theta],\[Phi]].{p,0,0});
(*------------------------*)
(*Power term corresponding to the scalar integral*)
pmndet[\[Theta]_,\[Phi]_,R0_]:= Evaluate[Simplify[Conjugate[Emndet[\[Theta],\[Phi],R0]]*Emndet[\[Theta],\[Phi],R0], Assumptions->simplfyGreenConditions]];
pcndet[\[Theta]_,\[Phi]_,R0_]:= Evaluate[Simplify[Conjugate[Ecndet[\[Theta],\[Phi],R0]]*Ecndet[\[Theta],\[Phi],R0], Assumptions->simplfyGreenConditions]];
(*Combining the contributions of the scalar part and the matrix part.
Integrals of the power in the measurement region and in the control region. Warning: different domains of integration.*)
Pmndet\[Theta][\[Theta]_,R0_]:= Evaluate[Integrate[Pmat[\[Theta],\[Phi]]*pmndet[\[Theta],\[Phi],R0]*Sin[\[Theta]],{\[Phi],0,2*\[Pi]}]];
Pmndet[\[Theta]D_,RR0_]:=Evaluate[Simplify[Integrate[Pmndet\[Theta][\[Theta],R0],{\[Theta],0,\[Theta]D}]+Integrate[Pmndet\[Theta][\[Theta],R0],{\[Theta],\[Pi]-\[Theta]D,\[Pi]}], Assumptions->simplfyGreenConditions]/.{R0->RR0}];
Pmndetout[\[Theta]D_,RR0_]:=Evaluate[Simplify[Integrate[Pmndet\[Theta][\[Theta],R0],{\[Theta],\[Theta]D,\[Pi]-\[Theta]D}], Assumptions->simplfyGreenConditions]/.{R0->RR0}];
Pcndet[\[Theta]D_,RR0_]:=Evaluate[Integrate[Integrate[Pmat[\[Theta],\[Phi]]*pcndet[\[Theta],\[Phi],R0]*Sin[\[Theta]],{\[Phi],0,2*\[Pi]}],{\[Theta],\[Theta]D,\[Pi]/2}]/.{R0->RR0}];
Pc2det[\[Theta]D_,R0_]:=Evaluate[Simplify[Normal[Series[Pcndet[\[Theta]D,R0],{R0,0,10}]], Assumptions->simplfyGreenConditions]];
Print["... radiated power calculated"]
Print["Power in the control region 4th order : ", Pcndet[\[Theta]D,R0]]
Print["Power in the control region 2nd order : ", Pc2det[\[Theta]D,R0]]
Print["Power in the measurement region 2nd order : ", Pmndet[\[Theta]D,R0]]
(*Estimate of how good is the series approximation as a function of the radius R0 the levitated sphere. This is done by comparing the calculation
in series expansion of higher and lower order.*)
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
plotconditions = {\[Lambda]->1, \[Mu]0->1, Rdet->1, \[Omega]->1};
Plot[relativeuncertainty[R0]/.plotconditions,{R0,0,0.3},PlotLegend->pltlgdxyz, LegendPosition->{-0.4,-0.2},LegendShadow->None,LegendSize->1,LegendBorder->Black,
	ImageSize->Large, AxesLabel -> xylabl , Axes -> True, Frame->frameflg, PlotStyle->xyzlinestyle, BaseStyle->{FontSize->16},
	 PlotTheme->"Monochrome", AspectRatio -> 3/4, Background->White]
xylabl = {Style["R0/\[Lambda]",Bold,Black,lgdfontsize],Style["Pc/Pm",Bold,Black,lgdfontsize]};
pltlgdxyz={Style["full NA - 2nd order",Black,lgdfontsize],Style["full NA - 4th order",Black,lgdfontsize],Style["NA=0.4 - 2nd order",Black,lgdfontsize], Style["NA=0.4 - 4th order",Black,lgdfontsize]};
LogPlot[{(Pc2det[ArcSin[0],R0]/Pm2det[ArcSin[1],R0])/.plotconditions,(Pc4det[ArcSin[0],R0]/Pm2det[ArcSin[1],R0])/.plotconditions,Pratio2[ArcSin[0.4],R0]/.{\[Lambda]->1},
     Pratio4[ArcSin[0.4],R0]/.plotconditions},{R0,0,0.3},PlotRange->{0.0001,1},PlotLegend->pltlgdxyz,
	LegendPosition->{-0.4,-0.2},LegendShadow->None,LegendSize->1,LegendBorder->Black,
	ImageSize->Large, AxesLabel -> xylabl , Axes -> True, Frame->frameflg, PlotStyle->xyzlinestyle, BaseStyle->{FontSize->16},
	 PlotTheme->"Monochrome", AspectRatio -> 3/4, Background->White]
xylabl = {Style["NA",Bold,Black,lgdfontsize],Style["Pc/Pm",Bold,Black,lgdfontsize]};
pltlgdxyz={Style["\[Lambda]/10 - 2nd order",Black,lgdfontsize],Style["\[Lambda]/10 - 4th order",Black,lgdfontsize],Style["\[Lambda]/5 - 2nd order",Black,lgdfontsize],Style["\[Lambda]/5 - 4th order",Black,lgdfontsize]};
LogPlot[{Pratio2[ArcSin[NA],\[Lambda]/10]/.plotconditions,Pratio4[ArcSin[NA],\[Lambda]/10]/.{\[Lambda]->1},Pratio2[ArcSin[NA],\[Lambda]/5]/.plotconditions,Pratio4[ArcSin[NA],\[Lambda]/5]/.plotconditions},{NA,0,1},
    PlotRange->{0.001,1},PlotLegend->pltlgdxyz,
	LegendPosition->{-0.4,-0.2},LegendShadow->None,LegendSize->1,LegendBorder->Black,
	ImageSize->Large, AxesLabel -> xylabl , Axes -> True, Frame->frameflg, PlotStyle->xyzlinestyle, BaseStyle->{FontSize->16},
	 PlotTheme->"Monochrome", AspectRatio -> 3/4, Background->White]
Print["scattered power ratio at conditions: with mirror, NA=0.4, R0=\[Lambda]/10 (2nd order) = ", N[Simplify[Pratio2[ArcSin[0.4],\[Lambda]/10]/.plotconditions]]]
Print["scattered power ratio at conditions: with mirror, NA=0.4, R0=\[Lambda]/10 (next to leading order) = ", N[Simplify[Pratio4[ArcSin[0.4],\[Lambda]/10]/.plotconditions]]]
Print["scattered power ratio at conditions: with mirror, NA=0.4, R0=\[Lambda]/5  (2nd order) = ", N[Simplify[Pratio2[ArcSin[0.4],\[Lambda]/5]/.plotconditions]]]
Print["scattered power ratio at conditions: with mirror, NA=0.4, R0=\[Lambda]/5  (next to leading order) = ", N[Simplify[Pratio4[ArcSin[0.4],\[Lambda]/5]/.plotconditions]]]
Print["scattered power ratio at conditions: without mirror, NA=0.4, R0=\[Lambda]/10 (2nd order) = ", N[Simplify[Pratio2out[ArcSin[0.4],\[Lambda]/10]/.plotconditions]]]
Print["scattered power ratio at conditions: without mirror, NA=0.4, R0=\[Lambda]/5  (2nd order) = ", N[Simplify[Pratio2out[ArcSin[0.4],\[Lambda]/5]/.plotconditions]]]


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
(*The factor of 2 in front account for both sides of the measurement region.*)
Sbar[\[Theta]_,\[Theta]0_,\[Phi]0_]:=Evaluate[Simplify[2*Sba[\[Theta], \[Theta]0, \[Phi]0]/Sba[\[Pi], 0, 0]]];
Sbarvet[\[Theta]_]:=Evaluate[Simplify[{Sbar[\[Theta],\[Pi]/2,0], Sbar[\[Theta],\[Pi]/2,\[Pi]/2], Sbar[\[Theta],0,0]}]];
(*Print*)
Print["power spectral density fluctuation of the radiation pressure force = ", dsba[\[Theta],\[Phi],\[Theta]p,\[Phi]p,\[Theta]0,\[Phi]0]]
Print["BackAction(x, free space) = ", Sbavet[\[Pi]][[1]]]
Print["BackAction(y, free space) = ", Sbavet[\[Pi]][[2]]]
Print["BackAction(z, free space) = ", Sbavet[\[Pi]][[3]]]
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
minS = FullSimplify[SSmirror[\[Theta]pmax,\[Theta]0max,\[Phi]0max]];
maxS = FullSimplify[SSmirror[\[Theta]pmax,\[Pi]/2,0]];
Print["Minimum imprecision (natural units): ", minS];
Print["Minimum imprecision (natural units): ", maxS];
Print[""]
Print["differential imprecision"]
Print["smirror = ", smirror[\[Theta],\[Phi],\[Pi]/2,0,\[Theta]0,\[Phi]0]]
Print[""]
Print["Total imprecision \[Theta]0=0"]
Print["Smirror(\[Theta]0=0)= ", FullSimplify[Smirror[\[Theta]D,0,0]]]
Print[""]
Print["Total imprecision \[Theta]0=\[Pi]/2"]
Print["Smirror(\[Theta]0=\[Pi]/2)= ", FullSimplify[Smirror[\[Theta]D,\[Pi]/2,\[Phi]0]]]
Print[""]
Print["Total imprecision \[Theta]D=\[Pi]/2"]
Print["Smirror(\[Theta]=\[Pi]/2)= ", FullSimplify[Smirror[\[Pi]/2,\[Theta]0,\[Phi]0]]]
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






