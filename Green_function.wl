(* ::Package:: *)

(*
Copyright: Copyright: Giovanni Cerchiari, Yannick Weiser, Tommaso Faorlin,
           Lorenz Panzl, Thomas Lafenthaler
e-mail: giovanni.cerchiari@uibk.ac.at
date : 08/2023
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

(*This calculation follows the derivation of the Green's function that can be found in this book:
Lukas Novotny, Bert Hecht - Principles of Nano-Optics - Cambridge University Press (2012) ISBN 978-1-107-00546-4 Section 2.12, Section 8.3, Appendix D. 
We will refer to this book with 'Novotny' in the rest of the text.
\.08
Here, we compute the Green function of the system in the presence of a hemipherical mirror in the far-field
*)

ClearAll["Global`*"]

(*-------------------------------------*)
(*COORDINATES*)
(*Cartesian distance*)
dist[x_,y_,z_,x0_,y0_,z0_]:=Sqrt[(x-x0)^2+(y-y0)^2+(z-z0)^2];

(*Cartesian to spherical coordinates transformation direct and inverse.
These are used often because distances and gradient for this problem have easier formulas in Cartesia coordinates,
but the Series expansion depends on the radii.*)

(*Direct change of coordinate*)
varchangextor[x_,y_,z_,r_,\[Theta]_,\[Phi]_]:={x->r*Sin[\[Theta]]*Cos[\[Phi]], y->r*Sin[\[Theta]]*Sin[\[Phi]], z->r*Cos[\[Theta]]};
(*Inverse changes of coordinate*)
varchange\[Theta]\[Phi]tox[r_,\[Theta]_,\[Phi]_,x_,y_,z_]:={Cos[\[Phi]]->x/(r*Sin[\[Theta]]), Sin[\[Phi]]->y/(r*Sin[\[Theta]]), Cos[\[Theta]]->z/r};
varchangeSin\[Theta]tox[\[Theta]_,x_,y_,z_]:={Sin[\[Theta]]->Sqrt[x^2+y^2]/Sqrt[x^2+y^2+z^2], Csc[\[Theta]]->Sqrt[x^2+y^2+z^2]/Sqrt[x^2+y^2]}
varchangertox[r_,x_,y_,z_]:={r->Sqrt[x^2+y^2+z^2]};
(*Distance, but expressed with spherical variables*)
distr[r_,\[Theta]_,\[Phi]_,r0_,\[Theta]0_,\[Phi]0_]:=(dist[x,y,z,x0,y0,z0]/.varchangextor[x,y,z,r,\[Theta],\[Phi]])/.varchangextor[x0,y0,z0,r0,\[Theta]0,\[Phi]0];
Eikror[k_,r_]:=\[Mu]*\[Mu]0*Exp[I*k*r]/(4*\[Pi]*r);


(*--------------------------------------------*)
(*Scalar Green's function for the vector potential A*)
(*r0 (x0, y0, z0) = location of the detector*)
(*r  (x, y, z)    = location of the object*)
G0f[k_,x_,y_,z_,x0_,y0_,z0_]:=(\[Mu]0*\[Mu]/(4*\[Pi]))*Exp[I*k*dist[x,y,z,x0,y0,z0]]/dist[x,y,z,x0,y0,z0];
Print[" G0f (no approximation) = ", G0f[k,x,y,z,x0,y0,z0]]

(*The far-field approximation leverages the expansion of the radial part assuming r0 >> \[Lambda] > r. 
Under this condition, it is possible to expand the distance in series of r and retain only the first terms.*)

distr0[r_,\[Theta]_,\[Phi]_,r0_,\[Theta]0_,\[Phi]0_]:=Evaluate[Normal[Simplify[Series[distr[rr,\[Theta]\[Theta],\[Phi]\[Phi],rr0,\[Theta]\[Theta]0,\[Phi]\[Phi]0],{rr,0,0}], Assumptions->{rr0>0}]]/.{rr->r,\[Theta]\[Theta]->\[Theta],\[Phi]\[Phi]->\[Phi],rr0->r0,\[Theta]\[Theta]0->\[Theta]0,\[Phi]\[Phi]0->\[Phi]0}];
distr1[r_,\[Theta]_,\[Phi]_,r0_,\[Theta]0_,\[Phi]0_]:=Evaluate[TrigExpand[Normal[Simplify[Series[distr[rr,\[Theta]\[Theta],\[Phi]\[Phi],rr0,\[Theta]\[Theta]0,\[Phi]\[Phi]0],{rr,0,1}], Assumptions->{rr0>0, rr>0}]]]/.{rr->r,\[Theta]\[Theta]->\[Theta],\[Phi]\[Phi]->\[Phi],rr0->r0,\[Theta]\[Theta]0->\[Theta]0,\[Phi]\[Phi]0->\[Phi]0}];
Print["Distance expanded (0-order) = ",distr0[r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0]]
Print["Distance expanded (1st-order) = ",distr1[r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0]]

(*At the exponent, the distance is expanded to 1st order to take into account that \[Lambda] > r, whereas at the denominator it is expanded to 0th order because r0 >> r.*)
G0fr\[Theta]\[Phi][k_,r_,\[Theta]_,\[Phi]_,r0_,\[Theta]0_,\[Phi]0_]:=(\[Mu]0*\[Mu]/(4*\[Pi]))*Exp[I*k*distr1[r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0]]/distr0[r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0];

(*Changing to Cartesian coordinates because the mirror is easier to introduce. The variable r0 is left for readability*)
G0fxyzr0tmp[kk_,xx_,yy_,zz_,xx0_,yy0_,zz0_,rr0_]:=Evaluate[FullSimplify[(G0fr\[Theta]\[Phi][k,r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0]/.varchange\[Theta]\[Phi]tox[r,\[Theta],\[Phi],x,y,z])/.varchange\[Theta]\[Phi]tox[r0,\[Theta]0,\[Phi]0,x0,y0,z0]]/.{k->kk,x->xx,y->yy,z->zz,x0->xx0,y0->yy0,z0->zz0,r0->rr0}];

(*Introducing the mirror*)
G0fxyzr0[k_,\[Rho]_,\[Psi]_,x_,y_,z_,x0_,y0_,z0_,r0_]:=Evaluate[Simplify[G0fxyzr0tmp[k,x,y,z,x0,y0,z0,r0]]-\[Rho]*G0fxyzr0tmp[k,x,y,z,-x0,-y0,-z0,r0]*Exp[-I*\[Psi]]];

(*Change of coordinates*)
G0fxyz[k_,\[Rho]_,\[Psi]_,x_,y_,z_,x0_,y0_,z0_]:=Evaluate[FullSimplify[G0fxyzr0[k,\[Rho],\[Psi],x,y,z,x0,y0,z0,r0]/.varchangertox[r0,x0,y0,z0]]];
G0fr\[Theta]\[Phi][k_,\[Rho]_,\[Psi]_,x_,y_,z_,r0_,\[Theta]0_,\[Phi]0_]:=Evaluate[Simplify[G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]/.varchangextor[x0,y0,z0,r0,\[Theta]0,\[Phi]0], Assumptions->{r0>0}]];
G0fr\[Theta]\[Phi]red[k_,\[Rho]_,\[Psi]_,x_,y_,z_,r0_,\[Theta]0_,\[Phi]0_]:=Evaluate[Simplify[G0fr\[Theta]\[Phi][k,\[Rho],\[Psi],x,y,z,r0,\[Theta]0,\[Phi]0]/Eikror[k,r0], Assumptions->{r0>0}]];
Print["-------------------"]
Print["Far-field scalar Green's function valid for the each component of A in Lorentz's gauge."]
Print["x,y,z = coordinate of the object (to be integrated to)"]
Print["x0,y0,z0 = coordinate at the detector"]
Print[" general = ", G0fxyzr0[k,\[Rho],\[Psi],x,y,z,x0,y0,z0,r0]]
Print[" free space = ", Eikror[k,r0], Simplify[G0fxyzr0[k,0,0,x,y,z,x0,y0,z0,r0]/Eikror[k,r0]]]
Print[" mirror suppression = ", Eikror[k,r0], "(",Simplify[ExpToTrig[G0fr\[Theta]\[Phi]red[k,1,0,x,y,z,r0,\[Theta]0,\[Phi]0]]/.varchange\[Theta]\[Phi]tox[r0,\[Theta]0,\[Phi]0,x0,y0,z0]], ")"]


(*Calculating the Green's function for the Electric field.*)
(* ----------------------------------------------------
The full and expanded expression is (Novotny pag. 30 Eq. 2.94):
Gt = i*\[Omega]*(I * G0 + (1/k^2)*Grad(Div(I*G0))),
where all the coefficients are already in the Green's function.
i*\[Omega]*(I * G0) is dA/dt and (1/k^2)*Grad(Div(I*G0)) originates from -Grad[\[Phi]]=-Grad[Div[A]] in Lorentz's gauge
The factor i*\[Omega] in front is due to the derivative of the temporal part.
------------------------------------------------------
There are two derivatives to do in the left term because we have to do Grad(Div()).
We do this term first. In each derivative step, only the leading term is retained.
To retain leading term, two actions are taken:
1) we use x=0, y=0, z=0 at evaluating the derivatives. This is equivalent of retaining the order 0 in the exapansion r/r0.
This can be done because the detector is far away from the sample r0>>r.
2) We keep only first term in the Series expansion of \[Alpha]=1/r0\[Rule]0 *)
GGt = G0fxyz[k,\[Rho],\[Psi],0,0,0,x0,y0,z0]*IdentityMatrix[3];
(*Perform the divergnence*)
GGt = FullSimplify[Div[GGt,{x0,y0,z0},"Cartesian"]];

(*Going to radial coordinates*)
GGt = Simplify[FullSimplify[(GGt/.varchangextor[x0,y0,z0,r0,\[Theta]0,\[Phi]0])/.{r0->1/\[Alpha]}],Assumptions->{\[Alpha]>0}];
(*Expanding assuming r0 big*)
GGt = Normal[Series[GGt,{\[Alpha],0,1}]]/.{\[Alpha]->1/r0};
(*Going back to cartesian coordinates*)
GGtxyz = (GGt/.varchange\[Theta]\[Phi]tox[r0,\[Theta]0,\[Phi]0,x0,y0,z0])/.varchangertox[r0,x0,y0,z0];

(*Reshaping for printing*)
GGtvet = FullSimplify[GGtxyz/((I*k)*G0fxyz[k,\[Rho],\[Psi],0,0,0,x0,y0,z0])];
GGtvetr = Simplify[GGtvet/.varchangextor[x0,y0,z0,r0,\[Theta]0,\[Phi]0], Assumptions->{r0>0}];
Print["GGt - first derivative only = ", I*k*Simplify[G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]/.varchangextor[x0,y0,z0,r0,\[Theta]0,\[Phi]0], Assumptions->{r0>0}], MatrixForm[GGtvetr]]


(*Performing the gradient*)
Gt = Grad[GGtxyz,{x0,y0,z0},"Cartesian"];
Gt = Simplify[(Gt/.varchangextor[x0,y0,z0,r0,\[Theta]0,\[Phi]0])/.{r0->1/\[Alpha]},Assumptions->{\[Alpha]>0}];
Gt = TrigExpand[Normal[Series[Gt,{\[Alpha],0,1}]]]/.{\[Alpha]->1/r0};
Gtxyz = FullSimplify[((Gt/.varchange\[Theta]\[Phi]tox[r0,\[Theta]0,\[Phi]0,x0,y0,z0])/.varchangertox[r0,x0,y0,z0])/.varchangeSin\[Theta]tox[\[Theta]0,x0,y0,z0]];

(*Divide by scalar to obtain the matrix part*)
Gtmat2 = FullSimplify[Gtxyz/G0fxyz[k,\[Rho],\[Psi],0,0,0,x0,y0,z0]/.varchangextor[x0,y0,z0,r0,\[Theta]0,\[Phi]0]];

(*The term with derivatives in far field is equal to
Grad(Div(I*G0)) = G0fxyz * Gtmat2,
but, before putting the matrix (Gtmat2) and the scalar (G0fxyz) togheter, we can consider that 
I * G0 = IdentityMatrix[3] * G0fxyz.
Thus, we can write the final Green function combining the matrices first:
Gf = G0 * (I + Gtmat2/k^2)
*)
Gtmat = FullSimplify[IdentityMatrix[3]+Simplify[Gtmat2/k^2, Assumptions->{k>0}]];
Gtr\[Theta]\[Phi][kk_,\[Rho]\[Rho]_,\[Psi]\[Psi]_,xx_,yy_,zz_,rr0_,\[Theta]\[Theta]0_,\[Phi]\[Phi]0_]:= Evaluate[I*\[Omega]*G0fr\[Theta]\[Phi][k,\[Rho],\[Psi],x,y,z,r0,\[Theta]0,\[Phi]0]*Gtmat/.{k->kk, \[Rho]->\[Rho]\[Rho], \[Psi]->\[Psi]\[Psi], x->xx,y->yy,z->zz,r0->rr0,\[Theta]0->\[Theta]\[Theta]0,\[Phi]0->\[Phi]\[Phi]0}];
Gtxyz = I*\[Omega]*G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]*Gtmat;

Print["-------------------"]
Print["Far-field Green's tensor valid for the Electric field"]
Print["x,y,z = coordinate of the object (to be integrated to)"]
Print["x0,y0,z0 = coordinate at the detector"]
Print["general = ", I*\[Omega], G0fxyzr0[k,\[Rho],\[Psi],x,y,z,x0,y0,z0,r0], MatrixForm[Gtmat]]
Print["free space = ", I*\[Omega], Eikror[k,r0], Simplify[G0fxyzr0[k,0,0,x,y,z,x0,y0,z0,r0]/Eikror[k,r0]], MatrixForm[Gtmat]]
Print["mirror suppression = ", I*\[Omega], Eikror[k,r0], "(",Simplify[ExpToTrig[G0fr\[Theta]\[Phi]red[k,1,0,x,y,z,r0,\[Theta]0,\[Phi]0]]/.varchange\[Theta]\[Phi]tox[r0,\[Theta]0,\[Phi]0,x0,y0,z0]], ")", MatrixForm[Gtmat]]


(*
Testing the expression with dipolar patterns emerging from a single source.
The current that is generating the field is j(r) = q*dx/dt = \[Minus]i\[Omega]p\[Delta][r], with p being the dipole moment.
The light intensity is given by I = (\[Epsilon]0*n*c/2) * |E|^2, where n is the index of refraction. The factor 1/2 comes from the time average.
*)
simconditions = {\[Mu]>0, \[Mu]0>0, \[Omega]>0, r0>0, \[Theta]0>0, \[Phi]0>0, k>0, p>0};
dPd\[CapitalOmega][jvet_,rr0_,\[Theta]\[Theta]0_,\[Phi]\[Phi]0_]:=Simplify[(Gtr\[Theta]\[Phi][k,0,0,0,0,0,r0,\[Theta]0,\[Phi]0] . jvet) . Conjugate[Gtr\[Theta]\[Phi][k,0,0,0,0,0,r0,\[Theta]0,\[Phi]0] . jvet], Assumptions->simconditions]/.{r0->rr0,\[Theta]0->\[Theta]\[Theta]0,\[Phi]0->\[Phi]\[Phi]0};
TotalIrradiatedPower[jvet_]:=(\[Epsilon]0*c/2)*Integrate[Integrate[r0^2*dPd\[CapitalOmega][jvet,r0,\[Theta]0,\[Phi]0]*Sin[\[Theta]0],{\[Theta]0,0,\[Pi]}],{\[Phi]0,0,2*\[Pi]}]/.{\[Mu]->1,\[Mu]0->1/(\[Epsilon]0*c^2)}
Print["Radiated power \[Pi] x-pol = ", TotalIrradiatedPower[-I*\[Omega]*p*{1,0,0}]]
Print["Radiated power \[Pi] y-pol = ", TotalIrradiatedPower[-I*\[Omega]*p*{0,1,0}]]
Print["Radiated power \[Pi] z-pol = ", TotalIrradiatedPower[-I*\[Omega]*p*{0,0,1}]]
Print["Radiated power \[Sigma] x-pol = ", TotalIrradiatedPower[-(I*\[Omega]*p/Sqrt[2])*{0,1,I}]]
Print["Radiated power \[Sigma] y-pol = ", TotalIrradiatedPower[-(I*\[Omega]*p/Sqrt[2])*{1,0,I}]]
Print["Radiated power \[Sigma] z-pol = ", TotalIrradiatedPower[-(I*\[Omega]*p/Sqrt[2])*{1,I,0}]]
(*plotting*)
plotconditions = {\[Omega]->1, \[Mu]0->1, \[Mu]->1, \[Epsilon]0->1, p->1, r0->1};

dPd\[CapitalOmega][-I*\[Omega]*p*{1,0,0},1,\[Theta],\[Phi]]
\[Pi]xpolplot = SphericalPlot3D[dPd\[CapitalOmega][-I*\[Omega]*p*{1,0,0},1,\[Theta],\[Phi]]/.plotconditions,{\[Theta],0,\[Pi]},{\[Phi],0,2*\[Pi]}];
\[Pi]ypolplot = SphericalPlot3D[dPd\[CapitalOmega][-I*\[Omega]*p*{0,1,0},1,\[Theta],\[Phi]]/.plotconditions,{\[Theta],0,\[Pi]},{\[Phi],0,2*\[Pi]}];
\[Pi]zpolplot = SphericalPlot3D[dPd\[CapitalOmega][-I*\[Omega]*p*{0,0,1},1,\[Theta],\[Phi]]/.plotconditions,{\[Theta],0,\[Pi]},{\[Phi],0,2*\[Pi]}];
\[Sigma]xpolplot = SphericalPlot3D[dPd\[CapitalOmega][-(I*\[Omega]*p/Sqrt[2])*{0,1,I},1,\[Theta],\[Phi]]/.plotconditions,{\[Theta],0,\[Pi]},{\[Phi],0,2*\[Pi]}];
\[Sigma]ypolplot = SphericalPlot3D[dPd\[CapitalOmega][-(I*\[Omega]*p/Sqrt[2])*{1,0,I},1,\[Theta],\[Phi]]/.plotconditions,{\[Theta],0,\[Pi]},{\[Phi],0,2*\[Pi]}];
\[Sigma]zpolplot = SphericalPlot3D[dPd\[CapitalOmega][-(I*\[Omega]*p/Sqrt[2])*{I,1,0},1,\[Theta],\[Phi]]/.plotconditions,{\[Theta],0,\[Pi]},{\[Phi],0,2*\[Pi]}];
Grid[{{Text[Style["dp/d\[CapitalOmega]", Larger, Bold]],Text[Style["x", Larger, Bold]],Text[Style["y", Larger, Bold]], Text[Style["z", Larger, Bold]]},
{Text[Style["\[Pi]", Larger, Bold]],\[Pi]xpolplot,\[Pi]ypolplot, \[Pi]zpolplot}, 
{Text[Style["\[Sigma]", Larger, Bold]],\[Sigma]xpolplot,\[Sigma]ypolplot, \[Sigma]zpolplot}}]



