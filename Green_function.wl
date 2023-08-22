(* ::Package:: *)

ClearAll["Global`*"]
(*-------------------------------------*)
(*COORDINATES*)
(*Cartesian distance*)
dist[x_,y_,z_,x0_,y0_,z0_]:=Sqrt[(x-x0)^2+(y-y0)^2+(z-z0)^2];
(*cartesian to spherical coordinates transformation direct and inverse.
These are used often because distances and gradient for this problem have easier formulas in Cartesia coordinates,
but the Series expansion depends on the radii.*)
varchangextor[x_,y_,z_,r_,\[Theta]_,\[Phi]_]:={x->r*Sin[\[Theta]]*Cos[\[Phi]], y->r*Sin[\[Theta]]*Sin[\[Phi]], z->r*Cos[\[Theta]]};
varchange\[Theta]\[Phi]tox[r_,\[Theta]_,\[Phi]_,x_,y_,z_]:={Cos[\[Phi]]->x/(r*Sin[\[Theta]]), Sin[\[Phi]]->y/(r*Sin[\[Theta]]), Cos[\[Theta]]->z/r};
varchangeSin\[Theta]tox[\[Theta]_,x_,y_,z_]:={Sin[\[Theta]]->Sqrt[x^2+y^2]/Sqrt[x^2+y^2+z^2], Csc[\[Theta]]->Sqrt[x^2+y^2+z^2]/Sqrt[x^2+y^2]}
varchangertox[r_,x_,y_,z_]:={r->Sqrt[x^2+y^2+z^2]};
(*Distance, but expressed with spherical variables*)
distr[r_,\[Theta]_,\[Phi]_,r0_,\[Theta]0_,\[Phi]0_]:=(dist[x,y,z,x0,y0,z0]/.varchangextor[x,y,z,r,\[Theta],\[Phi]])/.varchangextor[x0,y0,z0,r0,\[Theta]0,\[Phi]0];
Eikror[k_,r_]:=\[Mu]*\[Mu]0*\[Omega]*Exp[I*k*r]/(4*\[Pi]*r);
(*--------------------------------------------*)
(*Scalar Green's function for the vector potential A*)
(*r0 = location of the detector*)
(*r = location of the object*)
G0f[k_,x_,y_,z_,x0_,y0_,z0_]:=(\[Omega]*\[Mu]0*\[Mu]/(4*\[Pi]))*Exp[I*k*dist[x,y,z,x0,y0,z0]]/dist[x,y,z,x0,y0,z0];
Print[" G0f (no approximation) = ", G0f[k,x,y,z,x0,y0,z0]]
(*far-field approximation leverages the expansion of the radial part assuming r0 >> \[Lambda] > r.
Under this condition, it is possible to expand the distance in series of r and retain only the first terms.*)
distr0[r_,\[Theta]_,\[Phi]_,r0_,\[Theta]0_,\[Phi]0_]:=Evaluate[Normal[Simplify[Series[distr[rr,\[Theta]\[Theta],\[Phi]\[Phi],rr0,\[Theta]\[Theta]0,\[Phi]\[Phi]0],{rr,0,0}], Assumptions->{rr0>0}]]/.{rr->r,\[Theta]\[Theta]->\[Theta],\[Phi]\[Phi]->\[Phi],rr0->r0,\[Theta]\[Theta]0->\[Theta]0,\[Phi]\[Phi]0->\[Phi]0}];
distr1[r_,\[Theta]_,\[Phi]_,r0_,\[Theta]0_,\[Phi]0_]:=Evaluate[TrigExpand[Normal[Simplify[Series[distr[rr,\[Theta]\[Theta],\[Phi]\[Phi],rr0,\[Theta]\[Theta]0,\[Phi]\[Phi]0],{rr,0,1}], Assumptions->{rr0>0, rr>0}]]]/.{rr->r,\[Theta]\[Theta]->\[Theta],\[Phi]\[Phi]->\[Phi],rr0->r0,\[Theta]\[Theta]0->\[Theta]0,\[Phi]\[Phi]0->\[Phi]0}];
distr2[r_,\[Theta]_,\[Phi]_,r0_,\[Theta]0_,\[Phi]0_]:=Evaluate[TrigExpand[Normal[Simplify[Series[distr[rr,\[Theta]\[Theta],\[Phi]\[Phi],rr0,\[Theta]\[Theta]0,\[Phi]\[Phi]0],{rr,0,2}], Assumptions->{rr0>0, rr>0}]]]/.{rr->r,\[Theta]\[Theta]->\[Theta],\[Phi]\[Phi]->\[Phi],rr0->r0,\[Theta]\[Theta]0->\[Theta]0,\[Phi]\[Phi]0->\[Phi]0}];
Print["distance (0-order) = ",distr0[r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0]]
Print["distance (1st-order) = ",distr1[r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0]]
(*At the exponent, the distance is expanded to 1st order to take into account that \[Lambda] > r, whereas at the denominator it is expanded to 0th order because r0 >> r.*)
G0fr\[Theta]\[Phi][k_,r_,\[Theta]_,\[Phi]_,r0_,\[Theta]0_,\[Phi]0_]:=(\[Omega]*\[Mu]0*\[Mu]/(4*\[Pi]))*Exp[I*k*distr1[r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0]]/distr0[r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0];
(*Changing to Cartesian coordinates because the mirror is easier to introduce. The variable r0 is left for redeability*)
G0fxyzr0tmp[kk_,xx_,yy_,zz_,xx0_,yy0_,zz0_,rr0_]:=Evaluate[FullSimplify[(G0fr\[Theta]\[Phi][k,r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0]/.varchange\[Theta]\[Phi]tox[r,\[Theta],\[Phi],x,y,z])/.varchange\[Theta]\[Phi]tox[r0,\[Theta]0,\[Phi]0,x0,y0,z0]]/.{k->kk,x->xx,y->yy,z->zz,x0->xx0,y0->yy0,z0->zz0,r0->rr0}];
(*Introducing the mirror*)
G0fxyzr0[k_,\[Rho]_,\[Psi]_,x_,y_,z_,x0_,y0_,z0_,r0_]:=Evaluate[Simplify[G0fxyzr0tmp[k,x,y,z,x0,y0,z0,r0]]-\[Rho]*G0fxyzr0tmp[k,x,y,z,-x0,-y0,-z0,r0]*Exp[I*\[Psi]]];
G0fxyz[k_,\[Rho]_,\[Psi]_,x_,y_,z_,x0_,y0_,z0_]:=Evaluate[FullSimplify[G0fxyzr0[k,\[Rho],\[Psi],x,y,z,x0,y0,z0,r0]/.varchangertox[r0,x0,y0,z0]]];
G0fr\[Theta]\[Phi][k_,\[Rho]_,\[Psi]_,x_,y_,z_,r0_,\[Theta]0_,\[Phi]0_]:=Evaluate[Simplify[G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]/.varchangextor[x0,y0,z0,r0,\[Theta]0,\[Phi]0], Assumptions->{r0>0}]];
G0fr\[Theta]\[Phi]red[k_,\[Rho]_,\[Psi]_,x_,y_,z_,r0_,\[Theta]0_,\[Phi]0_]:=Evaluate[Simplify[G0fr\[Theta]\[Phi][k,\[Rho],\[Psi],x,y,z,r0,\[Theta]0,\[Phi]0]/Eikror[k,r0], Assumptions->{r0>0}]];
Print["-------------------"]
Print["Far-field scalar Green's function valid for the each component of A in Lorentz's gauge"]
Print["x,y,z = coordinate of the object (to be integrated to)"]
Print["x0,y0,z0 = coordinate at the detector"]
Print[" general = ", G0fxyzr0[k,\[Rho],\[Psi],x,y,z,x0,y0,z0,r0]]
Print[" free space = ", Eikror[k,r0], Simplify[G0fxyzr0[k,0,0,x,y,z,x0,y0,z0,r0]/Eikror[k,r0]]]
Print[" mirror suppression = ", Eikror[k,r0], "(",Simplify[ExpToTrig[G0fr\[Theta]\[Phi]red[k,1,0,x,y,z,r0,\[Theta]0,\[Phi]0]]/.varchange\[Theta]\[Phi]tox[r0,\[Theta]0,\[Phi]0,x0,y0,z0]], ")"]


derfactor[k_]:=I*k;
(*Calculating the Green's function for the Electric field. There are two derivatives to do.
In each derivative step, only the leading term is retained.
To retain leading term, two actions are taken:
1) x=0, y=0, z=0 at evaluating the derivatives. This is equivalent of retaining the order 0 in the exapansion r/r0.
This can be done because the detector is far away from the sample r0>>r.
2) Only first term in the Series expansion of \[Alpha]=1/r0\[Rule]0 is retained*)
GGt = G0fxyz[k,\[Rho],\[Psi],0,0,0,x0,y0,z0]*IdentityMatrix[3];
GGt = FullSimplify[Div[GGt,{x0,y0,z0},"Cartesian"]];
(*Going to radial coordinates, exapnd assuming r0 big and come back*)
GGt = Simplify[FullSimplify[(GGt/.varchangextor[x0,y0,z0,r0,\[Theta]0,\[Phi]0])/.{r0->1/\[Alpha]}],Assumptions->{\[Alpha]>0}];
GGt = Normal[Series[GGt,{\[Alpha],0,1}]]/.{\[Alpha]->1/r0};
GGtxyz = (GGt/.varchange\[Theta]\[Phi]tox[r0,\[Theta]0,\[Phi]0,x0,y0,z0])/.varchangertox[r0,x0,y0,z0];
GGtvet = FullSimplify[GGtxyz/(derfactor[k]*G0fxyz[k,\[Rho],\[Psi],0,0,0,x0,y0,z0])];
GGtvetr = Simplify[GGtvet/.varchangextor[x0,y0,z0,r0,\[Theta]0,\[Phi]0], Assumptions->{r0>0}];
Print["GGt - first derivative only = ", derfactor[k]*Simplify[G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]/.varchangextor[x0,y0,z0,r0,\[Theta]0,\[Phi]0], Assumptions->{r0>0}], MatrixForm[GGtvetr]]
Gt = Grad[GGtxyz,{x0,y0,z0},"Cartesian"];
Gt = Simplify[(Gt/.varchangextor[x0,y0,z0,r0,\[Theta]0,\[Phi]0])/.{r0->1/\[Alpha]},Assumptions->{\[Alpha]>0}];
Gt = TrigExpand[Normal[Series[Gt,{\[Alpha],0,1}]]]/.{\[Alpha]->1/r0};
Gtxyz = FullSimplify[((Gt/.varchange\[Theta]\[Phi]tox[r0,\[Theta]0,\[Phi]0,x0,y0,z0])/.varchangertox[r0,x0,y0,z0])/.varchangeSin\[Theta]tox[\[Theta]0,x0,y0,z0]];
(*Putting constant and derivative term togheter*)
Gtxyz = G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]*IdentityMatrix[3] + Gtxyz/derfactor[k]^2;
Gtmat = FullSimplify[Gtxyz/G0fxyz[k,\[Rho],\[Psi],0,0,0,x0,y0,z0]/.varchangextor[x0,y0,z0,r0,\[Theta]0,\[Phi]0]];
Print["-------------------"]
Print["Far-field Green's tensor valid for the Electric field"]
Print["x,y,z = coordinate of the object (to be integrated to)"]
Print["x0,y0,z0 = coordinate at the detector"]
Print["general = ", G0fxyzr0[k,\[Rho],\[Psi],x,y,z,x0,y0,z0,r0], MatrixForm[Gtmat]]
Print["free space = ", Eikror[k,r0], Simplify[G0fxyzr0[k,0,0,x,y,z,x0,y0,z0,r0]/Eikror[k,r0]], MatrixForm[Gtmat]]
Print["mirror suppression = ", Eikror[k,r0], "(",Simplify[ExpToTrig[G0fr\[Theta]\[Phi]red[k,1,0,x,y,z,r0,\[Theta]0,\[Phi]0]]/.varchange\[Theta]\[Phi]tox[r0,\[Theta]0,\[Phi]0,x0,y0,z0]], ")", MatrixForm[Gtmat]]
