(* ::Package:: *)

ClearAll["Global`*"]
(*-------------------------------------*)
(*COORDINATES*)
(*Cartesian distance*)
dist[x_,y_,z_,x0_,y0_,z0_]:=Sqrt[(x-x0)^2+(y-y0)^2+(z-z0)^2];
(*spherical coordinates transformation*)
rtx[r_,\[Theta]_,\[Phi]_]:=r*Sin[\[Theta]]*Cos[\[Phi]];
rty[r_,\[Theta]_,\[Phi]_]:=r*Sin[\[Theta]]*Sin[\[Phi]];
rtz[r_,\[Theta]_,\[Phi]_]:=r*Cos[\[Theta]];
(*Distance, but expressed with spherical variables*)
distr[r_,\[Theta]_,\[Phi]_,r0_,\[Theta]0_,\[Phi]0_]:=dist[rtx[r,\[Theta],\[Phi]],rty[r,\[Theta],\[Phi]],rtz[r,\[Theta],\[Phi]],rtx[r0,\[Theta]0,\[Phi]0],rty[r0,\[Theta]0,\[Phi]0],rtz[r0,\[Theta]0,\[Phi]0]];
(*--------------------------------------------*)
(*FREE SPACE*)
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
Print["dist (0-order) = ",distr0[r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0]]
Print["dist (1st-order) = ",distr1[r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0]]
(*At the exponent, the distance is expanded to 1st order to take into account that \[Lambda] > r, whereas at the denominator it is expanded to 0th order because r0 >> r.*)
G0fr\[Theta]\[Phi][k_,r_,\[Theta]_,\[Phi]_,r0_,\[Theta]0_,\[Phi]0_]:=(\[Omega]*\[Mu]0*\[Mu]/(4*\[Pi]))*Exp[I*k*distr1[r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0]]/distr0[r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0];
G0fxyzr0[kk_,xx_,yy_,zz_,xx0_,yy0_,zz0_,rr0_]:=Evaluate[FullSimplify[G0fr\[Theta]\[Phi][k,r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0]/.{Cos[\[Phi]]->x/(r*Sin[\[Theta]]), Sin[\[Phi]]->y/(r*Sin[\[Theta]]), Cos[\[Phi]0]->x0/(r0*Sin[\[Theta]0]), Sin[\[Phi]0]->y0/(r0*Sin[\[Theta]0]), Cos[\[Theta]]->z/r, Cos[\[Theta]0]->z0/r0}]/.{k->kk,x->xx,y->yy,z->zz,x0->xx0,y0->yy0,z0->zz0,r0->rr0}];
Print["G0fxyzr0 (far-field, r0) = ", G0fxyzr0[k,x,y,z,x0,y0,z0,r0]]
G0fxyz[k_,\[Rho]_,\[Psi]_,x_,y_,z_,x0_,y0_,z0_]:=Evaluate[FullSimplify[(G0fxyzr0[k,x,y,z,x0,y0,z0,r0]+\[Rho]*G0fxyzr0[k,x,y,z,-x0,-y0,-z0,r0]*Exp[I*\[Psi]])/.{r0->Sqrt[x0^2+y0^2+z0^2], r->Sqrt[x^2+y^2+z^2]}]];
G0fxr[kk_,xx_,yy_,zz_,rr0_,\[Theta]\[Theta]0_,\[Phi]\[Phi]0_]:=Evaluate[FullSimplify[G0fr\[Theta]\[Phi][k,r,\[Theta],\[Phi],r0,\[Theta]0,\[Phi]0]/.{Cos[\[Phi]]->x/(r*Sin[\[Theta]]), Sin[\[Phi]]->y/(r*Sin[\[Theta]]), Cos[\[Theta]]->z/r}]/.{k->kk,x->xx,y->yy,z->zz,r0->rr0,\[Theta]0->\[Theta]\[Theta]0,\[Phi]0->\[Phi]\[Phi]0}];
Print[" G0fxyz free space(far-field) = ", G0fxyz[k,0,\[Psi],x,y,z,x0,y0,z0]]
Print[" G0fxr (far-field) = ", G0fxr[k,x,y,z,r0,\[Theta]0,\[Phi]0]]
Print["-------------------"]
Print[" G0fxyz mirror (far-field) = ", G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]]
(*---------------------------------------------*)
(*Electric field*)



nxxyz = {1,0,0};
nyxyz = {0,1,0};
nzxyz = {0,0,1};

Gtx = FullSimplify[Div[G0fxyz[k,\[Rho],\[Psi],0,0,0,x0,y0,z0]*nxxyz,{x0,y0,z0},"Cartesian"]]
Gtx = Simplify[Gtx/.{(x0^2+y0^2+z0^2)->r0^2}, Assumptions->{r0>0}]/.{I+k*r0->k*r0}
Gtx = Grad[Gtx/.{r0->Sqrt[x0^2+y0^2+z0^2]},{x0,y0,z0},"Cartesian"]
Gtx = Simplify[Gtx/.{(x0^2+y0^2+z0^2)->r0^2}, Assumptions->{r0>0}]

Gtx = FullSimplify[Simplify[G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]*nxxyz+(1/k^2)*Grad[Div[G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]*nxxyz,{x,y,z},"Cartesian"],{x,y,z},"Cartesian"]]/G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]]
Gty = FullSimplify[Simplify[G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]*nyxyz+(1/k^2)*Grad[Div[G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]*nyxyz,{x,y,z},"Cartesian"],{x,y,z},"Cartesian"]]/G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]]
Gtz = FullSimplify[Simplify[G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]*nzxyz+(1/k^2)*Grad[Div[G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]*nzxyz,{x,y,z},"Cartesian"],{x,y,z},"Cartesian"]]/G0fxyz[k,\[Rho],\[Psi],x,y,z,x0,y0,z0]]

Gt = {Gtx, Gty, Gtz}

MatrixForm[FullSimplify[Gt/.{x0->r0*Sin[\[Theta]0]*Cos[\[Phi]0],y0->r0*Sin[\[Theta]0]*Sin[\[Phi]0],z0->r0*Cos[\[Theta]0]}]]


G0m[k_,x_,y_,z_,x0_,y0_,z0_,\[Phi]0_]:=G0f[k,x,y,z,x0,y0,z0]+G0f[k,x,y,z,-x0,-y0,-z0]*Exp[I*\[Phi]0];
Print[" G0m = ", G0m[k,x,y,z,x0,y0,z0,\[Phi]0]]

TrigExpand[Cos[2*\[Theta]]]

(*G0t[k_,x_,y_,z_]:=Div;*)



