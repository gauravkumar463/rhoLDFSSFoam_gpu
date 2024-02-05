// calculate square of speed of sound
// volScalarField c2 = (thermo.gamma()-1)*thermo.Cp()*T;

rpsi = 1/psi;

rpsi_pos = interpolate(rpsi, pos, T.name());
rpsi_neg = interpolate(rpsi, neg, T.name());

p_pos = interpolate(p, pos, p.name());
p_neg = interpolate(p, neg, p.name());

U_pos = interpolate(U, pos, U.name());
U_neg = interpolate(U, neg, U.name());

gammaSf = linearInterpolate(thermo.gamma());

// surfaceScalarField a_pos(sqrt(gammaSf*rpsi_pos));
// surfaceScalarField a_neg(sqrt(gammaSf*rpsi_pos));

rhoHalf = 0.5*((p_pos/rpsi_pos)+(p_neg/rpsi_neg));
aHalf = sqrt(0.5*gammaSf*(rpsi_pos+rpsi_neg));

// M_pos = ((U_pos & mesh.Sf()) - meshPhi)/(aHalf*mesh.magSf());
// M_neg = ((U_neg & mesh.Sf()) - meshPhi)/(aHalf*mesh.magSf());

M_pos = (U_pos & mesh.Sf())/(aHalf*mesh.magSf());
M_neg = (U_neg & mesh.Sf())/(aHalf*mesh.magSf());

Mavg = sqrt(0.5*(sqr(M_pos)+sqr(M_neg)));

Mr2Half = min(scalar(1.0), max(0.5*(magSqr(U_pos)+magSqr(U_neg))/sqr(aHalf),scalar(0.09)));

fHalf =  sqrt((sqr(1.0-Mr2Half)*sqr(Mavg))+(4.0*Mr2Half))/(1.0+Mr2Half);

aHalf *= fHalf;
M_pos /= fHalf;
M_neg /= fHalf;

Mavg = sqrt(0.5*(sqr(M_pos)+sqr(M_neg)));

amaxSf = 0.5*aHalf*(1+Mr2Half)*(Mavg+1)*mesh.magSf();

alpha_pos = 0.5*(1.0 + sign(M_pos));
alpha_neg = 0.5*(1.0 - sign(M_neg));

beta_pos = -0.5*(sign(1.0-mag(M_pos)) + mag(sign(1.0-mag(M_pos))));
beta_neg = -0.5*(sign(1.0-mag(M_neg)) + mag(sign(1.0-mag(M_neg))));
// beta_pos = max(scalar(0),1.0-mag(M_pos));
// beta_neg = max(scalar(0),1.0-mag(M_neg));

M4_pos = 0.25*sqr(M_pos+1.0) + (1.0/8.0)*sqr(sqr(M_pos)-1.0);
M4_neg= -0.25*sqr(M_neg-1.0) - (1.0/8.0)*sqr(sqr(M_neg)-1.0);

Mplus = (alpha_pos*(1+beta_pos)*M_pos) - (beta_pos*M4_pos);
Mminus = (alpha_neg*(1+beta_neg)*M_neg) - (beta_neg*M4_neg);

MHalf = 0.25*beta_pos*beta_pos*sqr(sqrt(0.5*(sqr(M_pos)+sqr(M_neg)))-1.0);

UPlus = aHalf*(Mplus - (MHalf*(1.0-((p_pos-p_neg)
                        /(2*p_pos*Mr2Half*sqr(aHalf/fHalf)/(rpsi_pos))))));
Uminus = aHalf*(Mminus + (MHalf*(1.0+((p_pos-p_neg)
                        /(2*p_neg*Mr2Half*sqr(aHalf/fHalf)/(rpsi_neg))))));

M5_pos = 0.25*sqr(M_pos+1)*(2.0-M_pos) 
                            + (3.0/16.0)*M_pos*sqr(sqr(M_pos)-1.0);
M5_neg = 0.25*sqr(M_neg-1)*(2.0+M_neg)
                            - (3.0/16.0)*M_neg*sqr(sqr(M_neg)-1.0);

p1_pos = alpha_pos*(1+beta_pos) - beta_pos*M5_pos;
p1_neg = alpha_neg*(1+beta_neg) - beta_neg*M5_neg;

pHalf = 0.5*(p_pos+p_neg) + 0.5*(p1_pos-p1_neg)*(p_pos-p_neg) 
                    + rhoHalf*(Mr2Half*sqr(aHalf/fHalf))*(p1_pos+p1_neg-1.0);

H_pos = interpolate(e, pos, e.name()) + interpolate(k_, pos, k_.name())
                + 0.5*magSqr(U_pos) + (pHalf*rpsi_pos/p_pos);

H_neg = interpolate(e, neg, e.name()) + interpolate(k_, neg, k_.name())
                + 0.5*magSqr(U_neg) + (pHalf*rpsi_neg/p_neg);


phi = mesh.magSf()*((p_pos*UPlus/rpsi_pos) + (p_neg*Uminus/rpsi_neg));

rhoFlux = - fvc::div(phi);

rhoUFlux = - fvc::div(mesh.magSf()*((p_pos*U_pos*UPlus/rpsi_pos) 
                                + (p_neg*U_neg*Uminus/rpsi_neg)) 
                                + pHalf*mesh.Sf());

rhoEFlux = -fvc::div(mesh.magSf()*((p_pos*H_pos*UPlus/rpsi_pos) 
                                + (p_neg*H_neg*Uminus/rpsi_neg))); 
                        // + pHalf*meshPhi);