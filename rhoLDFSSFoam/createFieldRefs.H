volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();

volScalarField rpsi
(
    IOobject
    (
        "rpsi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    1/psi
);

const Foam::Vector<label>& directions = mesh.solutionD();

surfaceScalarField rpsi_pos(interpolate(rpsi, pos, T.name()));
surfaceScalarField rpsi_neg(interpolate(rpsi, pos, T.name()));

surfaceScalarField p_pos(interpolate(p, pos, p.name()));
surfaceScalarField p_neg(interpolate(p, neg, p.name()));

surfaceVectorField U_pos(interpolate(U, pos, U.name()));
surfaceVectorField U_neg(interpolate(U, neg, U.name()));


surfaceScalarField gammaSf(linearInterpolate(thermo.gamma()));

surfaceScalarField rhoHalf = p_pos/rpsi_pos;
surfaceScalarField aHalf = sqrt(rpsi_pos);

surfaceScalarField M_pos = pos;
surfaceScalarField M_neg = pos;

surfaceScalarField Mavg = pos;

surfaceScalarField Mr2Half = pos;

surfaceScalarField fHalf = pos;

surfaceScalarField alpha_pos = pos;
surfaceScalarField alpha_neg = pos;

surfaceScalarField beta_pos = pos;
surfaceScalarField beta_neg = pos;

surfaceScalarField M4_pos = pos;
surfaceScalarField M4_neg= pos;

surfaceScalarField Mplus = pos;
surfaceScalarField Mminus = pos;

surfaceScalarField MHalf = pos;

surfaceScalarField UPlus = aHalf;
surfaceScalarField Uminus = aHalf;

surfaceScalarField M5_pos = pos;
surfaceScalarField M5_neg = pos;

surfaceScalarField p1_pos = pos;
surfaceScalarField p1_neg = neg;

surfaceScalarField pHalf = p_pos;

surfaceScalarField H_pos = rpsi_pos;
surfaceScalarField H_neg = rpsi_neg;

/*volVectorField fbody
(
    IOobject
    (
        "fbody",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    acceleration
);*/
