muEff = turbulence->muEff();
tauMC = muEff*dev2(Foam::T(fvc::grad(U)));

solve
(
	fvm::ddt(rho, U) - fvc::ddt(rho, U)
  - fvm::laplacian(muEff, U)
  - fvc::div(tauMC)
);
U.correctBoundaryConditions();

e.dimensionedInternalField() -= 
    0.5*magSqr(U.dimensionedInternalField()) 
    + k_.dimensionedInternalField() - q.dimensionedInternalField();

//e.ref() -= 0.5*magSqr(U.ref()) - q.ref();
e.correctBoundaryConditions();

solve
(
	fvm::ddt(rho,e) - fvc::ddt(rho, e)
  - fvm::laplacian(turbulence->alphaEff(), e)
  - fvc::div
    (
        (fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
          + (mesh.Sf() & fvc::interpolate(tauMC))
          // + fvc::dotInterpolate(mesh.Sf(), tauMC)
        ) & linearInterpolate(U)
    )
);

e.correctBoundaryConditions();

thermo.correct(); 

p.dimensionedInternalField() = rho.dimensionedInternalField()
            /psi.dimensionedInternalField();

p.correctBoundaryConditions();
rho.boundaryField() == psi.boundaryField()*p.boundaryField();
