p += runTime.deltaT()*pFlux;

U += runTime.deltaT()*UFlux;

e += runTime.deltaT()*eFlux;

if(stage>0){
	p = ((1-beta[stage])*tp)+(beta[stage]*p);
	U = ((1-beta[stage])*tU)+(beta[stage]*U);
	e = ((1-beta[stage])*te)+(beta[stage]*e);
}

for (Foam::direction dir=0; dir<directions.nComponents; dir++)
{
    if (directions[dir] == -1){
       U.replace(dir, 0.0*U.component(dir));
    }
}

p.correctBoundaryConditions();
U.correctBoundaryConditions();
e.correctBoundaryConditions();

thermo.correct();

rho.dimensionedInternalField() = 
	psi.dimensionedInternalField()*p.dimensionedInternalField();
	
rho.boundaryField() == psi.boundaryField()*p.boundaryField();

q.dimensionedInternalField() = 
	0.5*magSqr(U.dimensionedInternalField()) 
		+ k_.dimensionedInternalField();

q.boundaryField() == 0.5*magSqr(U.boundaryField()) 
						+ k_.boundaryField();
