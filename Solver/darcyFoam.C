/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Application
    darcyFoam

Description
    Solves for pressure & velocity fields and other variables 
    for drug delivery in brain tumor.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "IOMRFZoneList.H"
#include "IOporosityModelList.H"     

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"

simpleControl simple(mesh);
#include "CourantNo.H"

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

 // Finding Infusion Site
  
 scalar maxInfSite = -GREAT;
 label maxCellI;

 forAll(infSite, cellI)
 {
     if (maxInfSite < infSite[cellI])
     {
         maxInfSite = infSite[cellI];
         maxCellI = cellI;
     }

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

 Info << "\nCalculating IFP, IFV, Drug Bioavailability and Microvasculature Density Variation Scale\n" << endl;

    while (simple.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            solve
            (
                fvm::laplacian(K, p)
                - (1 - infSite) * fvm::Sp((ktrans / ktransavg) * (hc * sapuv), p)
                - (1 - infSite) * fvm::Sp((lfc), p)
                ==
                - (1 - infSite) * ((ktrans / ktransavg) * (hc * sapuv) * (Peff))
                - (1 - infSite) * (lfc * lp)
                - infSite * (Q / Vinf)
            );
        }

        U = -K * fvc::grad(p);

        //Cytotoxic Drug Transport Equations

        w = por * (1 + Kecs) + (Vcell / Vtissue) * Pics_ecs * (1 + Kics) + (1 - por - (Vcell / Vtissue)) * Pcm_ecs;
        kb = psi * ktrans;

        Ustar = (por / w) * U;
        phi = linearInterpolate(Ustar) & mesh.Sf();

        Dstar = (por / w) * Df;
        kstar = (por * kb + (por + Pics_ecs * (Vcell / Vtissue)) * ke) / w;

        solve
        (
              fvm::ddt(Cf)
            - fvm::laplacian(Dstar, Cf)
            + fvm::div(phi, Cf)
            + fvm::Sp(kstar, Cf)
        );

        Cf[maxCellI] = 7.0e-6 + Cf[maxCellI];   //For Continuous Infusion

        // Anti-angiogenic Drig Transport equations

       phi = linearInterpolate(U) & mesh.Sf();

        solve
        (
            fvm::ddt(Caa)           
          - fvm::laplacian(Daa * por, Caa)
          + fvm::div(phi, Caa)
          + fvm::Sp(kaae, Caa)
        );

        Caa[maxCellI] = 7.25e-5 + Caa[maxCellI];   //For Continuous Infusion

        // equation 7
        solve
        (
            fvm::ddt(psi)
          - fvm::Sp(alpha, psi)
          - fvm::Sp(beta * psi, psi)
          - fvm::Sp(gamma * psi * psi, psi)
          + fvm::Sp(kak * (Caa / C_infSite), psi)
         );

        // Plasma Volume Fraction and Porosity Update

        pvf = (psi * psi) * pvf;

        por = 1 - pvf - (Vcell / Vtissue);

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n" << endl;

    return 0;
}
  
// ************************************************************************* //
