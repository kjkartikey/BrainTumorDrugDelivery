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
    Solves for IFP, IFV, Drug Bioavailability and Microvasulature density
    variation scale.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "IOMRFZoneList.H"
#include "IOporosityModelList.H"     

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

int main(int argc, char* argv[])
{
#include "setRootCase.H"

#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"

simpleControl simple(mesh);

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

Info << "\nCalculating IFP, IFV, Drug Bioavailability and Microvasulature density variation scale.\n" << endl;

#include "CourantNo.H"

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

        phi = linearInterpolate(U) & mesh.Sf();

        /*forAll(Z, i)
        {
            Z[i] = mag(U[i]);
        }

        Cp.value() = cpValue->value(mesh.time().timeOutputValue());

        solve
        (
            fvm::ddt(C)
            + (1 / por) * fvm::div(phi, C)
            - fvm::laplacian(D, C)
            + fvm::Sp((ktrans / por), C)
            + fvm::Sp((lfc / por) * (p - lp), C)
            ==
            (ktrans * Cp)

        );*/

        //J = fvc::laplacian(D,C);
        //      I = (ktrans*(Cp-((1/por)*C)));
         //     L = (1/por)*fvc::div(phi, C);

        // Drug Bioavailabilty
        solve
        (
            fvm::ddt(Caa)           
          - (1 - infSite) * fvm::laplacian(Daa * por, Caa)
          + (1 - infSite) * fvm::div(phi, Caa)
          + (1 - infSite) * fvm::Sp(kaae, Caa)
        );

        // Drug Bioavailabilty
        solve
        (
            fvm::ddt(vsmd)
          - fvm::Sp(alpha, vsmd)
          - fvm::Sp(beta * vsmd, vsmd)
          - fvm::Sp(gamma * vsmd * vsmd, vsmd)
          + fvm::Sp(kak * (Caa / C_infSite), vsmd)
         );
      
      
      //* * * * * * * * * * * * * * * * * * * * * * * Properties Update * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
      
      pvf = (vsmd * vsmd) * pvf;
      
      por = 1 - pvf - (Vcell / Vtissue);
      
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
      
       runTime.write();

       Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info << "End\n" << endl;

    return 0;
}



// ************************************************************************* //
