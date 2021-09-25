/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    darcyFoam

Description

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

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info << "\nCalculating PRESSURE distribution\n" << endl;

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

              // equation 6
        solve
        (
            fvm::ddt(Caa)           
          - fvm::laplacian(Daa, Caa)
          - fvm::div(phi, Caa)
          + fvm::Sp(kaae, Caa)
        );

        // equation 7
        solve
        (
            fvm::ddt(pvf)
          - fvm::Sp(alpha, pvf)
          - fvm::Sp(beta * pvf, pvf)
          - fvm::Sp(gamma * pvf * pvf, pvf)
          + fvm::Sp(kak * Caa, pvf)
         );
      
        //porosity update
       por = por * ((Vtissue - Vcell - 2 * constant::mathematical::pi * pvf*pvf * r1*r1 * l) / (Vtissue - Vcell - constant::mathematical::pi * r1*r1 * l));

            runTime.write();

            Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
    }

    Info << "End\n" << endl;

    return 0;
}



// ************************************************************************* //
