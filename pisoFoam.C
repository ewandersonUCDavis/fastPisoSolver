/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    pisoFoam

Description
    Transient solver for incompressible flow.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "horizontalAxisWindTurbinesFAST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

extern "C" 
{
  void fastinit_( float& , int& );
  void fastread_( float*, float*, float*);
  void fastrun_( );
  void fastgetbldpos_( float*, float*, float*);
  void fastgetbldforce_(float*, float*, float*);
  void fastend_( );
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // initialize FAST
    Info << "Number of Turbs: " << turbfast.turbNum << endl;
    float tstep = runTime.deltaT().value();
    for(int turbNo=0; turbNo<turbfast.turbNum; turbNo++)
    {
      if(Pstream::myProcNo() == turbNo)
      {
        fastinit_(tstep, turbNo);
        fastgetbldpos_(turbfast.bldptx[turbNo], turbfast.bldpty[turbNo], turbfast.bldptz[turbNo]);
      }
      turbfast.getBldPos(turbNo);
    }


    while (runTime.loop())
    {
      Info<< "Time = " << runTime.timeName() << nl << endl;

      #include "readPISOControls.H"
      #include "CourantNo.H"

    // Pressure-velocity PISO corrector
    {

      for(int turbNo=0; turbNo<turbfast.turbNum; turbNo++)
      {
        turbfast.getWndVec(turbNo);
        if(Pstream::myProcNo() == turbNo) 
        {
          fastread_(turbfast.uin[turbNo], turbfast.vin[turbNo], turbfast.win[turbNo]);
          fastrun_();
          fastgetbldpos_(turbfast.bldptx[turbNo], turbfast.bldpty[turbNo], turbfast.bldptz[turbNo]); 
          fastgetbldforce_(turbfast.bldfx[turbNo], turbfast.bldfy[turbNo], turbfast.bldfz[turbNo]);
        }
        turbfast.computeBodyForce(turbNo);
      }

      // Momentum predictor
      fvVectorMatrix UEqn
      (
        fvm::ddt(U)
           + fvm::div(phi, U)
           + turbulence->divDevReff(U) - turbfast.force() 
      );

      UEqn.relax();

      if (momentumPredictor)
      {
        solve(UEqn == -fvc::grad(p));
      }

      // --- PISO loop
      for (int corr=0; corr<nCorr; corr++)
      {
        volScalarField rUA = 1.0/UEqn.A();

        U = rUA*UEqn.H();
        phi = (fvc::interpolate(U) & mesh.Sf())
            + fvc::ddtPhiCorr(rUA, U, phi);

        adjustPhi(phi, U, p);

        // Non-orthogonal pressure corrector loop
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
          // Pressure corrector
          fvScalarMatrix pEqn
          (
            fvm::laplacian(rUA, p) == fvc::div(phi)
          );

          pEqn.setReference(pRefCell, pRefValue);

          if(corr == nCorr-1 && nonOrth == nNonOrthCorr)
          {
            pEqn.solve(mesh.solver("pFinal"));
          }
          else
          {
            pEqn.solve();
          }

          if (nonOrth == nNonOrthCorr)
          {
            phi -= pEqn.flux();
          }
        }

        #include "continuityErrs.H"

        U -= rUA*fvc::grad(p);
        U.correctBoundaryConditions();
      }
    }

    turbulence->correct();

    runTime.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl; 
    }

    // terminate FAST
    for(int turbNo=0; turbNo<turbfast.turbNum; turbNo++)
    {
      if(Pstream::myProcNo() == turbNo)
      {
        fastend_();
      }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
