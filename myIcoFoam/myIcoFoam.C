/*---------------------------------------------------------------------------*\
   =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    icoFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}}
          + \div \left( \vec{U} \vec{U} \right)
          - \div \left(\nu \grad \vec{U} \right)
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

#include <stdio.h>
#include "time.h"

#define USE_DSPACES

#ifdef USE_DSPACES
    #include "ds_adaptor.h"

    // write data from veolscalar field into dspaces send buffer 
    // not sure if we need special buffer
    void mydump(volScalarField &p, float *tmp_buffer){
        //int i;
        float * tmp = tmp_buffer;
        forAll(p , i)
         {
             *(tmp++) = p[i];
         }
    }
    void mydump(volVectorField &U, float *tmp_buffer){
        float *tmp = tmp_buffer;
        forAll(U, i)
        {
            *tmp = U[i].x();
            *(tmp+1) = U[i].y();
            *(tmp+2) = U[i].z();

            //Info<<" point "<< i << ": "<< tmp[0] <<" " << tmp[1]<< " "<<tmp[2]<< endl;
            tmp+=3;
        }
    }

    void mydump(volScalarField &p, float *tmp_buffer,int points_per_dim, int bottom_most, int right_most){
        int ii, jj;
        ii = 0;
        jj = 0;

        float * tmp = tmp_buffer;
        forAll(p , i)
         {
             // if this isn't a block in the boundary, don't write right/bottom most 
             if((bottom_most !=1&&ii == points_per_dim -1) || (right_most !=1&& jj==points_per_dim-1)){
                     ;
              }
             else{
                *(tmp++) = p[i];
            }

             jj+=1;

             // control loop
             if(jj%points_per_dim == 0){
                 jj =0;
                 ii +=1;
             }
         }
    }
    void mydump(volVectorField &U, float *tmp_buffer,int points_per_dim, int bottom_most, int right_most){

        int ii, jj;
        ii = 0;
        jj = 0;
        float *tmp = tmp_buffer;
        forAll(U, i)
        {
            // don't writte right/bottom boundary for some block
            if((bottom_most !=1&&ii == points_per_dim -1) || (right_most !=1&& jj==points_per_dim-1)){
                     ;
              }
            else{
            *tmp = U[i].x();
            *(tmp+1) = U[i].y();
            *(tmp+2) = U[i].z();

            //Info<<" point "<< i << ": "<< tmp[0] <<" " << tmp[1]<< " "<<tmp[2]<< endl;
            tmp+=3;
            }

             jj+=1;

             // control loop
             if(jj%points_per_dim == 0){
                 jj =0;
                 ii +=1;
             }
        }
    }
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    

    // get size and version number
    int num_points = p.size();
    int timestep = 0;
    Info<< "size is \n" << num_points << endl;

    // how many steps
    int count;

#ifdef USE_DSPACES
    // init dspaces
    // dataspaces preparation
        int nprocs, rank, ret;
        MPI_Comm gcomm;


        // MPI communicator
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Barrier(MPI_COMM_WORLD);
        gcomm = MPI_COMM_WORLD;

    // check the size and get
        int procs_per_dim = (int)std::sqrt(nprocs);
        int points_per_dim = (int)std::sqrt(num_points);

        // cartisian position
        int proc_row_pos = rank/procs_per_dim;
        int proc_col_pos = rank%procs_per_dim;

        // whether its the last proc in each dim
        int flag_bottom_most = 0;
        int flag_right_most = 0;

        // icofoam size
        int case_length = CASE_LENGTH;

        //x_min,y_min,z_min,x_max_y_max_z_max
        int bounds[6] = {0};
        bounds[0]=(case_length)*proc_row_pos;
        bounds[1]=(case_length)*proc_col_pos;

        // don't overlap
        if(proc_col_pos == procs_per_dim-1){
            bounds[4]=(case_length)*(proc_col_pos+1);
            flag_right_most = 1;
        }
        else{
            bounds[4]=(case_length)*(proc_col_pos+1) - 1;
        }

        // we want 1 line overlap (consistent with dividing in consumer)
        if(proc_row_pos == procs_per_dim-1){
            bounds[3]=case_length*(proc_row_pos+1) ;
            flag_bottom_most = 1;
        }
        else{
            bounds[3]=case_length*(proc_row_pos+1)-1 ;
        }


        // Initalize DataSpaces
        // # of Peers, Application ID, ptr MPI comm, additional parameters
        // # Peers: Number of connecting clients to the DS server
        // Application ID: Unique idenitifier (integer) for application
        // Pointer to the MPI Communicator, allows DS Layer to use MPI barrier func
        // Addt'l parameters: Placeholder for future arguments, currently NULL.
        char msg[STRING_LENGTH];
        printf("trying init dspaces for %d process\n", nprocs);
        ret = dspaces_init(nprocs, 1, &gcomm, NULL);

        printf("dspaces init successfuly \n");

        if(ret == 0){
            sprintf(msg, "dataspaces init successfully");
            my_message(msg, rank, LOG_CRITICAL);
        }else{
            sprintf(msg, "dataspaces init error");
            my_message(msg, rank, LOG_CRITICAL);
            exit(-1);
        }


        // time 
        double time_comm_vel;
        
        char var_name_vel[STRING_LENGTH];
        char var_name_pres[STRING_LENGTH];
        sprintf(var_name_vel, "VEL");
        sprintf(var_name_pres, "PRES");


        // data layout
        uint64_t gdims_raw[3] = {POINTS_SIDE, POINTS_SIDE,1};
        dspaces_define_gdim(var_name_vel, 3,gdims_raw);
        dspaces_define_gdim(var_name_pres, 3,gdims_raw);

        /*
        char var_name_vel_2[STRING_LENGTH];
        char var_name_pres_2[STRING_LENGTH];
        sprintf(var_name_vel_2, "VEL_2");
        sprintf(var_name_pres_2, "PRES_2");
        */
        
        // prepare space
        float * vel_data = (float *)malloc(num_points* sizeof(float)*3);
        if(vel_data == NULL){
            perror("vel data allocated error");
            exit(-1);
        }
        float * pres_data = (float *)malloc(num_points* sizeof(float));

        if(pres_data == NULL){
            perror("pres data allocated error");
            exit(-1);
        }


// end of dspaces preparation
#endif

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        sleep(2);
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        // if the setting is 'no write', that field won't be written to files
        runTime.write();

        /*
         * also write pressure info
         */

#ifdef USE_DSPACES
    // if use  dataspces, write the correct 
        //put_vel_buffer(timestep, NULL,rank, &gcomm, &vel_data, &time_comm_vel);
        // dump data from U and P into buffer
        printf("points_per_dim %d, flag_bottom_most %d, flag_right_most %d\n", points_per_dim, flag_bottom_most, flag_right_most);
        mydump(p, pres_data, points_per_dim,flag_bottom_most, flag_right_most);
        mydump(U, vel_data, points_per_dim,flag_bottom_most, flag_right_most);

        /*
        Info<<" first data, address"<< vel_data << ": "<< vel_data[0] <<" " << vel_data[1]<< " "<<vel_data[2]<< endl;
       printf(" first data, address %p: %f %f %f\n", vel_data, vel_data[0], vel_data[1], vel_data[2]);
       */
        //if(count%interval == 0){

        put_raw_buffer(timestep, bounds, &num_points ,rank, &gcomm, var_name_vel,  &vel_data,var_name_pres, &pres_data, &time_comm_vel);

        MPI_Barrier(gcomm);

        // in case dspaces cannot read twice
        //put_raw_buffer(timestep, &num_points ,rank, &gcomm, var_name_vel_2,  &vel_data,var_name_pres_2, &pres_data, &time_comm_vel);
            timestep++;
        //}

        count++;


#else
    // use file writer
        OFstream ofs_p("p_"+ runTime.timeName() + ".txt");
        OFstream ofs_U("U_"+ runTime.timeName() + ".txt");
        
        // get dimensions here
        ofs_p << "dimensions "<< p.dimensions() << endl;

        // get name
        ofs_p << "internal field, size " << p.size() << endl;

        // Write contents
         forAll(p , i)
         {
             if (i) ofs_p << " ";
             ofs_p << p[i];
         }

        ofs_p << endl;

        // get dimension here
        ofs_U << "dimensions "<< p.dimensions() << endl;

        // get name
        ofs_U << "internal field, size " << p.size() << endl;

        // Write contents
         forAll(U , i)
         {

             // textfile output
             if (i) ofs_U << endl;
             ofs_U << U[i];
             ofs_U << endl << "splited"<< U[i].x() << U[i].y() << U[i].z();
         }

        ofs_U << endl;
#endif

      /*
        ofs_p << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        */

    }
#ifdef USE_DSPACES
// finalize dspaces
    // free                                                                                                                                                                                                                                   
    if(vel_data != NULL){      
        free(vel_data);
    }
    if(pres_data != NULL){     
        free(pres_data);
    }

    dspaces_finalize();

    MPI_Barrier(gcomm);
    MPI_Finalize();
    return 0;

#endif

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
