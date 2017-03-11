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

// file writer for simulation data
//#define debug_1

#define USE_DSPACES

#ifdef USE_DSPACES
    #include "ds_adaptor.h"
    #include "divide.h"

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
    // fake one, only select first num_elem_sample
    void prepare_sampled_buffer(float *buffer_region, float *buffer_sample, int num_elems_region, int num_elems_sample, int region_length){
        int i;
        int spacing = region_length*region_length*3;
        for(i = 0; i< num_elems_sample; i++){
            memcpy(buffer_sample + i*spacing, buffer_region+ i*spacing,spacing*sizeof(float));
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

        // cartisian position
        int proc_row_pos = rank/procs_per_dim;
        int proc_col_pos = rank%procs_per_dim;


        // icofoam size
        int case_length = CASE_LENGTH;


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
        double time_comm = 0;
        double time_comm_vel = 0;
        double time_comm_pres = 0;
        double time_comp = 0;
        double time_comm_sample = 0;
        double time_latency = 0;
        double t1, t2, t3;
        // timer for end-to-end
        double t_start, t_end;
        
        char var_name_vel[STRING_LENGTH];
        char var_name_pres[STRING_LENGTH];
        sprintf(var_name_vel, "VEL");
        sprintf(var_name_pres, "PRES");


        // data layout
//#ifdef FORCE_GDIM
        uint64_t gdims_raw[2] = {POINTS_SIDE, POINTS_SIDE};
        dspaces_define_gdim(var_name_vel, 2,gdims_raw);
        dspaces_define_gdim(var_name_pres, 2,gdims_raw);
//#endif

        /*
        char var_name_vel_2[STRING_LENGTH];
        char var_name_pres_2[STRING_LENGTH];
        sprintf(var_name_vel_2, "VEL_2");
        sprintf(var_name_pres_2, "PRES_2");
        */
        size_t elem_size_vel = sizeof(float)*3;
        size_t elem_size_pres = sizeof(float);
        unsigned int dims[3] = {case_length, case_length, 1};
        uint64_t num_elems = dims[0]*dims[1]*dims[2];
        
        int bounds[6] = {0};
        // xmin
        bounds[1]=(case_length)*proc_row_pos;
        // ymin
        bounds[0]=(case_length)*proc_col_pos;

        // xmax
        bounds[4]=case_length*(proc_row_pos+1)-1 ;
        // ymax
        bounds[3]=(case_length)*(proc_col_pos+1) - 1;

        // prepare space
        float * vel_data = (float *)malloc(num_elems*elem_size_vel);
        if(vel_data == NULL){
              perror("vel data allocated error");
              exit(-1);
          }
        // prepare space for pres
        float * pres_data = (float *)malloc(num_elems*elem_size_pres);
        if(pres_data == NULL){
              perror("pres data allocated error");
              exit(-1);
          }

        // regions stripped to this rank
        int region_length = REGION_LENGTH;
        int num_region_row = dims[0]/region_length;
        int num_region_col =  dims[1]/region_length;
        int num_region = num_region_row*num_region_col; // this will be (201-1)/10 = 20

        int num_elems_region = num_region;
        size_t elem_size_region = (region_length)*(region_length)*3*sizeof(float);

        float *buffer_region = (float *)malloc(num_elems_region*elem_size_region);
        if(buffer_region== NULL){
            perror("    allocate space for striped  regions");
            exit(-1);
        }

        /*
         * sampled vel data
         */
        // sampling related
        int sample_size = SAMPLE_SIZE;
        if(sample_size > num_region){
            printf("cannot sample more than number of regions\n");
            exit(-1);
        }

        char var_name_sample[STRING_LENGTH];
        sprintf(var_name_sample, "sample");

        int bounds_sample[6]={0};
        // x_min
        bounds_sample[0] = rank*sample_size;

        // x_max
        bounds_sample[3] = (rank+1)*(sample_size) -1; 

        int num_elems_sample = sample_size;
        size_t elem_size_sample = (region_length)*(region_length)*3*sizeof(float);

        float *buffer_sample = (float *)malloc(num_elems_sample*elem_size_sample);
        if(buffer_sample== NULL){
            perror("    allocate space for sampled  regions");
            exit(-1);
        }
        uint64_t gdims_sample[1] = {nprocs*sample_size};
        dspaces_define_gdim(var_name_sample, 1,gdims_sample);

        
// end of dspaces preparation
#endif
#ifdef debug_1

        OFstream ofs_p("p_.txt");
        OFstream ofs_U("U_.txt");
#endif

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        //sleep(2);
        Info<< "Time = " << runTime.timeName() << nl << endl;
        printf("********************timestep %d now start!\n",timestep);

        MPI_Barrier(gcomm);
        t_start = MPI_Wtime();

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
        // start timer
        
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

        t2 = MPI_Wtime();
        time_comp = t2-t_start;

        // if the setting is 'no write', that field won't be written to files
        runTime.write();

        /*
         * also write pressure info
         */

#ifdef USE_DSPACES
    // if use  dataspces, write the correct 
        // dump data from U and P into buffer
        mydump(p, pres_data);
        mydump(U, vel_data);

        /*
        Info<<" first data, address"<< vel_data << ": "<< vel_data[0] <<" " << vel_data[1]<< " "<<vel_data[2]<< endl;
       printf(" first data, address %p: %f %f %f\n", vel_data, vel_data[0], vel_data[1], vel_data[2]);
       */

        // 1. write raw buffer
        put_common_buffer(timestep,2, bounds,rank, &gcomm, var_name_vel, (void **)&vel_data, elem_size_vel, &time_comm_vel);
        put_common_buffer(timestep,2, bounds,rank, &gcomm, var_name_pres, (void **)&pres_data, elem_size_pres, &time_comm_pres);

        // 3. divide into regions
        int num_region_2;
        divide(vel_data, dims,region_length,&num_region_2, buffer_region);
        printf("divide completed %d regions generated\n", num_region_2);
#ifdef debug_1
        int spacing = elem_size_region/sizeof(float);
        printf(" first data of all regions: %f %f %f \n", buffer_region[0], buffer_region[1], buffer_region[2]);
        printf(" last data of all regions: %f %f %f \n", buffer_region[spacing*num_region -3], buffer_region[spacing*num_region -2], buffer_region[spacing*num_region-1]);
#endif
        
        // 3. sampling
        prepare_sampled_buffer(buffer_region, buffer_sample, num_elems_region, num_elems_sample, region_length);
        printf("local_sample generated\n");

        // 4. send own sampled regions(only velocity)
        put_common_buffer(timestep, 1, bounds_sample, rank, &gcomm, var_name_sample,(void **)&buffer_sample, elem_size_region,  &time_comm_sample);
        printf("local_sample sent\n");



#ifdef debug_1
        spacing = elem_size_region/sizeof(float);
        printf(" first data of local samples: %f %f %f \n", buffer_sample[0], buffer_sample[1], buffer_sample[2]);
        printf(" last data of local samples: %f %f %f \n", buffer_sample[spacing*num_elems_sample -3], buffer_sample[spacing*num_elems_sample -2], buffer_sample[spacing*num_elems_sample-1]);
#endif




        
        // write samples into dataspaces
        
        time_comm = time_comm_vel + time_comm_pres;


        
#endif


#ifdef debug_1
    // use file writer
        int spacing = elem_size_vel/sizeof(float);
        printf(" first data before send, address %p: %f %f %f\n", vel_data, vel_data[0], vel_data[1], vel_data[2]);
        printf(" last data before send, address %p: %f %f %f\n", vel_data[spacing*num_elems-3], vel_data[spacing*num_elems-2], vel_data[spacing*num_elems-1]);
        
        
        // get dimensions here
        ofs_p << "** time step " << timestep << endl <<"dimensions "<< p.dimensions() << endl;

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
        ofs_U << "** time step " << timestep << endl << "dimensions "<< U.dimensions() << endl;

        // get name
        ofs_U << "internal field, size " << U.size() << endl;

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



        MPI_Barrier(gcomm);
        t_end = MPI_Wtime();
        time_latency = t_end-t_start;

        double global_time_comm;
        double global_time_comp;
        double global_time_comm_sample;
        MPI_Reduce(&time_comm, &global_time_comm, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);
        MPI_Reduce(&time_comm_sample, &global_time_comm_sample, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);
        MPI_Reduce(&time_comp, &global_time_comp, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);

        // Print the result
        if (rank == 0) {
          printf("%d comp sim Total  %lf avg %lf\n",timestep,  global_time_comp , global_time_comp/ (nprocs));
          printf("%d comm raw Total %lf avg %lf\n",timestep,  global_time_comm , global_time_comm/ (nprocs));
          printf("%d comm sample Total %lf avg %lf\n",timestep,  global_time_comm_sample , global_time_comm_sample/ (nprocs));

          printf("%d all latency Total %lf\n",timestep,  time_latency );
        }
        timestep++;

        count++;


        if(timestep == MAX_VERSION){
            break;
        }
    }
    if(buffer_region != NULL){
        free(buffer_region);
        sprintf(msg,"-- buffer_region freed");
        my_message(msg, rank, LOG_CRITICAL);
    }
    if(buffer_sample != NULL){
        free(buffer_sample);
        sprintf(msg,"-- buffer_sample freed");
        my_message(msg, rank, LOG_CRITICAL);
    }

    MPI_Barrier(gcomm);
    // reduce all comm_time

 
    // get avg comm_time

    // reduce all timers
#ifdef USE_DSPACES
// finalize dspaces
    // free                                                                                                                                                                                                                                   
    if(vel_data != NULL){      
        free(vel_data);
    }
    if(pres_data != NULL){     
        free(pres_data);
    }
    printf("now finalize the dspaces and exit");

    dspaces_finalize();

    MPI_Barrier(gcomm);
    MPI_Finalize();
    return 0;

#endif

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
