//=================================================================
//  File Constants.h
//  Created:  March 17, 2017
//  Auhtor: Sathwik Bharadwaj
//  Modified: 
//  Last change: 
//  Function: 
//  PetscErrorCode main (int argc,char** argv)
//=================================================================

//=================================================================
static char help[] = "An Example code to solve Laplace problem in 2D using FEM within Parallel Computing Environment\n\n";
//header file with all external libraries and user defined classes.
#include "global_params.h"


//=================================================================
// argc is # of command line entries +1,
// argv is list of entries
// Note that the first entry in argv is the program name itself
// so argc[1] is the first input file.
//=================================================================
//PetscErrorCode will Enable the error tracing
//=================================================================

PetscErrorCode main (int argc,char** argv)
{

	//Obtain the current time

	PetscInitialize(&argc, &argv, (char*) 0, help);//initializing petsc
	PetscLogDouble _begintime;		
	PetscCall(PetscTime(&_begintime));
	//=====================================================
	// The program is called with one parameter, 
	// the FEMstructure.inp input file
	//=====================================================

	// Check if the user put the correct number of arguments 
	if(argc <=1)
	{
		std::cout << " Need input file name on the command line\n";
		exit(0);
	}

	// Open the file for reading
	std::ifstream input_structurefile(argv[1]);

	//Check if the file exists
	if(!input_structurefile)
	{
		std::cout << "ERROR: The file '" << 
			argv[1] << "' was not found in the working folder" << std::endl;
		exit(0);
	}

	//=====================================================
	// Declaring the structures:
	//=====================================================
	// The class data contains all the necessary info 
	// for FEM work
	//=====================================================

	data dat; // declaring data class defined in global_params.h

	//=====================================================
	// The struct global_matrices has the global matrices
	//=====================================================

	global_matrices gmat; //delcaring global_matrices class defined in global_params.h

	PetscErrorCode ierr;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"=========================================\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,":     Example FEM calculation\t:\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,":       in MPI environment\t:\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,":      2D Laplace equation \t:\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,":        Sathwik Bharadwaj\t:\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"=========================================\n");CHKERRQ(ierr);
	//=====================================================


	//=====================================================
	// Open log file for writing progress of the calculation
	// =================================================

	char buffer[200];
	sprintf(buffer, "../log.out");
	dat.log1 = fopen(buffer, "w"); 


	//================================================================
	// call input subroutines with file pointer input_structurefile, 
	// and mesh input.  
	// Subroutine input is in file "input_reader.cpp"
	//================================================================
	mesh_input(dat); //call mesh_input to read in node,elem, bnode and belem data
			 //and store in dat
	input_reader(input_structurefile, dat); //call input reader to read in FEMstruct
						//to store FEM parameters
	ierr = PetscPrintf(PETSC_COMM_WORLD,"=========================================\n");CHKERRQ(ierr);
	input_structurefile.close(); //close the data files 

	ierr = PetscPrintf(PETSC_COMM_WORLD,"=========================================\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Input files are read.\n");CHKERRQ(ierr);

	//====================================================
	// construct the global matrices from element matrices
	//====================================================
	make_global(gmat, dat); //construct the global matrices by first 
				//defining the element matrices
	ierr = PetscPrintf(PETSC_COMM_WORLD,"=========================================\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Global matrix and vector are constructed.\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"=========================================\n");CHKERRQ(ierr);
	//====================================================
	// Apply  boundary conditions on the global matrices
	//====================================================

	apply_bc(gmat, dat); //apply BCs - Assigning values at some nodes
			     //and carrying out benediction
	ierr = PetscPrintf(PETSC_COMM_WORLD,"=========================================\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Boundary conditions are applied.\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"=========================================\n");CHKERRQ(ierr);

	//====================================================
	// Call the solver routine 
	//=================================================
	//====================================================
	// The solve subroutine calls solution, where 
	// we output the interpolated solution is written
	// into a file.
	//====================================================

	ierr = PetscPrintf(PETSC_COMM_WORLD,"=========================================\n");CHKERRQ(ierr);
	petsc_diagonalizer(gmat, dat); 
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Performed Matrix Diagonalization.\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"=========================================\n");CHKERRQ(ierr);
	PetscLogDouble _endtime;		
	PetscCall(PetscTime(&_endtime));

	ierr = PetscPrintf(PETSC_COMM_WORLD,"=========================================\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Total time taken: %1.2g seconds\n", _endtime-_begintime);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"=========================================\n");CHKERRQ(ierr);

	PetscFinalize(); //closing petsc. This command should be at the end of the compilation. 

	// =========================================================
	// release the allocated memory and close open output files
	// =========================================================

	fclose (dat.log1); 

	return ierr;

}

//=============================================================
//=============================================================
//=============================================================
