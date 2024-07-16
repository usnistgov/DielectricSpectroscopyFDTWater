#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

/*

This program reads the static and induced polarization of a
TINKER trajectory as per atom values and produces the total box 
polarization as well as the total box static and induced 
polarizations separately.

Input:

This program takes the static polarization output file from
TINKER. This program will look for a file in the working 
directory with the extension ".ustc" for this.

This program takes the induced dipole moment output of a simulation
using the AMOEBA force field in TINKER. This program will look for 
a file in the working directory with the extension ".uind" for this.

This program also takes as input a file named "parameters". This plain 
text file should have the following format:

         name_of_files.xyz

Output:

This program will output the total dipole moment (i.e., the 
polarization) of the simulation box into a plain-text file named 
"polarization.dat". Each line of this file lists the value of the 
polarization at the time step corresponding to that line.

This program also determines the static and induced polarization.
These are recorded in separate files, labelled respectively 
"static_polarization.dat" and "induced_polarization.dat".

Usage in bash environments:

Step 1: Compile (for instance, using the GCC compiler, enter the 
command "gcc grab_tinker_polariz.c")
Step 2: Run (call to execute the output of the compilation) (e.g., 
after running the above command, enter the command "./a.out")

- Dr. Rebecca A. Bone

*/

// Subroutines---------------------------------------------------------
// --------------------------------------------------------------------

// File existence check
int chk_file(char name[]) {
    // Function takes string array containing file name and locations
    // Function outputs flag signalling if indicated file exists
    int pacc = 0;
    if ( access( name,F_OK ) == -1) {
        pacc = 0; // No File found
    } else {
        pacc = 1; // File found
    }
    return pacc;
}

// Main function-------------------------------------------------------
// --------------------------------------------------------------------

void main(int argc,char* argv[]) {

    // Initialization -------------------------------------------------

    // Variable declarations
    int x,t,i,j,comp; // Integer variables for loops and flags
    char fname[100]; // Character array for parsing
    char fline[1200]; // Character array for reading in lines of files
    char delim[] = " "; // Expected delimiter within lines of files
    char fxyz[100]; // Character array for the name of the tinker xyz file
    char fustc[100]; // Character array for the name of the arc file
    char fuind[100]; // Character array for the the name of the uind file

    // Read or request input parameters
    if (argc < 2) { // Number of input parameters is less than expected
        printf("Tinker xyz file: "); // Initial configuration
        scanf("%s",&fxyz);
	x = chk_file(fxyz);
	while (x == 0) {
	    if (strcmp(fxyz,"exit") == 0) {
	        exit(EXIT_FAILURE);
	    } else if (strcmp(fxyz,"Exit") == 0) {
		exit(EXIT_FAILURE);
	    } else if (strcmp(fxyz,"EXIT") == 0) {
		exit(EXIT_FAILURE);
	    }
	    printf("Tinker xyz file: ");
	    scanf("%s",&fxyz);
	    x = chk_file(fxyz);
	}
    } else if (argc == 2) { // Number of input parameters is as expected
	sprintf(fxyz,"%s",argv[1]); // Initial configuration
	x = chk_file(fxyz);
	while (x == 0) {
	    if (strcmp(fxyz,"exit") == 0) {
		exit(EXIT_FAILURE);
	    } else if (strcmp(fxyz,"Exit") == 0) {
		exit(EXIT_FAILURE);
	    } else if (strcmp(fxyz,"EXIT") == 0) {
		exit(EXIT_FAILURE);
	    }
	    printf("Tinker xyz file: ");
	    scanf("%s",&fxyz);
	    x = chk_file(fxyz);
	}
    } else { // Number of input parameters is greater than expected
	printf("Too many input parameters detected.\n");
        printf("Tinker xyz file: \n"); // Initial configuration
	scanf("%s",&fxyz);
	x = chk_file(fxyz);
	while (x == 0) {
            if (strcmp(fxyz,"exit") == 0) {
	        exit(EXIT_FAILURE);
	    } else if (strcmp(fxyz,"Exit") == 0) {
		exit(EXIT_FAILURE);
	    } else if (strcmp(fxyz,"EXIT") == 0) {
		exit(EXIT_FAILURE);
	    }
	    printf("Tinker xyz file: \n");
	    scanf("%s",&fxyz);
	    x = chk_file(fxyz);
	}
    }

    printf("All input collected.\n"); // Report to terminal

    // Check for file existence of all input files
    x = chk_file(fxyz);
    if (x == 0) { // Hard error and exit
	printf("ERROR:\n No tinker xyz file found.\n");
	printf("All files must be in the working directory.\n");
	exit(EXIT_FAILURE);
    }

    // Parse initial configuration file name to names of polarization files
    sprintf(fname,"%s",fxyz);
    char *ptr = strtok(fname,"."); // Retrieve name up to first period
    sprintf(fustc,"%s.ustc",ptr); // Static polarization
    sprintf(fuind,"%s.uind",ptr); // Induced polarization

    // Check for remaining input file existence
    x = chk_file(fustc); // Static polarization
    if (x == 0) { // Hard error and exit
        printf("ERROR:\n No ustc file found.\n");
	printf("All files must be in the working directory.\n");
	exit(EXIT_FAILURE);
    }
    x = chk_file(fuind); // Induced polarization
    if (x == 0) { // Hard error and exit
	printf("ERROR:\n No uind file found.\n");
	printf("All files must be in the working directory.\n");
	exit(EXIT_FAILURE);
    }

    printf("Beginning system initialization.\n"); // Report to terminal

    // Open tinker xyz file
    FILE * xyz;
    xyz = fopen(fxyz,"r");

    // Retrieve number of atoms in simulation from xyz file
    fgets(fline,1200,xyz); // Retrieve first line
    ptr = strtok(fline,delim); // Parse first line entry
    int natom = atoi(ptr); // Assign to variable

    printf("%d atoms in simulation.\n",natom); // Report to terminal

    // Close xyz file
    fclose(xyz);

    // Declare values for polarization collection
    double spx = 0; // Static polarization x component
    double spy = 0; // Staticpolarization y coponent
    double spz = 0; // Static polarization z component
    double ipx = 0; // Induced polarization x component
    double ipy = 0; // Induced polarization y component
    double ipz = 0; // Induced polarization z component
    double ptotx = 0; // Total polarization x component
    double ptoty = 0; // Total polarization y component
    double ptotz = 0; // Total polarization z component

    // Open ustc file
    FILE * ustc;
    ustc = fopen(fustc,"r");

    // Open uind file
    FILE * uind;
    uind = fopen(fuind,"r");

    // Open static polarization file for reporting
    FILE * spout;
    spout = fopen("static_polarization.dat","w");

    // Open induced polarization file for reporting
    FILE * ipout;
    ipout = fopen("induced_polarization.dat","w");

    // Open total polarization file for reporting
    FILE * totout;
    totout = fopen("total_polarization.dat","w");

    printf("Initialization complete.\n"); // Report to terminal

    // Analysis------------------------------------------------------

    printf("Beginning analysis.\n"); // Report to terminal

    t = 0;
    comp = 0;
    while (comp == 0) { // Loop through each frame

        t++;
	printf("t=%d\n",t);

	// Zero out polarization values
	spx = 0;
	spy = 0;
	spz = 0;
	ipx = 0;
	ipy = 0;
	ipz = 0;
	ptotx = 0;
	ptoty = 0;
	ptotz = 0;

        //printf("Determining static polarization.\n"); // Report to terminal -> DEBUG

        // Static polarization
	fgets(fline,1200,ustc); // Header line
	if (fline == NULL) { // Check for end of file
	    comp = 1;
	    break;
	}
	fgets(fline,1200,ustc); // Possibly header line
	if (fline == NULL) { // Check for end of file
	    comp = 1;
	    break;
	}
	ptr = strtok(fline,delim); // Parse first entry
	sprintf(fname,"%s",ptr);
	while (strcmp(fname,"1") != 0) {
            fgets(fline,1200,ustc); // Retrieve line
	    if (fline == NULL) { // Check for end of file
	        comp = 1;
		break;
	    }
            ptr = strtok(fline,delim); // Parse first entry
            sprintf(fname,"%s",ptr);
	    if (strcmp(fname," ") == 0) {
	        ptr = strtok(NULL,delim);
		sprintf(fname,"%s",ptr);
	    }
	}
	// Read per-atom values from file
	for (i = 0; i < natom; i++) { // Loop through each atom
            if (i != 0) { // On subsequent atoms, retrieve next line
		fgets(fline,1200,ustc); // Retrieve line
		if (fline == NULL) { // Check for end of file
	            comp = 1;
		    break;
		}
		ptr = strtok(fline,delim); // Parse first enry-> atom number
	    }
	    ptr = strtok(NULL,delim); // Atomic symbol
	    ptr = strtok(NULL,delim); // ustc(x,i)
            spx+=strtod(ptr,NULL); // Add to x component total
	    ptr = strtok(NULL,delim); // ustc(y,i)
	    spy+=strtod(ptr,NULL); // Add to y component total
	    ptr = strtok(NULL,delim); // stc(z,i)
	    spz+=strtod(ptr,NULL);
	}
        spx*=0.208; // Convert from Debye to e- * Ang
	spy*=0.208; // Convert from Debye to e- * Ang
	spz*=0.208; // Convert from Debye to e- * Ang

	// Induced polarization

        //printf("Determining induced polarization.\n"); // Report to terminal -> DEBUG

	fgets(fline,1200,uind); // Header line
	if (fline == NULL) { // Check for end of file
	    comp = 1;
	    break;
	}
        fgets(fline,1200,uind); // Possibly header line
	if (fline == NULL) { // Check for end of file
            comp = 1;
	    break;
	}
        ptr = strtok(fline,delim); // Parse first entry
	sprintf(fname,"%s",ptr);
	x = strcmp(fname,"1");
	while (x != 0) {
	    fgets(fline,1200,uind); // Retrieve line
	    if (fline == NULL) { // Check for end of file
		comp = 1;
		x = 0;
		break;
	    }
	    ptr = strtok(fline,delim);
            sprintf(fname,"%s",ptr);
	    x = strcmp(fname,"1");
        }
        // Read per-atom values from file
	for (i = 0; i < natom; i++) {
	    if (i != 0) { // On subsequent, read next line
                fgets(fline,1200,uind); // Retrieve line
		if (fline == NULL) { // Check for end of file
		    comp = 1;
		    break;
	        }
		ptr = strtok(fline,delim); // Atom number
	    }
            ptr = strtok(NULL,delim); // Atomic symbol
            ptr = strtok(NULL,delim); // uind(x)
            ipx+=strtod(ptr,NULL); // Add to x component total
            ptr = strtok(NULL,delim); // uind(y)
            ipy+=strtod(ptr,NULL); // Add to y component total
            ptr = strtok(NULL,delim); // uind(z)
            ipz+=strtod(ptr,NULL); // Add to z component total
	}
        ipx*=0.208; // Convert to e- * A from Debye
	ipy*=0.208; // Convert to e- * A
	ipz*=0.208; // Convert to e- * A

        printf("Printing.\n"); // Report to terminal -> DEBUG

	// Total polarization
	ptotx =  spx + ipx;
	ptoty =  spy + ipy;
	ptotz =  spz + ipz;

	// Report polarization
	fprintf(spout,"%lf %lf %lf\n",spx,spy,spz);
        fprintf(ipout,"%lf %lf %lf\n",ipx,ipy,ipz);
	fprintf(totout,"%lf %lf %lf\n",ptotx,ptoty,ptotz);

        //printf("Current frame complete.\n"); // Report ot terminal -> DEBUG

    }

    printf("Analysis complete.\n"); // Report to terminal

    // Clean up--------------------------------------------------------

    printf("Beginning clean up.\n"); // Report to terminal

    // Close files
    fclose(ustc);
    fclose(uind);
    fclose(spout);
    fclose(ipout);
    fclose(totout);

    printf("Complete.\n"); // Report to terminal
    exit(EXIT_SUCCESS);
}
