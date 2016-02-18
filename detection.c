// make a color map of sensitivity as a function scintillation time-scale and bandwidth
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "detection.h"
#include <mpi.h>

int main (int argc, char* argv[])
{
	int id;  //  process rank
	int p;   //  number of processes
	//int num;
	int nMax;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &p);

	//FILE *fin;
	controlStruct control;
	acfStruct acfStructure;
	noiseStruct noiseStructure;

	char fname[1024];   // read in parameter file
	char Tname[1024];   // read in scintillation timescale

	int nt;
	double *tn;

	//char oname[1024];   // output file name
	//char pname[1024];   // file to plot
	char dname[1024];   // graphics device

	int i;
	int n = 1;  // number of dynamic spectrum to simulate

	double flux0;
	double flux1;

	int plotMode = 0;  // plot existing dynamic spectrum, off by default
	control.noplot = 1;  // don't show dynamic spectrum while simulationg, on by default

	// read options
	for (i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
			strcpy(fname,argv[++i]);
			strcpy(Tname,argv[++i]);
			printf ("Parameters are in %s\n", fname);
		}
		//else if (strcmp(argv[i],"-o")==0)
		//{
		//	strcpy(oname,argv[++i]);
		//	//printf ("Dynamic spectrum is output into %s\n", oname);
		//}
		else if (strcmp(argv[i],"-n")==0)
		{
			n = atoi (argv[++i]);
			printf ("Number of dynamic spectrum to simulate %d\n", n);
		}
		//else if (strcmp(argv[i],"-p")==0) // just plot 
		//{
		//	strcpy(pname,argv[++i]);
		//	plotMode = 1;
		//	printf ("Plotting exiting dynamic spectrum.\n");
		//}
		//else if (strcmp(argv[i],"-noplot")==0)
		//{
		//	control.noplot = 1; // Don't show dynamic spetrum while simulationg
		//}
		//else if (strcmp(argv[i],"-dev")==0)
		//{
		//	strcpy(dname,argv[++i]);
		//	printf ("Graphic device: %s\n", dname);
		//}
	}

	nt = readDissNum (Tname);
	tn = (double *)malloc(sizeof(double)*nt);
	readDiss (Tname, tn);

	sprintf(dname, "%s", "1/xs");
	/*
	if ((fin=fopen(oname, "w"))==NULL)
	{
		printf ("Can't open file...\n");
		exit(1);
	}
	*/

	if (plotMode==0)
	{
		// Simulate dynamic spectrum
		// read parameters
		initialiseControl(&control);
		readParams (fname,dname,n,&control);
		//readParams (fname,oname,dname,n,&control);
		printf ("Finished reading parameters.\n");

		//preAllocateMemory (&acfStructure, &control);
		//allocateMemory (&acfStructure);

		//num = (int)((control.scint_ts1-control.scint_ts0)/control.scint_ts_step);
		//for (tdiff=control.scint_ts0; tdiff<control.scint_ts1; tdiff+=control.scint_ts_step)
		//for (i=id; i<=num; i+=p)
		for (i=id; i<nt; i+=p)
		{
			control.nsub = tn[i];
			control.tsub = control.T/control.nsub;
			control.scint_ts = control.tsub;
			control.whiteLevel = 0.1*sqrt(control.nsub*control.nchan);  // 0.1 gives 1 to a 10*10 dynamic spectrum
			//printf ("%d %lf %lf %lf\n", control.nsub, control.tsub, control.scint_ts, control.whiteLevel);

			calNoise (&noiseStructure, &control);

			flux0 = control.cFlux0;
			flux1 = control.cFlux1;
				
			acfStructure.probability = 1.0;

			nMax = 0;
			
			while (fabs(acfStructure.probability-0.8)>=0.01 && nMax <= 100)
			{
				control.cFlux = flux0+(flux1-flux0)/2.0;
				//printf ("%lf %lf %.8lf %.3f\n", tdiff, fdiff, control.cFlux, acfStructure.probability);
		
				// simulate dynamic spectra
				calculateScintScale (&acfStructure, &control);

				qualifyVar (&acfStructure, &noiseStructure, &control);

				if (acfStructure.probability>0.8)
				{
					flux1 = control.cFlux;
					
				}
				else 
				{
					flux0 = control.cFlux;
				}
				nMax++;
			}
				
			printf ("%d %lf %f %d\n", control.nsub, control.cFlux, acfStructure.probability, nMax);
			//fprintf (fin, "%lf %lf %lf %f\n", tdiff, fdiff, control.cFlux, acfStructure.probability);
		}

		fflush (stdout);
		MPI_Finalize ();

		printf ("Threshold: %f \n", noiseStructure.detection);
		// deallocate memory
		deallocateMemory (&acfStructure, &noiseStructure);
	}
	//else
	//{
	//	plotDynSpec(pname, dname);
	//}

	/*
	if (fclose(fin))
	{
		printf ("Can't close file...\n");
		exit(1);
	}
	*/

	free(tn);
	return 0;
}

