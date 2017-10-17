#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

// global variables
double *grid;
double min;
int *min_perm;
int *temp_perm;
double *distances;
char *ptype="TSP";
int use_greedy=0;

// functions used
void parse_args(int argc, char *argv[], int *n, char **input, char **output, int *timing, int* method, double *cooling, int *reject_base);
double euc_dist(int i, int j);
double travel_dist(int *a, int n);
void grid_set_up(char *input, int *n, double grid_limit);
void dist_setup(int n);

void swap(int *a, int i, int j);
void brute_force(int *a, int l, int r, int n);

void greedy(int n);
void randomise(int n);
void two_opt(int start, int end);
void sa(double temp, double cooling, int reject_base, int n);

double tsp(int n, int method, double temp, double cooling, int reject_base);
void grid_to_file(char *output, int n);
void print_sol(int n);







int main(int argc, char *argv[]){
	srand48(time(NULL));	

	//set default arguments
	int n=5;
	double grid_limit=2;
	char* input=NULL;
	char *output="plot.dat";
	int timing=0;
	double temp=10;
	double cooling=0.9999;
	int reject_base=1000;
	int method=2;
	parse_args(argc, argv, &n, &input, &output, &timing, &method, &cooling, &reject_base);

	// prints out type of tour
	#ifdef CLOSED
		printf("CLOSED LOOP\n");
	#else
		printf("OPEN LOOP\n");
	#endif
	
	// if just a single permutation is requested
	if(timing==0){
		grid_set_up(input, &n, grid_limit);		// sets up grid randomly or from file
		printf("Calculating %s for %d cities\n",ptype,n);
		double elapsed_time=tsp(n, method, temp, cooling, reject_base);		// times the tsp
		print_sol(n);					// prints the solution to screen
		printf("Time taken:\t%.2es\n",elapsed_time);	// and the time taken
		printf("Min Distance is %lf\n",min);
		grid_to_file(output, n);			// saves the solution to a file
		free(grid);free(min_perm);free(distances);
	}

	// otherwise does the timing for up to a specified number of cities and stores them in a file
	else{
		grid_set_up(NULL, &timing, grid_limit);		// sets up a grid of the max size
		printf("Creating data plot of timing for %s up to %d cities\n",ptype,timing);
		
		// sets up timing file
		FILE *time_data;
		if((time_data=fopen(output,"w"))==NULL){
			fprintf(stderr, "*** ERROR 1 ***\nCouldn't open timing.dat\n");
			exit(1);
		}
		int i;
		double old_time, new_time;

		// prints the first line for 1 city to file
		old_time=tsp(1, method, temp, cooling, reject_base);
		fprintf(time_data,"%d\t%lf\tN/A\n",1,old_time);

		printf("1 done\n");
		// loops through the remaining number of cities
		for(i=2;i<=timing;i++){
			new_time=tsp(i, method, temp, cooling, reject_base);
			printf("%d done\n",i);
			fprintf(time_data,"%d\t%lf\t%lf\n",i,new_time, new_time/old_time);	
			free(min_perm);
			old_time=new_time;	// keeps old time to see increase ratio
		}
		free(grid);free(distances);
		fclose(time_data);
	}


	return 0;
}

//parse command line arguments
void parse_args(int argc, char *argv[], int *n, char **input, char **output, int *timing, int *method, double *cooling, int *reject_base){
	int opt;
	while((opt=getopt(argc,argv,"n:f:o:t:m:bc:r:g"))!=-1){
		switch(opt){
			case 'n':
				*n=atoi(optarg);
				break;
			case 'f':
				*input=optarg;
				break;
			case 'o':
				*output=optarg;
				break;
			case 't':
				*timing=atoi(optarg);
				break;
			case 'm':
				*method=atoi(optarg);
				break;
			case 'b':
				*method=1;
				break;
			case 'c':
				*cooling=strtod(optarg, NULL);
				break;
			case 'r':
				*reject_base=atoi(optarg);
				break;
			case 'g':
				use_greedy=1;
				break;
			default:
				fprintf(stderr,"Usage: %s [-b] [-c double] [-f file] [-g] [-m int] [-n int] [-o file] [-r int] [-t int]\nSee README file for more details\n",argv[0]);
				exit(EXIT_FAILURE);
		}
	}
	return;
}

// returns Euclidean distance between two points in the grid
double euc_dist(int i, int j){
	return pow(pow(grid[2*i]-grid[2*j],2)+pow(grid[2*i+1]-grid[2*j+1],2),0.5);
}

// calculates the travel distance for a given permutation
double travel_dist(int *a, int n){
	double distance=0;
	int i;

	// adds up all the step distances
	for(i=0;i<n-1;i++){
		distance+=distances[a[i]*n+a[i+1]];
	}

	// adds the the distance back to the start if closed
	#ifdef CLOSED
		distance+=distances[a[n-1]*n+a[0]];
	#endif
	return distance;
}

// sets up the grid either randomly or from a file
void grid_set_up(char *input, int *n, double grid_limit){
	int i;
	// if we want a random grid
	if(input==NULL){
		grid=malloc(2*(*n)*sizeof(double));
		// randomly assigns values based on a grid_limit
		for(i=0;i<2*(*n);i++){
			grid[i]=2*grid_limit*(drand48()-.5);
		}
		dist_setup(*n);	
	}
	// else we set it up from a file
	else{
		FILE *data;
		// opens up the file		
		if((data=fopen(input, "r"))==NULL){
			fprintf(stderr, "*** ERROR 1 ***\nCouldn't open %s\n",input);
			exit(1);
		}
		char lineBuffer[50], *p, *end;

		fgets(lineBuffer,50,data);

		// if reading from a file on TSPLIB
		if(strncmp(lineBuffer,"NAME",4)==0){
			// finds dimension
			while(fgets(lineBuffer, 50, data)){
				if(strncmp(lineBuffer,"DIMENSION",9)==0){
					p=lineBuffer;
					while(strtod(p,&end)<1){
						p++;
					}
					*n=atoi(p);
					break;
				}
			}


			// finds problem type
			rewind(data);
			while(fgets(lineBuffer, 50, data)){
				if(strncmp(lineBuffer,"TYPE",4)==0){
					p=lineBuffer;
					while(strncmp(p,"ATSP",4)!=0 && strncmp(p,"TSP",3)){
						p++;
					}
					break;
				}
			}

			// if its a standard TSP with coordinates, it sets up grid and calculates distances
			if(strncmp(p,"TSP",3)==0){
				grid=malloc(2*(*n)*sizeof(double));	// sets up grid

				// looks for start position
				while(fgets(lineBuffer, 50, data)){
					// looks for header to coord data to start from
					if(strncmp(lineBuffer,"DISPLAY_DATA_SECTION",20)==0 || strncmp(lineBuffer,"NODE_COORD_SECTION",18)==0){
						break;
					}
				}

				// fills in all the values while error checking
				for(i=0;i<(*n);i++){	
					fgets(lineBuffer, 50, data);
					p=lineBuffer;
					strtod(p,&end);
					if(p==end){
						fprintf(stderr,"*** ERROR 2 ***\nToo few entries on line %d of %s\n",i+1, input);
							exit(2);
					}
					p=end;
					grid[2*i]=strtod(p, &end);
					if(p==end){
						fprintf(stderr,"*** ERROR 2 ***\nToo few entries on line %d of %s\n",i+1, input);
							exit(2);
					}
					p=end;
					grid[2*i+1]=strtod(p, &end);
					if(p==end){
						fprintf(stderr,"*** ERROR 2 ***\nToo few entries on line %d of %s\n",i+1, input);
						exit(2);
					}
				}
				dist_setup(*n);		// fill in distances matrix	
			}
			// Otherwise reading an Assymetric TSP file
			// reads in distances 
			else{
				ptype="Assymetric TSP";
				while(fgets(lineBuffer, 50, data)){
					// looks for header to matrix data to start from
					if(strncmp(lineBuffer,"EDGE_WEIGHT_SECTION",19)==0){
						break;
					}
				}
				int len=(*n)*(*n);
				distances=malloc(len*sizeof(double));
				char longlb[1000];
				double num;
				i=0;

				// loops through the matrix data storing them in distances
				while(fgets(longlb,1000,data)){
					p=longlb;
					while(strcmp(p,"\n")!=0){
						num=strtod(p,&end);

						// if it finds a value it stores it and bumps p forward the length of the number
						if(num!=0 || strncmp(p,"0",1)==0){
							distances[i]=num;
							i++;
							p=end;
						}
						// else just bumps it one place
						else{
							p++;
						}
					}
				}

				// if we didn't read in n*n values then there was an issue somewhere
				if(i!=len){
					fprintf(stderr,"*** ERROR 3 ***\nIssue reading in distances matrix from %s\n",input);
					exit(3);
				}

				// we don't have coord values for the grid so we set to NULL to avoid errors when freeing later
				grid=NULL;
			}
		}
		// else reading from a plain file like an output from previous run
		else{

			*n=0;
			// finds how many lines are in the file
			while(fgets(lineBuffer, 50, data)){
				(*n)++;
			}
			grid=malloc(2*(*n)*sizeof(double));	// sets up grid
			rewind(data);				// rewinds back to start of file
				
			// loops through each line storing the values while error checking for correct number of inputs
			for(i=0;i<(*n);i++){
				fgets(lineBuffer, 50, data);
				p=lineBuffer;
				grid[2*i]=strtod(p, &end);
				if(p==end){
					fprintf(stderr,"*** ERROR 2 ***\nToo few entries on line %d of %s\n",i+1, input);
						exit(2);
				}
				p=end;
				grid[2*i+1]=strtod(p, &end);
				if(p==end){
					fprintf(stderr,"*** ERROR 2 ***\nToo few entries on line %d of %s\n",i+1, input);
					exit(2);
				}
			}
			dist_setup(*n);		
		}
		fclose(data);
	}
}

// sets up matrix with distances where entry (i,j) represents the distance between i and j 
void dist_setup(int n){
	int i,j;
	distances=malloc(n*n*sizeof(double));
	
	// sets up the upper part
	for(i=0;i<n-1;i++){
		for(j=i+1;j<n;j++){
			distances[i*n+j]=euc_dist(i,j);
		}
	}

	// and copies into lower part
	for(j=0;j<n-1;j++){
		for(i=j+1;i<n;i++){
			distances[i*n+j]=distances[j*n+i];
		}
	}
}

// swaps two values in a given array
void swap(int *a, int i, int j){
	double temp;
	temp=a[i];
	a[i]=a[j];
	a[j]=temp;
}

// permutes the array recursively and finds the distance for each permuation updating min and min_perm if necessary
void brute_force(int *a, int l, int r, int n){
	int i;
	
	// if we've reached a final permutation
	if(l==r){
		// find the distance and updates the mins if necessary
		double dist=travel_dist(a,n);
		if(dist<min){
			min=dist;
			for(i=0;i<n;i++){
				min_perm[i]=a[i];
			}
		}
	}

	// else recursively permute the rest of the array
	else{
		for(i=l;i<=r;i++){
			swap(a,l,i);
			brute_force(a, l+1, r, n);
			swap(a,l,i);
		}
	}
}

// Greedy nearest neighbour algorithm
void greedy(int n){
	int i,j;

	// temp array indicates if a point is yet to be visited
	int *temp=malloc(n*sizeof(int));
	for(i=0;i<n;i++){temp[i]=1;}
	
	int index=0;

	// if it's open we try find a good start point if possible, finds the index with the lowest sum of coordinates (roughly in the bottom left corner)
	#ifndef CLOSED
		// only runs if we have coordinates (ie symmetric version)
		if(strcmp(ptype,"TSP")==0){
			double min=grid[0]+grid[1], min_temporary;	// initialises min point as the first one
			for(i=1;i<n;i++){
				min_temporary=grid[2*i]+grid[2*i+1];
		
				// updates min and index if new point found
				if(min_temporary<min){
					min=min_temporary;
					index=i;
				}
			}
		}
	#endif

	// uses this point as the start point of the tour
	temp[index]=0;
	min_perm[0]=index;
	int perm_index=1;
	int last;

	// at each point it looks for the closest point that hasn't already been visited and goes there next
	for(i=1;i<n;i++){
		j=0;
		last=min_perm[i-1];

		// finds the first point in the array that hasn't been visited yet to establish a min value
		while(j<n){
			if(temp[j]==1){
				min=distances[last*n+j];
				index=j;
				break;
			}
			j++;
		}
	
		// loops through the rest of the values and sees if there is a closer value
		for(j=j+1;j<n;j++){
			if(temp[j]==1 && distances[last*n+j]<min){
				min=distances[last*n+j];
				index=j;
			}
		}

		// sets it as the next value
		min_perm[perm_index]=index;
		temp[index]=0;
		perm_index++;
	}		
	free(temp);
}	

// randomises min_perm
void randomise(int n){
	int i, j, counter, rand;
	for(i=0;i<n;i++){temp_perm[i]=1;}	// sets the temp to all 1s to indicate unvisited

	// sets the n values in min_perm
	for(i=0;i<n;i++){
		rand=drand48()*(n-i);		// picks a random number in [0,n-i)
		counter=-1;
		j=0;

		// finds the "rand"th unvisited number and sets that 
		while(1){
			if(temp_perm[j]==1){counter++;}
			if(counter==rand){break;}
			j++;
		}
		min_perm[i]=j;
		temp_perm[j]=0;
	}
}

// reverses the path between some start and end point 
void two_opt(int start, int end){
	int gap=end-start;
	int i;

	// stores the current order in temp
	for(i=0;i<gap;i++){
		temp_perm[i]=min_perm[start+i];
	}

	// puts them back in min_perm in reverse order
	for(i=0;i<gap;i++){
		min_perm[start+i]=temp_perm[gap-i-1];
	}
}


// simulated annealing method
void sa(double temp, double cooling, int reject_base, int n){

	// if grid is small enough there is only one permutation (depends on open/closed)
	//	so it returns straight away to avoid unnecessary iterations
	#ifdef CLOSED
		if(n<=3){return;}
	#else
		if(n<=2){return;}
	#endif

	int rand, rand2;
	double old, new;
	old=travel_dist(min_perm, n);	// current distance
	int i;
	int rejects=0;
	
	int loop=1000;			// number of loops in between reject checks

	// if requested start by using the greedy nearest neighbour algorithm
	if(use_greedy==1){
		greedy(n);
	}
	// else randomise
	else{
		randomise(n);
	}

	// if we wish to plot the min over time then open up the file
	#ifdef PLOT
		FILE *mindata;
		if((mindata=fopen("mindata.dat","w"))==NULL){
			fprintf(stderr,"*** ERROR 1 ***\nCouldn't open mindata.dat\n");
			exit(1);
			}
		int counter=0;
	#endif

	int reject_limit=sqrt(n)*reject_base;		// sets the reject limit based on n

	// loops while the temp is above a certain threshold
	while(temp>.00001){
		// loops in increments of 'loop'
		for(i=0;i<loop;i++){

			// finds a random point and another point further on
			rand=drand48()*(n-1);
			rand2=rand+drand48()*(n-rand);

			// runs two_opt on these indexes
			two_opt(rand, rand2);

			// finds new tour distance
			new=travel_dist(min_perm,n);

			// if its less then keep that permutation
			if(new<old){
				old=new;		// update old distance
				rejects=0;		// reset rejects
				temp*=cooling;		// reduce temp by cooling factor
			}
			else{
				// otherwise potentially accept it anyway based on how much worse it is and temp
				if(exp((old-new)/(old*temp))>drand48()){
					old=new;
				}
				// else increment rejects and undo the change
				else{
					rejects++;
					two_opt(rand, rand2);
				}
			}
		}

		// if the number of rejects (in a row) is above the limit then break 
		if(rejects>reject_limit){
			printf("Broke at temp=%lf\n",temp);
			break;
		}
		// if we are plotting then insert the current min value in 
		#ifdef PLOT
		counter++;
		fprintf(mindata,"%d %lf\n",counter*loop,old);
		#endif
	}
	min=old;		// sets min to current value
	
	#ifdef PLOT
	fclose(mindata);	// close file if necessary
	#endif
}
	
// sets up and runs the tsp for a given n, returning the time taken
double tsp(int n, int method, double temp, double cooling, int reject_base){
	int i;

	// sets up a temp permuatation array and the min_perm holder for the solution so far
	min_perm=malloc(n*sizeof(int));
	temp_perm=malloc(n*sizeof(int));
	for(i=0;i<n;i++){temp_perm[i]=i;min_perm[i]=i;}		// sets both initially to 1234..


	struct timeval t1, t2;
	gettimeofday(&t1, NULL);	// starts timing


	min=travel_dist(min_perm,n);	// initialises min to the distance for 1234.. permutation

	// runs one of 3 methods
	if(method==1){
		printf("USING BRUTE FORCE\n");

		// finds all permutations and distances, storing solution in min_perm
		#ifdef CLOSED
			// if closed the starting point doesn't matter
			brute_force(temp_perm, 1, n-1, n);
		#else
			brute_force(temp_perm, 0, n-1, n);
		#endif
	}
	else if(method==2){
		printf("USING SIMULATED ANNEALING\n");
		sa(temp, cooling, reject_base, n);
	}
	else if(method==3){
		printf("USING NEAREST NEIGHBOUR\n");
		greedy(n);
		min=travel_dist(min_perm,n);
	}
	gettimeofday(&t2, NULL);	// ends timing

	free(temp_perm);
	return (t2.tv_sec-t1.tv_sec+(t2.tv_usec-t1.tv_usec)/1000000.0);
}

// stores the permuted grid in a file
void grid_to_file(char *output, int n){
	FILE *plot;

	// opens up file
	if((plot=fopen(output,"w"))==NULL){
		fprintf(stderr, "*** ERROR 1 ***\nCouldn't open %s\n",output);
		exit(1);
	}

	int i,index;
	// loops through finding the next values in the permutation
	if(strcmp(ptype,"TSP")==0){
		for(i=0;i<n;i++){
			index=min_perm[i];
			fprintf(plot,"%lf %lf\n",grid[2*index],grid[2*index+1]);
		}
		#ifdef CLOSED
			fprintf(plot,"%lf %lf\n",grid[2*min_perm[0]],grid[2*min_perm[0]+1]);
		#endif
	}
	// if it's an ATSP then just print the node order
	else{
		for(i=0;i<n;i++){
			fprintf(plot,"%d\n",min_perm[i]+1);
		}
		#ifdef CLOSED
			fprintf(plot,"%d\n",min_perm[0]+1);
		#endif
	}
	fclose(plot);	
}

// prints the solution to the screen based on the min_perm
void print_sol(int n){
	int i, index;
	// prints the coordinates of the tour to the screen
	if(strcmp(ptype,"TSP")==0){
		for(i=0;i<n;i++){
			index=min_perm[i];
			printf("%lf %lf\n",grid[2*index],grid[2*index+1]);
		}
		#ifdef CLOSED
			printf("%lf %lf\n",grid[2*min_perm[0]],grid[2*min_perm[0]+1]);
		#endif
	}
	// if it's an ATSP then just print the node order
	else{
		for(i=0;i<n;i++){
			printf("%d\n",min_perm[i]+1);
		}
		#ifdef CLOSED
			printf("%d\n",min_perm[0]+1);
		#endif
	}
		
}
