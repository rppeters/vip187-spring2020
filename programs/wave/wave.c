#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gd.h>
#include <assert.h>
#include <time.h>

//array index to km
#define XIND2KM(x) ((x)*(dimx/nx))
#define YIND2KM(y) ((y)*(dimy/ny))



#define FRAMEDELAY 0
#define PHEIGHT 1
#define PWIDTH 1
#define MAXCOLORS 256
#define PRESSURERANGE 2


//globals
double c = 1.5;		//speed of wave (km/s)
double lambda = 3; 	//wave length (km)
double sigma = 4;	//width of gaussian disturbance

double dimx = 100; //metric distance of x in km
double dimy = 100; //metric distance of y in km

int nx;
int ny;
int sx;
int sy;
int timesteps = 5;
double **p_arr;
int radius = 50;

FILE *giffile;	
gdImagePtr im, previm;
int *colors;
int framecount = 0;
double *buf;

int fps = 10;
int maxdist;


void printPressureArray();

void write_frame(double **p_arr, double time) {	
	im = gdImageCreate(nx*PWIDTH,ny*PHEIGHT);
	if (time == 0) {
		colors = (int *)malloc(MAXCOLORS*sizeof(int));
		for (int j = 0; j < MAXCOLORS; j++) {
			colors[j] = gdImageColorAllocate(im, j, 0, MAXCOLORS-j-1);
		}
		gdImageGifAnimBegin(im,giffile,1,-1);
	} else {
		gdImagePaletteCopy(im, previm);
	}
	
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			int color = (int)(((1+p_arr[i][j])*MAXCOLORS)/PRESSURERANGE);

			assert(color >= 0);
			if (color >= MAXCOLORS) color = MAXCOLORS-1;
			
			gdImageFilledRectangle(im, i*PWIDTH, j*PHEIGHT, (i+1)*PWIDTH-1, (j+1)*PHEIGHT-1, colors[color]);
		}

	}

	if (time == 0) {
		//use a large frame delay to give buffer time for eog to open .gif file
		gdImageGifAnimAdd(im,giffile, 0, 0, 0, 200, gdDisposalNone, NULL);	
	} else {	
		gdImageSetPixel(im, 0, 0, framecount%2);
		gdImageGifAnimAdd(im, giffile, 0, 0, 0, FRAMEDELAY , gdDisposalNone, previm);
		gdImageDestroy(previm);
	}

	previm = im;
	im = NULL;
	framecount++;

#ifdef DEBUG
	if (framecount < 10) printPressureArray();
#endif	
}

int calculateWaveProp(double time) {	
	for (int i = 0; i < nx; i++) {
		double rangex2 = pow(XIND2KM(i-sx),2);
		for (int j = 0; j < ny; j++) {
			double range = sqrt(rangex2 + pow(YIND2KM(j-sy),2)); 				
			
			//gaussian pulse generation based on range and time
			double p0 = exp(-0.5 * pow((((c*time)-range) / sigma ),2));
			double p = p0 * cos(2*M_PI*((c*time)-range)/lambda);			

			p_arr[i][j] = p;
		}
	}
}

void printPressureArray() {
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			printf("%6.2f ",p_arr[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void init(int argc, char *argv[]) {
	nx = atoi(argv[1]);
	ny = atoi(argv[2]);
	sx = atoi(argv[3]);
	sy = atoi(argv[4]);
	timesteps = atoi(argv[5]);

	maxdist = atoi(argv[6]);
	if (nx >=  ny) {
		dimx = maxdist; 
		dimy = ny*maxdist/nx;
	} else {
		dimy = maxdist;
		dimx = nx*maxdist/ny;
	}
	long int totalPoints = nx*ny*timesteps*fps;
	printf("Total Points: %ld\n",totalPoints);

	//ignore warnings for gif opening size if gif shall not be opened
	char *filename = argv[7];
	if (filename != "DONOTOPEN.gif") assert(totalPoints < 500000000);

	giffile = fopen(argv[7],"wb");

	if (argc > 8) {
		fps = atoi(argv[8]);
		lambda = atoi(argv[9]);
		sigma = lambda * 1.25;
		printf("l=%f,s=%f\n",lambda,sigma);
	}
}

int main(int argc, char *argv[]) {
	init(argc, argv);
	//allocate array
	p_arr = (double **)malloc(nx*sizeof(double));
	for (int i = 0; i < nx; i++) {
		p_arr[i] = (double *)malloc(ny*sizeof(double));
	}

	double timeinc = (double)(1.0/fps);

	clock_t start,end;
	start = clock();
	for (double t = 0; t < timesteps; t += timeinc) { 	
		calculateWaveProp(t);
		write_frame(p_arr, t);
	}
	end = clock();
	printf("Time: %3.5f\n",((double)(end-start)/CLOCKS_PER_SEC));

	printf("Gif of size %dx%d created:\n\tNumber of Frames: %d \n\tFPS: %d \n\tTotal Time: %d\n",nx,ny,timesteps*fps,fps,timesteps);
	printf("\tLength of X Dim: %4.2f km\n\tLength of Y Dim: %4.2f km\n", dimx, dimy);
	
	//free
	for (int i = 0; i < nx; i++) {
		free(p_arr[i]);
	}

	gdImageGifAnimEnd(giffile);
	fclose(giffile);

	free(p_arr);
	free(colors);
	gdImageDestroy(previm);
}
