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
#define MAXBLOCKS 512
#define MAXTHREAD 512

//globals
double c = 1.5;		//speed of wave (km/s)
double lambda = 3; 	//wave length (km)
double sigma = 4;	//width of gaussian disturbance

int dimx = 100; //metric distance of x in km
int dimy = 100; //metric distance of y in km

int nx, ny, sx, sy; 
int timesteps = 5;
double *p_arr;
int radius = 50;

FILE *giffile;	
gdImagePtr im, previm;
int *colors;
int framecount = 0;
double *buf;

int fps = 10;
int maxdist;

int nthreads = 512;

__device__ int d_nx, d_ny, d_sx, d_sy, d_dimx, d_dimy;
__device__ double d_c, d_lambda, d_sigma;

__global__ void calculateWaveProp(double *arr, double time) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;

	double range = sqrt(pow((x-d_sx)*d_dimx/ d_nx, 2) + pow((y-d_sy)*d_dimy/d_ny, 2)); 				
	

	//gaussian pulse generation based on range and time
	double p0 = exp(-0.5 * pow((((d_c * time)-range) / d_sigma ),2));
	double p = p0 * cos(2 * M_PI * ((d_c * time)-range)/ d_lambda);			

	//printf("p[%d][%d]=%3.3f, range=%3.3f, p0=%3.3f,xdim=%d,ydim=%d\n",x,y,p,sqrt(pow(x-d_sx,2) + pow(y-d_sy,2)),p0,d_dimx,d_dimy);
	arr[x*d_nx + y] = p;
}



void printPressureArray();

void write_frame(double *p_arr, double time) {	
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
			int color = (int)(((1+p_arr[i*nx + j])*MAXCOLORS)/PRESSURERANGE);

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

void printArray(double* arr) {
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			printf("%6.2f ",arr[i*nx + j]);
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
}

void calculateEvenSquareDistribution(int *nxb, int *nyb) {
	int xb = (int)(sqrt(nx));
	while (xb > 0) {
		if (nx % xb == 0) break;
		xb--;
	}

	int yb = (int)(sqrt(ny));
	while (yb > 0) {
		if (ny % yb == 0) break;
		yb--;
	}
	*nxb = xb; *nyb = yb;
}

void calculateEvenMaxDistribution(int *nxb, int *nyb) {
	int xb = MAXBLOCKS;
	while (xb > 0) {
		if (nx % xb == 0 && xb < MAXBLOCKS) break;
		xb--;
	}

	int yb = MAXBLOCKS;
	while (yb > 0) {
		if (ny % yb == 0 && yb < MAXBLOCKS) break;
		yb--;
	}
	*nxb = xb; *nyb = yb;
}

void calculateMaxCacheHitsDistribution(int *nxb, int *nyb) {
	//maximize cache  (ie blocks are 1 row but split maximally)
	*nxb = nx;
	int yb;
	if (MAXBLOCKS > nx) {
		yb = (int) (sqrt(nx));
	} else {yb = nx;}
	
	while (yb > 0) {
		if (ny % yb == 0 && yb < MAXBLOCKS) break;
		yb--;
	}
	*nyb = yb;
}

int main(int argc, char *argv[]) {
	init(argc, argv);
	double *d_p_arr, *p_arr;		
	//allocate array on CPU 
	p_arr = (double *)malloc(ny*nx*sizeof(double));	

	//allocate and copy to GPU
	cudaMalloc((void**)&d_p_arr, nx*ny*sizeof(double));
	

	//copy CPU globals to GPU variables
	cudaMemcpyToSymbol(d_c, &c, sizeof(double), 0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_nx, &nx, sizeof(int), 0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_ny, &ny, sizeof(int), 0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_sx, &sx, sizeof(int), 0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_sy, &sy, sizeof(int), 0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_dimx, &dimx, sizeof(int), 0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_dimy, &dimy, sizeof(int), 0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_lambda, &lambda, sizeof(double),0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_sigma, &sigma, sizeof(double), 0,cudaMemcpyHostToDevice);
	
	//run wave propagation
	double timeinc = (double)(1.0/fps); //time increments

	int nxblocks = 0;
	int nyblocks = 0;

	//determine which block distribution scheme to use 
	int dist_type = atoi(argv[8]);
	if (dist_type == 0) calculateEvenSquareDistribution(&nxblocks,&nyblocks);
	else if (dist_type == 1) calculateEvenMaxDistribution(&nxblocks,&nyblocks);
	else if (dist_type == 2) calculateMaxCacheHitsDistribution(&nxblocks, &nyblocks);
	else assert(0);
	
	dim3 numBlocks(nxblocks,nyblocks,1);
	dim3 numThreadsPerBlock(nx/nxblocks,ny/nyblocks,1);

	clock_t start,end;
	start = clock();
	for (double t = 0; t < timesteps; t += timeinc) { 	
		calculateWaveProp<<<numBlocks,numThreadsPerBlock>>>(d_p_arr, t);
		cudaDeviceSynchronize();

		cudaMemcpy(p_arr, d_p_arr, nx*ny*sizeof(double),cudaMemcpyDeviceToHost);
		
		write_frame(p_arr, t);
	}

	
	//print timings and results
	end = clock();
	printf("Time: %3.5f\n",((double)(end-start)/CLOCKS_PER_SEC));

	printf("Gif of size %dx%d created:\n\tNumber of Frames: %d \n\tFPS: %d \n\tTotal Time: %d\n",nx,ny,timesteps*fps,fps,timesteps);
	printf("\tLength of X Dim: %d km\n\tLength of Y Dim: %d km\n", dimx, dimy);
	printf("\tX-Dim Blocks: %d\n\tY-Dim Blocks: %d\n\tX-Dim of Threads: %d\n\tY-Dim of Threads: %d\n\n\n",nxblocks, nyblocks, nx/nxblocks, ny/nyblocks);
	
	
	//free
	cudaFree(d_p_arr);
	gdImageGifAnimEnd(giffile);
	fclose(giffile);

	free(p_arr);
	free(colors);
	gdImageDestroy(previm);

	return 0;
}
