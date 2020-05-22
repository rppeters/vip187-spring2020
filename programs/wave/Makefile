LINKS = -lgd -lm
CC = cc -pedantic -std=c11
CCCUDA = nvcc
CUDARUN = srun -n 1 --gres=gpu:1

#	       nx ny sx sy sec maxdist outputfile
CONFIG_SMALL = 20 20 1 1 10 100 small.gif
CONFIG_LARGE = 200 200 50 50 10 100 large.gif
CONFIG_PACIFIC = 200 200 0 500 1000 7500 pacific.gif 1 200 
CONFIG_TEST = 20 20 1 1 10 10 test.gif 
CONFIG_CUDASMALL = 20 20 0 0 10 100 small.gif
CONFIG_HUGE = 500 500 0 370 10 100 huge.gif
CONFIG_MASSIVE = 2000 2000 250 250 10 100 DONOTOPEN.gif
CONFIG_VERYHUGE = 1000 1000 0 0 10 100 DONOTOPEN.gif
CONFIG_BIGGEST = 5000 5000 0 0 10 100 DONOTOPEN.gif

wave_debug : wave.c
	$(CC) -DDEBUG -o wave wave.c $(LINKS)

wave : wave.c
	$(CC) -o wave wave.c $(LINKS)

wave_small : wave.c
	./wave $(CONFIG_SMALL)

wave_large : wave.c
	./wave $(CONFIG_LARGE)

wave_huge : wave.c
	./wave $(CONFIG_HUGE)

#If you open this generated gif, you may/will crash from VM :)
#This make is purely for testing large runtimes
wave_massive: wave.c
	./wave $(CONFIG_MASSIVE)

wave_pacific : wave.c
	./wave $(CONFIG_PACIFIC)

#CUUUUUUUDA
wave_cuda.exec : wave.cu
	$(CCCUDA) -o wave_cuda.exec wave.cu $(LINKS) 

wave_cuda_test : wave_cuda.exec
	$(CUDARUN) wave_cuda.exec $(CONFIG_TEST)

wave_cuda_small : wave_cuda.exec
	$(CUDARUN) wave_cuda.exec $(CONFIG_CUDASMALL) 0 
	$(CUDARUN) wave_cuda.exec $(CONFIG_CUDASMALL) 1
	$(CUDARUN) wave_cuda.exec $(CONFIG_CUDASMALL) 2 

wave_cuda_large : wave_cuda.exec
	$(CUDARUN) wave_cuda.exec $(CONFIG_LARGE) 0
	$(CUDARUN) wave_cuda.exec $(CONFIG_LARGE) 1
	$(CUDARUN) wave_cuda.exec $(CONFIG_LARGE) 2

wave_cuda_huge : wave_cuda.exec
	$(CUDARUN) wave_cuda.exec $(CONFIG_HUGE) 0 
	$(CUDARUN) wave_cuda.exec $(CONFIG_HUGE) 1
	$(CUDARUN) wave_cuda.exec $(CONFIG_HUGE) 2 

wave_cuda_massive : wave_cuda.exec
	$(CUDARUN) wave_cuda.exec $(CONFIG_MASSIVE) 0
	$(CUDARUN) wave_cuda.exec $(CONFIG_MASSIVE) 1
	$(CUDARUN) wave_cuda.exec $(CONFIG_MASSIVE) 2

wave_cuda_veryhuge : wave_cuda.exec
	$(CUDARUN) wave_cuda.exec $(CONFIG_VERYHUGE) 0
	$(CUDARUN) wave_cuda.exec $(CONFIG_VERYHUGE) 1
	$(CUDARUN) wave_cuda.exec $(CONFIG_VERYHUGE) 2


wave_cuda_biggest : wave_cuda.exec
	$(CUDARUN) wave_cuda.exec $(CONFIG_BIGGEST) 0
	$(CUDARUN) wave_cuda.exec $(CONFIG_BIGGEST) 1
	$(CUDARUN) wave_cuda.exec $(CONFIG_BIGGEST) 2

clean:
	rm -f *.gif wave_cuda.exec wave *.mp4 
