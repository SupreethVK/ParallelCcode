#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <omp.h>

#define PI acos(-1.0)  //to be precise

typedef double complex dcom;

double calc_time(struct timespec start, struct timespec stop) 
{
	double t;
	t = (stop.tv_sec - start.tv_sec) * 1000; 
	t += (stop.tv_nsec - start.tv_nsec) * 0.000001; 
	return t;
}


void display(dcom * array, int len)
{
	for(int i=0; i<len; i++)
	{
		printf("%g, %g\n", creal(array[i]), cimag(array[i])); 
	}
}

// if flag is 0, it takes only the even elements, else if flag is 1 it takes only the odd ones 
dcom *sub_array(dcom * array, int len, int flag)
{
	int c = 0;
	dcom * result = malloc((len/2)*sizeof(dcom));

	for(int i=0+flag; i<len; i=i+2)
	{
		result[c] = array[i];
		c++;
	}
	return result;
}

dcom * fastFourierTrans(dcom * inputArray, int len)
{
	
	
	if(len == 1)
  {
		return inputArray;
	}
	
	dcom omega = 1;
	dcom omegaN = cexp(2*PI*I/len); //the fourier coefficient

	dcom * evenSubArray;
	dcom * oddSubArray;

	#pragma omp parallel sections
	{
		#pragma omp section
		{
			evenSubArray = fastFourierTrans(sub_array(inputArray, len, 0), len/2);	
		}
		#pragma omp section
		{
			oddSubArray = fastFourierTrans(sub_array(inputArray, len, 1), len/2);	
		}
	}
	
	
	
	dcom * outputArray = malloc(len*sizeof(dcom));

	#pragma omp parallel for
	for(int i=0; i<(len/2); i++)
	{
		outputArray[i] = evenSubArray[i] + omega * oddSubArray[i];
		outputArray[i+(len/2)] = evenSubArray[i] - omega * oddSubArray[i];
		omega = omegaN*omega;
	}
	free(evenSubArray);
	free(oddSubArray);
	return outputArray;
}

int main()
{
		struct timespec start, stop;
		int arrayLength, temp;
		scanf("%d", &arrayLength);
		
		if(arrayLength%2==0)
		{
			printf("%d\n", arrayLength);
		
			//int * inpInt = (int *) malloc(arrayLength*sizeof(int));

			dcom * inputArray = (dcom *) malloc(arrayLength*sizeof(dcom));
			
			for(int i=0; i<arrayLength; i++)
			{	
				scanf("%d", &temp);
				inputArray[i] = temp;
				//printf("%g, %g\n", creal(inputArray[i]), cimag(inputArray[i]));
			}
			
			//printf("Input Array: \n");
			//display(inputArray, arrayLength);

			clock_t start_t,end_t;
			start_t=clock();
			clock_gettime(CLOCK_REALTIME,&start);
			dcom * outputArray = fastFourierTrans(inputArray, arrayLength);
			clock_gettime(CLOCK_REALTIME,&stop);
			end_t=clock();
			printf("%lf milliseconds\n",calc_time(start,stop));
			printf("Clock Ticks: %ld\n",end_t-start_t);

			//printf("Output Array: \n");
			//display(outputArray, arrayLength);
			

			free(inputArray);
			free(outputArray);

		}
		else
		{
			printf("Input array length must be a multiple of 2\n");
		}
		return 0;
}
