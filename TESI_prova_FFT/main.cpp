#include "mbed.h"
#include <complex>
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define MAX 1024

#define M_PI 3.1415926535897932384

using namespace std;

Serial pc(USBTX,USBRX);

int log2(int N)    /*function to calculate the log2(.) of int numbers*/
{
  int k = N, i = 0;
  while(k) {
    k >>= 1;
    i++;
  }
  return i - 1;
}

int check(int n)    //checking if the number of element is a power of 2
{
  return n > 0 && (n & (n - 1)) == 0;
}

int reverse(int N, int n)    //calculating revers number
{
  int j, p = 0;
  for(j = 1; j <= log2(N); j++) {
    if(n & (1 << (log2(N) - j)))
      p |= 1 << (j - 1);
  }
  return p;
}

void ordina(complex<double>* f1, int N) //using the reverse order in the array
{
  complex<double> f2[MAX];
  for(int i = 0; i < N; i++)
    f2[i] = f1[reverse(N, i)];
  for(int j = 0; j < N; j++)
    f1[j] = f2[j];
}

void transform(complex<double>* f, int N) //
{
  ordina(f, N);    //first: reverse order
  complex<double> *W;
  W = (complex<double> *)malloc(N / 2 * sizeof(complex<double>));
  W[1] = polar(1., -2. * M_PI / N);
  W[0] = 1;
  for(int i = 2; i < N / 2; i++)
    W[i] = pow(W[1], i);
  int n = 1;
  int a = N / 2;
  for(int j = 0; j < log2(N); j++) {
    for(int i = 0; i < N; i++) {
      if(!(i & n)) {
        complex<double> temp = f[i];
        complex<double> Temp = W[(i * a) % (n * a)] * f[i + n];
        f[i] = temp + Temp;
        f[i + n] = temp - Temp;
      }
    }
    n *= 2;
    a = a / 2;
  }
}

void FFT(complex<double>* f, int N, double d)
{
  transform(f, N);
  for(int i = 0; i < N; i++)
    f[i] *= d; //multiplying by step
}

int main()
{
  srand(time(NULL));
  int n = MAX;
  /*do {
    pc.printf("Specifiy array dimension (MUST be a power of 2)\n\r");
    pc.scanf("%d",&n);
  } while(!check(n)); */
  double d = 1;
  //pc.printf("Specify the sampling step\n\r");
  //pc.scanf("%lf",&d);
  
  complex<double> vec[MAX];
  
  double Fs = 1000;
  double T = 1/Fs;
  double t;  
  //pc.printf("Specify the array\n\r");
  for(int i = 0; i < n; i++) {
    t = i*T;
    vec[i] = complex<double>(0.7*sin(2*M_PI*50*t) + sin(2*M_PI*120*t) /*+ 2*(rand()%RAND_MAX)*/,0);
    //vecc[i] = 0.7*sin(2*M_PI*50*t) + sin(2*M_PI*120*t) + 2*(1+rand()%100); //funzione di esempio
    wait_ms(10); //prende un campione ogni 10 ms - 100 Hz
    //S[i] = vecc[i] + 2*AWGN_Generator();
  }
  
  
  /*for(int i = 0; i < n; i++) {  
    //pc.printf("vec[%d] = %f+i%f\n\r",i,vec[i].real(),vec[i].imag());
    pc.printf("vec[%d] = %f+i%f\n\r",i,vec[i].real(),vec[i].imag());
  }*/
  
  clock_t begin = clock();
  
  FFT(vec, n, d);
  
  //clock_t end = clock();
  //double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  
  //pc.printf("Elapsed time (seconds): %f for %d samples \n\r",time_spent,MAX);
  
  
  //pc.printf("Printing the fft of the array specified.\n\r");
  /*for(int j = 0; j < n; j++){
    //pc.printf("FFT[%d]: %f+i%f\n\r",j,vec[j].real(),vec[j].imag());
    pc.printf("FFT[%d]: %f+i%f\n\r",j,vec[j].real(),vec[j].imag());
  }*/
  
  complex<double> vecc[MAX];
  
  for(int j = 0; j < n; j++){
      vecc[j] = complex<double>(vec[j].real()/MAX,vec[j].imag()/MAX); 
  }
  
  double S[MAX];
  for(int j = 0; j < n; j++){
      S[j] = abs(vecc[j]); 
  }
  
  /*for(int j = 0; j < n; j++){
      pc.printf("S[%d]: %f\n\r",j,S[j]); 
  }*/
  double positiveband[(MAX/2)+1];
  for(int j = 1; j < (MAX/2)+1; j++){
      positiveband[j] = S[j]; 
  }
  
  for(int j = 2; j < (MAX/2)+1; j++){
      positiveband[j] = 2*positiveband[j]; 
  }
  
  /*for(int j = 0; j < (MAX/2)+1; j++){
      pc.printf("positiveband[%d]: %f\n\r",j,positiveband[j]); 
      //pc.printf("%f\n\r",positiveband[j]); 
  }*/
  
  double frequencies[(MAX/2)];
  for(int j = 0; j < (MAX/2); j++){
      frequencies[j] = (Fs*j)/MAX; 
  }
  //pc.printf("\n\n");
  
  /*for(int j = 0; j < (MAX/2); j++){
      pc.printf("frequencies[%d]: %f\n\r",j,frequencies[j]); 
      //pc.printf("%f\n\r",frequencies[j]);
  }*/
  
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  
  pc.printf("Elapsed time (seconds): %f for %d samples \n\r",time_spent,MAX);
  
  return 0;
}

