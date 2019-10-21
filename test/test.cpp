#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <chrono>
#include "dinterpl.h"

using namespace std;

#define N 6000
#define M 200000

/*
MAIN PROGRAM
*/

int main(int argc, char *argv[])
{

double ax[N]; //array being read in
double ay1[N]; //array being read in
double ay2[N]; //array being read in
ifstream myfile; myfile.open("it_10");
if (myfile.is_open())
	{
        int i=0;
		while ( true)
		{
			double x,y1,y2;
			myfile >> x;
			myfile >> y1;
			myfile >> y2;
			if (myfile.eof()) break;
			ax[i] = x;
			ay1[i] = y1;
  		ay2[i] = y2;
  		i++;
		}
	}
myfile.close();
/*
{
std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
dinterpl interpl(ax,ay2,N);
std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Avg Time CPU (Cubic init) = " << std::chrono::duration_cast<std::chrono::microseconds>((end - begin)).count() << "[µs]" << std::endl;
  std::cout << "Avg Time CPU (Cubic init) = " << std::chrono::duration_cast<std::chrono::nanoseconds> ((end - begin)).count() << "[ns]" << std::endl;
}
*/

dinterpl interpl(ax,ay2,N);
//Use linear interpolation at discontinuities
const double threshold = 500.0;
interpl.cspline_filter(threshold);


gsl_set_error_handler_off();

//Linear Test and benchmark


{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  ofstream file;
  file.open ("lin_interp.txt");

  //#pragma omp parallel for
  for (int i=0;i<=M;i++){
    double x = -6.0+i*12.0/M;
  	double test = interpl.linear_cached_eval(x);
    file << x << ' ' << test << '\n';
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Tot Time CPU (linear internal cache) = " << std::chrono::duration_cast<std::chrono::microseconds>((end - begin)).count() << "[µs]" << std::endl;
  std::cout << "Tot Time CPU (linear internal cache) = " << std::chrono::duration_cast<std::chrono::nanoseconds> ((end - begin)).count() << "[ns]" << std::endl;
  std::cout << "Avg Time CPU (linear internal cache) = " << std::chrono::duration_cast<std::chrono::microseconds>((end - begin)/(M+1)).count() << "[µs]" << std::endl;
  std::cout << "Avg Time CPU (linear internal cache) = " << std::chrono::duration_cast<std::chrono::nanoseconds> ((end - begin)/(M+1)).count() << "[ns]" << std::endl;

  file.close();
}

{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  //ofstream file;
  //file.open ("lin_interp.txt");

  //#pragma omp parallel for
  for (int i=0;i<=M;i++){
    double x = -6.0+i*12.0/M;
  	double test = interpl.linear_eval(x);
    //file << x << ' ' << test << '\n';
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Tot Time CPU (linear no cache) = " << std::chrono::duration_cast<std::chrono::microseconds>((end - begin)).count() << "[µs]" << std::endl;
  std::cout << "Tot Time CPU (linear no cache) = " << std::chrono::duration_cast<std::chrono::nanoseconds> ((end - begin)).count() << "[ns]" << std::endl;
  std::cout << "Avg Time CPU (linear no cache) = " << std::chrono::duration_cast<std::chrono::microseconds>((end - begin)/(M+1)).count() << "[µs]" << std::endl;
  std::cout << "Avg Time CPU (linear no cache) = " << std::chrono::duration_cast<std::chrono::nanoseconds> ((end - begin)/(M+1)).count() << "[ns]" << std::endl;

  //file.close();
}

{
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, N);
  gsl_spline_init(spline,ax,ay2,N);
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (int i=0;i<=M;i++){
    double x = -6.0+i*12.0/M;
  	double test = gsl_spline_eval(spline,x,acc);
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Avg Time CPU (gsl linear w/ cache) = " << std::chrono::duration_cast<std::chrono::microseconds>((end - begin)/(M+1)).count() << "[µs]" << std::endl;
  std::cout << "Avg Time CPU (gsl linear w/ cache) = " << std::chrono::duration_cast<std::chrono::nanoseconds> ((end - begin)/(M+1)).count() << "[ns]" << std::endl;

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
}

{
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, N);
  gsl_spline_init(spline,ax,ay2,N);
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (int i=0;i<=M;i++){
    double x = -6.0+i*12.0/M;
  	double test = gsl_spline_eval(spline,x,0);
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Avg Time CPU (gsl linear w/o cache) = " << std::chrono::duration_cast<std::chrono::microseconds>((end - begin)/(M+1)).count() << "[µs]" << std::endl;
  std::cout << "Avg Time CPU (gsl linear w/o cache) = " << std::chrono::duration_cast<std::chrono::nanoseconds> ((end - begin)/(M+1)).count() << "[ns]" << std::endl;

  gsl_spline_free (spline);
}


//Cubic Test

{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  
  //ofstream file;
  //file.open ("cubic_interp.txt");

  for (int i=0;i<=M;i++){
    double x = -6.0+i*12.0/M;
  	double test = interpl.cspline_cached_eval(x);
    //file << x << ' ' << test << '\n';
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Tot Time CPU (linear internal cache) = " << std::chrono::duration_cast<std::chrono::microseconds>((end - begin)).count() << "[µs]" << std::endl;
  std::cout << "Tot Time CPU (linear internal cache) = " << std::chrono::duration_cast<std::chrono::nanoseconds> ((end - begin)).count() << "[ns]" << std::endl;

  std::cout << "Avg Time CPU (cubic internal cache) = " << std::chrono::duration_cast<std::chrono::microseconds>((end - begin)/(M+1)).count() << "[µs]" << std::endl;
  std::cout << "Avg Time CPU (cubic internal cache) = " << std::chrono::duration_cast<std::chrono::nanoseconds> ((end - begin)/(M+1)).count() << "[ns]" << std::endl;

  //file.close();
}

{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  
  ofstream file;
  file.open ("cubic_interp.txt");

  for (int i=0;i<=M;i++){
    double x = -6.0+i*12.0/M;
  	double test = interpl.cspline_eval(x);
    file << x << ' ' << test << '\n';
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Tot Time CPU (cubic no cache) = " << std::chrono::duration_cast<std::chrono::microseconds>((end - begin)).count() << "[µs]" << std::endl;
  std::cout << "Tot Time CPU (cubic no cache) = " << std::chrono::duration_cast<std::chrono::nanoseconds> ((end - begin)).count() << "[ns]" << std::endl;

  std::cout << "Avg Time CPU (cubic no cache) = " << std::chrono::duration_cast<std::chrono::microseconds>((end - begin)/(M+1)).count() << "[µs]" << std::endl;
  std::cout << "Avg Time CPU (cubic no cache) = " << std::chrono::duration_cast<std::chrono::nanoseconds> ((end - begin)/(M+1)).count() << "[ns]" << std::endl;

  file.close();
}

{
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, N);
  gsl_spline_init(spline,ax,ay2,N);
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (int i=0;i<=M;i++){
    double x = -6.0+i*12.0/M;
  	double test = gsl_spline_eval(spline,x,acc);
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Avg Time CPU (gsl cspline w/ cache) = " << std::chrono::duration_cast<std::chrono::microseconds>((end - begin)/(M+1)).count() << "[µs]" << std::endl;
  std::cout << "Avg Time CPU (gsl cspline w/ cache) = " << std::chrono::duration_cast<std::chrono::nanoseconds> ((end - begin)/(M+1)).count() << "[ns]" << std::endl;

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
}

{
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, N);
  gsl_spline_init(spline,ax,ay2,N);
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (int i=0;i<=M;i++){
    double x = -6.0+i*12.0/M;
  	double test = gsl_spline_eval(spline,x,0);
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Avg Time CPU (gsl cspline w/o cache) = " << std::chrono::duration_cast<std::chrono::microseconds>((end - begin)/(M+1)).count() << "[µs]" << std::endl;
  std::cout << "Avg Time CPU (gsl cspline w/o cache) = " << std::chrono::duration_cast<std::chrono::nanoseconds> ((end - begin)/(M+1)).count() << "[ns]" << std::endl;

  gsl_spline_free (spline);
}

return 0;
}
