#include "stdlib.h"
#include <math.h>
#include "complex.h"
#include "constants.h"


Complex::Complex(){
	real = 0;
	im = 0;
}

Complex::Complex(double a, double b){
	real = a;
	im = b;
}

Complex::Complex(double a){
	real = a;
	im = 0;
}

Complex& Complex::operator=(Complex value){
	real = value.real;
	im = value.im;
	return *this;
}

Complex& Complex::operator=(double value){
	real = value;
	im = 0;
	return *this;
}

Complex Complex::operator+(Complex a){
	Complex result = Complex(real + a.real, a.im + im);
	return result;
}

Complex Complex::operator-(Complex a){
	Complex result = Complex(real - a.real, im - a.im );
	return result;
}

Complex Complex::operator*(Complex a){
	Complex result = Complex(real*a.real - im*a.im, real*a.im + im*a.real);
	return result;
}

Complex Complex::operator/(Complex a){
	double module2 = a.real*a.real + a.im*a.im;
	Complex result = Complex((real*a.real + im*a.im)/module2, (-real*a.im + im*a.real)/module2);
	return result;
}

Complex Complex::operator+(double a){
	Complex result = Complex(real + a, im);
	return result;
}

Complex Complex::operator-(double a){
	Complex result = Complex(real - a, im);
	return result;
}

Complex Complex::operator*(double a){
	Complex result = Complex(real*a, im*a);
	return result;
}

Complex Complex::operator/(double a){
	Complex result = Complex(real/a, im/a);
	return result;
}

Complex Complex::conjugate(){
	return Complex(real, -im);
}

Complex csqrt(Complex a){
	double module = sqrt(sqrt(a.real*a.real + a.im*a.im));
	double phi = atan2(a.im, a.real);
	if (phi < 0){
		phi += 2*pi;
	}

	phi = phi/2;

	Complex result = Complex(module*cos(phi), module*sin(phi));
	return result;
}