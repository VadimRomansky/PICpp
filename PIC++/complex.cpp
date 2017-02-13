#include "math.h"
#include "stdlib.h"
#include "stdio.h"
//#include <crtdbg.h>

//#include "memory_debug.h"

#include "complex.h"

Complex::Complex(){
	re = 0;
	im = 0;
}

Complex::Complex(double a, double b){
	re = a;
	im = b;
}

Complex::Complex(const Complex& a){
	re = a.re;
	im = a.im;
}

Complex& Complex::operator=(const Complex& a) {
	re = a.re;
	im = a.im;

	return *this;
}

double Complex::module(){
	return sqrt(re*re + im*im);
}

double Complex::phase(){
	return atan2(im, re);
}

Complex Complex::conjugate(){
	return Complex(re, -im);
}

Complex Complex::operator-(const Complex& a){
	return Complex(re - a.re, im - a.im);
}

Complex Complex::operator+(const Complex& a){
	return Complex(re + a.re, im + a.im);
}

Complex& Complex::operator+=(const Complex& a){
	re += a.re;
	im += a.im;
	return *this;
}

Complex& Complex::operator-=(const Complex& a){
	re -= a.re;
	im -= a.im;
	return *this;
}

Complex Complex::operator*(const double& a){
	return Complex(re*a, im*a);
}

Complex Complex::operator/(const double& a){
	return Complex(re/a, im/a);
}

Complex Complex::operator*(const Complex& a){
	double newRe = re*a.re - im*a.im;
	double newIm = re*a.im + im*a.re;
	return Complex(newRe, newIm);
}

Complex Complex::operator/(const Complex& a){
	double mod2 = a.re*a.re + a.im*a.im;
	double newRe = (re*a.re + im*a.im)/mod2;
	double newIm = (im*a.re - re*a.im)/mod2;
	return Complex(newRe, newIm);
}

Complex complexExp(double phase){
	return Complex(cos(phase), sin(phase));
}