#ifndef COMPLEX_H
#define COMPLEX_H


class Complex{
public:
	double real;
	double im;

	Complex();
	Complex(double a, double b);
	Complex(double a);


	Complex& operator=(Complex value);
	Complex& operator=(double value);

	Complex operator+(Complex value);
	Complex operator-(Complex value);
	Complex operator*(Complex value);
	Complex operator/(Complex value);

	Complex operator+(double value);
	Complex operator-(double value);
	Complex operator*(double value);
	Complex operator/(double value);

	Complex conjugate();
};

Complex csqrt(Complex v);

#endif
