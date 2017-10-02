#ifndef COMPLEX_H
#define COMPLEX_H

class Complex{
public:
	double re;
	double im;

	Complex();
	Complex(double a, double b);
	Complex(const Complex& a);

	Complex& operator=(const Complex& a);

	double module();
	double phase();

	Complex conjugate();

	Complex operator-(const Complex& vector);
	Complex operator+(const Complex& vector);
	Complex& operator+=(const Complex& vector);
	Complex& operator-=(const Complex& vector);
	Complex operator*(const double& value);
	Complex operator/(const double& value);
	Complex operator*(const Complex& value);
	Complex operator/(const Complex& value);
};

Complex complexExp(double phase);

#endif