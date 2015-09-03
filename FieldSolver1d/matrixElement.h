#ifndef _MATRIX_ELEMENT_H_
#define _MATRIX_ELEMENT_H_

struct MatrixElement
{
	double value;
	int i;
	int l;

	MatrixElement();
	MatrixElement(double v, int iv, int lv);

	bool equalsIndex(MatrixElement& element);
};

#endif