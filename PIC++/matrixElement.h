#ifndef _MATRIX_ELEMENT_H_
#define _MATRIX_ELEMENT_H_

struct MatrixElement
{
	double value;
	int i;
	int j;
	int k;
	int l;

	MatrixElement();
	MatrixElement(double v, int iv, int jv, int kv, int lv);

	bool equalsIndex(MatrixElement& element);
};

#endif