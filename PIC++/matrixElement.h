#ifndef _MATRIX_ELEMENT_H_
#define _MATRIX_ELEMENT_H_

class MatrixElement
{
public:
	double value;
	int i;
	int j;
	int k;
	int l;

	MatrixElement();
	MatrixElement(const double& v, int iv, int jv, int kv, int lv);

	bool equalsIndex(const MatrixElement& element);
};

#endif