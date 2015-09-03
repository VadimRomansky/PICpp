#include "matrixElement.h"

MatrixElement::MatrixElement() {
	value = 0;
	i = 0;
	l = 0;
}

MatrixElement::MatrixElement(double v, int iv, int lv) {
	value = v;
	i = iv;
	l = lv;
}

bool MatrixElement::equalsIndex(MatrixElement& element) {
	if(i != element.i) return false;
	if(l != element.l) return false;
	
	return true;
}