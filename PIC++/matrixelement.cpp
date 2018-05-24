//#include <crtdbg.h>

//#include "memory_debug.h"
#include "matrixElement.h"
#include "paths.h"

MatrixElement::MatrixElement() {
	value = 0;
	i = 0;
	l = 0;
}

MatrixElement::MatrixElement(const double& v, int iv, int jv, int kv, int lv) {
	value = v;
	i = iv;
	j = jv;
	k = kv;
	l = lv;
}

bool MatrixElement::equalsIndex(const MatrixElement& element) {
	if (i != element.i) return false;
	if (j != element.j) return false;
	if (k != element.k) return false;
	if (l != element.l) return false;

	return true;
}
