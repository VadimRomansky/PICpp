#include "particle.h"
#include "constants.h"

CorrelationMap::CorrelationMap() {
	for(int i = 0; i < splineOrder + 2; ++i) {
		xindex[i] = - splineOrder;
		yindex[i] = - splineOrder;
		zindex[i] = - splineOrder;
		xcorrelation[i] = 0;
		ycorrelation[i] = 0;
		zcorrelation[i] = 0;
	}
}
