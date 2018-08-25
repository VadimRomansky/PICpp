//#include <crtdbg.h>

//#include "memory_debug.h"
#include "particle.h"
#include "constants.h"
#include "paths.h"

CorrelationMap::CorrelationMap() {
	for (int i = 0; i < splineOrder + 2; ++i) {
		xindex[i] = i;
		yindex[i] = i;
		zindex[i] = i;
		xcorrelation[i] = 0;
		ycorrelation[i] = 0;
		zcorrelation[i] = 0;
	}
	xcorrelation[0] = 1.0;
	ycorrelation[0] = 1.0;
	zcorrelation[0] = 1.0;
}
