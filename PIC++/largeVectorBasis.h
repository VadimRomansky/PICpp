#ifndef LARGEVECTORBASiS_H
#define LARGEVECTORBASIS_H

class LargeVectorBasis {
public:
	int size;
	int xnumber;
	int ynumber;
	int znumber;
	int lnumber;
	int capacity;
	double***** array;

	LargeVectorBasis(int sizev, int xnumberv, int ynumberv, int znumberv, int lnumberv);
	~LargeVectorBasis();

	void resize(int capacityv);
	void clear();
};

#endif