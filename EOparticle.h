#pragma once
class EOparticle {
public:
	double* value;
	void setDim(int d){
		value = new double[d];
	}
};
