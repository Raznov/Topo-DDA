#include "definition.h"

SiCi::SiCi() {
	ifstream SiInput;
	ifstream CiInput;
	SiInput.open("./SiCi/Si.txt");
	CiInput.open("./SiCi/Ci.txt");
	SiInput >> numberSi;
	CiInput >> numberCi;
	SiInput >> disSi;
	CiInput >> disCi;
	Si = VectorXd::Zero(numberSi);
	Ci = VectorXd::Zero(numberCi);
	for (int i = 0; i <= numberSi - 1; i++) {
		SiInput >> Si(i);
		CiInput >> Ci(i);
	}
	SiInput.close();
	CiInput.close();
}

double SiCi::get_Si(double x) {
	int pos1 = floor(x / disSi) - 1;
	int pos2 = ceil(x / disSi) - 1;
	double val1 = Si(pos1);
	double val2 = Si(pos2);
	double axis1 = (pos1 + 1) * disSi;
	double axis2 = (pos2 + 1) * disSi;
	return val1 + ((val2 - val1) / (axis2 - axis1)) * (x - axis1);

}

double SiCi::get_Ci(double x) {
	int pos1 = floor(x / disCi) - 1;
	int pos2 = ceil(x / disCi) - 1;
	double val1 = Ci(pos1);
	double val2 = Ci(pos2);
	double axis1 = (pos1 + 1) * disCi;
	double axis2 = (pos2 + 1) * disCi;
	return val1 + ((val2 - val1) / (axis2 - axis1)) * (x - axis1);

}