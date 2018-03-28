#include "SparseMatrix_.h"

using namespace std;

void pause() {
	char i;
	cin >> i;
}

int main(void) {
	sparseMatrix<int> mat(2, 2);
	mat.insert(1, 0, 0);
	mat.insert(2, 0, 1);
	mat.insert(3, 1, 0);
	mat.insert(4, 1, 1);
	mat.insert(22, 0, 1);
	//mat.printValCols();
	cout << mat.at(0, 0) << endl;
	cout << mat.at(0, 1) << endl;
	cout << mat.at(1, 0) << endl;
	cout << mat.at(1, 1) << endl;
	pause();
	return 0;
}