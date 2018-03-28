#include "SparseMatrix_.h"

using namespace std;

void pause() {
	char i;
	cin >> i;
}

int main(void) {
	sparseMatrix<int> mat(20, 20);
	vector<int> col;
	vector<int> row;
	vector<int> value;
	for (int i = 0; i < 20; i++) {
		col.push_back(i);
		row.push_back(i);
		value.push_back(i);
	}
	mat.initializeFromVector(row, col, value);

	for (int i = 19; i >= 0; i--) {
		cout << mat.at(i, i) << endl;
	}

	pause();
	return 0;
}