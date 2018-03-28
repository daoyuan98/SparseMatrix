#include "SparseMatrix_.h"

using namespace std;

void pause() {
	char i;
	cin >> i;
}

template <class T>
void printVector(vector<T> t) {
	cout << "[";
	for (int i = 0; i < t.size()-1; i++) {
		cout << t[i] << "  ,";
	}
	cout << *(--t.end()) << "]" << endl;
}

int main(void) {
	sparseMatrix<double> mat(10, 10);
	vector<int> col;
	vector<int> row;
	vector<double> value;
	for (int i = 0; i < 10; i++) {
		col.push_back(i);
		row.push_back(i);
		value.push_back(i+1);
	}
	mat.initializeFromVector(row, col, value);
	for (int i = 0; i < 10; i++)
		cout << mat.at(i, i) << endl;


	vector<double> B(10, 1.0);
	vector<double> res = mat.Gauss_Seidel_Iter(B);
	printVector<double>(res);
	pause();
	return 0;
}