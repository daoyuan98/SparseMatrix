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
	sparseMatrix<double> mat(4, 4);
	vector<int> col;
	vector<int> row;
	vector<double> value;
	mat.insert(10, 0, 0);
	mat.insert(-1, 0, 1);
	mat.insert(2, 0, 2);
	mat.insert(-1, 1, 0);
	mat.insert(11, 1, 1);
	mat.insert(-1, 1, 2);
	mat.insert(3,  1, 3);
	mat.insert(2, 2, 0);
	mat.insert(-1, 2, 1);
	mat.insert(10, 2, 2);
	mat.insert(-1, 2, 3);
	mat.insert(3, 3, 1);
	mat.insert(-1, 3, 2);
	mat.insert(8, 3, 3);

	mat.initializeFromVector(row, col, value);

	vector<double> B(4);

	B[0] = 6; B[1] = 25; B[2] = -11; B[3] = 15;

	vector<double> res = mat.Gauss_Seidel_Iter(B);
	printVector<double>(res);
	pause();
	return 0;
}