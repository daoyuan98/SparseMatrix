#include "SparseMatrix_.h"

using namespace std;

void pause() {
	char i;
	cin >> i;
}

template <class T>
void printVector(vector<T> t) {
	cout << "[";
	if (t.size() == 0) { 
		cout << " ]" << endl;
		return;
	}
	for (int i = 0; i < t.size()-1; i++) {
		cout << t[i] << ", ";
	}
	cout << *(--t.end()) << "]" << endl;
}

int main(void) {
	sparseMatrix<double> mat(4, 4);

	//Gauss-Seidel Example
	//test correct!
	//vector<int> col;
	//vector<int> row;
	//vector<double> value;
	//mat.insert(10, 0, 0);
	//mat.insert(-1, 0, 1);
	//mat.insert(2, 0, 2);
	//mat.insert(-1, 1, 0);
	//mat.insert(11, 1, 1);
	//mat.insert(-1, 1, 2);
	//mat.insert(3,  1, 3);
	//mat.insert(2, 2, 0);
	//mat.insert(-1, 2, 1);
	//mat.insert(10, 2, 2);
	//mat.insert(-1, 2, 3);
	//mat.insert(3, 3, 1);
	//mat.insert(-1, 3, 2);
	//mat.insert(8, 3, 3);

	//mat.initializeFromVector(row, col, value);

	vector<double> B(4);

	//B[0] = 6; B[1] = 25; B[2] = -11; B[3] = 15;
	mat.insert(1, 0, 1);
	mat.insert(1, 0, 3);
	mat.insert(1, 1, 0);
	mat.insert(1, 1, 2);
	mat.insert(1, 2, 1);
	mat.insert(1, 3, 0);
	B[0] = 2; B[1] = 2; B[2] = 1; B[3] = 1;


	/*sparseMatrix<int> t(2, 2);
	t.insert(1, 0, 0);
	t.insert(2, 0, 1);
	t.insert(3, 1, 0);
	t.insert(4, 1, 1);
	sparseMatrix<int> r(2, 1);
	t.insert(1, 0, 0);
	t.insert(2, 0, 1);
	cout << "!!!" << endl;
	cout << t.dot(r).at(0, 1) << endl;;*/

	sparseMatrix<double> mat2(4, 4);
	mat2 = mat;
	/*cout << "mat2: 3,3  ";
	cout << mat2.at(3, 3) << endl;
	mat2.insert(3, 3, 3);
	cout << "mat2: 3,3  ";
	cout << mat2.at(3, 3) << endl;

	cout << "mat: 3,3  ";
	cout << mat.at(3, 3) << endl;*/

	//vector<double> res = mat.Gauss_Seidel_Iter(B);
	vector<double> res = mat.ConGrad(B);
	printVector<double>(res);
	pause();
	return 0;
}