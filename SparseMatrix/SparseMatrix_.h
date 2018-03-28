#ifndef sparseMatrix_h
#define sparseMatrix_h
#include <list>
#include <vector>
#include <iostream>

using namespace std;

#define ConvergeLimit 0.001
#define MaxIter 1000

double ZERO = 0;

template <class T> class sparseMatrix
{
	class val_col {
	public:
		T val;
		int col;
		val_col() {
			col = 0;
		}
		val_col(const T& val, int col) {
			this->val = val;
			this->col = col;
		}

		val_col(T& t) {
			this->val = t.val;
			this->col = t.col;
		}

		T get() {
			return val;
		}
		void set(T val, int col) {
			this->val = val;
			this->col = col;
		}
		void set(T val) {
			this->val = val;
		}
	};

private:
	vector<val_col> val_cols_list;
	int* row_list;
	int Nrow;
	int Ncol;
	int Nele;
public:
	sparseMatrix();
	sparseMatrix(int nrow, int ncol);
	~sparseMatrix();
	T at(int row, int col);
	int cols();
	bool insert(const T& val, int row, int col);
	bool initializeFromVector(vector<int>& rows,
		vector<int>& cols, vector<T>& vals);

	void updateRow(int row) {
		for (int i = row + 1; i < Nrow; i++) {
			row_list[i]++;
		}
	}

	void printValCols() {
		for (int i = 0; i < val_cols_list.size(); i++) {
			cout << val_cols_list[i].col << "  " << val_cols_list[i].val << endl;
		}
	}

	vector<double> Gauss_Seidel_Iter(vector<double> B);
};

template<class T>
sparseMatrix<T>::sparseMatrix()
{
	row_list = new int[10];
	Nrow = 100;
	Ncol = 0;
	Nele = 0;
}

template<class T>
sparseMatrix<T>::sparseMatrix(int nrow, int ncol)
{
	row_list = new int[nrow + 1];
	for (int i = 0; i <= nrow; i++)
		row_list[i] = 0;
	Nrow = nrow;
	Ncol = ncol;
	Nele = 0;
}

template<class T>
sparseMatrix<T>::~sparseMatrix()
{

}

template<class T>
int sparseMatrix<T>::cols() {
	return Ncol;
}


template<class T>
bool sparseMatrix<T>::insert(const T& val, int row, int col) {
	try {
		int row_idx = row_list[row];
		int next_row_idx = row_list[row + 1];
		int i;
		vector<val_col>::iterator iter = val_cols_list.begin() + row_idx;
		for (i = row_idx; i < next_row_idx && i<val_cols_list.size(); i++, iter++) {
			if (val_cols_list[i].col >= col)
				break;
		}



		if (i == val_cols_list.size()) { // should push_back
			val_col t(val, col);
			val_cols_list.push_back(t);
			//printf("val: %d  push back\n",val);
			row_list[Nrow]++;
			updateRow(row);
			return true;
		}
		else if (col == val_cols_list[i].col) {//already exist
			val_cols_list[i].val = val;
			//printf("val: %d already exit\n",val);
			return true;
		}
		else {// common insertion
			val_col t(val, col);
			val_cols_list.insert(++iter, t);
			row_list[Nrow]++;
			//printf("val: %d  insertion\n",val);
			updateRow(row);
			return true;
		}
	}
	catch (exception e) {
		cout << e.what() << endl;
		return false;
	}

}

//The L1 distance is small, then it is converged.
static bool converged(vector<double>& v1, vector<double>& v2) {
	double res = 0.0;
	_ASSERT(v1.size() == v2.size());
	for (int i = 0; i < v1.size(); i++) {
		res += abs(v1[i] - v2[i]);
	}
	if (res <= ConvergeLimit)
		return true;
	else
		return false;
}

template<class T>
vector<double> sparseMatrix<T>::Gauss_Seidel_Iter(vector<double> B) {
	_ASSERT(Nrow == Ncol && Nrow == B.size());
	vector<double> result;
	vector<double> prev_result;
	//Nrow elements init with value of 1.
	result.assign(Nrow, 5.0);

	int iterTimes = 0;
	do {

		iterTimes++;
		//element-wise copy
		prev_result = result;
		if (iterTimes == 1)
			prev_result.assign(Nrow,5.0);

		for (int i = 0; i < Nrow; i++) {
			if (at(i, i)==ZERO)  
				continue;
			double bi = B[i];
			double sigma1 = 0.0;
			double sigma2 = 0.0;
			for (int j = 0; j < i; j++) {
				sigma1 += at(i, j)*result[j];
			}
			for (int j = i + 1; j < Nrow; j++) {
				sigma2 += at(i, j)*prev_result[j];
			}
			result[i] = (bi - sigma1 - sigma2) / at(i, i);
			cout << result[i] << endl;
		}
	} while(!converged(result, prev_result) && iterTimes <= MaxIter);

	
	cout << iterTimes << endl;

	return result;
}

template<class T> 
T sparseMatrix<T>::at(int row, int col) {
	try {
		int row_idx = row_list[row];
		int next_row_idx = row_list[row + 1];
		for (int i = row_idx; i < next_row_idx; i++) {
			if (val_cols_list[i].col == col) {
				T t = (val_cols_list[i].get());
				return t;
			}
		}
		return ZERO;
	}
	catch (exception e) {
		std::cout << "invalid access" << std::endl;
		return ZERO;
	}
}

template<class T>
bool sparseMatrix<T>::initializeFromVector(vector<int>& rows,
	vector<int>& cols, vector<T>& vals) {
	_ASSERT(rows.size() == cols.size() && cols.size() == vals.size());
	int n_error = 0;
	for (unsigned int i = 0; i < rows.size(); i++) {
		if (!insert(vals[i], rows[i], cols[i]))
			n_error++;
	}
	if (n_error == 0)
		return true;
	else {
		cout << "n_error: " << n_error << endl;
		return false;
	}
}

#endif