#ifndef sparseMatrix_h
#define sparseMatrix_h
#include <list>
#include <vector>
#include <iostream>

using namespace std;
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
	bool initializeFromVector(vector<int> rows,
		vector<int> cols, vector<T> vals);

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

template<class T> T sparseMatrix<T>::at(int row, int col) {
	try {
		int row_idx = row_list[row];
		int next_row_idx = row_list[row + 1];
		for (int i = row_idx; i < next_row_idx; i++) {
			if (val_cols_list[i].col == col) {
				T t = (val_cols_list[i].get());
				return t;
			}
		}
		return 0;
	}
	catch (exception e) {
		std::cout << "invalid access" << std::endl;
		return 0;
	}
}

template<class T>
bool sparseMatrix<T>::initializeFromVector(vector<int> rows,
	vector<int> cols, vector<T> vals) {
	_ASSERT(rows.size() == cols.size() && cols.size() == vals.size());
	int n_success = 0;
	int n_error = 0;
	for (unsigned int i = 0; i < rows.size(); i++) {
		if (insert(vals[i], rows[i], cols[i]))
			n_success++;
		else
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