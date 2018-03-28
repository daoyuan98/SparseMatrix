#include "SparseMatrix_.h"
template<class T>
sparseMatrix<T>::sparseMatrix()
{
	row_list = new int[100];
	Nrow = 0;
	Ncol = 0;
	Nele = 0;
}

template<class T>
sparseMatrix<T>::sparseMatrix(int nrow, int ncol)
{
	row_list = new int[nrow];
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
		vector<val_col>::iterator iter = val_cols_list.begin();
		for (i = row_idx; i < next_row_idx; i++, iter++) {
			if (val_cols_list[i].col >= col)
				break;
		}
		if (i == next_row_idx)
			return false;
		if (val_cols_list[i].col == col) {
			val_cols_list[i].val = val;
			return true;
		}
		else {
			val_col t(val, col);
			val_cols_list.insert(iter, t);
			return true;
		}
	}
	catch (exception e) {
		cout << e.what() << endl;
		return false;
	}

}

template<class T>
T* sparseMatrix<T>::at(int row, int col) {
	try {
		int row_idx = row_list[row];
		int next_row_idx = row_list[row + 1];
		for (int i = row_idx; i < next_row_idx; i++) {
			if (val_cols_list[i].col == col) {
				T* t = (val_cols_list[i].get());
				return t;
			}
		}
		return NULL;
	}
	catch (exception e) {
		std::cout << "invalid access" << std::endl;
		return NULL;
	}
}

template<class T>
bool sparseMatrix<T>::initializeFromVector(vector<int> rows,
	vector<int> cols, vector<T> vals) {
	_ASSERT(rows.size() == cols.size() && col.size() == vals.size());
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