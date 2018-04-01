#ifndef sparseMatrix_h
#define sparseMatrix_h
#include <list>
#include <vector>
#include <iostream>
using namespace std;

#define ConvergeLimit 10
#define MaxIter 100

int ZERO = 0;

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
	
	int Nrow;
	int Ncol;
	int Nele;
public:
	int* row_list;
	vector<val_col> val_cols_list;
	sparseMatrix();
	sparseMatrix(int nrow, int ncol);
	sparseMatrix(vector<T> v);
	~sparseMatrix();
	T at(int row, int col);
	int cols() {
		return Ncol;
	}
	int rows() {
		return Nrow;
	}
	bool insert(const T& val, int row, int col);
	bool initializeFromVector(vector<int>& rows,
		vector<int>& cols, vector<T>& vals);

	bool isSymmetric(void);
	bool isPosDefinite(void);

	sparseMatrix dot(sparseMatrix<T> mat2);

	vector<T> dot_v(vector<T> v);
	vector<T> v_dot(vector<T> v);

	void updateRow(int row) {
		for (int i = row + 1; i < Nrow; i++) {
			row_list[i]++;
		}
		row_list[Nrow] = 0x7fffffff;
	}

	void printValCols() {
		for (int i = 0; i < val_cols_list.size(); i++) {
			cout << val_cols_list[i].col << "  " << val_cols_list[i].val << endl;
		}
	}

	void clear() {
		Nrow = Ncol = 0;
		memset(row_list, 0, sizeof(row_list));
		for (int i = 0; i < val_cols_list.size(); i++)
			val_cols_list.pop_back();
	}

	vector<double> Gauss_Seidel_Iter(vector<double> B);
	vector<double> ConGrad(vector<double> B);

	vector<T> toVec();

	sparseMatrix<T> operator-(sparseMatrix<T> m2);
	sparseMatrix<T> operator+(sparseMatrix<T> m2);
	sparseMatrix<T> operator*(double fac);
	sparseMatrix<T> operator=(sparseMatrix<T> m2);
	sparseMatrix<T> Transpose();
};

template<class T>
sparseMatrix<T> sparseMatrix<T>::operator=(sparseMatrix<T> m) {
	Nrow = m.rows();
	Ncol = m.cols();
	row_list = new int[Nrow + 1];
	for (int i = 0; i < Nrow; i++)
		row_list[i] = m.row_list[i];
	int n = val_cols_list.size();
	val_cols_list.clear();
	n = m.val_cols_list.size();
	for (int i = 0; i < n; i++) {
		val_cols_list.push_back(m.val_cols_list[i]);
	}
	/*for (int i = 0; i < Nrow; i++) {
		for (int j = 0; j < Ncol; j++) {
				this->insert(m.at(i, j), i, j);
		}
	}*/
	return *this;
}

//矩阵相乘
template<class T>
sparseMatrix<T> sparseMatrix<T>::dot(sparseMatrix<T> mat) {
	_ASSERT(Ncol = mat.rows());
	int row = Nrow;
	int col = mat.cols();
	/*cout << Nrow << Ncol << endl;
	cout << mat.rows() << mat.cols() << endl;*/
	sparseMatrix<T> r_m(row,col);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			double res = 0;
			for (int k = 0; k < Ncol; k++) {
				if(at(i,k)&&mat.at(k,j))
					res += at(i, k)*mat.at(k, j);
			}
			r_m.insert(res, i, j);
		}
	}
	return r_m;
}

//向量相乘
//template<class T>
//static double inner_dot<T>(vector<T> v1, vector<T> v2) {
//	double res = 0;
//	_ASSERT(v1.size() == v2.size());
//	for (int i = 0; i < v1.size(); i++)
//		res += v1[i] * v2[i];
//	return res;
//}

//矩阵×向量
template<class T>
vector<T> sparseMatrix<T>::dot_v(vector<T> v) {
	_ASSERT(Ncol == v.size());
	vector<T> r_v(Nrow, 0);
	for (int i = 0; i < Nrow; i++) {
		double res = 0.0;
		for (int j = 0; j < Ncol; j++) {
			if (!at(i, j))
				res += at(i, j)*v[j];
		}
		r_v[i] = res;
	}
	return r_v;
}

template<class T>
vector<T> sparseMatrix<T>::v_dot(vector<T> v) {
	_ASSERT(v.size() == Nrow);
	vector<T> r_v(Ncol, 0);
	for (int i = 0; i < Ncol; i++) {
		double val = 0.0;
		for (int j = 0; j < v.size(); j++) {
			if (!at(j, i))
				val += at(j, i)*v[j];
		}
		r_v[i] = val;
	}
	return r_v;
}

//矩阵转置
template<class T>
sparseMatrix<T> sparseMatrix<T>::Transpose() {
	sparseMatrix<T> res(this->cols(), this->rows());
	for (int i = 0; i < Nrow; i++) {
		for (int j = 0; j < Ncol; j++) {
			double val = at(i, j);
			if(val)
				res.insert(val, j, i);
		}
	}
	return res;
}

//矩阵减法
template<class T>
sparseMatrix<T> sparseMatrix<T>::operator-(sparseMatrix<T> m) {
	//printf("%d %d\n%d %d\n", Ncol, Nrow, m.cols(), m.rows());
	_ASSERT_EXPR(Ncol == m.cols() && Nrow == m.rows(),
		"");
	sparseMatrix<T> res(m.rows(), m.cols());
	for (int i = 0; i < m.rows(); i++) {
		for (int j = 0; j < m.cols(); j++) {
			res.insert(this->at(i, j) - m.at(i, j), i, j);
		}
	}
	return res;
}

template<class T>
sparseMatrix<T> sparseMatrix<T>::operator*(double d) {
	for (int i = 0; i < Nrow; i++) {
		for (int j = 0; j < Ncol; j++) {
			if (at(i, j))
				insert(d*at(i, j), i, j);
		}
	}
	return *this;
}

sparseMatrix<double> res(4,1);
//矩阵加法
template<class T>
sparseMatrix<T> sparseMatrix<T>::operator+(sparseMatrix<T> m) {
	_ASSERT(this->Ncol == m.cols() && this->Nrow == m.rows());
	res.clear();
	for (int i = 0; i < m.rows(); i++) {
		for (int j = 0; j < m.cols(); j++) {
			T v1 = at(i, j);
			T v2 = m.at(i, j);
			//cout << v1 << " " << v2 << endl;
			T r = v1 + v2;
			if (r) {
				//cout << "insert val: " << r << endl;
				this->insert(r, i, j);
			}
		}
	}
	//printf("ans vec:\n");
	//printVector(this->toVec());
	return *this;
}

template<class T>
vector<T> sparseMatrix<T>::toVec() {
	if (Nrow == 1) {
		vector<T> v(Ncol);
		for (int i = 0; i < Ncol; i++ )
			v[i] = at(0, i);
		return v;
	}
	else if (Ncol == 1) {
		vector<T> v(Nrow);
		for (int i = 0; i < Nrow; i++)
			v[i] = at(i, 0);
		return v;
	}
}

template<class T>
bool sparseMatrix<T>::isSymmetric(void) {
	
	if (Nrow != Ncol)
		return false;

	int i,j, n;
	n = Nrow;
	for (i = 0; i < Nrow; i++) {
		for (j = i + 1; j < Nrow; j++) {
			if (at(i, j) != at(j, i))
				return false;
		}
	}
	return true;
}

template<class T>
bool sparseMatrix<T>::isPosDefinite(void) {
	if (!isSymmetric())
		return false;
	else {
		//....
		return true;
	}
}


template<class T>
sparseMatrix<T>::sparseMatrix()
{
	sparseMatrix(10, 10);
}

template<class T>
sparseMatrix<T>::sparseMatrix(int nrow, int ncol)
{
	row_list = new int[nrow + 1];
	for (int i = 0; i <= nrow; i++)
		row_list[i] = 0;
	val_cols_list.clear();
	Nrow = nrow;
	Ncol = ncol;
	Nele = 0;
}


template<class T>
sparseMatrix<T>::sparseMatrix(vector<T> v) {
	row_list = new int[2];
	row_list[0] = row_list[1] = 0;
	//val_cols_list.push_back(val_col(0, 0));
	Nrow = 1;
	Ncol = v.size();
	for (int i = 0; i < v.size(); i++)
		if (v[i])
			insert(v[i], 0, i);
}

template<class T>
sparseMatrix<T>::~sparseMatrix()
{

}

template<class T>
bool sparseMatrix<T>::insert(const T& val, int row, int col) {
	try {
		int row_idx = row_list[row];
		int next_row_idx = row_list[row + 1];
		int i = 0;
		//vector<val_col>::iterator iter = val_cols_list.begin() + row_idx;
		vector<val_col>::iterator iter = val_cols_list.begin();
		
		for (i = row_idx; i < next_row_idx && i<val_cols_list.size()&&iter!=val_cols_list.end(); i++, iter++) {
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
			iter = iter + row_idx;
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
static bool converged(vector<T> v) {
	double res = 0.0;
	for (int i = 0; i < v.size(); i++)
		res += v[i];
	if (res <= ConvergeLimit)
		return true;
	else
		return false;
}

template<class T>
vector<double> sparseMatrix<T>::ConGrad(vector<double> B) {
	/*if (!this->isSymmetric())
		return vector<double>();*/
	
	//real calculation
	sparseMatrix b(B);
	sparseMatrix<double> x(Nrow, 1);
	for (int i = 0; i < Nrow; i++) {
		x.insert(0, i, 0);
	}
	sparseMatrix<double> r(1, B.size());
	sparseMatrix<double> p(1, B.size());
	sparseMatrix<double> new_r(1, B.size());
	r = b - dot(x).Transpose();
	p = r;
	int niter = 0;
	do {
		niter++;
		double rkTrk = (r.dot(r.Transpose()).at(0, 0));	
		//cout << "rkTrk： "<<rkTrk << endl;
		double div = (p.Transpose().dot(*this).dot(p)).at(0, 0);
		double alpha;
		if(div==0)	alpha =  1 ;
		else alpha = rkTrk / div;
		//cout << alpha << endl;
		//printVector(x.toVec());
		/*cout << "inter" << endl;
		printVector(inter.toVec());*/
		x = x + p.Transpose()*alpha;
		//cout << "x" << endl;
		//printVector(x.toVec());

		//cout << "x" << endl;
		new_r = r - (dot(p.Transpose())).Transpose()*alpha;
		//printVector(new_r.toVec());
		//cout << "after new r" << endl;
		double beta = new_r.dot(new_r.Transpose()).at(0, 0) / rkTrk;
		//cout << beta << endl;
		p = new_r + p*beta;
		r = new_r;
	} while (!converged<T>(r.toVec())&&niter<=MaxIter);
	return x.toVec();
}

template<class T>
vector<double> sparseMatrix<T>::Gauss_Seidel_Iter(vector<double> B) {
	_ASSERT(Nrow == Ncol && Nrow == B.size());
	vector<double> result;
	vector<double> prev_result;

	//Nrow elements init with value of 1.
	result.assign(Nrow, 1.0);

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
		if (row == Nrow - 1) next_row_idx = val_cols_list.size();
		if (row_idx < 0) return ZERO;
		for (int i = row_idx; i < next_row_idx; i++) {
			if (val_cols_list[i].col == col) {
				T t = (val_cols_list[i].get());
				if (t < 0) t = -t;
				//cout <<"at: "<< t << endl;
				return t;
			}
		}
		//cout << "at: " << ZERO << endl;
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