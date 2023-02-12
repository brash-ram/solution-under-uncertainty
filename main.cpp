#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>
using namespace std;

typedef vector<vector<float>> Matrix;

struct MatrixWithVector
{
	Matrix matrix;
	vector<float> vectorPriority;
};



long numberOfDigits(double n) {//для красивого вывода матриц
	std::ostringstream strs;

	strs << n;
	return strs.str().size();
}

void printMatrix(const Matrix& M) {//красивый вывод матрицы
	long max_len_per_column[M.size()];
	long n = M.size(), m = M[0].size();

	for (long j = 0; j < m; ++j) {
		long max_len {};

		for (long i = 0; i < n; ++i){
			long num_length {numberOfDigits(M[i][j])};
			if (num_length > max_len)
				max_len = num_length;
		}

		max_len_per_column[j] = max_len;
	}

	for (long i = 0; i < n; ++i)
		for (long j = 0; j < m; ++j)
			std::cout << (j == 0 ? "\n| " : "") << std::setw(max_len_per_column[j]) << M[i][j] << (j == m - 1 ? " |" : " ");

	std::cout << '\n';
}

void initializeMatrix(Matrix& mat) {//заполняем матрицу еденицами
	for (int i=0; i<mat.size(); i++){
		for (int j=0; j<mat.size(); j++){
			mat[i][j] = 1;
		}
	}
}

float parseNumber(string line) {
	return (float) (line[0]-'0') / (float) (line[2]-'0');//получаем  числа из строки типа 1\2
}

Matrix enterMatrix(int numberParams) {//ввод матрицы особым способом
	Matrix matrix(numberParams, vector<float>(numberParams));
	string number;
	initializeMatrix(matrix);
	for (int i=0; i<numberParams-1; i++) {//вводятся построчно элементы верхней половины матрцы
		for (int j=i+1; j<numberParams; j++) {//т.е. всех элментов выше главное диагонали
			cout << "Relation " << i+1 << "\\" << j+1 << endl;
			cin >> number;
			if (number.size() == 1){
				matrix[i][j] = stof(number);//выше главное дигонали
				matrix[j][i] = 1/stof(number);//симметричное ему ниже гл диагонали обратное число
			}
			else {// если число введено в виде 1\2
				float numberFromInput = parseNumber(number);
				matrix[i][j] = numberFromInput;
				matrix[j][i] = 1/numberFromInput;
			}
		}
	}
	return matrix;
}

void createVectorPriority(MatrixWithVector& matrix) {//рассчитывается вектор приоритетов
	vector<float> vectorGeometricMean(matrix.matrix.size());
	
	for (int i = 0; i < matrix.matrix.size(); i++) {
		float geometricMean = 1.0;
		int size =  matrix.matrix[0].size();

		for (int j = 0; j < size; j++) {
			geometricMean *= matrix.matrix[i][j];
		}
		vectorGeometricMean[i] = pow(geometricMean, (1/ (float) size));
	}
	float normalizingFactor = 0.0;
	for (int i = 0; i < vectorGeometricMean.size(); ++i) {
		normalizingFactor += vectorGeometricMean[i];
	}
	matrix.vectorPriority = vector<float>(matrix.matrix.size());
	for (int i = 0; i < matrix.matrix.size(); i++) {
		matrix.vectorPriority[i] = vectorGeometricMean[i]/normalizingFactor;
	}
}

void printCR(float CR) {
    if (CR < 0.1) {
        cout << "CR = " << CR << " - Ok (" << CR << "<0.1)" << endl;
    } else {
        cout << "CR = " << CR << " - Bad (" << CR << ">0.1)" << endl;
    }
}

float calculateCR(const MatrixWithVector& M) {//тут считается отношение согласованности как в контрольное примере
    Matrix matrix = M.matrix;
    
    if (matrix.size() <= 2) {
        return -1;
    }

    float l = 0;
    for (int i = 0; i < matrix.size(); i++) {
        float s = 0;
        for (int j = 0; j < matrix[i].size(); j++) {
            s += 1 / matrix[i][j];
        }
        // cout << "s" << i << " = " << s << endl;
        l += s * M.vectorPriority[i];
    }
    // cout << "l" << " = " << l << endl;
    float CI = (l - matrix.size())/(matrix.size() - 1);
    cout << "CI" << " = " << CI << endl;
    float RI = (1.98*(matrix.size()-2))/matrix.size();
    cout << "RI" << " = " << RI << endl;
    return CI / RI;
}


void printMatrixWithVector(const MatrixWithVector& M) {//выводит матрицу с вектором 
    printMatrix(M.matrix);
	for (int j = 0; j < M.vectorPriority.size(); j++) {
		cout << "q" << j+1 << " - " << M.vectorPriority[j] << "\t";
	}
    cout << endl;
    printCR(calculateCR(M));
    cout << endl;
}

int main() {
	cout << "Enter number criteria" << endl;
    int numberK; cin >> numberK;
	cout << "Enter number of alternative" << endl;
	int numberA; cin >> numberA;

	vector<MatrixWithVector> listAllMatrix(numberK+1);

	cout << "Enter matrix of criteria" << endl;
	listAllMatrix[0].matrix = enterMatrix(numberK);
	createVectorPriority(listAllMatrix[0]);
	printMatrixWithVector(listAllMatrix[0]);

	
	for (int i = 1; i < numberK+1; i++) {
		cout << "Enter matrix of alternative on criteria "  << i << endl;
		listAllMatrix[i].matrix = enterMatrix(numberA);
		createVectorPriority(listAllMatrix[i]);
		printMatrixWithVector(listAllMatrix[i]);
	}

	cout << "Full matrixes" << endl;
	for (int i = 0; i < numberK+1; ++i) {
		printMatrixWithVector(listAllMatrix[i]);
		cout << endl;
	}

    return 0;
}