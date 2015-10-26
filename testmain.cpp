#include <cstdio>
#include <cstdlib>
#include "minigl.h"

#include <iostream>
#include <stack>
#include <vector>
#include "math.h"

using namespace std;

// classes I created
class Matrix
{
	// variables
	public:
	MGLfloat matrix[4][4];
	
	// by default, it's the identity matrix
	Matrix()
	{
		clearMatrix();
		initMatrix (0, 0, 0);
	}
	
	Matrix(MGLfloat X, MGLfloat Y, MGLfloat Z)
	{
		clearMatrix();
		initMatrix (X, Y, Z);
	}
	
	Matrix& operator= (const Matrix& rhs)
	{
		// this is for avoiding assignment of the same object
		if (this != &rhs) 
		{
			for (int row = 0; row < 4; ++row)
				for (int col = 0; col < 4; ++col)
					matrix[row][col] = rhs.matrix[row][col];
		}
		
		return *this;
	}
	
	friend ostream& operator<< (ostream& os, const Matrix& m)
	{
		for (int row = 0; row < 4; ++row)
		{
			for (int col = 0; col < 4; ++col)
				os << m.matrix[row][col] << " ";
			os << "\n";
		}
		return os;
	}
	
	void clearMatrix()
	{
		for (int row = 0; row < 4; ++row)
			for (int col = 0; col < 4; ++col)
				matrix[row][col] = 0;
	}		
	
	private:
	void initMatrix (MGLfloat X, MGLfloat Y, MGLfloat Z)
	{	
		// set diagonals
		matrix[0][0] = 1;
		matrix[1][1] = 1;
		matrix[2][2] = 1;
		matrix[3][3] = 1;
		
		// set the x, y, z coordinate
		matrix[0][3] = X;
		matrix[1][3] = Y; 
		matrix[2][3] = Z;
	}

};


class Vertex
{
	public:
	MGLfloat x, y, z;
	
	Vertex() :x(0), y(0), z(0) {}
	Vertex(MGLfloat X, MGLfloat Y, MGLfloat Z) :x(X), y(Y), z(Z) {}
};


// pixel contains the coordinates, 
class Pixel
{
	MGLfloat x, y;
	
	Pixel(MGLfloat X, MGLfloat Y)
		: x(X), y(Y)
	{}
};

Matrix currMatrix;


void mglMultMatrix(const Matrix& m)
{
	Matrix result;
	result.clearMatrix();
	
	// [ matrix ] [currMatrix]
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			for (int k = 0; k < 4; ++k)
				result.matrix[i][j] += currMatrix.matrix[i][k] * m.matrix[k][j];
				
	currMatrix = result;
}

void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
	Matrix t(x, y, z); 
	mglMultMatrix(t);
	
}

// test that matrix is set up properly
void test1()
{
	Matrix m1(2,3,4);
	cout << m1 << endl;
}

// test that assignment operator works
void test2()
{
	Matrix m1 (22, 11, 4);
	currMatrix = m1; 
	
	cout << currMatrix;
}

// test that matrix multiplication works
void test3()
{
	Matrix m1(3,4,4);
	Matrix m2(1,1,1);
	
	currMatrix = m1;
	
	m2.matrix[0][0] = 5;
	m2.matrix[1][1] = 3;
	m2.matrix[2][2] = 9;
	
	mglMultMatrix(m2);
	
	cout << currMatrix;
}

void test4()
{
	Matrix m1(3, 1, 3);
	currMatrix = m1;
	
	cout << "Before translation:\n" << currMatrix << endl;
	
	mglTranslate(5, 0, 9);
	
	cout << "After translation:\n" << currMatrix << endl;
	
}

int main()
{
	//test1();
	//test2();
	//test3();
	test4();
	return 0;
}
