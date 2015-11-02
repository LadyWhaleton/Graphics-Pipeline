#include <cstdio>
#include <cstdlib>
#include "minigl.h"

#include <iostream>
#include <stack>
#include <vector>
#include "math.h"

using namespace std;

int mgl_ShapeMode;
int mgl_MatrixMode;

int HALF_WIDTH = 160;
int HALF_HEIGHT = 120;

class Matrix;
class Vertex;
class Pixel;

stack <Matrix> ModelMatrixStack;
stack <Matrix> ProjMatrixStack;

vector<Pixel> frameBuffer;
vector<Vertex> vertexList;

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
	
	Matrix operator* (const Matrix& rhs)
	{
		Matrix result;
		result.clearMatrix();
		
		// [ currMatrix ] [ m ]
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				for (int k = 0; k < 4; ++k)
					result.matrix[i][j] += matrix[i][k] * rhs.matrix[k][j];
		return result;
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
	
	void createScaler(float x, float y, float z)
	{
		clearMatrix();
		matrix[0][0] = x;
		matrix[1][1] = y;
		matrix[2][2] = z;
	}
		
	
	private:
	void initMatrix (MGLfloat X, MGLfloat Y, MGLfloat Z)
	{	
		// set diagonals
		matrix[3][3] = 1; // w
		
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

void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
	Vertex Vec3(x, y, z);
	vertexList.push_back(Vec3);
	//cout << "pushing vertex 3D" << endl;

}

// pixel contains the coordinates, 
class Pixel
{
	MGLfloat x, y;
	
	Pixel(MGLfloat X, MGLfloat Y)
		: x(X), y(Y)
	{}
};

Matrix topMatrix()
{
	if (mgl_MatrixMode == MGL_MODELVIEW && !ModelMatrixStack.empty())
		return ModelMatrixStack.top();
	
	else if (mgl_MatrixMode == MGL_PROJECTION && !ProjMatrixStack.empty())
		return ProjMatrixStack.top();
	
}

void mglMultMatrix(Matrix& left, const Matrix& m)
{
	Matrix result;
	result.clearMatrix();
	
	// [ currMatrix ] [ m ]

    MGLfloat s;
    for (int i=0;i<4;i++)
      for (int j=0;j<4;j++){
        s=0;
        for (int e=0;e<4;e++)
          s+=left.matrix[i][e]* m.matrix[e][j];

        
        result.matrix[i][j]=s;
	}
	
	left = result;
}

void mglLoadIdentity()
{
	// sets the matrix currently on the top of the stack to the identity matrix
	// this function seems to appear after mglModelView calls
	
	Matrix Identity;
	Identity.matrix[0][0] = 1;
	Identity.matrix[1][1] = 1;
	Identity.matrix[2][2] = 1;
	
	if (mgl_MatrixMode == MGL_PROJECTION)
	{

		if (!ProjMatrixStack.empty())
			ProjMatrixStack.pop();
			
		ProjMatrixStack.push(Identity);
		
	}
	else if (mgl_MatrixMode == MGL_MODELVIEW)
	{
		if (!ModelMatrixStack.empty())
			ModelMatrixStack.pop();
			
		ModelMatrixStack.push(Identity);
	}
		
}

void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
	Matrix t(x, y, z); 
	if (mgl_MatrixMode == MGL_PROJECTION)
		mglMultMatrix(ProjMatrixStack.top(), t); 
	else if (mgl_MatrixMode == MGL_MODELVIEW)
		mglMultMatrix(ModelMatrixStack.top(), t);
	
}

void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
	Matrix s;
	
	s.matrix[0][0] = x;
	s.matrix[1][1] = y;
	s.matrix[2][2] = z;
	
	if (mgl_MatrixMode == MGL_PROJECTION)
		mglMultMatrix(ProjMatrixStack.top(), s); 
	else if (mgl_MatrixMode == MGL_MODELVIEW)
		mglMultMatrix(ModelMatrixStack.top(), s);
}

void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
	MGLfloat x, y, A, B, C, D;
	
	x = (2 * near)/(right - left);
	y = (2 * near)/(top - bottom);
	A = (right + left);
	B = (top + bottom)/(top - bottom);
	C = -(far + near)/(far - near);
	D = -(2 * far)/(far - near);
	
	Matrix frustrum(0, 0, D);
	
	frustrum.matrix[0][0] = x;
	frustrum.matrix[1][1] = y;
	
	frustrum.matrix[0][2] = A;
	frustrum.matrix[1][2] = B;
	frustrum.matrix[2][2] = C;
	frustrum.matrix[3][2] = -1;
	frustrum.matrix[3][3] = 0;
	
	if (mgl_MatrixMode == MGL_PROJECTION)
		mglMultMatrix(ProjMatrixStack.top(), frustrum); 
	else if (mgl_MatrixMode == MGL_MODELVIEW)
		mglMultMatrix(ModelMatrixStack.top(), frustrum);
	
}

void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
	MGLfloat t_x, t_y, t_z, x, y, z;
	
	x = 2/(right - left);
	y = 2/(top - bottom);
	z = -2/(far - near);
	
	t_x = -(right + left)/(right - left);
	t_y = -(top + bottom)/(top-bottom);
	t_z = -(far + near)/(far - near);
	
	Matrix ortho;
	ortho.matrix[0][0] = x;
	ortho.matrix[1][1] = y;
	ortho.matrix[2][2] = z;
	
	ortho.matrix[0][3] = t_x;
	ortho.matrix[1][3] = t_y;
	ortho.matrix[2][3] = t_z;
	
	if (mgl_MatrixMode == MGL_PROJECTION)
		mglMultMatrix(ProjMatrixStack.top(), ortho); 
	else if (mgl_MatrixMode == MGL_MODELVIEW)
		mglMultMatrix(ModelMatrixStack.top(), ortho);
}

void mglMatrixMode(MGLmatrix_mode mode)
{
	if (mode == MGL_MODELVIEW || mode == MGL_PROJECTION)
		mgl_MatrixMode = mode;
	else
		mgl_MatrixMode = -1;
}

void mglPushMatrix()
{
	// depending on which mode you're on, you want to push the matrix onto
	// a specific stack. So you don't accidentally modify the projection matrix
	// when you wanted to modify your shape/object
	
	// push a copy of the top matrix
	if (mgl_MatrixMode == MGL_MODELVIEW && !ModelMatrixStack.empty())
		ModelMatrixStack.push(ModelMatrixStack.top());
	
	else if (mgl_MatrixMode == MGL_PROJECTION && !ProjMatrixStack.empty())
		ProjMatrixStack.push(ModelMatrixStack.top());
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
	Matrix m2;
	
	m2 = m1;
	
	cout << m2;
}

// test that matrix multiplication works
void test3()
{
	Matrix m1(3,4,4);
	Matrix m2(1,1,1);
	
	m2.matrix[0][0] = 5;
	m2.matrix[1][1] = 3;
	m2.matrix[2][2] = 9;
	
	mglMultMatrix(m1, m2);
	
	cout << m1;
}

// testing that translations work
void test4()
{
	mglMatrixMode(MGL_MODELVIEW);
	mglLoadIdentity();
	
	cout << "Before translation:\n" << topMatrix() << endl;
	
	mglTranslate(5, 0, 9);
	
	cout << "After translation:\n" << topMatrix() << endl;
	
}

// [model][proj][vertices]
// result = proj * vert
// then model * result


// test that vertices are being translated to screen coordinates
void test5()
{
     mglMatrixMode(MGL_PROJECTION);
     mglLoadIdentity();
     mglOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
     mglMatrixMode(MGL_MODELVIEW);
     mglLoadIdentity();
     
     Matrix m(.25, .25, 0);
     
     cout << "Proj Matrix:\n" << ProjMatrixStack.top() << endl;
	 mglMultMatrix(m, ModelMatrixStack.top());
	 mglMultMatrix(m, ProjMatrixStack.top());
	 
	 cout << m << endl;
     
   
}

void test6()
{
	mglMatrixMode(MGL_MODELVIEW);
    mglLoadIdentity();
    
	mglTranslate(2,1,1);
	
	cout << "\nModel Matrix:\n" << ModelMatrixStack.top() << endl;
}

int main()
{
	//test1();
	//test2();
	//test3();
	//test4();
	test5();
	//test6();
	return 0;
}
