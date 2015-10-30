/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 * Do not use any additional files
 * -------------------------------
 * Name: Stephanie Tong
 * Project 1: Simplified Rendering Pipeline
 * Email: stong002@ucr.edu
 * CS130: Computer Graphics
 * Fall 2015
 */


// http://www.falloutsoftware.com/tutorials/gl/gl0.htm
// http://www.songho.ca/opengl/gl_transform.html
// http://www.learnopengles.com/understanding-opengls-matrices/

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

class Matrix;
class Vertex;
class Pixel;
class BoundingBox;

stack <Matrix> ModelMatrixStack;
stack <Matrix> ProjMatrixStack;

vector<Pixel> frameBuffer;
vector<Vertex> vertexList;

MGLpixel color; // MGLpixel is unsigned int

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


void mglMultMatrix(Matrix& left, const Matrix& m);

class Vertex
{
	public:
	MGLfloat x, y, z;
	int x_screen, y_screen, z_screen;
	
	Vertex() :x(0), y(0), z(0) {}
	Vertex(MGLfloat X, MGLfloat Y, MGLfloat Z) :x(X), y(Y), z(Z) {}
	
	// convert world to screen coordinates
	// by default, points are relative to the world space
	// http://webglfactory.blogspot.com/2011/05/how-to-convert-world-to-screen.html
	// might want to deal with negatives here...
	
	/*
	void convert2ScreenCoord(MGLsize width, MGLsize height)
	{
		cout << "converting" << endl;
		
		// create the matrix representing the vertex
		Matrix m(x, y, z);
		
		cout << m << endl;
		
		// create the matrix that will scale according to screen size
		Matrix s;
		s.createScaler(width, height, 1);
		cout << "s:\n" << s << endl;
		
		// create the matrix that will translate according to screen size
		Matrix t(1, 1, 1);
		
		cout << "t:\n" << t << endl;
		
		// obtain the projection matrix from the stack
		Matrix proj;
		proj = ProjMatrixStack.top();
		cout << proj << endl;
		
		Matrix model;
		model = ModelMatrixStack.top();
		cout << model << endl;
		
		// cout << "matrixMode: " << mgl_MatrixMode << endl;

		
		//mglMultMatrix(m, t);
		mglMultMatrix(m, s);
		
		

		
		cout << "m:\n" << m << endl;
		
		
	}
	*/
	
	
	
	void convert2ScreenCoord(MGLsize width, MGLsize height)
	{
		if (x < 0)
			x_screen = width + (int) (x*width);
			
		else
			x_screen = (int) (x*width);
		
		if (y < 0)
			y_screen =  height + (int) (y*height);
		else
			y_screen = (int) (y*height);	
		
		cout <<  x << endl;
		
		cout << "x screen: " << x_screen << ", y screen: " << y_screen << endl;
		
	}
	
	
	
};


// pixel contains the coordinates, 
class Pixel
{
	public:
	int x, y;
	MGLpixel pcolor;
	
	Pixel(int X, int Y, MGLpixel c)
		: x(X), y(Y), pcolor(color)
	{}
};


class BoundingBox
{
	public:
		float min_x, max_x, min_y, max_y;
		
	BoundingBox()
		: min_x(0), max_x(0), min_y(0), max_y(0)
	{}
	
	void initBB(float width, float height)
	{
		min_x = getMin_X();
		min_y = getMin_Y();
		max_x = getMax_X();
		max_y = getMax_Y();
	}
	
	private:
	// returns the minimum value 
	float getMin_X()
	{
		int numVertices = vertexList.size();
		
		float min = vertexList[0].x_screen;
		
		for (int i = 1; i < numVertices; ++i)
		{
			if (vertexList[i].x_screen < min)
				min = vertexList[i].x_screen;
		}
		
		return min;
		
	}

	float getMin_Y()
	{
		int numVertices = vertexList.size();
		
		float min = vertexList[0].y_screen;
		
		for (int i = 1; i < numVertices; ++i)
		{
			if (vertexList[i].y_screen < min)
				min = vertexList[i].y_screen;
		}
		
		return min;
	}

	float getMax_X()
	{
		int numVertices = vertexList.size();
		
		float max = vertexList[0].x_screen;
		
		for (int i = 1; i < numVertices; ++i)
		{
			if (vertexList[i].x_screen > max)
				max = vertexList[i].x_screen;
		}
		
		return max;
	}

	float getMax_Y()
	{
		int numVertices = vertexList.size();
		
		float max = vertexList[0].y_screen;
		
		for (int i = 1; i < numVertices; ++i)
		{
			if (vertexList[i].y_screen > max)
				max = vertexList[i].y_screen;
		}
		
		return max;
	}
		
};

BoundingBox mgl_BoundingBox;

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

Matrix topMatrix()
{
	if (mgl_MatrixMode == MGL_MODELVIEW && !ModelMatrixStack.empty())
		return ModelMatrixStack.top();
	
	else if (mgl_MatrixMode == MGL_PROJECTION && !ProjMatrixStack.empty())
		return ProjMatrixStack.top();
	
}


// from lab 1

// set pixel (x,y) in framebuffer to color col, where
// col is a float array of three values between 0 and 1
// which specify the amount of red, green, and blue to mix (e.g.
// RED: (1,0,0) GREEN: (0,1,0) BLUE: (0,0,1) 
// YELLOW: (1,1,0) MAGENTA: (1,0,1) CYAN: (0,1,1)
// )
//
void set_pixel(int x, int y, MGLpixel c)
{
    Pixel pixy(x, y, c);
    frameBuffer.push_back(pixy);
}


// convert each vertex to screen coordinates
void convert2ScreenCoord(MGLsize width, MGLsize height)
{
	int numVertices = vertexList.size();
	cout << "numVert: " << numVertices << endl;
	
	for (int i = 0; i < numVertices; ++i)
		vertexList[i].convert2ScreenCoord(width, height);
}

// barycentric coordinate stuff
float f (float x, float y, float x_b, float y_b, float x_c, float y_c)
{
	return (y_b - y_c)*x + (x_c - x_b)*y + (x_b*y_c) - (x_c*y_b);
}

// a is vertexList[0];
// b is vertexList[1];
// c is vertexList[2];


// determines if points are inside the triangle
void drawTriangle(float x, float y)
{
	// vertex A
	float x_a = vertexList[0].x_screen;
	float y_a = vertexList[0].y_screen;
	
	// vertex B
	float x_b = vertexList[1].x_screen;
	float y_b = vertexList[1].y_screen;
	
	// vertex C 
	float x_c = vertexList[2].x_screen;
	float y_c = vertexList[2].y_screen;
	
	float alpha = f(x, y, x_b, y_b, x_c, y_c) / f(x_a, y_a, x_b, y_b, x_c, y_c);
	float beta = f(x, y, x_c, y_c, x_a, y_a) / f(x_b, y_b, x_c, y_c, x_a, y_a);
	float gamma = f(x, y, x_a, y_a, x_b, y_b) / f(x_c, y_c, x_a, y_a, x_b, y_b);
	
	// if it's inside the triangle
	if (alpha >= 0 && beta >= 0 && gamma >= 0)
	{
		//cout << 'a: ' << alpha << ' b: ' << beta << ' c: ' << gamma << endl;
		// MGLpixel c = alpha *color + beta*color + gamma*color;
		
		set_pixel(x, y, color);
		
		/*
		int position = ((int) y) * width + ((int) x);
		data[position] = color;
		*/
	}
	
}

void rasterizeTriangle(MGLsize width, MGLsize height)
{
	// convert vertices to screen coordinates
	convert2ScreenCoord(width, height);
	
	// obtain the bounding box in screen coordinates
	mgl_BoundingBox.initBB(width, height);
	
	float x_min = mgl_BoundingBox.min_x;
	float x_max = mgl_BoundingBox.max_x;
	float y_min = mgl_BoundingBox.min_y;
	float y_max = mgl_BoundingBox.max_y;
	
	//cout << "rasterizing triangle" << endl;
	
	for (float x = x_min; x <= x_max; ++x)
		for (float y = y_min; y <= y_max; ++y)
		{
			//cout << "x: " << x << " y: " << y << endl;
			drawTriangle(x, y);
		}
}

// for quads
// determines if points are inside the triangle
void drawTriangle1(float x, float y)
{
	// vertex A
	float x_a = vertexList[0].x_screen;
	float y_a = vertexList[0].y_screen;
	
	// vertex B
	float x_b = vertexList[1].x_screen;
	float y_b = vertexList[1].y_screen;
	
	// vertex C 
	float x_c = vertexList[2].x_screen;
	float y_c = vertexList[2].y_screen;
	
	float alpha = f(x, y, x_b, y_b, x_c, y_c) / f(x_a, y_a, x_b, y_b, x_c, y_c);
	float beta = f(x, y, x_c, y_c, x_a, y_a) / f(x_b, y_b, x_c, y_c, x_a, y_a);
	float gamma = f(x, y, x_a, y_a, x_b, y_b) / f(x_c, y_c, x_a, y_a, x_b, y_b);
	
	// if it's inside the triangle
	if (alpha >= 0 && beta >= 0 && gamma >= 0)
	{
		//cout << 'a: ' << alpha << ' b: ' << beta << ' c: ' << gamma << endl;
		// MGLpixel c = alpha *color + beta*color + gamma*color;
		
		set_pixel(x, y, color);
		
		/*
		int position = ((int) y) * width + ((int) x);
		data[position] = color;
		*/
	}

	
	
}

void drawTriangle2(float x, float y)
{
	// vertex A
	float x_a = vertexList[0].x_screen;
	float y_a = vertexList[0].y_screen;
	
	// vertex B
	float x_b = vertexList[3].x_screen;
	float y_b = vertexList[3].y_screen;
	
	// vertex C 
	float x_c = vertexList[2].x_screen;
	float y_c = vertexList[2].y_screen;
	
	float alpha = f(x, y, x_b, y_b, x_c, y_c) / f(x_a, y_a, x_b, y_b, x_c, y_c);
	float beta = f(x, y, x_c, y_c, x_a, y_a) / f(x_b, y_b, x_c, y_c, x_a, y_a);
	float gamma = f(x, y, x_a, y_a, x_b, y_b) / f(x_c, y_c, x_a, y_a, x_b, y_b);
	
	// if it's inside the triangle
	if (alpha >= 0 && beta >= 0 && gamma >= 0)
	{
		//cout << 'a: ' << alpha << ' b: ' << beta << ' c: ' << gamma << endl;
		// MGLpixel c = alpha *color + beta*color + gamma*color;
		
		set_pixel(x, y, color);
		
		/*
		int position = ((int) y) * width + ((int) x);
		data[position] = color;
		*/
	}

	
	
}

void rasterizeQuad(MGLsize width, MGLsize height)
{
	// convert vertices to screen coordinates
	convert2ScreenCoord(width, height);
	
	// obtain the bounding box in screen coordinates
	mgl_BoundingBox.initBB(width, height);
	
	float x_min = mgl_BoundingBox.min_x;
	float x_max = mgl_BoundingBox.max_x;
	float y_min = mgl_BoundingBox.min_y;
	float y_max = mgl_BoundingBox.max_y;
	
	cout << "rasterizing quad" << endl;
	
	for (float x = x_min; x <= x_max; ++x)
		for (float y = y_min; y <= y_max; ++y)
		{
			//cout << "x: " << x << " y: " << y << endl;
			drawTriangle1(x, y);
			drawTriangle2(x, y);
		}
}

/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
	// this function is called in main. it colors all of the pixels
	// that you want it to all at once. Need to store these pixels with the 
	// corresponding color somehow
	
	
	// data contains buffer of coordinates for each pixel on screen
	// data is an 1D array containing x,y coordinates of the screen
	// for ex: starting from bottom left. 
	// data[0] contains the point (0,0),
	// data[1] contains the point (1,0),
	// data[2] contains the point (2,0), etc.

	
	 int size = frameBuffer.size();
	 for (int i = 0; i < size; ++i)
	 {
		int x = frameBuffer[i].x;
		int y = frameBuffer[i].y;
		// MGLpixel c = frameBuffer[i].pcolor;
		data[y*width + x] = color;
		
	}
	

	
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
	// check what the mode is
	if (mode == MGL_TRIANGLES || mode == MGL_QUADS)
		mgl_ShapeMode = mode;

	else
		mgl_ShapeMode = -1;
}

/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
	if (mgl_ShapeMode == MGL_TRIANGLES)
	{
		rasterizeTriangle(320, 240);
		cout << "Done creating triangle" << endl;
	}
	
	else if (mgl_ShapeMode == MGL_QUADS)
	{
		rasterizeQuad(320, 240);
		cout << "Done creating rectangle" << endl;
	}
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
	Vertex Vec2(x, y, 0);
	vertexList.push_back(Vec2);
	//cout << "pushing vertex 2D" << endl;

}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
	Vertex Vec3(x, y, z);
	vertexList.push_back(Vec3);
	//cout << "pushing vertex 3D" << endl;

}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
	if (mode == MGL_MODELVIEW || mode == MGL_PROJECTION)
		mgl_MatrixMode = mode;
	else
		mgl_MatrixMode = -1;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 * 
 * glPushMatrix(); and glPopMatrix(); are used when you want to transform
 * an object, such as scaling, translating, or rotating. 
 * glPushMatrix(); is placed before the code to transform the object and
 * glPopMatrix(); is placed after the code to draw the object.
 */
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

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
	// depending on which mode you're on, you want to pop the matrix from
	// a specific stack. So you don't accidentally pop the projection matrix
	// when you wanted to pop your shape/object
	if (mgl_MatrixMode == MGL_MODELVIEW && !ModelMatrixStack.empty())
	{
		if (!ModelMatrixStack.empty())
			ModelMatrixStack.pop();
	}
	
	else if (mgl_MatrixMode == MGL_PROJECTION && !ProjMatrixStack.empty())
	{
		if (!ProjMatrixStack.empty())
			ProjMatrixStack.pop();
	}

}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
	// sets the matrix currently on the top of the stack to the identity matrix
	// this function seems to appear after mglModelView calls
	
	Matrix Identity;
	
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

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const Matrix& m)
{
	// a specific stack. So you don't accidentally modify the projection matrix
	// when you wanted to modify your shape/object
	
	// push a copy of the top matrix
	if (mgl_MatrixMode == MGL_MODELVIEW && !ModelMatrixStack.empty())
		ModelMatrixStack.push(m);
	
	else if (mgl_MatrixMode == MGL_PROJECTION && !ProjMatrixStack.empty())
		ProjMatrixStack.push(m);
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(Matrix& left, const Matrix& m)
{
	Matrix result;
	result.clearMatrix();
	
	// [ currMatrix ] [ m ]
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			for (int k = 0; k < 4; ++k)
				result.matrix[i][j] += left.matrix[i][k] * m.matrix[k][j];
				
	left = result;
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
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

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
 // http://www.dei.isep.ipp.pt/~matos/cg/docs/manual/glRotate.3G.html
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
	MGLfloat s = sin(angle * M_PI / 180);
	MGLfloat c = cos(angle * M_PI / 180);	
	
	Matrix r(0, 0, 0);
	
	r.matrix[0][0] = x*x * (1 - c) + c;
	r.matrix[0][1] = x*y * (1 - c) - z*s;
	r.matrix[0][2] = x*z * (1 - c) + y*s;
	
	r.matrix[1][0] = y*x * (1 - c) + z*s;
	r.matrix[1][1] = y*y * (1 - c) + c;
	r.matrix[1][2] = y*z * (1 - c) - x*s;
	
	r.matrix[2][0] = x*z * (1 - c) - y*s;
	r.matrix[2][1] = y*z * (1 - c) + x*s;
	r.matrix[2][2] = z*z * (1 - c) + c;
	
	if (mgl_MatrixMode == MGL_PROJECTION)
		mglMultMatrix(ProjMatrixStack.top(), r); 
	else if (mgl_MatrixMode == MGL_MODELVIEW)
		mglMultMatrix(ModelMatrixStack.top(), r);
	
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
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

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
 // this handles the 3D stuff for a camera so objects in a scene
 // apear to tilt 
 // https://www.opengl.org/sdk/docs/man2/xhtml/glFrustum.xml
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

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
// ortho positions your camera relative to the center of the screen
// https://www.opengl.org/sdk/docs/man2/xhtml/glOrtho.xml
// http://stackoverflow.com/questions/2571402/explain-the-usage-of-glortho
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
	z = 2/(far - near);
	
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

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLbyte red,
              MGLbyte green,
              MGLbyte blue)
{
	
	MGL_SET_RED(color, red);
	MGL_SET_GREEN(color, green);
	MGL_SET_BLUE(color, blue);
}
