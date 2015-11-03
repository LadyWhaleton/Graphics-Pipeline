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

#include <cstdio>
#include <cstdlib>
#include "minigl.h"

#include <iostream>
#include <stack>
#include <vector>
#include <algorithm>
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
vector<Pixel> zBuffer;

vector< vector<Vertex> > shapeList;
vector<Vertex> vertexList;

void mglMultMatrix(Matrix& left, const Matrix& m);

MGLpixel RGB[3]; // MGLpixel is unsigned int

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

// classes I created
class Matrix
{
	// variables
	public:
	MGLfloat matrix[4][4];
	
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
	
	// Currently unused. mglMultMatrix is used instead.
	Matrix operator* (const Matrix& rhs)
	{
		Matrix result;
		result.clearMatrix();
		
		// need to reset something to 0.
		MGLfloat sum = 0;
		
		// [ currMatrix ] [ m ]
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				sum = 0;
				for (int k = 0; k < 4; ++k)
					sum += matrix[i][k] * rhs.matrix[k][j];
				result.matrix[i][j] = sum;
			}
		}
			
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
		matrix[3][3] = 1;
	}
	
	void createTranslater(MGLfloat X, MGLfloat Y, MGLfloat Z)
	{
		matrix[0][0] = 1;
		matrix[1][1] = 1;
		matrix[2][2] = 1;
		matrix[3][3] = 1;
		
		matrix[0][3] = X;
		matrix[1][3] = Y;
		matrix[2][3] = Z;
	}	
	
	private:
	void initMatrix (MGLfloat X, MGLfloat Y, MGLfloat Z)
	{	
		// set w
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
	MGLfloat x, y, z, w;
	MGLfloat x_screen, y_screen, z_screen, w_screen;
	MGLpixel vColor[3];
	
	Vertex() :x(0), y(0), z(0), w(1)
	{
		vColor[0] = RGB[0];
		vColor[1] = RGB[1];
		vColor[2] = RGB[2];
	}
	
	Vertex(MGLfloat X, MGLfloat Y, MGLfloat Z, MGLfloat W) :x(X), y(Y), z(Z), w(W)
	{
		vColor[0] = RGB[0];
		vColor[1] = RGB[1];
		vColor[2] = RGB[2];
	}
	
	friend ostream& operator<< (ostream& os, const Vertex& v)
	{
		os << "( " << v.x << ", " << v.y << ", " << v.z << ", " << v.w << " ) " << endl;
		return os;
		
	}
	
	// This is the assignment operator I have overloaded.
	Vertex& operator= (const Vertex& rhs)
	{
		// this is for avoiding assignment of the same object
		if (this != &rhs) 
		{
			x = rhs.x;
			y = rhs.y;
			z = rhs.z;
			w = rhs.w;
			vColor[0] = rhs.vColor[0];
			vColor[1] = rhs.vColor[1];
			vColor[2] = rhs.vColor[2];
			
			x_screen = rhs.x_screen;
			y_screen = rhs.y_screen;
			z_screen = rhs.y_screen;
			w_screen = rhs.w_screen;
		}
		
		return *this;
	}
	
	/* This operation performs matrix and vector multiplication.
	 * The vector must be on the left of *, and the matrix must be on
	 * the right.
	 * 
	 * It follows the matrix * vector order of multiplication.
	 */
	Vertex operator* (const Matrix& rhs)
	{
		MGLfloat currVertex[4] = {x, y, z, w};
		MGLfloat newVertex[4] = {0, 0, 0, 0};
		
		MGLfloat sum, val = 0;
		
		// [ Vertex ] [ Matrix ]
		for (int row = 0; row < 4; ++row)
		{
			sum = 0;
			for (int col = 0; col < 4; ++col)
			{
				val = currVertex[col];
				sum += rhs.matrix[row][col] * val;	
			}
			
			newVertex[row] = sum;
			//cout << "row: " << row << " sum: " << sum << endl;
		}
		
		Vertex v (newVertex[0], newVertex[1], newVertex[2], newVertex[3]);
		return v;
	}
	
	/* This function converts world to screen coordinates.
	 * It should be called when you're about to rasterize the
	 * triangle or quadrilateral.
	 * 
	 * By now, I should know the screen width and height because
	 * it has been passed through mglReadPixels. 
	 */	
	void scaleToScreen(MGLsize width, MGLsize height)
	{			
		x = (x*width)/2;
		y = (y*height)/2;
		
		// convert to NDC, aka divide by w
		x = x/w;
		y = y/w;
		z = z/w;
		w = w/w;
		
		if (x > width)
			x = width;
		
		x_screen = x;
		y_screen = y;
		z_screen = z;
		w_screen = w;
		
		cout << "x: "  << x_screen << ", " << "y: " << y << endl;
	}
	
	/* This function, applyTransformation, applys transformations to
	* the vertex (x, y, z, w).
	* 
	* For my implementation, the order of multiplication is
	* 1. vector * model matrix
	* 2. step 1 * projection matrix
	* 3. step 2 * translation matrix ( this part handles negative values )
	* 
	* Each step of my multiplication returns a vector.
	* 
	* After that, I update the x, y, and z values.
	* 
	* I don't scale to the screen yet because I don't know the screen
	* resolution until mglReadPixels (which is called at the very very
	* end).
	*/ 
	void applyTransformations()
	{
		Matrix model = ModelMatrixStack.top();
		Matrix proj = ProjMatrixStack.top();
		
		Matrix trans;
		trans.createTranslater(1, 1, 1);
			
		Vertex v(x,y,z,w);			
				
		v = v * model;
		v = v * proj;
		
		// you multiply the vector by the translation matrix 
		// to handle negatives
		v = v * trans;
		
		// update the x, y, z, w values.
		x = v.x;
		y = v.y;
		z = v.z;
		w = v.w;
	}
};

class Pixel
{
	public:
	int x, y;
	MGLpixel pColor;
	MGLfloat z;
	
	Pixel(int X, int Y, MGLpixel c, MGLfloat Z)
		: x(X), y(Y), pColor(c), z(Z)
	{}
};


class BoundingBox
{
	public:
		float min_x, max_x, min_y, max_y;
		
	BoundingBox()
		: min_x(0), max_x(0), min_y(0), max_y(0)
	{}
	
	void initBB(const vector<Vertex>& vl)
	{
		min_x = getMin_X(vl);
		min_y = getMin_Y(vl);
		max_x = getMax_X(vl);
		max_y = getMax_Y(vl);
	}
	
	private:
	// returns the minimum value 
	float getMin_X(const vector<Vertex>& vl)
	{
		int numVertices = vl.size();
		float min = vl[0].x_screen;
		
		for (int i = 1; i < numVertices; ++i)
		{
			if (vl[i].x_screen < min)
				min = vl[i].x_screen;
		}
		
		return min;
	}

	float getMin_Y(const vector<Vertex>& vl)
	{
		int numVertices = vl.size();
		float min = vl[0].y_screen;
		
		for (int i = 1; i < numVertices; ++i)
		{
			if (vl[i].y_screen < min)
				min = vl[i].y_screen;
		}
		
		return min;
	}

	float getMax_X(const vector<Vertex>& vl)
	{
		int numVertices = vl.size();	
		float max = vl[0].x_screen;
		
		for (int i = 1; i < numVertices; ++i)
		{
			if (vl[i].x_screen > max)
				max = vl[i].x_screen;
		}
		
		return max;
	}

	float getMax_Y(const vector<Vertex>& vl)
	{
		int numVertices = vl.size();
		float max = vl[0].y_screen;
		
		for (int i = 1; i < numVertices; ++i)
		{
			if (vl[i].y_screen > max)
				max = vl[i].y_screen;
		}
		
		return max;
	}
		
};

/* set_pixel takes in (x, y) which determines where on the screen
 * we want to draw the pixel. 
 * 
 * z is used for the zBuffer, which will 
 * later be sorted to give the illusion of layers / clipping.
 * 
 * c is the color we will use for the pixel at (x, y).
 */ 
void set_pixel(int x, int y, MGLpixel c, MGLfloat z)
{
    Pixel pixy(x, y, c, z);
    frameBuffer.push_back(pixy);
    zBuffer.push_back (pixy);
}


// convert each vertex to screen coordinates
void convert2ScreenCoord(MGLsize width, MGLsize height, vector <Vertex>& v)
{
	int numVertices = v.size();
	
	for (int i = 0; i < numVertices; ++i)
		v[i].scaleToScreen(width, height);
}

// barycentric coordinate stuff
float f (float x, float y, float x_b, float y_b, float x_c, float y_c)
{
	return (y_b - y_c)*x + (x_c - x_b)*y + (x_b*y_c) - (x_c*y_b);
}

// returns the new color 
MGLpixel mixColors(float alpha, float beta, float gamma, const MGLpixel* vColor1, const MGLpixel* vColor2, const MGLpixel* vColor3)
{
	float newRed = alpha*vColor1[0] + beta*vColor2[0]+ gamma*vColor3[0];
	float newGreen = alpha*vColor1[1] + beta*vColor2[1] + gamma*vColor3[1]; 
	float newBlue = alpha*vColor1[2] + beta*vColor2[2] + gamma*vColor3[2]; 
	
	MGLpixel newColor = 0;
	MGL_SET_RED(newColor, (MGLpixel) newRed);
	MGL_SET_GREEN(newColor, (MGLpixel) newGreen);
	MGL_SET_BLUE(newColor, (MGLpixel) newBlue);
	
	return newColor;	
}

// determines if points are inside the triangle
void drawTriangle(float x, float y, const Vertex& a, const Vertex &b, const Vertex &c)
{
	// vertex A
	float x_a = a.x_screen;
	float y_a = a.y_screen;
	
	// vertex B
	float x_b = b.x_screen;
	float y_b = b.y_screen;
	
	// vertex C 
	float x_c = c.x_screen;
	float y_c = c.y_screen;
	
	float alpha = f(x, y, x_b, y_b, x_c, y_c) / f(x_a, y_a, x_b, y_b, x_c, y_c);
	float beta = f(x, y, x_c, y_c, x_a, y_a) / f(x_b, y_b, x_c, y_c, x_a, y_a);
	float gamma = f(x, y, x_a, y_a, x_b, y_b) / f(x_c, y_c, x_a, y_a, x_b, y_b);
	
	// If it's inside the triangle
	if (alpha >= 0 && beta >= 0 && gamma >= 0)
	{
		MGLpixel pixelColor = mixColors(alpha, beta, gamma, a.vColor, b.vColor, c.vColor);
		MGLfloat newZ = alpha*a.z_screen + beta*b.z_screen + gamma*c.z_screen;
		
		set_pixel( (int) x, (int) y, pixelColor, newZ);
	}
	
}

/* vl contains the list of vertices.
 * width is the screen width, and height is the screen height.
 * 
 * After we set up the bounding box, we look at all of the pixels between
 * x_min, x_max, y_min, and y_max via a for loop. Inside the loop, I call
 * drawTriangle which preps up the values to be used in barycentric coordinates.
 */
void rasterizeTriangle(MGLsize width, MGLsize height, const vector<Vertex>& vl)
{
	BoundingBox mgl_BoundingBox;
	
	// obtain the bounding box in screen coordinates
	mgl_BoundingBox.initBB(vl);
	
	float x_min = mgl_BoundingBox.min_x;
	float x_max = mgl_BoundingBox.max_x;
	float y_min = mgl_BoundingBox.min_y;
	float y_max = mgl_BoundingBox.max_y;
	
	cout << "Rasterizing Triangle" << endl;
	
	for (float x = x_min; x <= x_max; ++x)
		for (float y = y_min; y <= y_max; ++y)
			drawTriangle(x, y, vl[0], vl[1], vl[2]);
}

/* vl contains the list of vertices.
 * width is the screen width, and height is the screen height.
 * 
 * After we set up the bounding box, we look at all of the pixels between
 * x_min, x_max, y_min, and y_max via a for loop. Inside the loop, I call
 * drawTriangle which preps up the values to be used in barycentric coordinates.
 * 
 * I call drawTriangle twice because you can represent a quadrilateral
 * by creating two triangles that share two vertices.
 */ 
void rasterizeQuad(MGLsize width, MGLsize height, const vector<Vertex>& vl)
{
	BoundingBox mgl_BoundingBox;
	
	// obtain the bounding box in screen coordinates
	mgl_BoundingBox.initBB(vl);
	
	float x_min = mgl_BoundingBox.min_x;
	float x_max = mgl_BoundingBox.max_x;
	float y_min = mgl_BoundingBox.min_y;
	float y_max = mgl_BoundingBox.max_y;
	
	cout << "Rasterizing Quad" << endl;
	
	for (float x = x_min; x <= x_max; ++x)
		for (float y = y_min; y <= y_max; ++y)
		{
			drawTriangle(x, y, vl[0], vl[1], vl[2]);
			drawTriangle(x, y, vl[0], vl[3], vl[2]);
		}
}


// want to sort Z by descending order
bool sortByZ(Pixel a, Pixel b) { return a.z > b.z; }

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
	/* This function is called in main. It colors all of the pixels
	 * that you want it to all at once.
	 * data is a 1D of coordinates for each pixel starting from the 
	 * bottom lefthand corner of the screen.
	 * for example:
	 * data[0] contains the point (0,0),
	 * data[1] contains the point (1,0),
	 * data[2] contains the point (2,0), etc.
	 * data[y*width] contains row y.
	 * data[y*width + x] indicates the column x.
	 * data[width*height] is the top right corner of the screen.
	 */
	
	int numShapes = shapeList.size();
	
	cout << "Num of Shapes: " << numShapes << endl << endl;
	
	// for each shape
	for (int i = 0; i < numShapes; ++i)
	{
		// convert vertices of shape i to screen coordinates
		convert2ScreenCoord(width, height, shapeList[i]);
		
		// get the number of vertices of the current shape
		int numVertices = shapeList[i].size();
	
		cout << "Num of Vertices of Shape " << i+1 << ": " << numVertices << endl;
		
		// Rasterize the shape according to the number of vertices.
		if (numVertices == 3)
			rasterizeTriangle(width, height, shapeList[i]);
		else if (numVertices == 4)
			rasterizeQuad(width, height, shapeList[i]);
			
		cout << endl;
	}	

	// sort the zBuffer by descending order.
	sort(zBuffer.begin(), zBuffer.end(), sortByZ);

	int size = zBuffer.size();
	 
	 for (int i = 0; i < size; ++i)
	 {
		int x = zBuffer[i].x;
		int y = zBuffer[i].y;
		
		MGLpixel color = zBuffer[i].pColor;
		data[y*width + x] = color;
		
	}

	/*
	 int size = frameBuffer.size();
	 
	 for (int i = 0; i < size; ++i)
	 {
		int x = frameBuffer[i].x;
		int y = frameBuffer[i].y;
		
		MGLpixel c = frameBuffer[i].pcolor;
		data[y*width + x] = c;
		
	}
	*/
	
	shapeList.clear();
		
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
	// Add all of these vertices to a Shape
	if(!vertexList.empty())
	{
		shapeList.push_back(vertexList);
	
		vertexList.clear();
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
	if (mgl_ShapeMode == -1)
	{
		MGL_ERROR("Error: Missing mglBegin! Aborting mission.\n");
		exit(1);
	}
	
	
	if (mgl_ShapeMode == MGL_TRIANGLES && vertexList.size() == 3)
	{
		shapeList.push_back(vertexList);
		vertexList.clear();
	}
	
	if (mgl_ShapeMode == MGL_QUADS && vertexList.size() == 4)
	{
		shapeList.push_back(vertexList);
		vertexList.clear();
	}
	
	Vertex Vec2(x, y, 0, 1);
	
	// Apply the transformations before you push to vertexList
	Vec2.applyTransformations();
	vertexList.push_back(Vec2);

}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
	if (mgl_ShapeMode == -1)
	{
		MGL_ERROR("Missing mglBegin! Aborting mission.\n");
		exit(1);
	}
	
	if (mgl_ShapeMode == MGL_TRIANGLES && vertexList.size() == 3)
	{
		shapeList.push_back(vertexList);
		vertexList.clear();
	}
	
	if (mgl_ShapeMode == MGL_QUADS && vertexList.size() == 4)
	{
		shapeList.push_back(vertexList);
		vertexList.clear();
	}
	
	Vertex Vec3(x, y, z, 1);
	
	// Apply the transformations before you push to vertexList
	Vec3.applyTransformations();
	vertexList.push_back(Vec3);
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
 */
void mglPushMatrix()
{
	// Depending on which mode you're on, you want to push the matrix onto
	// a specific stack. So you don't accidentally modify the projection matrix
	// when you wanted to modify your model matrix.
	
	// Push a copy of the top matrix
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
	// Sets the matrix currently on the top of the stack to the identity matrix
	// This function seems to appear after mglMatrixMode calls.
	
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
	// Push a copy of the matrix m onto the stack
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
	Matrix t;
	t.createTranslater(x, y, z);
	
	if (mgl_MatrixMode == MGL_PROJECTION)
		mglMultMatrix(ProjMatrixStack.top(), t); 
		
	else if (mgl_MatrixMode == MGL_MODELVIEW)
		mglMultMatrix(ModelMatrixStack.top(), t);	
}

// This function is used for mglRotation.
void normalize(MGLfloat &x, MGLfloat &y, MGLfloat &z)
{
	MGLfloat mag = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2) );
		
	if (mag == 0)
		return;
		
	x = x/mag;
	y = y/mag;
	z = z/mag;	
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
	
	normalize(x, y, z);
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
 // This handles the 3D stuff for a camera so objects in a scene
 // apear to be tilted. 
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

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLbyte red,
              MGLbyte green,
              MGLbyte blue)
{
	RGB[0] = red;
	RGB[1] = green;
	RGB[2] = blue;
}
