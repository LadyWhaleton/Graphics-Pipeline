/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 * Do not use any additional files
 */

#include <cstdio>
#include <cstdlib>
#include "minigl.h"

using namespace std;


/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

int mgl_ShapeMode = -1;
int mgl_MatrixMode = -1;

// from lab 1
void readFromFile(vector<float> &monkeyVector)
{
    ifstream file;
    file.open("monkey.raw");
    float data;
    
    while (!file.eof())
    {
        file >> data;
        monkeyVector.push_back(data/1.5);
    }
    // some reason it adds the last vertex twice
    monkeyVector.pop_back();
    
    file.close();
}

// ==================from lab2, the DDA algorithm=======================
void draw_line_shallow_pos(int x0, int x1, int y0, int y1, float m)
{
    float y = y0;
    for(int x = x0; x < x1; ++x)
    {
        // compute successive y values
        // cout << "shallow slope pos: " << m << endl;
        set_pixel(x, y);
        y = y + m;
    }
}

void draw_line_shallow_neg(int x0, int x1, int y0, int y1, float m)
{
    float y = y0;
    for(int x = x0; x > x1; --x)
    {
        // compute successive y values
        // cout << "shallow slope neg: " << m << endl;
        set_pixel(x, y);
        y = y - m;
    }
}

void draw_line_steep_pos(int x0, int x1, int y0, int y1, float m)
{
    float x = x0; 
    // sample at dy = 1
    for(int y = y0; y < y1; ++y)
    {
        // << "steep sleep pos: " << m << endl;
        // compute success x values
        set_pixel(x, y);
        x = x + 1/m;            
    }
}

void draw_line_steep_neg(int x0, int x1, int y0, int y1, float m)
{
    float x = x0; 
    // sample at dy = 1
    for(int y = y0; y > y1; --y)
    {
        //cout << "steep sleep neg: " << m << endl;
        // compute success x values
        set_pixel(x, y);
        x = x - 1/m;            
    }
}

void draw_line(int x0, int y0, int x1, int y1)
{
    //NOT WORKING CODE(PUT BETTER CODE HERE!!)
    float dx = x1 - x0;
    float dy = y1 - y0;
    
    if (dx == 0)
    {
        cout << "Error: Undefined slope!\n";
        return;
    }
    
    // calculate slope
    float m = dy/dx;

    // shallow slope
    if (abs(m) <= 1)
    {
        if (dx < 0)
            draw_line_shallow_neg(x0, x1, y0, y1, m);
        else
            draw_line_shallow_pos(x0, x1, y0, y1, m);
    }
    
    // steep slope
    else if (abs(m) > 1)
    {
        if (dy < 0)
            draw_line_steep_neg(x0, x1, y0, y1, m);
        else
            draw_line_steep_pos(x0, x1, y0, y1, m);
    }
        
    return;
}
// ===================== end of DDA stuff ============================

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
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
	// check what the mode is
	if (mode == MGL_TRIANGLES || mode == MGL_QUADS)
		mgl_shapeMode = mode;
	
	else
		mgl_ShapeMode = -1;
}

/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
	mgl_ShapeMode = -1;
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
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
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
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
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
void mglLoadMatrix(const MGLfloat *matrix)
{
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
void mglMultMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLbyte red,
              MGLbyte green,
              MGLbyte blue)
{
}
