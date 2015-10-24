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

using namespace std;

// classes 
class Matrix
{
	vector <MGLfloat*> matrix;
	
	// by default, it's the identity matrix
	Matrix()
	{
		setRows(0,0,0);
	}
	
	Matrix(MGLfloat x, MGLfloat y, MGLfloat z)
	{
		setRows(x,y,z);
	}
	
	void setRows(MGLfloat x, MGLfloat y, MGLfloat z)
	{
		MGLfloat Row1[4] = {1, 0, 0, x};
		MGLfloat Row2[4] = {0, 1, 0, y};
		MGLfloat Row3[4] = {0, 0, 1, z};
		MGLfloat Row4[4] = {0, 0, 0, 1};
		
		matrix.push_back(Row1);
		matrix.push_back(Row2);
		matrix.push_back(Row3);
		matrix.push_back(Row4);
	}
	
	
};


// Global variables I added

int mgl_ShapeMode = -1;
int mgl_MatrixMode = -1;

MGLpixel color;
stack <Matrix> MatrixStack;



/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

/*
// from lab 1

// set pixel (x,y) in framebuffer to color col, where
// col is a float array of three values between 0 and 1
// which specify the amount of red, green, and blue to mix (e.g.
// RED: (1,0,0) GREEN: (0,1,0) BLUE: (0,0,1) 
// YELLOW: (1,1,0) MAGENTA: (1,0,1) CYAN: (0,1,1)
// )
//
void set_pixel(int x, int y, float col[3])
{
    // write a 1x1 block of pixels of color col to framebuffer
    // coordinates (x, y)
    glRasterPos2i(x, y);
    glDrawPixels(1, 1, GL_RGB, GL_FLOAT, col);
}

void set_pixel(int x, int y)
{
    float col[] = { 1.0, 1.0, 1.0 };
    set_pixel(x,y,col);
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
*/

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
	// data contains buffer of coordinates for each pixel on screen
	// data is an 1D array containing x,y coordinates of the screen
	// for ex: starting from bottom left. 
	// data[0] contains the point (0,0),
	// data[1] contains the point (1,0),
	// data[2] contains the point (2,0), etc.
	
	
	// MGLpixel color
	// set red,blue,green
	// data[y*width+x] = color
	MGL_SET_RED(color, 255);
	MGL_SET_BLUE(color, 255);
	MGL_SET_GREEN(color, 255);
	
	data[650] = color;
	
	// double for loop.
	// for y < height
	// for x < width

	// you times y by width to go up rows and add x to indicate which
	// row & column of the screen
	
	// this function is called in main. it colors all of the pixels
	// that you want it to. Need to store these pixels with the 
	// corresponding color somehow
	
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
	// check what the mode is
	if (mgl_ShapeMode == MGL_TRIANGLES || mgl_ShapeMode == MGL_QUADS)
		mgl_ShapeMode = mode;
	
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
