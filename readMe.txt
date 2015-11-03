Name: Stephanie Tong
Email: stong002@ucr.edu
Course: CS130 - Computer Graphics
Project 1: Simplified Rendering Pipeline

Functions:

1. mglReadPixels (Line 557)
	Iterates through a vector of shapes. Each shape
	contains a certain number of vertices (depending whether it's a
	triangle or quadrilateral). I convert each shape's vertices to 
	screen coordinates. After that, I check whether the shape is a 
	triangle or quad based on the number of vertices and all the 
	corresponding rasterization function.
	
	My rasterization functions utilizes bounding boxes and barycentric 
	coordinates to determine which pixels to be added to the 
	z buffer. After that process is complete, I sort
	the z buffer by descending order. My z buffer is just a vector of
	pixels, each of which contain x, y, z and a color. After the sort,
	I iterate through all pixels in the z buffer and set
	data[y*height + x] to the current z buffer pixel.
	
2. mglBegin (Line 635)
	Checks the parameter mode for a valid value.
	If the mode is either MGL_TRIANGLES or MGL_QUADS, I set my global
	variable called mgl_ShapeMode to mode. Otherwise, mgl_ShapeMode = -1.
	
3. mglEnd (Line 648)
	Checks vertexList if it's empty or not. If it isn't, push the vertexList 
	onto a shape. This means that the shape now is specified by those vertices 
	in the vertexList. For my implementation, vertexList is a vector of class Vertex. 
	I do not have a class named shape, but I use a 
	vector < vector<Vertex> > shapeList which essentially contains all 
	of the vertices for each shape.

4. mglVertex2 (Line 666)
	First, check if mgl_shapeMode is -1. If it is, print out an error
	message and exit the program. After that, check what mgl_shapeMode is.
	If mgl_shapeMode == MGL_TRIANGLES and the vertexList is currently size 3,
	then push the vertexList onto shapeList. That way, I can
	start specifying vertices for another triangle. I also do a similiar 
	check if mgl_shapeMode == MGL_QUADS and vertexList.size() == 4.
	
	After that check, I simply create a Vertex with x & y, and I set z to 0
	and w to 1. Within my Vertex constructor, I store the current RGB color 
	values into an size three array of MGLpixels. After I create that vertex,
	I apply all projection and modelview projections onto the vertex. 
	(I don't scale to screen size yet.) Then, I push that vertex onto
	vertexList.
	
5. mglVertex3 (Line 700)
	For mglVertex3, I perform the same checks and operations as mglVertex2 
	except that z has an actual value now.
	
6. mglMatrixMode (Line 732)
	Check the parameter mode for a valid value.
	If the mode is either MGL_MODELVIEW or MGL_PROJECTION, set the global
	variable called mgl_MatrixMode to mode. Otherwise, mgl_MatrixMode = -1.
	
7. mglPushMatrix (Line 744)
	Depending on which mode I'm on, I push a copy of the current matrix 
	specified by mgl_MatrixMode onto a specific stack. For example,
	if mgl_MatrixMode == MGL_MODELVIEW, I push a copy of the current matrix
	onto ModelMatrixStack. For my implementation, the current matrix is
	the matrix is the one currently on top of the stack specified by the
	matrix mode.
	
8. mglPopMatrix (Line 762)
	This function works similiarly to mglPushMatrix except I'm
	popping off the top of the stack.
	
9. mglLoadIdentity (Line 780)
	Sets the current matrix to the identity matrix. This actually replaces
	the matrix, not create a copy. If the stack specified by mgl_MatrixMode
	is empty, then just push the identity matrix.
	
10. mglLoadMatrix (Line 822)
	Similiar to mglLoadIdentity, but replace with the parameter matrix
	instead of the identity. Likewise, if the stack is empty, just
	push the matrix.

11. mglMultMatrix (Line 844)
	I modified this function so that it takes two parameters instead of one.
	It multiplies the two matricies from the parameters together, and modifies
	the left matrix parameter. See line 844 for the function algorithm.
	
12. mglTranslate (Line 862)
	Creates a translation matrix. Then calls 
	mglMultMatrix (current Matrix, translation matrix) based on the current
	matrix mode. 
	
13. mglRotate (Line 896)
	First, normalizes the parameters x, y, z by their magnitude. Then
	creates the rotation matrix after various calculations. Then calls 
	mglMultMatrix (current Matrix, rotation matrix) based on the current
	matrix mode. 
	
14. mglScale (Line 930)
	Creates a scaling matrix. Then calls 
	mglMultMatrix (current Matrix, scaling matrix) based on the current
	matrix mode. 
	
15. mglFrustrum (Line 953)
	Creates a matrix representing the frustrum. Then calls 
	mglMultMatrix (current Matrix, frustrum matrix) based on the current
	matrix mode. See function for calculations.
	
16. mglOrtho (Line 994)
	Creates the orthographic matrix. Then calls 
	mglMultMatrix (current Matrix, orthographic matrix) based on the current
	matrix mode. See function for calculations.
	
17. mglColor (Line 1029)
	Sets RGB[0] to the parameter red, RGB[1] to green, and RGB[2] to blue.
	
