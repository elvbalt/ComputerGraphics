//------------ Assignment ----------------
// Add vertices for the cube as well as colors for corresponding vertices.
// Remember that:
// (1) The vertices for each triangle has to be provided in a counter-clockwise order to make the face culling work properly. If the order of vertices won't be correct the triangles won't be visible.
// (2) The current coordinate system has: x-axis pointing right, y-axis pointing up, and z-axis pointing towards the screen
// (3) Everything which is within [-1..1]x[-1..1]x[-1..1] volume will be projected onto the screen along z-axis
//-----------------------------------------

//var cube_vertices = [...];

var cube_vertices = [
    /**front */
    0.5,0.5,-0.5,
    -0.5,0.5,-0.5,
    -0.5,-0.5,-0.5,
    0.5,0.5,-0.5,
    -0.5,-0.5,-0.5,
    0.5,-0.5,-0.5,
    /**back*/
    0.5,-0.5,0.5,
    -0.5,-0.5,0.5,
    0.5,0.5,0.5,
    -0.5,-0.5,0.5,
    -0.5,0.5,0.5,
    0.5,0.5,0.5,
    /**left*/
    -0.5,0.5,-0.5,
    -0.5,0.5,0.5,
    -0.5,-0.5,0.5,
    -0.5,0.5,-0.5,
    -0.5,-0.5,0.5,
    -0.5,-0.5,-0.5,
    /**right */
    0.5,-0.5,-0.5,
    0.5,-0.5,0.5,
    0.5,0.5,-0.5,
    0.5,-0.5,0.5,
    0.5,0.5,0.5,
    0.5,0.5,-0.5,
    /**up */
    0.5,0.5,0.5,
    -0.5,0.5,0.5,
    -0.5,0.5,-0.5,
    0.5,0.5,0.5,
    -0.5,0.5,-0.5,
    0.5,0.5,-0.5,
    /**down */
    0.5,-0.5,-0.5,
    -0.5,-0.5,-0.5,
    0.5,-0.5,0.5,
    -0.5,-0.5,-0.5,
    -0.5,-0.5,0.5,
    0.5,-0.5,0.5
   
];

//var cube_colors = [...];

var cube_colors = [
    /*CYAN front face*/
    0.0, 1.0, 1.0, 
    0.0, 1.0, 1.0, 
    0.0, 1.0, 1.0, 
    0.0, 1.0, 1.0, 
    0.0, 1.0, 1.0, 
    0.0, 1.0, 1.0,
    /*blue back face*/
    0.0, 0.0, 1.0,
    0.0, 0.0, 1.0,
    0.0, 0.0, 1.0,
    0.0, 0.0, 1.0,
    0.0, 0.0, 1.0,
    0.0, 0.0, 1.0,
    /*red left*/
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    /*orange right*/
    1.0, 0.5, 0.0,
    1.0, 0.5, 0.0,
    1.0, 0.5, 0.0,
    1.0, 0.5, 0.0,
    1.0, 0.5, 0.0,
    1.0, 0.5, 0.0,
    /*green up */
    0.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    /* */
    0.0, 1.0, 0.5,
    0.0, 1.0, 0.5,
    0.0, 1.0, 0.5,
    0.0, 1.0, 0.5,
    0.0, 1.0, 0.5,
    0.0, 1.0, 0.5,
];

