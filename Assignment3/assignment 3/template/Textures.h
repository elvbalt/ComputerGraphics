//
//  Textures.h
//  Raytracer
//
//  Created by Piotr Didyk on 19.08.21.
//

#ifndef Textures_h
#define Textures_h


#include "glm/glm.hpp"

glm::vec3 checkerboardTexture(glm::vec2 uv){

    /*
     
     
        Exercise 2 (3 points)
     
     
    */
    int u = static_cast<int>(floor(uv.x * 50));
    int v = static_cast<int>(floor(uv.y * 25));

    
    if ((u + v) % 2 == 0) {
        return glm::vec3(1.0f, 1.0f, 1.0f);  // white
    } else {
        return glm::vec3(0.0f, 0.0f, 0.0f);  // black
    }
    return glm::vec3(0.0);
}
glm::vec3 rainbowTexture(glm::vec2 uv){
    /*
     
     
        Exercise 2 (5 points)
     
     
    */
    float diagonal = 2*uv.x + uv.y;

        // Determina si la coordenada está en una franja roja, verde o azul
        int stripe_index = int(diagonal * 25) % 3;

        glm::vec3 color;
        if (stripe_index == 0) {
            color = glm::vec3(1.0, 0.0, 0.0);  // Red
        } else if (stripe_index == 1) {
            color = glm::vec3(0.0, 1.0, 0.0);  // green
        } else {
            color = glm::vec3(0.0, 0.0, 1.0);  // blue
        }

        return color;
}

#endif /* Textures_h */
