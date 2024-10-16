/**
@file main.cpp
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include "glm/glm.hpp"

#include "Image.h"

#define PI 3.141592653589793

using namespace std;

/**
 Class representing a single ray.
 */
class Ray{
public:
    glm::vec3 origin; ///< Origin of the ray
    glm::vec3 direction; ///< Direction of the ray
	/**
	 Contructor of the ray
	 @param origin Origin of the ray
	 @param direction Direction of the ray
	 */
    Ray(glm::vec3 origin, glm::vec3 direction) : origin(origin), direction(direction){
    }
};


class Object;

/**
 Structure representing the even of hitting an object
 */
struct Hit{
    bool hit; ///< Boolean indicating whether there was or there was no intersection with an object
    glm::vec3 normal; ///< Normal vector of the intersected object at the intersection point
    glm::vec3 intersection; ///< Point of Intersection
    float distance; ///< Distance from the origin of the ray to the intersection point
    Object *object; ///< A pointer to the intersected object
};

/**
 General class for the object
 */
class Object{
public:
	glm::vec3 color; ///< Color of the object
	/** A function computing an intersection, which returns the structure Hit */
    virtual Hit intersect(Ray ray) = 0;
};

/**
 Implementation of the class Object for sphere shape.
 */
class Sphere : public Object{
private:
    float radius; ///< Radius of the sphere
    glm::vec3 center; ///< Center of the sphere

public:
	/**
	 The constructor of the sphere
	 @param radius Radius of the sphere
	 @param center Center of the sphere
	 @param color Color of the sphere
	 */
    Sphere(float radius, glm::vec3 center, glm::vec3 color) : radius(radius), center(center){
		this->color = color;
    }
	/** Implementation of the intersection function*/
    Hit intersect(Ray ray){
        
        Hit hit;
		hit.hit = false;
		
        float r=this->radius;
        //calc. in Ray-sphere intersection.pdf method 2 (not seen class)
        float a=glm::dot(this->center, ray.direction);  //origin to intersection point
        float normpow2=glm::dot(this->center,this->center); //since <a,a>=||a^2|| , calculates the norm of a vec3
        float D=sqrt(normpow2 - glm::pow(a,2));  //distance intersection point to center of sphere
        /*cases: if D:
         * <r 2sol
         * =r 1sol
         * >r no sol
         * */
        if (D<=r){    //!(D>r)
            float sqrtdelta=sqrt(glm::pow(this->radius,2)-glm::pow(D,2));
            float t1=a+sqrtdelta;
            float t2=a-sqrtdelta;
            if (D==r){
                //1 collision
                //I can use t1 or t2 ...are the same
                hit.intersection = ray.origin+ (ray.direction* t1);
            } else{
                //2 collisions
                //calc. the closest intersection point
                hit.intersection = ray.origin+ (ray.direction* min(t1,t2));
            }
            hit.hit = true;
            hit.normal =glm::normalize(hit.intersection-this->center);  //only works for a sphere
            hit.distance = distance(ray.origin,hit.intersection);
            hit.object = this;
        }
        //no collision
		return hit;
    }
};


vector<Object *> objects; ///< A list of all objects in the scene

/**
 Functions that computes a color along the ray
 @param ray Ray that should be traced through the scene
 @return Color at the intersection point
 */
glm::vec3 trace_ray(Ray ray){
		
	Hit closest_hit;
	
	closest_hit.hit = false;
	closest_hit.distance = INFINITY;

	//Loop over all objects to find the closest intersection
	for(int k = 0; k<objects.size(); k++){
		Hit hit = objects[k]->intersect(ray);
		if(hit.hit == true && hit.distance < closest_hit.distance)
			closest_hit = hit;
	}
	
	glm::vec3 color(0.0);
		
	if(closest_hit.hit){
		color = closest_hit.object->color;
	}else{
		color = glm::vec3(0.0, 0.0, 0.0); //black
	}
	return color;
}
/**
 Function defining the scene
 */
void sceneDefinition (){
	
	//first sphere
	objects.push_back(new Sphere(1.0, glm::vec3(-0,-2,8), glm::vec3(0.6, 0.9, 0.6)));
	//2
	objects.push_back(new Sphere(1.0, glm::vec3(1.0, -2.0, 8.0), glm::vec3(0.6, 0.6, 0.9)));
}

int main(int argc, const char * argv[]) {
	
    clock_t t = clock(); // variable for keeping the time of the rendering
    
    int width = 1024; //width of the image
    int height = 768; // height of the image
    float fov = 90; // field of view
	
	sceneDefinition(); // Let's define a scene
	
	Image image(width,height); // Create an image where we will store the result
    
    float s = 2*tan(0.5*fov/180*M_PI)/width;
    float X = -s * width / 2; 
    float Y = s * height / 2;
    
    for(int i = 0; i < width ; i++)
        for(int j = 0; j < height ; j++){
            
			float dx = X + i*s + s/2;
            float dy = Y - j*s - s/2;
            float dz = 1;
			
			glm::vec3 origin(0, 0, 0);
            glm::vec3 direction(dx, dy, dz);
            direction = glm::normalize(direction);
            
            Ray ray(origin, direction);
			
			image.setPixel(i, j, trace_ray(ray));
        }
    
    t = clock() - t;
    cout<<"It took " << ((float)t)/CLOCKS_PER_SEC<< " seconds to render the image."<< endl;
    cout<<"I could render at "<< (float)CLOCKS_PER_SEC/((float)t) << " frames per second."<<endl;
    
	// Writing the final results of the rendering
	if (argc == 2){
		image.writeImage(argv[1]);
	}else{
		image.writeImage("./result.ppm");
	}
	
    return 0;
}
