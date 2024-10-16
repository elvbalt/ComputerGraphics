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
#include "Material.h"

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
	glm::vec2 uv; ///< Coordinates for computing the texture (texture coordinates)
};

/**
 General class for the object
 */
class Object{
public:
	glm::vec3 color; ///< Color of the object
	Material material; ///< Structure describing the material of the object
	/** A function computing an intersection, which returns the structure Hit */
    virtual Hit intersect(Ray ray) = 0;

	/** Function that returns the material struct of the object*/
	Material getMaterial(){
		return material;
	}
	/** Function that set the material
	 @param material A structure describing the material of the object
	*/
	void setMaterial(Material material){
		this->material = material;
	}
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
	Sphere(float radius, glm::vec3 center, Material material) : radius(radius), center(center){
		this->material = material;
	}
	/** Implementation of the intersection function*/
    Hit intersect(Ray ray){

        glm::vec3 c = center - ray.origin;

        float cdotc = glm::dot(c,c);
        float cdotd = glm::dot(c, ray.direction);

        Hit hit;

        float D = 0;
		if (cdotc > cdotd*cdotd){
			D =  sqrt(cdotc - cdotd*cdotd);
		}
        if(D<=radius){
            hit.hit = true;
            float t1 = cdotd - sqrt(radius*radius - D*D);
            float t2 = cdotd + sqrt(radius*radius - D*D);

            float t = t1;
            if(t<0) t = t2;
            if(t<0){
                hit.hit = false;
                return hit;
            }

			hit.intersection = ray.origin + t * ray.direction;
			hit.normal = glm::normalize(hit.intersection - center);
			hit.distance = glm::distance(ray.origin, hit.intersection);
			hit.object = this;
			
			
			/*
			 
			 
			 Excercise 2 - computing texture coordintes for the sphere. You can refer to them as s and t, or x and y, whatever you prefer ;)
			 
			 hit.uv.s =
			 hit.uv.t =
			 
			 
			 
			 */

			glm::vec3 d = hit.normal;
			
			//en los apuntes

			float u = (atan2(d.z, d.x) + 3.14159265f)/(2.0f * 3.14159265f);
			float v = (asin(d.y) + 3.14159265f/2)/3.14159265f;

			hit.uv = glm::vec2(u, v);
			

        }
		else{
            hit.hit = false;
		}
		return hit;
    }
};


class Plane : public Object{

private:
	glm::vec3 normal;
	glm::vec3 point;

public:
	Plane(glm::vec3 point, glm::vec3 normal) : point(point), normal(normal){
	}
	Plane(glm::vec3 point, glm::vec3 normal, Material material) : point(point), normal(normal){
		this->material = material;
	}
	Hit intersect(Ray ray){
		
		Hit hit;
		hit.hit = false;
		/*
		 
		 
		 
		 Excercise 1 - Plane-ray intersection
		 
		 
		 
		 
		 */

		glm::vec3 a = this->point - ray.origin;
		float arg1 = glm::dot(a, normal);
		float arg2 = glm::dot(ray.direction, normal);

		if (arg2 != 0){
			
			float t = arg1/arg2;

			if (t < 0){
				return hit;
			}

			hit.hit = true;
			hit.intersection = ray.origin + (ray.direction*t);

			hit.normal = -normal;
			hit.distance = glm::distance(ray.origin, hit.intersection);
			hit.object = this;
		}
		
		return hit;
	}
};


/**
 Light class
 */
class Light{
public:
	glm::vec3 position; ///< Position of the light source
	glm::vec3 color; ///< Color/intentisty of the light source
	Light(glm::vec3 position): position(position){
		color = glm::vec3(1.0);
	}
	Light(glm::vec3 position, glm::vec3 color): position(position), color(color){
	}
};

vector<Light *> lights; ///< A list of lights in the scene
glm::vec3 ambient_light(0.5,0.5,0.5);
vector<Object *> objects; ///< A list of all objects in the scene


/** Function for computing color of an object according to the Phong Model
 @param point A point belonging to the object for which the color is computer
 @param normal A normal vector the the point
 @param uv Texture coordinates
 @param view_direction A normalized direction from the point to the viewer/camera
 @param material A material structure representing the material of the object
*/
glm::vec3 PhongModel(glm::vec3 point, glm::vec3 normal, glm::vec2 uv, glm::vec3 view_direction, Material material){

	glm::vec3 color(0.0);
	for(int light_num = 0; light_num < lights.size(); light_num++){

		glm::vec3 light_direction = glm::normalize(lights[light_num]->position - point);
		glm::vec3 reflected_direction = glm::reflect(-light_direction, normal);

		float NdotL = glm::clamp(glm::dot(normal, light_direction), 0.0f, 1.0f);
		float VdotR = glm::clamp(glm::dot(view_direction, reflected_direction), 0.0f, 1.0f);

		
		/*
		 
		 
		 Excercise 2 - Modify the code by adding texturing, i.e., 
		 the diffuse color should be computed using one 
		 of the texture functions according to the texture
		  coordinates stored in the uv variable. 
		  Make sure that the code works also for objects 
		  that should not have texture.
		 
		 
		 */
		

		glm::vec3 diffuse = material.diffuse;
        if (material.texture != nullptr) {
            // Si hay una textura definida en el material, utiliza las coordenadas uv para muestrearla
			
            diffuse = material.texture(uv);;
        }
        
        diffuse *= glm::vec3(NdotL);

		glm::vec3 specular = material.specular * glm::vec3(pow(VdotR, material.shininess));
		
		float distance = glm::distance(lights[light_num]->position, point);
        float light_attenuation = 1 / distance;
        
		/*
		 
		 
		 Excercise 3 - Modify the code by adding 
		 attenuation of the light due to distance 
		 from the intersection point to the light 
		 source
		 
		 
		 
		 */
		

		color += lights[light_num]->color * (diffuse + specular) * light_attenuation;
	}
	color += ambient_light * material.ambient;

	color = glm::clamp(color, glm::vec3(0.0), glm::vec3(1.0));
	return color;
}

/**
 Functions that computes a color along the ray
 @param ray Ray that should be traced through the scene
 @return Color at the intersection point
 */
glm::vec3 trace_ray(Ray ray){

	Hit closest_hit;

	closest_hit.hit = false;
	closest_hit.distance = INFINITY;

	for(int k = 0; k<objects.size(); k++){
		Hit hit = objects[k]->intersect(ray);
		if(hit.hit == true && hit.distance < closest_hit.distance)
			closest_hit = hit;
	}

	glm::vec3 color(0.0);

	if(closest_hit.hit){
		color = PhongModel(closest_hit.intersection, closest_hit.normal, closest_hit.uv, glm::normalize(-ray.direction), closest_hit.object->getMaterial());
	}else{
		color = glm::vec3(0.0, 0.0, 0.0);
	}
	return color;
}
/**
 Function defining the scene
 */
void sceneDefinition (){

	Material green_diffuse;
	green_diffuse.ambient = glm::vec3(0.07f, 0.09f, 0.07f);
	green_diffuse.diffuse = glm::vec3(0.7f, 0.9f, 0.7f);

	Material red_specular;
	
	red_specular.ambient = glm::vec3(0.01f, 0.03f, 0.03f);
	red_specular.diffuse = glm::vec3(1.0f, 0.3f, 0.3f);
	red_specular.specular = glm::vec3(0.5);
	red_specular.shininess = 10.0;

	Material blue_specular;
	blue_specular.ambient = glm::vec3(0.07f, 0.07f, 0.1f);
	blue_specular.diffuse = glm::vec3(0.7f, 0.7f, 1.0f);
	blue_specular.specular = glm::vec3(0.6);
	blue_specular.shininess = 100.0;


	objects.push_back(new Sphere(1.0, glm::vec3(1,-2,8), blue_specular));
	objects.push_back(new Sphere(0.5, glm::vec3(-1,-2.5,6), red_specular));
	objects.push_back(new Sphere(1.0, glm::vec3(3,-2,6), green_diffuse));
	
	
	// Excercise 2 - Textured sphere
	Material textured;
	//textured.texture = &checkerboardTexture;
	textured.texture = &rainbowTexture;
	objects.push_back(new Sphere(7.0, glm::vec3(-6,4,23), textured));
	
	Material textured2;
	textured2.texture = &checkerboardTexture;
	objects.push_back(new Sphere(2.0, glm::vec3(4,1,8), textured2));
	
	
	
	/*
	 
	 
	 Excercise 1 - Definition of planes and the materials
	 
	 
	 
	 */
	Material white;
	white.ambient = glm::vec3(0.06f, 0.06f, 0.06f);
	white.diffuse = glm::vec3(0.8f, 0.8f, 0.8f);
	white.specular = glm::vec3(0.1);
	white.shininess = 0.0;

	Material blue;
	blue.ambient = glm::vec3(0.001f, 0.001f, 0.004f);
	blue.diffuse = glm::vec3(1.0f, 1.0f, 1.3f);
	blue.shininess = 0.0;
	blue.specular = glm::vec3(0.0);

	Material red;
	red.ambient = glm::vec3(0.004f, 0.001f, 0.001f);
	red.diffuse = glm::vec3(1.5f, 0.8f, 0.8f);
	red.shininess = 0.0;
	red.specular = glm::vec3(0.0);

	Material green;
	green.ambient = glm::vec3(0.06f, 0.09f, 0.06f);
	green.diffuse=  glm::vec3(0.7f, 0.9f, 0.7f);
	green.shininess = 0.0;
	green.specular = glm::vec3(0.0);

	objects.push_back(new Plane(glm::vec3(0,0,30.0), glm::vec3(0,0,1.0), green));
	objects.push_back(new Plane(glm::vec3(0,0,-0.01), glm::vec3(0,0,-1), green));
	
	objects.push_back(new Plane(glm::vec3(15.0,0,0), glm::vec3(1.0,0,0), blue));
	objects.push_back(new Plane(glm::vec3(-15.0,0,0), glm::vec3(-1.0,0,0), red));
	
	objects.push_back(new Plane(glm::vec3(0,27.0,0), glm::vec3(0,1.0,0), white));
	objects.push_back(new Plane(glm::vec3(0,-3.0,0), glm::vec3(0,-1.0,0),white));
	

	
	lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(4)));
	lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(4)));
	lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(4)));
}

/**
 Function performing tonemapping of the intensities computed using the raytracer
 @param intensity Input intensity
 @return Tonemapped intensity in range (0,1)
 */
glm::vec3 toneMapping(glm::vec3 intensity){
	/*
	 
	 
	 Excercise 3 - Tone mapping
	 
	 
	 
	 */
	
	

	glm::vec3 alpha(2.0);
	glm::vec3 beta(2.0);
	glm::vec3 gamma(2.2);

	glm::vec3 tone_mapped = glm::pow(alpha * glm::pow(intensity, beta), glm::vec3(1.0) / gamma);
	
	return glm::clamp(tone_mapped, glm::vec3(0.0), glm::vec3(1.0));
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

			image.setPixel(i, j, toneMapping(trace_ray(ray)));
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
