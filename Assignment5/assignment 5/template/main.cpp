/**
@file main.cpp
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"

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

glm::mat4 genTRMat(glm::vec3 t, // Translation
                   glm::vec3 r, // Rotation
                   glm::vec3 s  // Scaling
) {
    glm::mat4 trans = glm::translate(glm::mat4(1.0f),t);
    glm::mat4 rotat = glm::rotate(trans,glm::radians(r.x), glm::vec3(1,0,0));
    rotat = glm::rotate(rotat,glm::radians(r.y), glm::vec3(0,1,0));
    rotat = glm::rotate(rotat,glm::radians(r.z), glm::vec3(0,0,1));
    glm::mat4 scale = glm::scale(rotat, s);
    return scale;
}

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
	
protected:
	glm::mat4 transformationMatrix; ///< Matrix representing the transformation from the local to the global coordinate system
	glm::mat4 inverseTransformationMatrix; ///< Matrix representing the transformation from the global to the local coordinate system
	glm::mat4 normalMatrix; ///< Matrix for transforming normal vectors from the local to the global coordinate system
	
public:
	glm::vec3 color; ///< Color of the object
	Material material; ///< Structure describing the material of the object
	bool trans = false;
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
	/** Functions for setting up all the transformation matrices
	 @param matrix The matrix representing the transformation of the object in the global coordinates */
	void setTransformation(glm::mat4 matrix){
		
		transformationMatrix = matrix;
		inverseTransformationMatrix = glm::inverse(transformationMatrix);
		normalMatrix = glm::transpose(inverseTransformationMatrix);
		
		trans = true;
		/* ----- Assignment 5 ---------
		 Set the two remaining matrices
		
		inverseTransformationMatrix =
		normalMatrix =
		 
		 */
	}

	Ray toLocal(Ray ray) {
        ray.origin = glm::vec3(inverseTransformationMatrix * glm::vec4(ray.origin, 1.0f ));
        ray.direction = glm::normalize(glm::vec3(inverseTransformationMatrix * glm::vec4(ray.direction, 0.0f)));
        return ray;
    }

	Hit toGlobal(Hit localHit, Ray ray) {
        if(!localHit.hit) { return localHit; }

		Hit globalHit = localHit;
        globalHit.normal = glm::normalize(glm::vec3(normalMatrix * glm::vec4(globalHit.normal, 0.0f)));
        globalHit.intersection = glm::vec3(transformationMatrix * glm::vec4(globalHit.intersection, 1.0f));

        globalHit.distance = glm::length(globalHit.intersection - ray.origin);

        return globalHit;
    }
};

/**
 Implementation of the class Object for sphere shape.
 */
class Sphere : public Object{
private:
    float radius = 1; ///< Radius of the sphere
    glm::vec3 center = glm::vec3(0.0f); ///< Center of the sphere

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
	Sphere(glm::vec3 color) {
		this->color = color;
    }
	Sphere(Material material) {
		this->material = material;
    }
	/** Implementation of the intersection function*/
    Hit intersect(Ray ray){

		Ray local = toLocal(ray);

        glm::vec3 c = center - local.origin;

        float cdotc = glm::dot(c,c);
        float cdotd = glm::dot(c, local.direction);

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

			hit.intersection = local.origin + t * local.direction;
			hit.normal = glm::normalize(hit.intersection - center);
			hit.distance = glm::distance(local.origin, hit.intersection);
			hit.object = this;
			
			hit.uv.s = (asin(hit.normal.y) + M_PI/2)/M_PI;
			hit.uv.t = (atan2(hit.normal.z,hit.normal.x) + M_PI) / (2*M_PI);
        }
		else{
            hit.hit = false;
		}
		return toGlobal(hit, ray);
    }
};


class Plane : public Object{

private:
	glm::vec3 normal = glm::vec3(0.0f, 1.0f, 0.0f);
	glm::vec3 point = glm::vec3(0.0f);

public:
	Plane(glm::vec3 point, glm::vec3 normal) : point(point), normal(normal){
	}
	Plane(glm::vec3 point, glm::vec3 normal, Material material) : point(point), normal(normal){
		this->material = material;
	}
	Plane(Material material, glm::mat4 transformation){
		this->material = material;
		setTransformation(transformation);
	}
	Hit intersect(Ray ray){
		
		Hit hit;
		hit.hit = false;

		Ray local = ray;

		if(trans){
			local = toLocal(ray);
		}

		float DdotN = glm::dot(local.direction, normal);
		if(DdotN < 0){
			
			float PdotN = glm::dot (point-local.origin, normal);
			float t = PdotN/DdotN;
			
			if(t > 0){
				hit.hit = true;
				hit.normal = normal;
				hit.distance = t;
				hit.object = this;
				hit.intersection = t * local.direction + local.origin;
			}
		}
		

		if (!trans){
			return hit;
		}
		else{
			return toGlobal(hit, ray);
		}
	}
};

class Cone : public Object{
private:
	const float height = 1.0; ///< Local height of the cone
    const glm::vec3 H = glm::vec3(0.0f, height, 0.0f); ///< Local height of the cone
    const glm::vec3 O = glm::vec3(0.0f); ///< origin of the cone
    const float radius = 1.0; ///< Local radius of the cone

public:
	Cone(Material material){
		this->material = material;
	}
	Hit intersect(Ray ray){
		
		Hit hit;
		hit.hit = false;

		Ray local = toLocal(ray);
		
		float tan = (radius / height) * (radius / height);

        float a = (local.direction.x * local.direction.x) +
                  (local.direction.z * local.direction.z) -
                  (tan*(local.direction.y * local.direction.y));
        float b = (2*local.origin.x*local.direction.x) +
                  (2*local.origin.z*local.direction.z) +
                  (2*tan*(height-local.origin.y)*local.direction.y);
        float c = (local.origin.x*local.origin.x) +
                  (local.origin.z*local.origin.z) -
                  (tan*((height-local.origin.y)*(height-local.origin.y)));

        float delta = b*b - 4*a*c;
        if (delta < 0) {
            return hit;
        }

        float t1 = (-b - sqrt(delta)) / (2*a);
        float t2 = (-b + sqrt(delta)) / (2*a);

        float closest_t;
        if (t1 < 0 && t2 < 0) {
            return hit;
        } else if (t1 < 0) {
            closest_t = t2;
        } else if (t2 < 0) {
            closest_t = t1;
        } else {
            closest_t = t1 < t2 ? t1 : t2;
        }
        hit.distance = closest_t;
        hit.intersection = local.origin + local.direction * hit.distance;

        // create a plain at the base of the cone
        
		glm::mat4 plane_t = genTRMat(glm::vec3(0),glm::vec3(180.0, 0, 0),glm::vec3(0.5f));
		Plane base = Plane(material, plane_t);


        Hit baseHit = base.intersect(local);

        if (baseHit.hit && baseHit.distance < closest_t || (hit.intersection.y < 0 || hit.intersection.y > height)) {
            float distance = glm::distance(baseHit.intersection, O);
            if (distance <= radius) {
                return toGlobal(baseHit, ray);
            }
        }

        if (hit.intersection.y < 0 || hit.intersection.y > height) {
            return hit;
        }
        hit.hit = true;
        hit.normal = glm::normalize(hit.intersection);
        hit.object = this;

        // Texture mapping
        hit.uv.x = 0.5 - asin(hit.normal.y) / M_PI;
        hit.uv.y = 2*(0.5 + atan2(hit.normal.z, hit.normal.x) / (2 * M_PI));

        return toGlobal(hit, ray);
	}
		/*  ---- Assignemnt 5 -----
		
		 Implement the ray-cone intersection. Before intersecting the ray with the cone,
		 make sure that you transform the ray into the local coordinate system.
		 Remember about normalizing all the directions after transformations.
		 
		*/
	
		/* If the intersection is found, you have to set all the critical fields in the Hit strucutre
		 Remember that the final information about intersection point, normal vector and distance have to be given
		 in the global coordinate system.
		 
		hit.hit = true;
		hit.object = this;
		hit.intersection =
		hit.normal =
		hit.distance =
		
		 */
		
		//return hit;
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
glm::vec3 ambient_light(0.001,0.001,0.001);
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

		
		glm::vec3 diffuse_color = material.diffuse;
		if(material.texture){
			diffuse_color = material.texture(uv);
		}
		
		glm::vec3 diffuse = diffuse_color * glm::vec3(NdotL);
		glm::vec3 specular = material.specular * glm::vec3(pow(VdotR, material.shininess));
		
		
		// distance to the light
		float r = glm::distance(point,lights[light_num]->position);
		r = max(r, 0.1f);
		

		color += lights[light_num]->color * (diffuse + specular) / r/r;
	}
	color += ambient_light * material.ambient;
	
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
void sceneDefinition (float x){

	Material green_diffuse;
	green_diffuse.ambient = glm::vec3(0.03f, 0.1f, 0.03f);
	green_diffuse.diffuse = glm::vec3(0.3f, 1.0f, 0.3f);

	Material red_specular;
	red_specular.diffuse = glm::vec3(1.0f, 0.2f, 0.2f);
	red_specular.ambient = glm::vec3(0.01f, 0.02f, 0.02f);
	red_specular.specular = glm::vec3(0.5);
	red_specular.shininess = 10.0;

	Material blue_specular;
	blue_specular.ambient = glm::vec3(0.02f, 0.02f, 0.1f);
	blue_specular.diffuse = glm::vec3(0.2f, 0.2f, 1.0f);
	blue_specular.specular = glm::vec3(0.6);
	blue_specular.shininess = 100.0;

	//Textured sphere
	Material textured;
	textured.texture = &rainbowTexture;
	Sphere *s3 = new Sphere(textured);

	if (x == -1){	
		s3->setTransformation(genTRMat(glm::vec3(-6,2.0,16),glm::vec3(0.0),glm::vec3(4.5f)));
	}else{
		objects.clear();
		float a = 2.0 + (2*x);
		s3->setTransformation(genTRMat(glm::vec3(-6,a,16),glm::vec3(0.0)+glm::vec3(0,180,0)*x,glm::vec3(4.5f)));
	}
	
	objects.push_back(s3);


	Sphere *s1 = new Sphere(blue_specular);
	s1->setTransformation(genTRMat(glm::vec3(1.0, -2.0, 9.0),glm::vec3(0.0),glm::vec3(1.0f)));
	objects.push_back(s1);

	Sphere *s2 = new Sphere(red_specular);
	s2->setTransformation(genTRMat(glm::vec3(-1.0, -2.5, 6.0),glm::vec3(0.0),glm::vec3(0.5f)));
	objects.push_back(s2);
	
	
	// ------ Assignment 5 -------
	
	// You can remove the green sphere as it should be replaced with a green cone
	//objects.push_back(new Sphere(1.0, glm::vec3(3,-2,6), green_diffuse));

	//Planes
	Material red_diffuse;
	red_diffuse.ambient = glm::vec3(0.09f, 0.06f, 0.06f);
	red_diffuse.diffuse = glm::vec3(0.9f, 0.6f, 0.6f);
		
	Material blue_diffuse;
	blue_diffuse.ambient = glm::vec3(0.06f, 0.06f, 0.09f);
	blue_diffuse.diffuse = glm::vec3(0.6f, 0.6f, 0.9f);
	objects.push_back(new Plane(glm::vec3(0,-3,0), glm::vec3(0.0,1,0)));
	objects.push_back(new Plane(glm::vec3(0,1,30), glm::vec3(0.0,0.0,-1.0), green_diffuse));
	objects.push_back(new Plane(glm::vec3(-15,1,0), glm::vec3(1.0,0.0,0.0), red_diffuse));
	objects.push_back(new Plane(glm::vec3(15,1,0), glm::vec3(-1.0,0.0,0.0), blue_diffuse));
	objects.push_back(new Plane(glm::vec3(0,27,0), glm::vec3(0.0,-1,0)));
	objects.push_back(new Plane(glm::vec3(0,1,-0.01), glm::vec3(0.0,0.0,1.0), green_diffuse));
	
	
	/* ----- Assignment 5 -------
	Create two conse and add them to the collection of our objects.
	Remember to create them with corresponding materials and transformation matrices
	
	
	Cone *cone1 = new Cone(...);
	cone1->setTransformation(...);
	objects.push_back(cone1);
	
	Cone *cone2 = new Cone(...);
	cone2->setTransformation(...);
	objects.push_back(cone2);
	
	*/

	Material yellow;
    yellow.ambient = glm::vec3(0.07f, 0.07f, 0.07f);
    yellow.diffuse = glm::vec3(1.0f, 1.0f, 0.0f);
    yellow.specular = glm::vec3(1.0);
    yellow.shininess = 100.0;
	
	Cone *cone1 = new Cone(yellow);
	glm::mat4 cone1_t = genTRMat(glm::vec3(5.0, -3.0, 14.0),glm::vec3(0.0),glm::vec3(3.0f, 12.0f, 3.0f));
	cone1->setTransformation(cone1_t);
	objects.push_back(cone1);
	
	Cone *cone2 = new Cone(green_diffuse);
	glm::mat4 cone2_t = genTRMat(glm::vec3(3.5, -2.0, 8.0),glm::vec3(0.0, 0.0, -110.0),glm::vec3(1.0f, 3.0f, 1.0f));
	cone2->setTransformation(cone2_t);
	objects.push_back(cone2);

	lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(1.0, 1.0, 1.0)));
	lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(0.1)));
	lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.4)));

}

/**
 Function performing tonemapping of the intensities computed using the raytracer
 @param intensity Input intensity
 @return Tonemapped intensity in range (0,1)
 */
glm::vec3 toneMapping(glm::vec3 intensity){
	float gamma = 1.0/2.0;
	float alpha = 12.0f;
	return glm::clamp(alpha * glm::pow(intensity, glm::vec3(gamma)), glm::vec3(0.0), glm::vec3(1.0));
}

int main(int argc, const char * argv[]) {

	if(argc < 2){

    clock_t t = clock(); // variable for keeping the time of the rendering

    int width = 1024; //width of the image
    int height = 768; // height of the image
    float fov = 90; // field of view

	sceneDefinition(-1); // Let's define a scene

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
}else{
		int total_frames = atoi(argv[1]); // Número total de fotogramas en la animación

		float animation = 0;

        for (int frame = 0; frame < total_frames/2; frame++) {
            clock_t t = clock(); // variable para medir el tiempo de renderización

            int width = 1024; // ancho de la imagen
            int height = 768; // alto de la imagen
            float fov = 90; // campo de visión

			
			animation = (frame/(float)total_frames)*2.0;

            sceneDefinition(animation); // Definir la escena para el fotograma actual


            Image image(width, height); // Crear una imagen donde se almacenará el resultado

            float s = 2 * tan(0.5 * fov / 180 * M_PI) / width;
            float X = -s * width / 2;
            float Y = s * height / 2;

            for (int i = 0; i < width; i++){
                for (int j = 0; j < height; j++) {

                    float dx = X + i * s + s / 2;
                    float dy = Y - j * s - s / 2;
                    float dz = 1;

                    glm::vec3 origin(0, 0, 0);
                    glm::vec3 direction(dx, dy, dz);
                    direction = glm::normalize(direction);

                    Ray ray(origin, direction);

					image.setPixel(i, j, toneMapping(trace_ray(ray)));

                }
			}
			t = clock() - t;
            cout << "Frame " << frame << ": Rendering took " << ((float)t) / CLOCKS_PER_SEC << " seconds." << endl;

            // Escribir el resultado de la renderización en un archivo con nombre único
            char filename[100];
            snprintf(filename, sizeof(filename), "result_%d.ppm", frame);
            image.writeImage(filename);
		}

		for (int frame = total_frames/2; frame >= 0 ; frame--) {
            clock_t t = clock(); // variable para medir el tiempo de renderización

            int width = 1024; // ancho de la imagen
            int height = 768; // alto de la imagen
            float fov = 90; // campo de visión

			
			animation = (frame/(float)total_frames)*2.0;

            sceneDefinition(animation); // Definir la escena para el fotograma actual

            Image image(width, height); // Crear una imagen donde se almacenará el resultado

            float s = 2 * tan(0.5 * fov / 180 * M_PI) / width;
            float X = -s * width / 2;
            float Y = s * height / 2;

            for (int i = 0; i < width; i++){
                for (int j = 0; j < height; j++) {
                    float dx = X + i * s + s / 2;
                    float dy = Y - j * s - s / 2;
                    float dz = 1;

                    glm::vec3 origin(0, 0, 0);
                    glm::vec3 direction(dx, dy, dz);
                    direction = glm::normalize(direction);

                    Ray ray(origin, direction);

					image.setPixel(i, j, toneMapping(trace_ray(ray)));

                }
		}
            t = clock() - t;
            cout << "Frame " << frame << ": Rendering took " << ((float)t) / CLOCKS_PER_SEC << " seconds." << endl;

            // Escribir el resultado de la renderización en un archivo con nombre único
            char filename[100];
            snprintf(filename, sizeof(filename), "result_%d_0.ppm", frame);
            image.writeImage(filename);
		}
}

    return 0;
}
