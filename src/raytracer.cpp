#include "raytracer.h"
#include "../extern/linalg/linalg.h"  // Inclusion de la librairie linalg.h
#include <cmath>  // Pour tan()
#include <cstdlib> // Pour rand() et srand()
#include "basic.h"



using namespace linalg::aliases; // Simplifie les types (double3, etc.)


void Raytracer::render(const Scene& scene, Frame* output)
{       
    // Crée le z_buffer.
    double *z_buffer = new double[scene.resolution[0] * scene.resolution[1]];
    for(int i = 0; i < scene.resolution[0] * scene.resolution[1]; i++) {
        z_buffer[i] = scene.camera.z_far; //Anciennement DBL_MAX. À remplacer avec la valeur de scene.camera.z_far
    }


	//---------------------------------------------------------------------------------------------------------------
	// Nous vous fournissons ci-bas du code pour une caméra orthographique très simple. Cette caméra peut être utilisée pour tester l’intersection avec la sphère.
	// Vous devez utiliser la scène de test portho.ray pour utiliser cette caméra. 
	// Notez que votre code de caméra ne doit pas être basé sur ce code de caméra. Ce code n’est là que pour prendre en compte le développement initial du test d’intersection.
	// Pour utiliser cette caméra, vous devez supprimer les commentaires qui rendent inactive cette partie du code, et mettre en commentaires la boucle d’image originale.

	// CameraOrthographic camOrth;
	// double3 uVec{ 0,1,0 };
	// double3 vVec{ 0,0,1 };
	// double y_shift = 2.0 / scene.resolution[1];
	// double x_shift = 2.0 / scene.resolution[0];
	//
	// for (int y = 0; y < scene.resolution[1]; y++) {
	// 	if (y % 40) {
	// 		std::cout << "\rScanlines completed: " << y << "/" << scene.resolution[1] << '\r';
	// 	}
	//
	// 	for (int x = 0; x < scene.resolution[0]; x++) {
	// 		double3 color{ 0,0,0 };
	//
	// 		Intersection hit;
	// 		double3 rayOrigin = camOrth.minPosition + uVec * x_shift * x + vVec * y_shift * y;
	// 		double3 rayDirection{ 1,0,0 };
	// 		Ray ray = Ray(rayOrigin, rayDirection);
	// 		double itHits = 0;
	//
	// 		double z_depth = scene.camera.z_far;
	// 		if (scene.container->intersect(ray, EPSILON, z_depth, &hit)) {
	// 			Material& material = ResourceManager::Instance()->materials[hit.key_material];
	// 			color = material.color_albedo;
	// 			itHits = 1.0f;
	// 		}
	//
	// 		output->set_color_pixel(x, y, color);
	// 		output->set_depth_pixel(x, y, itHits);
	// 	}
	// }

	//---------------------------------------------------------------------------------------------------------------


	// @@@@@@ VOTRE CODE ICI
	// Calculez les paramètres de la caméra pour les rayons.

	double aspect_ratio = scene.camera.aspect;  // Utilisation du ratio d'aspect fourni
	double fov_rad = deg2rad(scene.camera.fovy); // Conversion du champ de vision (fovy) en radians
	double scale = tan(fov_rad / 2.0);  // Calcul du facteur de mise à l'échelle

	// Vecteurs de base pour la caméra
	double3 forward = normalize(scene.camera.center - scene.camera.position);  // Vecteur avant
	double3 right = normalize(cross(forward, scene.camera.up));  // Vecteur droite
	// double3 up = scene.camera.up;  // Vecteur haut déjà fourni
	double3 up= normalize(cross(right, forward));



	// debut commente pour tester portho.ray
    // Itère sur tous les pixels de l'image.
    for(int y = 0; y < scene.resolution[1]; y++) {
		if (y % 40){
			std::cout << "\rScanlines completed: " << y << "/" << scene.resolution[1] << '\r';
		}

        for(int x = 0; x < scene.resolution[0]; x++) {

			int avg_z_depth = 0;
			double3 avg_ray_color{0,0,0};
        	// Itérations pour chaque échantillon de pixel
			for(int iray = 0; iray < scene.samples_per_pixel; iray++) {
				// Génère le rayon approprié pour ce pixel.
				Ray ray;
				// Initialise la profondeur de récursivité du rayon.
				int ray_depth = 0;
				// Initialize la couleur du rayon
				double3 ray_color{0,0,0};

				// @@@@@@ VOTRE CODE ICI
				// Mettez en place le rayon primaire en utilisant les paramètres de la caméra.
				// Lancez le rayon de manière uniformément aléatoire à l'intérieur du pixel dans la zone délimité par jitter_radius.
				//Faites la moyenne des différentes couleurs obtenues suite à la récursion.

				// Ajout de jitter pour l'anti-aliasing
				double jitter_x = (rand_double() - 0.5) * scene.jitter_radius;
				double jitter_y = (rand_double() - 0.5) * scene.jitter_radius;

				// Calcul des coordonnées dans l'écran
				double pixel_ndc_x = (x + jitter_x + 0.5) / scene.resolution[0];
				double pixel_ndc_y = (y + jitter_y + 0.5) / scene.resolution[1];

				double pixel_screen_x = (2 * pixel_ndc_x - 1) * aspect_ratio * scale;
				// double pixel_screen_y = (1 - 2 * pixel_ndc_y) * scale;
				double pixel_screen_y = (2 * pixel_ndc_y - 1) * scale;

				// Calcul de la direction du rayon
				double3 ray_direction = normalize(
					forward + pixel_screen_x * right + pixel_screen_y * up
				);

				// Création du rayon
				ray = Ray (scene.camera.position, ray_direction);

				double z_depth;

				trace(scene, ray, ray_depth, &ray_color, &z_depth);
				avg_ray_color += ray_color;
				avg_z_depth += z_depth;
			}

        	// Moyenne des profondeurs et couleurs après les itérations
			avg_z_depth = avg_z_depth / scene.samples_per_pixel;
			avg_ray_color = avg_ray_color / scene.samples_per_pixel;

			// Test de profondeur
			if(avg_z_depth >= scene.camera.z_near && avg_z_depth <= scene.camera.z_far &&
				avg_z_depth < z_buffer[x + y*scene.resolution[0]]) {
				z_buffer[x + y*scene.resolution[0]] = avg_z_depth;

				// Met à jour la couleur de l'image (et sa profondeur)
				output->set_color_pixel(x, y, avg_ray_color);
				output->set_depth_pixel(x, y, (avg_z_depth - scene.camera.z_near) /
										(scene.camera.z_far-scene.camera.z_near));
			}
        }
    }

	// fin commente pour tester portho.ray
	// Libération du z_buffer
    delete[] z_buffer;
}

// @@@@@@ VOTRE CODE ICI
// Veuillez remplir les objectifs suivants:
// 		- Détermine si le rayon intersecte la géométrie.
//      	- Calculer la contribution associée à la réflexion.
//			- Calculer la contribution associée à la réfraction.
//			- Mettre à jour la couleur avec le shading + 
//			  Ajouter réflexion selon material.reflection +
//			  Ajouter réfraction selon material.refraction 
//            pour la couleur de sortie.
//          - Mettre à jour la nouvelle profondeure.
void Raytracer::trace(const Scene& scene,
					  Ray ray, int ray_depth,
					  double3* out_color, double* out_z_depth)
{
	Intersection hit;
	// Fait appel à l'un des containers spécifiées.
	if(scene.container->intersect(ray,EPSILON,scene.camera.z_far,&hit)) {
		Material& material = ResourceManager::Instance()->materials[hit.key_material];

		double3 shade_color = shade(scene,hit);

		// @@@@@@ VOTRE CODE ICI
		// Déterminer la couleur associée à la réflection d'un rayon de manière récursive.
		double3 reflection_color = {0, 0, 0};
		if (material.k_reflection > 0 && ray_depth < 16 ) {
			double3 reflected_dir = normalize(ray.direction - 2 * dot(ray.direction, hit.normal) * hit.normal);
			Ray reflected_ray(hit.position, reflected_dir);

			double3 temp_color = {0, 0, 0};
			double temp_z_depth = std::numeric_limits<double>::max();

			trace(scene, reflected_ray, ray_depth + 1, &temp_color, &temp_z_depth);

			reflection_color = material.k_reflection * temp_color;
		}
		// @@@@@@ VOTRE CODE ICI
		// Déterminer la couleur associée à la réfraction d'un rayon de manière récursive.
		//
		// Assumez que l'extérieur/l'air a un indice de réfraction de 1.
		//
		// Toutes les géométries sont des surfaces et non pas de volumes.

		double3 refraction_color = {0, 0, 0};
		if (material.k_refraction > 0 && ray_depth < 16) {
			double eta = 1.0 / material.refractive_index;
			double cos_i = dot(-ray.direction, hit.normal);
			if (cos_i < 0) {
				eta = material.refractive_index;
				cos_i = fabs(cos_i);
			}
			double sin_t2 = eta * eta * (1 - cos_i * cos_i);

			if (sin_t2 <= 1) {
				double3 refracted_dir = normalize(eta * ray.direction +
												  (eta * cos_i - sqrt(1 - sin_t2)) * hit.normal);
				Ray refracted_ray(hit.position, refracted_dir);

				double3 temp_color = {0, 0, 0};
				double temp_z_depth = std::numeric_limits<double>::max();

				trace(scene, refracted_ray, ray_depth + 1, &temp_color, &temp_z_depth);

				refraction_color = material.k_refraction * temp_color;
			}
		}
		// *out_color =
		*out_color = shade_color + reflection_color + refraction_color;
		// *out_z_depth =
		*out_z_depth = hit.depth;

	}else {
		// Si aucune intersection, la couleur de sortie est noire
		*out_color = {0, 0, 0};
		*out_z_depth = std::numeric_limits<double>::max();
	}
}


// @@@@@@ VOTRE CODE ICI
// Veuillez remplir les objectifs suivants:
// 		* Calculer la contribution des lumières dans la scène.
//			- Itérer sur toutes les lumières.
//				- Inclure la contribution spéculaire selon le modèle de Blinn en incluant la composante métallique.
//	          	- Inclure la contribution diffuse. (Faites attention au produit scalare. >= 0)
//   	  	- Inclure la contribution ambiante
//      * Calculer si le point est dans l'ombre
//			- Itérer sur tous les objets et détecter si le rayon entre l'intersection et la lumière est occludé.
//				- Ne pas considérer les points plus loins que la lumière.
//			- Par la suite, intégrer la pénombre dans votre calcul
//		* Déterminer la couleur du point d'intersection.
//        	- Si texture est présente, prende la couleur à la coordonnées uv
//			- Si aucune texture, prendre la couleur associé au matériel.

double3 Raytracer::shade(const Scene& scene, Intersection hit) {
    Material& material = ResourceManager::Instance()->materials[hit.key_material];
    bool has_texture = false;
    double3 texture;

    // Ambiante
    double3 couleur = scene.ambient_light * material.k_ambient;

	if(material.texture_albedo.width() > 0 && material.texture_albedo.height() > 0) {
		has_texture = true;
		int u = static_cast<int>(hit.uv[0] * material.texture_albedo.width());
		int v = static_cast<int>(hit.uv[1] * material.texture_albedo.height());
		rgb_t texture_point = material.texture_albedo.get_pixel(u, v);
		double normalized_r = texture_point.red / 255.0;
		double normalized_g = texture_point.green / 255.0;
		double normalized_b = texture_point.blue / 255.0;
		texture = {normalized_r, normalized_g, normalized_b};
		couleur *= texture;
	}else {
		couleur *= material.color_albedo;
	}

    // Itère sur chaque lumière pour diffuse et spéculaire
    for (const SphericalLight& spherical_light : scene.lights) {
        double3 lightDir = normalize(spherical_light.position - hit.position);
        double3 viewDir = normalize(scene.camera.position - hit.position); // Direction vers la caméra
        double3 halfDir = normalize(lightDir + viewDir); // Direction moyenne pour Blinn
        double3 difference = (hit.position - spherical_light.position);
        double distance_squared = dot(difference, difference); // Distance entre point et lumière
        double distance = std::sqrt(distance_squared);

        // Ombres
        double terme_eclairage = 1.0;
        Intersection shadow_hit;



                int num_samples = 16;
                double occlusion = 0.0;

                for (int i = 0; i < num_samples; ++i) {
                    double2 disk_offset_2d = random_in_unit_disk();
                    double3 disk_offset = {disk_offset_2d.x * spherical_light.radius, disk_offset_2d.y * spherical_light.radius, 0.0};
                    double3 sample_pos = spherical_light.position + disk_offset;
                    double3 sample_dir = normalize(sample_pos - hit.position);
                    double sample_distance = length(sample_pos - hit.position);

                    Ray sample_ray(hit.position + hit.normal * EPSILON, sample_dir);
                    Intersection shadow_sample_hit;

                    if (!scene.container->intersect(sample_ray, EPSILON, sample_distance, &shadow_sample_hit)) {
                        occlusion += 1.0;
                    }
                }

                terme_eclairage = occlusion / num_samples;


        // Calcul de la composante diffuse et spéculaire
        double d = std::max(dot(hit.normal, lightDir), 0.0);
        double diffuse = d * material.k_diffuse * terme_eclairage;

        double spec = std::max(dot(hit.normal, halfDir), 0.0);
        double speculaire = pow(spec, material.shininess) * material.k_specular * terme_eclairage;

        // Couleur du matériau ou de la texture
        double3 base_color = has_texture ? texture : material.color_albedo;
		double3 metallic = material.metallic * base_color + (1 - material.metallic);
		double ponctuelle = 1.0 / distance_squared;

        // Application des contributions avec la pénombre et ombre complète
        couleur += ponctuelle * (diffuse * spherical_light.emission) * base_color; // Ajoute la contribution diffuse
		couleur += metallic * (ponctuelle * (speculaire * spherical_light.emission));  // Ajoute la contribution spéculaire
    }

    return couleur;
}
