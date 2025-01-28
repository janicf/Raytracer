#include "object.h"
#include "../extern/linalg/linalg.h"  // Inclusion de la librairie linalg.h

#include <cmath>  // Pour sqrt() et atan2()
#include "basic.h"



using namespace linalg::aliases; // Simplifie les types (double3, etc.)



// Fonction retournant soit la valeur v0 ou v1 selon le signe.
int rsign(double value, double v0, double v1) {
	return (int(std::signbit(value)) * (v1-v0)) + v0;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection d'une sphère.
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Sphere::local_intersect(Ray ray,
							 double t_min, double t_max,
							 Intersection *hit)
{
	// Centre de la sphère à l'origine
	double3 center = {0, 0, 0};

	// Calcul du vecteur de l'origine du rayon au centre de la sphère
	double3 oc = ray.origin - center;

	// Coefficients de l'équation quadratique
	double a = dot(ray.direction, ray.direction);
	double b = 2.0 * dot(oc, ray.direction);
	double c = dot(oc, oc) - radius * radius;

	// Calcul du discriminant
	double discriminant = b * b - 4 * a * c;

	// Calcul du discriminant (Utilisation de EPSILON pour la gestion des erreurs numériques)
	if (discriminant < EPSILON) return false;

	// Calcul des racines de l'équation quadratique
	double sqrt_discriminant = sqrt(discriminant);

	// Trouve la plus petite racine dans l'intervalle [t_min, t_max]
	double root = (-b - sqrt_discriminant) / (2.0 * a);

	if (root < t_min || root > t_max) {
		root = (-b + sqrt_discriminant) / (2.0 * a);
		if (root < t_min || root > t_max) return false;
	}

	// Si nous avons une intersection valide, mettons à jour les informations
	hit->depth = root;
	hit->position = ray.origin + root * ray.direction;
	hit->normal = normalize(hit->position - center);

	// Calcul des coordonnées UV (paramétrisation de la sphère)

	double u = atan2(hit->normal.x, hit->normal.z) / (2 * PI);
	double v = 0.5 - (asin(hit->normal.y) / PI);
	hit->uv = {u, v};

	return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour la sphère.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Sphere::compute_aabb() {
	// Calcul du AABB en espace local
	double3 local_min = {-radius, -radius, -radius};
	double3 local_max = {radius, radius, radius};

	// Transforme les coins dans l'espace global
	std::vector<double3> corners = retrieve_corners(AABB{local_min, local_max});
	for (double3& corner : corners) {
		corner = mul(transform, {corner, 1}).xyz(); // Passage dans l'espace global
	}

	// Construit et retourne le AABB englobant global
	return construct_aabb(corners);
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un quad (rectangle).
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Quad::local_intersect(Ray ray,
							double t_min, double t_max,
							Intersection *hit)
{
	// Plan forme (A,B,C,D) avec normale Z+
	double A = 0, B = 0, C = 1, D = 0;

	// Origine du rayon et direction
	double x0 = ray.origin.x;
	double y0 = ray.origin.y;
	double z0 = ray.origin.z;
	double dx = ray.direction.x;
	double dy = ray.direction.y;
	double dz = ray.direction.z;

	// Calcul dénominateur
	double denom = A * dx + B * dy + C * dz; // AΔx + BΔy + CΔz

	// Si dénominateur == 0, rayon parallèle au plan et il n'y a pas d'intersection
	if (std::abs(denom) == 0) {
		return false;
	}

	// Calcule t
	double t = -(A * x0 + B * y0 + C * z0 + D) / denom;

	// Vérifie si t dans les bornes
	if (t < t_min || t > t_max) {
		return false;
	}

	// Calcul point d'intersection
	double3 intersect_point = ray.origin + t * ray.direction;

	// Vérifie point d'intersection dans les bornes
	if (std::abs(intersect_point.x) > half_size || std::abs(intersect_point.y) > half_size) {
		return false;
	}

	if(denom > 0) {
		hit->normal = {0,0,-1};
	}else {
		hit->normal = {0,0,1};
	}

	// Si intersection, update information
	hit->depth = t;
	hit->position = intersect_point;

	// Calcul des coordonnées UV
	hit->uv = {
		(intersect_point.x + half_size) / (2 * half_size),
		 (-intersect_point.y + half_size) / (2 * half_size)
	};

	return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le quad (rectangle).
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Quad::compute_aabb() {
	// Calcul du AABB en espace local
	double3 local_min = {-half_size, -half_size, 0};
	double3 local_max = {half_size, half_size, 0};

	// Transforme les coins dans l'espace global
	std::vector<double3> corners = retrieve_corners(AABB{local_min, local_max});
	for (double3& corner : corners) {
		corner = mul(transform, {corner, 1}).xyz(); // Passage dans l'espace global
	}

	// Construit et retourne le AABB englobant global
	return construct_aabb(corners);;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un cylindre.
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Cylinder::local_intersect(Ray ray,
							   double t_min, double t_max,
							   Intersection *hit)
{
	// Origine et direction du rayon
    double x0 = ray.origin.x;
	double y0 = ray.origin.y;
	double z0 = ray.origin.z;
	double dx = ray.direction.x;
	double dy = ray.direction.y;
	double dz = ray.direction.z;

	// Vérifie si le rayon est parallèle au cylindre
	if (std::abs(dx) == 0 && std::abs(dz) == 0) return false;

	// Coefficients de l'équation quadratique
	double a = dx * dx + dz * dz;
	double b = 2 * (dx * x0 + dz * z0);
	double c = x0 * x0 + z0 * z0 - radius * radius;

	// Calcul du discriminant
	double discriminant = b * b - 4 * a * c;
	if (discriminant < 0) return false; // Pas d'intersection réelle

	// Calcul des racines (deux intersections possibles)
	double sqrt_discriminant = sqrt(discriminant);
	double root1 = (-b - sqrt_discriminant) / (2.0 * a);
	double root2 = (-b + sqrt_discriminant) / (2.0 * a);

	// Pour vérifier que l'intersection est dans la hauteur du cylindre
	double y1 = y0 + root1 * dy;
	double y2 = y0 + root2 * dy;

	bool hit_found = false;
	double3 normal;

	// Vérifie la première intersection
	if (root1 >= t_min && root1 <= t_max && y1 >= -half_height && y1 <= half_height) {
		hit->depth = root1;
		hit->position = ray.origin + root1 * ray.direction;
		normal = normalize(double3(hit->position.x, 0, hit->position.z));

		// Vérifie si le rayon touche l'intérieur du cylindre
		if (dot(ray.direction, normal) > 0) {
			hit->normal = -normal;
		} else {
			hit->normal = normal;
		}

		// Calcul des coordonnées UV pour root1
		double theta = atan2(hit->position.z, hit->position.x);
		double u = 1.0 - (theta / (2 * PI));
		if (u < 0) u += 1.0;
		double v = (-hit->position.y + half_height) / (2 * half_height);
		hit->uv = {u, v};

		hit_found = true;
	}

	// Vérifie la deuxième intersection
	if (root2 >= t_min && root2 <= t_max && y2 >= -half_height && y2 <= half_height) {
		if (root2 < hit->depth) {
			hit->depth = root2;
			hit->position = ray.origin + root2 * ray.direction;
			normal = normalize(double3(hit->position.x, 0, hit->position.z));

			// Vérifie si le rayon touche l'intérieur du cylindre
			if (dot(ray.direction, normal) > 0) {
				hit->normal = -normal;
			} else {
				hit->normal = normal;
			}

			// Calcul des coordonnées UV pour root2
			double theta = atan2(hit->position.z, hit->position.x);
			double u = 1.0 -(theta / (2 * PI));
			if (u < 0) u += 1.0;
			double v = (-hit->position.y + half_height) / (2 * half_height);
			hit->uv = {u, v};

			hit_found = true;
		}
	}

	// Retourne vrai si une intersection valide a été trouvée
	return hit_found;

}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le cylindre.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Cylinder::compute_aabb() {
	// Calcul du AABB en espace local
	double3 local_min = {-radius, -half_height, -radius};
	double3 local_max = {radius, half_height, radius};

	// Transforme les coins dans l'espace global
	std::vector<double3> corners = retrieve_corners(AABB{local_min, local_max});
	for (double3& corner : corners) {
		corner = mul(transform, {corner, 1}).xyz(); // Passage dans l'espace global
	}

	// Construit et retourne le AABB englobant global
	return construct_aabb(corners);
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un mesh.
//
// Référez-vous au PDF pour la paramétrisation pour les coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
//
bool Mesh::local_intersect(Ray ray,
						   double t_min, double t_max,
						   Intersection* hit)
{

	bool intersected = false;
	double closest_t = t_max;

	// Parcourt tous les triangles du maillage
	for (const auto& tri : triangles) {
		Intersection temp_hit;

		// Test d'intersection avec chaque triangle du maillage
		if (intersect_triangle(ray, t_min, closest_t, tri, &temp_hit)) {
			intersected = true;

			// Mise à jour du point d'intersection le plus proche trouvé
			closest_t = temp_hit.depth;
			*hit = temp_hit;  // Copie de l'information d'intersection
		}
	}
	return intersected;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un triangle.
// S'il y a intersection, remplissez hit avec l'information sur la normale et les coordonnées texture.
bool Mesh::intersect_triangle(Ray  ray,
							  double t_min, double t_max,
							  Triangle const tri,
							  Intersection *hit)
{
	// Extrait chaque position de sommet des données du maillage.
	double3 const &p0 = positions[tri[0].pi]; // ou Sommet A (Pour faciliter les explications)
	double3 const &p1 = positions[tri[1].pi]; // ou Sommet B
	double3 const &p2 = positions[tri[2].pi]; // ou Sommet C

	// Triangle en question. Respectez la convention suivante pour vos variables.
	//
	//     A
	//    / \
	//   /   \
	//  B --> C
	//
	// Respectez la règle de la main droite pour la normale.

	// @@@@@@ VOTRE CODE ICI
	// Décidez si le rayon intersecte le triangle (p0,p1,p2).
	// Si c'est le cas, remplissez la structure hit avec les informations
	// de l'intersection et renvoyez true.
	// Pour plus de d'informations sur la géométrie, référez-vous à la classe dans object.hpp.
	//
	// NOTE : hit.depth est la profondeur de l'intersection actuellement la plus proche,
	// donc n'acceptez pas les intersections qui occurent plus loin que cette valeur.

	// Calcul des arêtes du triangle
	double3 edge1 = p1 - p0;
	double3 edge2 = p2 - p0;

	// Calcul du déterminant pour détecter le parallélisme
	double3 pvec = cross(ray.direction, edge2);
	double det = dot(edge1, pvec);

	// Si le déterminant est proche de 0, le rayon est parallèle au triangle
	if (fabs(det) < EPSILON) return false;

	double inv_det = 1.0 / det;

	// Calcul du vecteur tvec
	double3 tvec = ray.origin - p0;

	// Calcul des coordonnées barycentriques u et v
	double u = dot(tvec, pvec) * inv_det;
	if (u < 0 || u > 1) return false;

	double3 qvec = cross(tvec, edge1);
	double v = dot(ray.direction, qvec) * inv_det;
	if (v < 0 || u + v > 1) return false;

	// Calcul de t (la profondeur)
	double t = dot(edge2, qvec) * inv_det;

	// Vérification si t est dans les bornes
	if (t < t_min || t > t_max) return false;

	// Mise à jour des informations d'intersection
	hit->depth = t;
	hit->position = ray.origin + t * ray.direction;

	// Interpolation des normales en fonction des coordonnées barycentriques
	double3 n0 = normals[tri[0].ni];
	double3 n1 = normals[tri[1].ni];
	double3 n2 = normals[tri[2].ni];
	hit->normal = normalize((1 - u - v) * n0 + u * n1 + v * n2);

	// Interpolation des coordonnées UV
	double2 uv0 = tex_coords[tri[0].ti];
	double2 uv1 = tex_coords[tri[1].ti];
	double2 uv2 = tex_coords[tri[2].ti];
	hit->uv = (1 - u - v) * uv0 + u * uv1 + v * uv2;

	return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le Mesh.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Mesh::compute_aabb() {
	// Initialise les limites du AABB
	double3 min_point = {DBL_MAX, DBL_MAX, DBL_MAX};
	double3 max_point = {-DBL_MAX, -DBL_MAX, -DBL_MAX};

	// Parcourt toutes les positions du mesh pour trouver le AABB minimal
	for (const double3& position : positions) {
		double3 global_position = mul(transform, {position, 1}).xyz();
		min_point = min(min_point, global_position);
		max_point = max(max_point, global_position);
	}

	return AABB{min_point, max_point};
}