#include "aabb.h"

#include <iostream>

// @@@@@@ VOTRE CODE ICI
// Implémenter l'intersection d'un rayon avec un AABB dans l'intervalle décrit.
bool AABB::intersect(Ray ray, double t_min, double t_max)  {

    // Boucle sur les 3 axes (x, y, z) pour vérifier l'intersection avec les faces de la boîte AABB
    for (int i = 0; i < 3; i++) {
        // Inverse de la direction du rayon sur l'axe actuel
        const double invD = 1.0 / ray.direction[i];

        // Calcul des distances d'entrée et de sortie du rayon sur l'axe
        double t0 = (min[i] - ray.origin[i]) * invD;
        double t1 = (max[i] - ray.origin[i]) * invD;

        // Si la direction est négative, on inverse t0 et t1 pour avoir les bons points d'intersection
        if (invD < 0.0) std::swap(t0, t1);

        // Ajustement de t_min et t_max pour inclure EPSILON comme marge de tolérance
        t_min = std::max(t_min, t0 - EPSILON);
        t_max = std::min(t_max, t1 + EPSILON);

        // Vérifie si les intervalles d'entrée et de sortie sont toujours valides
        // Si t_max est strictement inférieur à t_min, alors le rayon n'intersecte pas la boîte AABB
        if (t_max + EPSILON <= t_min) return false;
    }
    // Si toutes les conditions sont remplies, il y a une intersection
    return true;
};


// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction qui permet de trouver les 8 coins de notre AABB.
std::vector<double3> retrieve_corners(AABB aabb) {
    std::vector<double3> corners;
    corners.push_back(aabb.min);
    corners.push_back(double3{aabb.min.x, aabb.min.y, aabb.max.z});
    corners.push_back(double3{aabb.min.x, aabb.max.y, aabb.min.z});
    corners.push_back(double3{aabb.max.x, aabb.min.y, aabb.min.z});
    corners.push_back(double3{aabb.max.x, aabb.max.y, aabb.min.z});
    corners.push_back(double3{aabb.max.x, aabb.min.y, aabb.max.z});
    corners.push_back(double3{aabb.min.x, aabb.max.y, aabb.max.z});
    corners.push_back(aabb.max);
    return corners;
};

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction afin de créer un AABB qui englobe tous les points.
AABB construct_aabb(std::vector<double3> points) {
    // Initialise les coins minimum et maximum de la boîte AABB
    double3 min_point{DBL_MAX, DBL_MAX, DBL_MAX};
    double3 max_point{-DBL_MAX, -DBL_MAX, -DBL_MAX};

    // Boucle sur tous les points donnés pour trouver les limites de la boîte AABB
    for (const auto& point : points) {
        // Met à jour min et max
        min_point = min(min_point, point);
        max_point = max(max_point, point);
    }

    // Retourne la boîte AABB définie par les points minimum et maximum trouvés
    return AABB{min_point, max_point};
};

AABB combine(AABB a, AABB b) {
    return AABB{min(a.min,b.min),max(a.max,b.max)};
};

bool compare(AABB a, AABB b, int axis){
    return a.min[axis] < b.min[axis];
};