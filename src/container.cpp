#include "container.h"
#include <limits>  // Pour std::numeric_limits

// @@@@@@ VOTRE CODE ICI
// - Parcourir l'arbre DEPTH FIRST SEARCH selon les conditions suivantes:
// 		- S'il s'agit d'une feuille, faites l'intersection avec la géométrie.
//		- Sinon, il s'agit d'un noeud altérieur.
//			- Faites l'intersection du rayon avec le AABB gauche et droite. 
//				- S'il y a intersection, ajouter le noeud à ceux à visiter. 
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool BVH::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
    // Initialisation de la pile pour stocker les noeuds à traiter
    std::vector<BVHNode*> stack;
    stack.push_back(root);

    bool hit_found = false;      // Indicateur pour vérifier si une intersection a été trouvée
    double closest_t = t_max;    // Distance de la plus proche intersection trouvée

    // Boucle noeuds à traiter dans la pile
    while (!stack.empty()) {
        // Récupère le noeud en haut de la pile et le retire
        BVHNode* node = stack.back();
        stack.pop_back();

        // Vérifie si le rayon intersecte la boîte englobante (AABB) du noeud
        if (node->aabb.intersect(ray, t_min, closest_t)) {
            // Cas feuille : vérifier l'intersection avec la géométrie de l'objet contenu
            if (node->left == nullptr && node->right == nullptr) {
                Object* obj = objects[node->idx];  // Récupère l'objet dans la feuille
                Intersection temp_hit;

                // Test d'intersection avec l'objet de la feuille
                if (obj->intersect(ray, t_min, closest_t, &temp_hit)) {
                    // Mise à jour si une intersection plus proche est trouvée
                    if (temp_hit.depth < closest_t) {
                        hit_found = true;
                        closest_t = temp_hit.depth;  // Mise à jour de la profondeur minimale trouvée
                        *hit = temp_hit;             // Mise à jour des détails de l'intersection
                    }
                }
            }
            // Cas nœud interne
            else {
                if (node->left) stack.push_back(node->left);    // Empile le nœud enfant gauche
                if (node->right) stack.push_back(node->right);  // Empile le nœud enfant droit
            }
        }
    }
    return hit_found;
}


// @@@@@@ VOTRE CODE ICI
// - Parcourir tous les objets
// 		- Détecter l'intersection avec l'AABB
//			- Si intersection, détecter l'intersection avec la géométrie.
//				- Si intersection, mettre à jour les paramètres.
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool Naive::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {

	Intersection closest_hit;  // Stocke l'intersection la plus proche
	closest_hit.depth = std::numeric_limits<double>::max();  // Initialise à une valeur maximale
	bool has_hit = false;  // Indicateur pour vérifier s’il y a eu une intersection

	// Parcours de tous les objets du container
	for (const auto& object : objects) {
		Intersection temp_hit;

		// Test d’intersection entre le rayon et l'AABB de l'objet
		//if (object->compute_aabb().intersect(ray, t_min, t_max)) {
			// Si le rayon intersecte l'AABB, alors on teste l'intersection avec la géométrie
			if (object->intersect(ray, t_min, t_max, &temp_hit)) {
				// Vérifie si l'intersection trouvée est plus proche que la précédente
				if (temp_hit.depth < closest_hit.depth) {
					closest_hit = temp_hit;  // Met à jour l'intersection la plus proche
					has_hit = true;  // Marque qu'une intersection a été trouvée
				}
			}
		//}
	}

	// Si une intersection a été trouvée, met à jour le pointeur 'hit'
	if (has_hit) {
		*hit = closest_hit;
		return true;  // Retourne 'true' si une intersection a été détectée
	}

	return false;  // Retourne 'false' si aucune intersection n'a été trouvée
}
