#ifndef VECTOR_ALGEBRA_H
#define VECTOR_ALGEBRA_H

#include "cgal_typedefs.h"

#include <CGAL/Aff_transformation_3.h>

using Affine = CGAL::Aff_transformation_3<Kernel>;

/* Returns the unit vector
 *
 * if length of input is 0 it returns the input
 *
 */
Vector_3 unitVector(const Vector_3& t);

/* Returns the centroids of the points
 *
 */
Point_3 getCentroid(const std::vector<Point_3>& points);

/* computes the cos angle between vectors t1 and t2
 *
 * if any of t1 or t2 is a zero vector it returns std::nullopt
 *
 */
std::optional<double> cosAngle(const Vector_3& t1, const Vector_3& t2, bool sign = true);

/* Computes the cos of angle between vectors (p2-p1) and (p3-p2)
 *
 * if nay of (p2-p1) or (p3-p2) is a zero vector it returns std::nullopt
 *
 */
std::optional<double> cosAngle(Point_3 p1, Point_3 p2, Point_3 p3, bool sign = true);

/* Projects a vector onto a plane specified by the normal n
 *
 */
Vector_3 projectVector2Plane(const Vector_3& n, const Vector_3& v);

/* returns a rotation that takes vector "a" (direction only) to vector "b" (direction only)
 *
 * if either of a or b is a zero vector, returns std::nullopt
 * if a is exactly in the direction of -b, returns std::nullopt
 */
std::optional<Affine> getRotation(const Vector_3& a, const Vector_3& b);

/* Returns the Affine Transformation corresponding to aligning the "ref" direction to
 * the "target" direction, and then translating with "translate".
 *
 * returns std::nullopt in following cases
 * 1. reference, target or both are zero vectors
 * 2. reference and target are parallel but pointing in opposite directions
 *
 */
std::optional<Affine> alignAndTranslate(const Vector_3& reference,   // reference direction
                                        const Vector_3& target,      // target direction
                                        const Vector_3& translate);  // translations

/* Intersects the line connecting q1 to q2 with the line connecting p1 to p2, allowing a fold around p1-p2
 *
 * a) all points are in the same plane --> normal segment (q1-q2) to segment (p1-p2) intersection.
 * b) all points are not in the same plane -->
 *      1. fold p1-p2-q2 (about p1p2) to get p1-p2-q2u in the same plane as p1-q1-p2 (orders are important)
 *      2. perform a regular (q1-q2u) to (p1-p2) intersection
 *
 * Outcome:
 * a parametric value "u" the specifies the intersection point p(u) = (1-u).p1 + u.p2
 * 1. u<0   --> intersection point outside interval on the p1 side
 * 2. u=0   --> intersection point exactly on p1
 * 3. 0<u<1 --> intersection point between p1 and p2
 * 4. u=1   --> intersection point exactly at p2
 * 5. u>0   --> intersection point outside interval on the p2 side
 *
 * Special cases:
 * 1. either q1 or q2 is on p1-p2 --> parametric value (u) of that q is return
 * 2. both q1 and q2 are on p1-p2 --> std::nullopt
 * 3. both p1 == p2 --> std::nullopt
 * 4. both q1 == q2 --> std::nullopt
 *
 */
std::optional<double> foldingIntersection(const Point_3& p1,  // source of segment (used for unfolding)
                                          const Point_3& p2,  // target of segmtet (used for unfolding)
                                          const Point_3& q1, const Point_3& q2, double eps = 1.0e-6);

#endif  // VECTOR_ALGEBRA_H
