#include "vector_algebra_utils.h"

namespace
{
using Matrix3 = std::array<std::array<double, 3>, 3>;

Matrix3 getIdentityMatrix()
{
  Matrix3 res;
  for (std::size_t i = 0; i < 3; i++)
    for (std::size_t j = 0; j < 3; j++) res[i][j] = i == j ? 1.0 : 0.0;
  return res;
}

Matrix3 getZeroMatrix()
{
  Matrix3 res;
  for (std::size_t i = 0; i < 3; i++)
    for (std::size_t j = 0; j < 3; j++) res[i][j] = 0.0;
  return res;
}

Affine toAffine(const Matrix3& rotation)
{
  return Affine(rotation[0][0], rotation[0][1], rotation[0][2], rotation[1][0], rotation[1][1], rotation[1][2],
                rotation[2][0], rotation[2][1], rotation[2][2], 1.);
}

Affine toAffine(const Matrix3& rotation, const Vector_3& translation)
{
  return Affine(rotation[0][0], rotation[0][1], rotation[0][2], translation[0], rotation[1][0], rotation[1][1],
                rotation[1][2], translation[1], rotation[2][0], rotation[2][1], rotation[2][2], translation[2], 1.);
}

Matrix3 toRotation(const Affine& affine)
{
  Matrix3 res;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) res[i][j] = affine.m(i, j);
  return res;
}

Matrix3 asymmetricTensorProduct(const Vector_3& x, const Vector_3& y)
{
  Matrix3 res;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) res[i][j] = x[i] * y[j] - x[j] * y[i];
  return res;
}

Matrix3 operator+(const Matrix3& lhs, const Matrix3& rhs)
{
  Matrix3 res;
  for (std::size_t i = 0; i < 3; i++)
    for (std::size_t j = 0; j < 3; j++) res[i][j] = lhs[i][j] + rhs[i][j];
  return res;
}

Matrix3 operator*(const Matrix3& lhs, const Matrix3& rhs)
{
  Matrix3 res = getZeroMatrix();
  for (std::size_t i = 0; i < 3; i++)
    for (std::size_t j = 0; j < 3; j++)
      for (std::size_t k = 0; k < 3; k++) res[i][j] += lhs[i][k] * rhs[k][j];
  return res;
}

Matrix3 operator*(const Matrix3& lhs, double rhs)
{
  Matrix3 res;
  for (std::size_t i = 0; i < 3; i++)
    for (std::size_t j = 0; j < 3; j++) res[i][j] = lhs[i][j] * rhs;
  return res;
}

// fuzzy comparison ops using length threshold of eps

bool isZero(double a, double eps) { return std::fabs(a) < std::fabs(eps); }

bool isOne(double a, double eps) { return isZero(1. - std::fabs(a), eps); }

bool pointsEqual(const Point_3& p, const Point_3& q, double eps) { return (p - q).squared_length() < eps * eps; }

bool pointsCollinear(const Point_3& p, const Point_3& q, const Point_3& r, double eps)
{
  return CGAL::cross_product(q - p, r - p).squared_length() < eps * eps * eps * eps;
}

}  // namespace

/* Returns the unit vector
 *
 * if input is 0 vector it returns the input
 *
 */
Vector_3 unitVector(const Vector_3& t)
{
  if (t == Vector_3(0, 0, 0)) return t;
  return t / std::sqrt(t.squared_length());
}

/* Returns the centroid of the points
 *
 */
Point_3 getCentroid(const std::vector<Point_3>& points)
{
  if (points.empty()) return Point_3(0., 0., 0.);

  if (points.size() == 1) return points[0];

  Point_3 origin(0., 0., 0.);
  Vector_3 center(0., 0., 0.);
  for (const auto& p : points) center += (p - origin);
  return origin + center / static_cast<double>(points.size());
}

/* computes the cos angle between vectors t1 and t2
 *
 */
std::optional<double> cosAngle(const Vector_3& t1, const Vector_3& t2, bool sign)
{
  Vector_3 zero(0, 0, 0);
  if (t1 == zero || t2 == zero) return std::nullopt;
  double ca = scalar_product(unitVector(t1), unitVector(t2));
  return sign ? ca : std::fabs(ca);
}

/* Computes the cos of angle between vectors (p2-p1) and (p3-p2)
 *
 */
std::optional<double> cosAngle(Point_3 p1, Point_3 p2, Point_3 p3, bool sign)
{
  return cosAngle(unitVector(p2 - p1), unitVector(p3 - p2), sign);
}

/* Computes the difference vector
 *
 * Note that size of the output would be 1 less size of the input
 *
 */
template <class T, class S>
void getDifference(const T& in, S& out)
{
  out.clear();
  for (size_t i = 0; i < in.size() - 1; i++) out.push_back(in[i + 1] - in[i]);
}

/* Projects a vector onto a plane specified by the normal n
 *
 */
Vector_3 projectVector2Plane(const Vector_3& n, const Vector_3& v) { return v - scalar_product(v, n) * n; }

/* returns a rotation that takes vector "a" (direction only) a to vector "b" (direction only)
 *
 *
 */
std::optional<Affine> getRotation(const Vector_3& a, const Vector_3& b)
{
  if (a == Vector_3(0, 0, 0) || b == Vector_3(0, 0, 0)) return std::nullopt;

  Vector_3 m = unitVector(a);
  Vector_3 n = unitVector(b);

  if (m == -n) return std::nullopt;

  if (m == n) return toAffine(getIdentityMatrix());

  Vector_3 mxn = CGAL::cross_product(m, n);
  double s2 = mxn.squared_length();
  if (s2 > 0)
  {
    Vector_3 x = m;
    Vector_3 y = n;

    double factor = 1.0 / (1.0 + CGAL::scalar_product(x, y));
    const auto S = asymmetricTensorProduct(y, x);
    const auto SS = S * S;
    const auto R = getIdentityMatrix() + S + SS * factor;

    return toAffine(R);
  }

  return toAffine(getIdentityMatrix());
}

std::optional<Affine> alignAndTranslate(const Vector_3& reference, const Vector_3& target, const Vector_3& translate)
{
  // make sure reference and target are unit
  Vector_3 x = unitVector(reference);
  Vector_3 y = unitVector(target);

  const auto aff = getRotation(x, y);
  if (aff == std::nullopt) return std::nullopt;

  const auto rotation = toRotation(*aff);
  return toAffine(rotation, translate);
}

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
 * 5. the set {p1,p2} is the same as the set {q1, q2} --> std::nulloptA
 * 6. normal for p1-q2-p2 (n) and normal p1-p2-q2 make an angle that is close to 180 degrees
 *
 */
std::optional<double> foldingIntersection(const Point_3& p1,  // source of segment (used for unfolding)
                                          const Point_3& p2,  // target of segment (used for unfolding)
                                          const Point_3& q1, const Point_3& q2, double eps)
{
  // Note the following checks should not be exact checks. For example to check if a point p is the same as a point q
  // instead of using CGAL's p==q we use (p-q).squared_length() < eps*eps

  // check if p1 and p2 are the same
  if (pointsEqual(p1, p2, eps)) return std::nullopt;

  // check if q1 and q2 are the same
  if (pointsEqual(q1, q2, eps)) return std::nullopt;

  // check if both q1 and q2 are collinear with p1-p2
  if (pointsCollinear(p1, p2, q1, eps) && pointsCollinear(p1, p2, q2, eps)) return std::nullopt;

  // check if the line segments p1p2 and q1p2 are the same
  if ((pointsEqual(p1, q1, eps) && pointsEqual(p2, q2, eps)) || (pointsEqual(p1, q2, eps) && pointsEqual(p2, q1, eps)))
    return std::nullopt;

  // check if either q1 or q2 is the same as p1, and if so return 0.0
  if (pointsEqual(p1, q1, eps) || pointsEqual(p1, q2, eps)) return 0.0;

  // check if either q1 or q2 is the same as p2, and if so return 1.0
  if (pointsEqual(p2, q1, eps) || pointsEqual(p2, q2, eps)) return 1.0;

  // check if q1 is collinear with p1-p2, and if so use the line p1-p2 and the q1 to compute u.
  if (pointsCollinear(p1, p2, q1, eps)) return CGAL::scalar_product((q1 - p1), (p2 - p1)) / (p2 - p1).squared_length();

  // check if q2 is collinear with p1-p2, and if so use the line p1-p2 and the q2 to compute u.
  if (pointsCollinear(p1, p2, q2, eps)) return CGAL::scalar_product((q2 - p1), (p2 - p1)) / (p2 - p1).squared_length();

  // compute the normal for triangle p1-q1-p2
  Vector_3 n = unitVector(CGAL::cross_product(q1 - p1, p2 - p1));

  // compute the normal for triangle p1-p2-q2
  Vector_3 m = unitVector(CGAL::cross_product(p2 - p1, q2 - p1));

  if (CGAL::scalar_product(n, m) < -0.9) return std::nullopt;

  // find a rotation that transforms m to n
  auto rot = getRotation(m, n);
  if (rot == std::nullopt) return std::nullopt;
  Affine rotation = *rot;

  // right-handed coordinate triad t, s, n
  Vector_3 s = unitVector(p2 - p1);
  Vector_3 t = unitVector(CGAL::cross_product(s, n));

  Point_3 q2u = p1 + rotation.transform(q2 - p1);

  // now find the intersection of the line segments p1 -> p2  and q1 -> q2u
  // Note that this is technically a 2D intersection in the plane of triangle 1 as specified by
  // unit vectors s and t that satisfy s.t = 0. Assuming p1 is the origin the intersection can be found
  // by solving

  // p1 (1-u) + p2 u = q1 (1-v) + q2u v -> rearrange
  // p1 + (p2-p1) u = q1 + (q2u-q1) v   -> rearrange
  // (p2-p1) u - (q2u-q1) v = (q1-p1)  -> project to t and s directions to get 2 equations
  //
  // |(p2-p1).t  -(q2u-q1).t| [u] = [(q1-p1).t]
  // |(p2-p1).s  -(q2u-q1).s| [v] = [(q1-p1).s]
  //
  // |a  b| [u] = [e]
  // |c  d| [v] = [f]
  //
  // u = (de-bf)/(ad-bc)
  // v = (af-ce)/(ad-bc)

  // compute the coefficients
  double a = scalar_product(p2 - p1, t);
  double c = scalar_product(p2 - p1, s);
  double b = -scalar_product(q2u - q1, t);
  double d = -scalar_product(q2u - q1, s);
  double e = scalar_product(q1 - p1, t);
  double f = scalar_product(q1 - p1, s);

  if (std::fabs(a * d - b * c) > 0.0)
  {
    double u = (d * e - b * f) / (a * d - b * c);
    double v = (a * f - c * e) / (a * d - b * c);

    // if close to 0 (or 1), snap to 0 (or 1)
    if (isZero(u, eps)) u = 0.;
    if (isZero(v, eps)) v = 0.;
    if (isOne(u, eps)) u = 1.;
    if (isOne(v, eps)) v = 1.;

    return u;
  }
  // should never reach this point, just in case it did return nullopt
  return std::nullopt;
}
