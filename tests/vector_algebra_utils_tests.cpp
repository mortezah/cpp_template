#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "vector_algebra_utils.h"

typedef Point_3 Point_3;
typedef Vector_3 Vector_3;

TEST(VectorAlgebraUtilities, CanComputeUnitVectorForZeroVector)
{
  Vector_3 zeroVector(0.0, 0.0, 0.0);
  EXPECT_EQ(unitVector(zeroVector), zeroVector);
}

TEST(VectorAlgebraUtilities, CanComputeUnitVectorForNonZeroVector)
{
  Vector_3 vec(1.0, 1.0, 1.0);
  Vector_3 expectedVec = (1. / std::sqrt(3.0)) * vec;
  EXPECT_EQ(unitVector(vec), expectedVec);
}

TEST(VectorAlgebraUtilities, CannotComputeCentroidForEmptySet)
{
  std::vector<Point_3> points;
  EXPECT_EQ(getCentroid(points), Point_3(0, 0, 0));
}

TEST(VectorAlgebraUtilities, CanComputeCentroidForPointsAtSquareCorners)
{
  std::vector<Point_3> points = {Point_3(0, 0, 0), Point_3(1, 0, 0), Point_3(0, 1, 0), Point_3(1, 1, 0)};
  EXPECT_EQ(getCentroid(points), Point_3(0.5, 0.5, 0.0));
}

TEST(VectorAlgebraUtilities, CanComputeCentroidForPointsAtCubeCorners)
{
  std::vector<Point_3> points = {Point_3(0, 0, 0), Point_3(1, 0, 0), Point_3(0, 1, 0), Point_3(1, 1, 0),
                                 Point_3(0, 0, 1), Point_3(1, 0, 1), Point_3(0, 1, 1), Point_3(1, 1, 1)};
  EXPECT_EQ(getCentroid(points), Point_3(0.5, 0.5, 0.5));
}

TEST(VectorAlgebraUtilities, CannotComputeCosAngleForZeroVectors)
{
  Vector_3 t1(0, 0, 0);
  Vector_3 t2(0, 0, 0);
  Vector_3 t3(1, 0, 0);
  EXPECT_EQ(cosAngle(t1, t2), std::nullopt);
  EXPECT_EQ(cosAngle(t2, t1), std::nullopt);
  EXPECT_EQ(cosAngle(t1, t3), std::nullopt);
  EXPECT_EQ(cosAngle(t3, t1), std::nullopt);
  EXPECT_EQ(cosAngle(t2, t3), std::nullopt);
  EXPECT_EQ(cosAngle(t3, t2), std::nullopt);
}

TEST(VectorAlgebraUtilities, CanComputeCosAngleForCloseToZeroVectors)
{
  Vector_3 t1(0.001, 0.000, 0.000);
  Vector_3 t2(0.000, 0.001, 0.000);
  Vector_3 t3(0.000, 0.002, 0.000);
  Vector_3 t4(0.001, 0.000, 0.001);
  Vector_3 t5(-0.001, 0.000, 0.000);

  // t1 and t2 at 90 degrees -- expected outcome 0
  auto t1t2 = cosAngle(t1, t2);
  EXPECT_NE(t1t2, std::nullopt);
  EXPECT_NEAR(*t1t2, 0.0, 1e-12);

  // t2 and t3 at 0 degrees -- expected outcome 1
  auto t2t3 = cosAngle(t2, t3);
  EXPECT_NE(t2t3, std::nullopt);
  EXPECT_NEAR(*t2t3, 1.0, 1e-12);

  // t4 and t1 at 45 degrees -- expected outcome 1/sqrt(2)
  auto t4t1 = cosAngle(t4, t1);
  EXPECT_NE(t4t1, std::nullopt);
  EXPECT_NEAR(*t4t1, 1.0 / std::sqrt(2.0), 1e-12);

  // t1 and t5 at 180 degrees -- expected outcome -1
  auto t1t5 = cosAngle(t1, t5);
  EXPECT_NE(t1t5, std::nullopt);
  EXPECT_NEAR(*t1t5, -1.0, 1e-12);
}

TEST(VectorAlgebraUtilities, CanComputeCosAngleForClockVectors)
{
  double theta = 2 * M_PI / 12;
  std::vector<Vector_3> hours = {
      Vector_3(std::sin(0. * theta), std::cos(0. * theta), 0),    // 0
      Vector_3(std::sin(1. * theta), std::cos(1. * theta), 0),    // 1
      Vector_3(std::sin(2. * theta), std::cos(2. * theta), 0),    // 2
      Vector_3(std::sin(3. * theta), std::cos(3. * theta), 0),    // 3
      Vector_3(std::sin(4. * theta), std::cos(4. * theta), 0),    // 4
      Vector_3(std::sin(5. * theta), std::cos(5. * theta), 0),    // 5
      Vector_3(std::sin(6. * theta), std::cos(6. * theta), 0),    // 6
      Vector_3(std::sin(7. * theta), std::cos(7. * theta), 0),    // 7
      Vector_3(std::sin(8. * theta), std::cos(8. * theta), 0),    // 8
      Vector_3(std::sin(9. * theta), std::cos(9. * theta), 0),    // 9
      Vector_3(std::sin(10. * theta), std::cos(10. * theta), 0),  // 10
      Vector_3(std::sin(11. * theta), std::cos(11. * theta), 0),  // 11
      Vector_3(std::sin(12. * theta), std::cos(12. * theta), 0)   // 12
  };

  Vector_3 h0 = hours[0];
  for (std::size_t i = 1; i < hours.size() - 1; i++)
  {
    auto ca = cosAngle(h0, hours[i]);
    EXPECT_NE(ca, std::nullopt);
    EXPECT_NEAR(*ca, std::cos(i * theta), 1.e-12);
  }
}

TEST(VectorAlgebraUtilities, CanComputeCosAngleForClockPoints)
{
  double theta = 2 * M_PI / 12;
  std::vector<Vector_3> hours = {
      Vector_3(std::sin(0. * theta), std::cos(0. * theta), 0),    // 0
      Vector_3(std::sin(1. * theta), std::cos(1. * theta), 0),    // 1
      Vector_3(std::sin(2. * theta), std::cos(2. * theta), 0),    // 2
      Vector_3(std::sin(3. * theta), std::cos(3. * theta), 0),    // 3
      Vector_3(std::sin(4. * theta), std::cos(4. * theta), 0),    // 4
      Vector_3(std::sin(5. * theta), std::cos(5. * theta), 0),    // 5
      Vector_3(std::sin(6. * theta), std::cos(6. * theta), 0),    // 6
      Vector_3(std::sin(7. * theta), std::cos(7. * theta), 0),    // 7
      Vector_3(std::sin(8. * theta), std::cos(8. * theta), 0),    // 8
      Vector_3(std::sin(9. * theta), std::cos(9. * theta), 0),    // 9
      Vector_3(std::sin(10. * theta), std::cos(10. * theta), 0),  // 10
      Vector_3(std::sin(11. * theta), std::cos(11. * theta), 0),  // 11
      Vector_3(std::sin(12. * theta), std::cos(12. * theta), 0)   // 12
  };

  // make the third points (o, p0, and pi will be passed to cosAngle, where pi is the ith point in the list)
  Point_3 o(0, 0, 0);
  Point_3 p0(0, 1, 0);
  std::vector<Point_3> points;
  for (const auto& h : hours) points.push_back(p0 + h);

  for (std::size_t i = 0; i < points.size() - 1; i++)
  {
    auto ca = cosAngle(o, p0, points[i]);
    EXPECT_NE(ca, std::nullopt);
    EXPECT_NEAR(*ca, std::cos(i * theta), 1.e-12);
  }
}

TEST(VectorAlgebraUtilities, CanProjectVectorsToPlaneWithZeroNormal)
{
  Vector_3 normal(0, 0, 0);

  std::vector<Vector_3> vectors = {Vector_3(0, 0, 0), Vector_3(1, 0, 0), Vector_3(0, 1, 0), Vector_3(0, 0, 1),
                                   Vector_3(1, 1, 1)};

  for (const auto& v : vectors) EXPECT_EQ(projectVector2Plane(normal, v), v);
}

TEST(VectorAlgebraUtilities, CanProjectVectorsToXYPlane)
{
  Vector_3 normal(0, 0, 1);

  std::vector<Vector_3> vectors = {Vector_3(0, 0, 0),  Vector_3(1, 0, 0),  Vector_3(0, 1, 0),
                                   Vector_3(0, 0, 1),  Vector_3(1, 1, 1),  Vector_3(-1, 0, 0),
                                   Vector_3(0, -1, 0), Vector_3(0, 0, -1), Vector_3(-1, -1, -1)};

  std::vector<Vector_3> expectedProjections = {Vector_3(0, 0, 0),  Vector_3(1, 0, 0), Vector_3(0, 1, 0),
                                               Vector_3(0, 0, 0),  Vector_3(1, 1, 0), Vector_3(-1, 0, 0),
                                               Vector_3(0, -1, 0), Vector_3(0, 0, 0), Vector_3(-1, -1, 0)};

  for (std::size_t i = 0; i < vectors.size(); i++)
    EXPECT_EQ(projectVector2Plane(normal, vectors[i]), expectedProjections[i]);
}

TEST(VectorAlgebraUtilities, CanProjectVectorsToYZPlane)
{
  Vector_3 normal(1, 0, 0);

  std::vector<Vector_3> vectors = {Vector_3(0, 0, 0),  Vector_3(1, 0, 0),  Vector_3(0, 1, 0),
                                   Vector_3(0, 0, 1),  Vector_3(1, 1, 1),  Vector_3(-1, 0, 0),
                                   Vector_3(0, -1, 0), Vector_3(0, 0, -1), Vector_3(-1, -1, -1)};

  std::vector<Vector_3> expectedProjections = {Vector_3(0, 0, 0),  Vector_3(0, 0, 0),  Vector_3(0, 1, 0),
                                               Vector_3(0, 0, 1),  Vector_3(0, 1, 1),  Vector_3(0, 0, 0),
                                               Vector_3(0, -1, 0), Vector_3(0, 0, -1), Vector_3(0, -1, -1)};

  for (std::size_t i = 0; i < vectors.size(); i++)
    EXPECT_EQ(projectVector2Plane(normal, vectors[i]), expectedProjections[i]);
}

TEST(VectorAlgebraUtilities, CanProjectVectorsToZXPlane)
{
  Vector_3 normal(0, 1, 0);

  std::vector<Vector_3> vectors = {Vector_3(0, 0, 0),  Vector_3(1, 0, 0),  Vector_3(0, 1, 0),
                                   Vector_3(0, 0, 1),  Vector_3(1, 1, 1),  Vector_3(-1, 0, 0),
                                   Vector_3(0, -1, 0), Vector_3(0, 0, -1), Vector_3(-1, -1, -1)};

  std::vector<Vector_3> expectedProjections = {Vector_3(0, 0, 0), Vector_3(1, 0, 0),  Vector_3(0, 0, 0),
                                               Vector_3(0, 0, 1), Vector_3(1, 0, 1),  Vector_3(-1, 0, 0),
                                               Vector_3(0, 0, 0), Vector_3(0, 0, -1), Vector_3(-1, 0, -1)};

  for (std::size_t i = 0; i < vectors.size(); i++)
    EXPECT_EQ(projectVector2Plane(normal, vectors[i]), expectedProjections[i]);
}

TEST(VectorAlgebraUtilities, CannotGetRotationForZeroVectors)
{
  EXPECT_EQ(getRotation(Vector_3(0, 0, 0), Vector_3(1, 0, 0)), std::nullopt);
  EXPECT_EQ(getRotation(Vector_3(1, 0, 0), Vector_3(0, 0, 0)), std::nullopt);
  EXPECT_EQ(getRotation(Vector_3(0, 0, 0), Vector_3(0, 0, 0)), std::nullopt);
}

TEST(VectorAlgebraUtilities, CannotGetRotationForParallelVectorsInOppositeDirection)
{
  EXPECT_EQ(getRotation(Vector_3(1, 0, 0), Vector_3(-3, 0, 0)), std::nullopt);
}

TEST(VectorAlgebraUtilities, CanGetRotationForParallelVectors)
{
  EXPECT_EQ(getRotation(Vector_3(1, 0, 0), Vector_3(3, 0, 0)), Affine(1, 0, 0, 0, 1, 0, 0, 0, 1));
  EXPECT_EQ(getRotation(Vector_3(1, 3, 4), Vector_3(0.1, 0.3, 0.4)), Affine(1, 0, 0, 0, 1, 0, 0, 0, 1));
}

TEST(VectorAlgebraUtilities, CannotAlignZeroVectors)
{
  EXPECT_EQ(alignAndTranslate(Vector_3(0, 0, 0), Vector_3(1, 0, 0), Vector_3(1, 1, 1)), std::nullopt);
  EXPECT_EQ(alignAndTranslate(Vector_3(1, 0, 0), Vector_3(0, 0, 0), Vector_3(1, 1, 1)), std::nullopt);
  EXPECT_EQ(alignAndTranslate(Vector_3(0, 0, 0), Vector_3(0, 0, 0), Vector_3(1, 1, 1)), std::nullopt);
}

TEST(VectorAlgebraUtilities, CannotAlignParallelVectorsInOppositeDirection)
{
  EXPECT_EQ(alignAndTranslate(Vector_3(1, 0, 0), Vector_3(-3, 0, 0), Vector_3(1, 1, 1)), std::nullopt);
}

TEST(VectorAlgebraUtilities, CanAlignParallelVectors)
{
  EXPECT_EQ(alignAndTranslate(Vector_3(1, 0, 0), Vector_3(3, 0, 0), Vector_3(1, 2, 3)),
            Affine(1, 0, 0, 1, 0, 1, 0, 2, 0, 0, 1, 3));
  EXPECT_EQ(alignAndTranslate(Vector_3(1, 3, 4), Vector_3(0.1, 0.3, 0.4), Vector_3(10, 20, 30)),
            Affine(1, 0, 0, 10, 0, 1, 0, 20, 0, 0, 1, 30));
}

TEST(VectorAlgebraUtilities, CannotFindIntersectionForFourCollinearPoints)
{
  Point_3 p1(0, 0, 0);
  Point_3 p2(0, 0, 10);
  Point_3 q1(0, 0, 3);
  Point_3 q2(0, 0, 7);
  EXPECT_EQ(foldingIntersection(p1, p2, q1, q2), std::nullopt);
}

TEST(VectorAlgebraUtilities, CannotFindIntersectionForSamePairOfPoints)
{
  Point_3 p1(0, 0, 0);
  Point_3 p2(0, 0, 10);
  Point_3 q1(0, 0, 10);
  Point_3 q2(0, 0, 0);
  EXPECT_EQ(foldingIntersection(p1, p2, q1, q2), std::nullopt);
}

TEST(VectorAlgebraUtilities, CannotFindIntersectionWhenLineIsDegenerate)
{
  Point_3 p1(0, 0, 0);
  Point_3 p2(0, 0, 0);  // p2 = p1 -> line is a point
  Point_3 q1(0, 0, 10);
  Point_3 q2(0, 0, 0);
  EXPECT_EQ(foldingIntersection(p1, p2, q1, q2), std::nullopt);
}

TEST(VectorAlgebraUtilities, CannotFindIntersectionForPositiveLargeFolds)
{
  Point_3 p1(0, 0, 0);
  Point_3 p2(10, 0, 0);
  Point_3 q1(0, 10, 0);
  Point_3 q2(0, 10, 0.05);
  EXPECT_EQ(foldingIntersection(p1, p2, q1, q2), std::nullopt);
}

TEST(VectorAlgebraUtilities, CannotFindIntersectionForNegativeLargeFolds)
{
  Point_3 p1(0, 0, 0);
  Point_3 p2(10, 0, 0);
  Point_3 q1(0, 10, 0);
  Point_3 q2(0, 10, -0.05);
  EXPECT_EQ(foldingIntersection(p1, p2, q1, q2), std::nullopt);
}

TEST(VectorAlgebraUtilities, CanFindIntersectionForThreeCollinearPointsInBound)
{
  Point_3 p1(0, 0, 0);
  Point_3 p2(0, 0, 10);
  Point_3 q1(0, 0, 3);
  Point_3 q2(1, 1, 7);
  auto res = foldingIntersection(p1, p2, q1, q2);
  EXPECT_NE(res, std::nullopt);
  EXPECT_NEAR(*res, 0.3, 1.e-12);
}

TEST(VectorAlgebraUtilities, CanFindIntersectionForThreeCollinearPointsOutOfBound)
{
  Point_3 p1(0, 0, 0);
  Point_3 p2(0, 0, 10);
  Point_3 q1_1(0, 0, 13);  // out of bound on p2 side
  Point_3 q1_2(0, 0, -3);  // out of bound on p1 side
  Point_3 q2(1, 1, 7);
  auto res1 = foldingIntersection(p1, p2, q1_1, q2);
  EXPECT_NE(res1, std::nullopt);
  EXPECT_NEAR(*res1, 1.3, 1.e-12);

  auto res2 = foldingIntersection(p1, p2, q1_2, q2);
  EXPECT_NE(res2, std::nullopt);
  EXPECT_NEAR(*res2, -0.3, 1.e-12);
}

TEST(VectorAlgebraUtilities, CanFindIntersectionZereExtremeOnPlane)
{
  /*
   *  q1
   *  |
   *  |
   *  |
   *  |
   *  p1----p2
   *  |
   *  |
   *  |
   *  |
   *  q2
   */
  Point_3 p1(0, 0, 0);
  Point_3 p2(5, 0, 0);
  Point_3 q1(0, 10, 0);
  Point_3 q2(0, -10, 0);
  auto res = foldingIntersection(p1, p2, q1, q2);
  EXPECT_NE(res, std::nullopt);
  EXPECT_NEAR(*res, 0.0, 1.e-12);
}

TEST(VectorAlgebraUtilities, CanFindIntersectionOneExtremeOnPlane)
{
  /*
   *  q1
   *    \
   *     \
   *      \
   *       \
   *  p1----p2
   *         \
   *          \
   *           \
   *            \
   *            q2
   */
  Point_3 p1(0, 0, 0);
  Point_3 p2(5, 0, 0);
  Point_3 q1(0, 10, 0);
  Point_3 q2(10, -10, 0);
  auto res = foldingIntersection(p1, p2, q1, q2);
  EXPECT_NE(res, std::nullopt);
  EXPECT_NEAR(*res, 1.0, 1.e-12);
}

TEST(VectorAlgebraUtilities, CanFindIntersectionHalfWayPointOnPlane)
{
  /*
   *  q1
   *    \
   *     \
   *      \
   *       \
   *  p1-----------p2
   *         \
   *          \
   *           \
   *            \
   *             q2
   */
  Point_3 p1(0, 0, 0);
  Point_3 p2(10, 0, 0);
  Point_3 q1(0, 10, 0);
  Point_3 q2(10, -10, 0);
  auto res = foldingIntersection(p1, p2, q1, q2);
  EXPECT_NE(res, std::nullopt);
  EXPECT_NEAR(*res, 0.5, 1.e-12);
}

TEST(VectorAlgebraUtilities, CanFindIntersectionHalfWayPoinForFoldedCases)
{
  Point_3 p1(0, -10, 0);
  Point_3 p2(0, 10, 0);
  Point_3 q1(10, 0, 0);

  double rootTwo = std::sqrt(2.);
  std::vector<Point_3> q2s = {
      Point_3(10 / rootTwo, 0, 10 / rootTwo),   Point_3(0, 0, 10),
      Point_3(-10 / rootTwo, 0, 10 / rootTwo),  Point_3(-10, 0, 0),
      Point_3(-10 / rootTwo, 0, -10 / rootTwo), Point_3(0, 0, -10),
      Point_3(10 / rootTwo, 0, -10 / rootTwo),
  };

  for (const auto& q2 : q2s)
  {
    auto res = foldingIntersection(p1, p2, q1, q2);
    EXPECT_NE(res, std::nullopt);
    EXPECT_NEAR(*res, 0.5, 1.e-12);
  }
}
