#ifndef _MEDIA_GEOMETRY2D_H_
#define _MEDIA_GEOMETRY2D_H_

#include <iostream>
#include <vector>
#include <cfloat>
#include <cassert>

#define EPS 1e-10

struct Point2{
	double x, z;
	Point2(double x = 0, double z = 0):x(x), z(z) {}
    // just used for set, no actual meaning
    bool operator < (const Point2 &A) const {
        return x < A.x;
    }
};

typedef Point2 Vector2;

int dcmp(double x);

Vector2 operator + (const Vector2 &A, const Vector2 &B);
Vector2 operator - (const Vector2 &A, const Vector2 &B);
Vector2 operator * (const Vector2 &A, double p);
Vector2 operator / (const Vector2 &A, double p); 

bool operator == (const Point2 &a, const Point2 &b);

double Dot(const Vector2 &A, const Vector2 &B);
double Length(const Vector2 &A);
double Angle(const Vector2 &A, const Vector2 &B);
double Cross(const Vector2 &A, const Vector2 &B);
double Area2(const Point2 &A, const Point2 &B, const Point2 &C);

double DistanceP2P(const Point2 &A, const Point2 &B);

Vector2 Normal(const Vector2 &A);

double DistanceToLine(const Point2 &P, const Point2 &A, const Point2 &B);
double DistanceToSegment(const Point2 P, const Point2 &A, const Point2 &B);

bool isPointOnSegment(const Point2 &p, const Point2 &a1, const Point2 &a2);


bool SegmentProperIntersection(
	const Point2 &a1, const Point2 &a2, 
	const Point2 &b1, const Point2 &b2);

Point2 GetLineIntersection(
	const Point2 &P, const Vector2 &v, 
	const Point2 &Q, const Vector2 &w);

/*
 *    ↑ +z  
 *    |        2
 *         3-------2 
 *         |       |
 *        3|       |1      
 *         |       |
 *         0-------1
 *             0
 */
struct Mesh2 {
    Point2 v[4]; 
    Mesh2 (Point2 A, Point2 B, Point2 C, Point2 D) {
        this->v[0] = A;
        this->v[1] = B;
        this->v[2] = C;
        this->v[3] = D;
    }
};

#endif /*media_geometry*/
