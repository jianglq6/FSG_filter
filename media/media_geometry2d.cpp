/***************************************************************************
 *
 * Geometries funcs for medium parameterization
 *
 * Authors: Luqian Jiang <jianglq@mail.ustc.edu.cn>
 *          Wei Zhang <zhangwei@sustech.edu.cn>
 *
 * Copyright (c) 2021 zwlab
 *
 * ChangeLog: 
 *    09/2021: Created by Luqian Jiang 
 *
 ***************************************************************************/

#include <iostream>
#include <vector>
#include <cfloat>
#include <cmath>
#include "media_geometry2d.hpp"
//using namespace std;

// ignore 5 and 10
int edgeTable[16] = {
    0b0000, 0b1001, 0b0011, 0b1010,  
    0b0110, 0b0000, 0b0101, 0b1100,
    0b1100, 0b0101, 0b0000, 0b0110,
    0b1010, 0b0011, 0b1001, 0b0000 };


const double eps = 1e-10;
int dcmp(double x) {
	if (fabs(x) < EPS) 
		return 0;
	else
		return x < 0?-1:1;
}

Vector2 operator + (const Vector2 &A, const Vector2 &B) {
	return Vector2(A.x+B.x, A.z+B.z);
}

Vector2 operator - (const Vector2 &A, const Vector2 &B) {
	return Vector2(A.x-B.x, A.z-B.z);
}

Vector2 operator * (const Vector2 &A, double p) {
	return Vector2(A.x*p, A.z*p);
}

Vector2 operator / (const Vector2 &A, double p) {
	return Vector2(A.x/p, A.z/p);
}

bool operator == (const Point2 &a, const Point2 &b) {
    return dcmp(a.x-b.x) == 0 && (dcmp(a.z-b.z) == 0);
}

double Dot(const Vector2 &A, const Vector2 &B) {
	return A.x*B.x + A.z*B.z;
}

double Length(const Vector2 &A) {
	return sqrt(Dot(A, A));
}

double Angle(const Vector2 &A, const Vector2 &B) {
	return acos(Dot(A, B)/Length(A)/Length(B));
}

double Cross(const Vector2 &A, const Vector2 &B) {
	return A.x*B.z - A.z*B.x;
}

double Area2(const Point2 &A, const Point2 &B, const Point2 &C) {
	return Length(Cross(B-A, C-A));
}

double DistanceP2P(const Point2 &A, const Point2 &B) {
	return Length(B-A);	
}

// A must not be a zero vector
Vector2 Normal(const Vector2 &A) {
	double L = Length(A);
	return Vector2(-A.z/L, A.x/L);
}

// point P to line AB
double DistanceToLine(const Point2 &P, const Point2 &A, const Point2 &B) {
    Vector2 v1 = B-A, v2 = P-A;
    return fabs(Cross(v1, v2)) / Length(v1);
}

double DistanceToSegment(const Point2 P, const Point2 &A, const Point2 &B) {
    if (A == B) return Length(P-A);
    Vector2 v1 = B-A, v2 = P-A, v3 = P-B;
    if (dcmp(Dot(v1, v2)) < 0) 
        return Length(v2);
    else if (dcmp(Dot(v1, v3)) > 0)
        return Length(v3);
    else
        return fabs(Cross(v1, v2))/Length(v1);
}

// If Seg(a1, a2) and Seg(b1, b2) intersect 
bool SegmentProperIntersection(
	const Point2 &a1, const Point2 &a2, 
	const Point2 &b1, const Point2 &b2)
{
	double c1 = Cross(a2-a1,b1-a1), c2 = Cross(a2-a1,b2-a1);
	double c3 = Cross(b2-b1,a1-b1), c4 = Cross(b2-b1,a2-b1);
	return dcmp(c1)*dcmp(c2)<0 && dcmp(c3)*dcmp(c4)<0;
}

// Find the intersection of two lines, 
// first make sure that the two lines have a unique intersection
Point2 GetLineIntersection(
	const Point2 &P, const Vector2 &v, 
	const Point2 &Q, const Vector2 &w)
{
	Vector2 u = P-Q;
	double t = Cross(w, u)/Cross(v, w);
	return P+v*t;
}

bool isPointOnSegment(const Point2 &p, const Point2 &a1, const Point2 &a2) {
	return dcmp(Cross(a1-p, a2-p)) == 0 && dcmp(Dot(a1-p, a2-p)) < 0;
}
