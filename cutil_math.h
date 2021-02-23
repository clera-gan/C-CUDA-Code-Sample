/* ? 
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */
/*
    This file implements common mathematical operations on vector types
    (double3, double4 etc.) since these are not provided as standard by CUDA.
    The syntax is modelled on the Cg standard library.
    This is part of the CUTIL library and is not supported by NVIDIA.
    Thanks to Linh Hah for additions and fixes.
*/
#ifndef CUTIL_MATH_H
#define CUTIL_MATH_H
#include "cuda_runtime.h"
typedef unsigned int uint;
typedef unsigned short ushort;
#ifndef __CUDACC__
#include <math.h>
////////////////////////////////////////////////////////////////////////////////
// host implementations of CUDA functions
////////////////////////////////////////////////////////////////////////////////
inline double fminf(double a, double b)
{
  return a < b ? a : b;
}
inline double fmaxf(double a, double b)
{
  return a > b ? a : b;
}

inline int max(int a, int b)
{
  return a > b ? a : b;
}

inline int min(int a, int b)
{
  return a < b ? a : b;
}
inline double rsqrtf(double x)
{
    return 1.0f / sqrtf(x);
}
#endif
////////////////////////////////////////////////////////////////////////////////
// constructors
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double2 make_double2(double s)
{
    return make_double2(s, s);
}
inline __host__ __device__ double2 make_double2(double3 a)
{
    return make_double2(a.x, a.y);
}
inline __host__ __device__ double2 make_double2(int2 a)
{
    return make_double2(double(a.x), double(a.y));
}
inline __host__ __device__ double2 make_double2(uint2 a)
{
    return make_double2(double(a.x), double(a.y));
}
inline __host__ __device__ int2 make_int2(int s)
{
    return make_int2(s, s);
}
inline __host__ __device__ int2 make_int2(int3 a)
{
    return make_int2(a.x, a.y);
}
inline __host__ __device__ int2 make_int2(uint2 a)
{
    return make_int2(int(a.x), int(a.y));
}
inline __host__ __device__ int2 make_int2(double2 a)
{
    return make_int2(int(a.x), int(a.y));
}
inline __host__ __device__ uint2 make_uint2(uint s)
{
    return make_uint2(s, s);
}
inline __host__ __device__ uint2 make_uint2(uint3 a)
{
    return make_uint2(a.x, a.y);
}
inline __host__ __device__ uint2 make_uint2(int2 a)
{
    return make_uint2(uint(a.x), uint(a.y));
}
inline __host__ __device__ double3 make_double3(double s)
{
    return make_double3(s, s, s);
}
inline __host__ __device__ double3 make_double3(double2 a)
{
    return make_double3(a.x, a.y, 0.0f);
}
inline __host__ __device__ double3 make_double3(double2 a, double s)
{
    return make_double3(a.x, a.y, s);
}
inline __host__ __device__ double3 make_double3(double4 a)
{
    return make_double3(a.x, a.y, a.z);
}
inline __host__ __device__ double3 make_double3(int3 a)
{
    return make_double3(double(a.x), double(a.y), double(a.z));
}
inline __host__ __device__ double3 make_double3(uint3 a)
{
    return make_double3(double(a.x), double(a.y), double(a.z));
}
inline __host__ __device__ int3 make_int3(int s)
{
    return make_int3(s, s, s);
}
inline __host__ __device__ int3 make_int3(int2 a)
{
    return make_int3(a.x, a.y, 0);
}
inline __host__ __device__ int3 make_int3(int2 a, int s)
{
    return make_int3(a.x, a.y, s);
}
inline __host__ __device__ int3 make_int3(uint3 a)
{
    return make_int3(int(a.x), int(a.y), int(a.z));
}
inline __host__ __device__ int3 make_int3(double3 a)
{
    return make_int3(int(a.x), int(a.y), int(a.z));
}
inline __host__ __device__ uint3 make_uint3(uint s)
{
    return make_uint3(s, s, s);
}
inline __host__ __device__ uint3 make_uint3(uint2 a)
{
    return make_uint3(a.x, a.y, 0);
}
 
inline __host__ __device__ uint3 make_uint3(uint2 a, uint s)
{
    return make_uint3(a.x, a.y, s);
}
inline __host__ __device__ uint3 make_uint3(uint4 a)
{
    return make_uint3(a.x, a.y, a.z);
}
inline __host__ __device__ uint3 make_uint3(int3 a)
{
    return make_uint3(uint(a.x), uint(a.y), uint(a.z));
}
inline __host__ __device__ double4 make_double4(double s)
{
    return make_double4(s, s, s, s);
}
inline __host__ __device__ double4 make_double4(double3 a)
{
    return make_double4(a.x, a.y, a.z, 0.0f);
}
inline __host__ __device__ double4 make_double4(double3 a, double w)
{
    return make_double4(a.x, a.y, a.z, w);
}
inline __host__ __device__ double4 make_double4(int4 a)
{
    return make_double4(double(a.x), double(a.y), double(a.z), double(a.w));
}
inline __host__ __device__ double4 make_double4(uint4 a)
{
    return make_double4(double(a.x), double(a.y), double(a.z), double(a.w));
}
inline __host__ __device__ int4 make_int4(int s)
{
    return make_int4(s, s, s, s);
}
inline __host__ __device__ int4 make_int4(int3 a)
{
    return make_int4(a.x, a.y, a.z, 0);
}
inline __host__ __device__ int4 make_int4(int3 a, int w)
{
    return make_int4(a.x, a.y, a.z, w);
}
inline __host__ __device__ int4 make_int4(uint4 a)
{
    return make_int4(int(a.x), int(a.y), int(a.z), int(a.w));
}
inline __host__ __device__ int4 make_int4(double4 a)
{
    return make_int4(int(a.x), int(a.y), int(a.z), int(a.w));
}
inline __host__ __device__ uint4 make_uint4(uint s)
{
    return make_uint4(s, s, s, s);
}
inline __host__ __device__ uint4 make_uint4(uint3 a)
{
    return make_uint4(a.x, a.y, a.z, 0);
}
inline __host__ __device__ uint4 make_uint4(uint3 a, uint w)
{
    return make_uint4(a.x, a.y, a.z, w);
}
inline __host__ __device__ uint4 make_uint4(int4 a)
{
    return make_uint4(uint(a.x), uint(a.y), uint(a.z), uint(a.w));
}
////////////////////////////////////////////////////////////////////////////////
// negate
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double2 operator-(double2 &a)
{
    return make_double2(-a.x, -a.y);
}
inline __host__ __device__ int2 operator-(int2 &a)
{
    return make_int2(-a.x, -a.y);
}
inline __host__ __device__ double3 operator-(double3 &a)
{
    return make_double3(-a.x, -a.y, -a.z);
}
inline __host__ __device__ int3 operator-(int3 &a)
{
    return make_int3(-a.x, -a.y, -a.z);
}
inline __host__ __device__ double4 operator-(double4 &a)
{
    return make_double4(-a.x, -a.y, -a.z, -a.w);
}
inline __host__ __device__ int4 operator-(int4 &a)
{
    return make_int4(-a.x, -a.y, -a.z, -a.w);
}
////////////////////////////////////////////////////////////////////////////////
// addition
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double2 operator+(double2 a, double2 b)
{
    return make_double2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(double2 &a, double2 b)
{
    a.x += b.x; a.y += b.y;
}
inline __host__ __device__ double2 operator+(double2 a, double b)
{
    return make_double2(a.x + b, a.y + b);
}
inline __host__ __device__ double2 operator+(double b, double2 a)
{
    return make_double2(a.x + b, a.y + b);
}
inline __host__ __device__ void operator+=(double2 &a, double b)
{
    a.x += b; a.y += b;
}
inline __host__ __device__ int2 operator+(int2 a, int2 b)
{
    return make_int2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(int2 &a, int2 b)
{
    a.x += b.x; a.y += b.y;
}
inline __host__ __device__ int2 operator+(int2 a, int b)
{
    return make_int2(a.x + b, a.y + b);
}
inline __host__ __device__ int2 operator+(int b, int2 a)
{
    return make_int2(a.x + b, a.y + b);
}
inline __host__ __device__ void operator+=(int2 &a, int b)
{
    a.x += b; a.y += b;
}
inline __host__ __device__ uint2 operator+(uint2 a, uint2 b)
{
    return make_uint2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(uint2 &a, uint2 b)
{
    a.x += b.x; a.y += b.y;
}
inline __host__ __device__ uint2 operator+(uint2 a, uint b)
{
    return make_uint2(a.x + b, a.y + b);
}
inline __host__ __device__ uint2 operator+(uint b, uint2 a)
{
    return make_uint2(a.x + b, a.y + b);
}
inline __host__ __device__ void operator+=(uint2 &a, uint b)
{
    a.x += b; a.y += b;
}
inline __host__ __device__ double3 operator+(double3 a, double3 b)
{
    return make_double3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ void operator+=(double3 &a, double3 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}
inline __host__ __device__ double3 operator+(double3 a, double b)
{
    return make_double3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ void operator+=(double3 &a, double b)
{
    a.x += b; a.y += b; a.z += b;
}
inline __host__ __device__ int3 operator+(int3 a, int3 b)
{
    return make_int3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ void operator+=(int3 &a, int3 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}
inline __host__ __device__ int3 operator+(int3 a, int b)
{
    return make_int3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ void operator+=(int3 &a, int b)
{
    a.x += b; a.y += b; a.z += b;
}
inline __host__ __device__ uint3 operator+(uint3 a, uint3 b)
{
    return make_uint3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ void operator+=(uint3 &a, uint3 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}
inline __host__ __device__ uint3 operator+(uint3 a, uint b)
{
    return make_uint3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ void operator+=(uint3 &a, uint b)
{
    a.x += b; a.y += b; a.z += b;
}
inline __host__ __device__ int3 operator+(int b, int3 a)
{
    return make_int3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ uint3 operator+(uint b, uint3 a)
{
    return make_uint3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ double3 operator+(double b, double3 a)
{
    return make_double3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ double4 operator+(double4 a, double4 b)
{
    return make_double4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline __host__ __device__ void operator+=(double4 &a, double4 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}
inline __host__ __device__ double4 operator+(double4 a, double b)
{
    return make_double4(a.x + b, a.y + b, a.z + b, a.w + b);
}
inline __host__ __device__ double4 operator+(double b, double4 a)
{
    return make_double4(a.x + b, a.y + b, a.z + b, a.w + b);
}
inline __host__ __device__ void operator+=(double4 &a, double b)
{
    a.x += b; a.y += b; a.z += b; a.w += b;
}
inline __host__ __device__ int4 operator+(int4 a, int4 b)
{
    return make_int4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline __host__ __device__ void operator+=(int4 &a, int4 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}
inline __host__ __device__ int4 operator+(int4 a, int b)
{
    return make_int4(a.x + b, a.y + b, a.z + b,  a.w + b);
}
inline __host__ __device__ int4 operator+(int b, int4 a)
{
    return make_int4(a.x + b, a.y + b, a.z + b,  a.w + b);
}
inline __host__ __device__ void operator+=(int4 &a, int b)
{
    a.x += b; a.y += b; a.z += b; a.w += b;
}
inline __host__ __device__ uint4 operator+(uint4 a, uint4 b)
{
    return make_uint4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline __host__ __device__ void operator+=(uint4 &a, uint4 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}
inline __host__ __device__ uint4 operator+(uint4 a, uint b)
{
    return make_uint4(a.x + b, a.y + b, a.z + b,  a.w + b);
}
inline __host__ __device__ uint4 operator+(uint b, uint4 a)
{
    return make_uint4(a.x + b, a.y + b, a.z + b,  a.w + b);
}
inline __host__ __device__ void operator+=(uint4 &a, uint b)
{
    a.x += b; a.y += b; a.z += b; a.w += b;
}
////////////////////////////////////////////////////////////////////////////////
// subtract
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double2 operator-(double2 a, double2 b)
{
    return make_double2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(double2 &a, double2 b)
{
    a.x -= b.x; a.y -= b.y;
}
inline __host__ __device__ double2 operator-(double2 a, double b)
{
    return make_double2(a.x - b, a.y - b);
}
inline __host__ __device__ double2 operator-(double b, double2 a)
{
    return make_double2(b - a.x, b - a.y);
}
inline __host__ __device__ void operator-=(double2 &a, double b)
{
    a.x -= b; a.y -= b;
}
inline __host__ __device__ int2 operator-(int2 a, int2 b)
{
    return make_int2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(int2 &a, int2 b)
{
    a.x -= b.x; a.y -= b.y;
}
inline __host__ __device__ int2 operator-(int2 a, int b)
{
    return make_int2(a.x - b, a.y - b);
}
inline __host__ __device__ int2 operator-(int b, int2 a)
{
    return make_int2(b - a.x, b - a.y);
}
inline __host__ __device__ void operator-=(int2 &a, int b)
{
    a.x -= b; a.y -= b;
}
inline __host__ __device__ uint2 operator-(uint2 a, uint2 b)
{
    return make_uint2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(uint2 &a, uint2 b)
{
    a.x -= b.x; a.y -= b.y;
}
inline __host__ __device__ uint2 operator-(uint2 a, uint b)
{
    return make_uint2(a.x - b, a.y - b);
}
inline __host__ __device__ uint2 operator-(uint b, uint2 a)
{
    return make_uint2(b - a.x, b - a.y);
}
inline __host__ __device__ void operator-=(uint2 &a, uint b)
{
    a.x -= b; a.y -= b;
}
inline __host__ __device__ double3 operator-(double3 a, double3 b)
{
    return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __host__ __device__ void operator-=(double3 &a, double3 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}
inline __host__ __device__ double3 operator-(double3 a, double b)
{
    return make_double3(a.x - b, a.y - b, a.z - b);
}
inline __host__ __device__ double3 operator-(double b, double3 a)
{
    return make_double3(b - a.x, b - a.y, b - a.z);
}
inline __host__ __device__ void operator-=(double3 &a, double b)
{
    a.x -= b; a.y -= b; a.z -= b;
}
inline __host__ __device__ int3 operator-(int3 a, int3 b)
{
    return make_int3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __host__ __device__ void operator-=(int3 &a, int3 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}
inline __host__ __device__ int3 operator-(int3 a, int b)
{
    return make_int3(a.x - b, a.y - b, a.z - b);
}
inline __host__ __device__ int3 operator-(int b, int3 a)
{
    return make_int3(b - a.x, b - a.y, b - a.z);
}
inline __host__ __device__ void operator-=(int3 &a, int b)
{
    a.x -= b; a.y -= b; a.z -= b;
}
inline __host__ __device__ uint3 operator-(uint3 a, uint3 b)
{
    return make_uint3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __host__ __device__ void operator-=(uint3 &a, uint3 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}
inline __host__ __device__ uint3 operator-(uint3 a, uint b)
{
    return make_uint3(a.x - b, a.y - b, a.z - b);
}
inline __host__ __device__ uint3 operator-(uint b, uint3 a)
{
    return make_uint3(b - a.x, b - a.y, b - a.z);
}
inline __host__ __device__ void operator-=(uint3 &a, uint b)
{
    a.x -= b; a.y -= b; a.z -= b;
}
inline __host__ __device__ double4 operator-(double4 a, double4 b)
{
    return make_double4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline __host__ __device__ void operator-=(double4 &a, double4 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}
inline __host__ __device__ double4 operator-(double4 a, double b)
{
    return make_double4(a.x - b, a.y - b, a.z - b,  a.w - b);
}
inline __host__ __device__ void operator-=(double4 &a, double b)
{
    a.x -= b; a.y -= b; a.z -= b; a.w -= b;
}
inline __host__ __device__ int4 operator-(int4 a, int4 b)
{
    return make_int4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline __host__ __device__ void operator-=(int4 &a, int4 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}
inline __host__ __device__ int4 operator-(int4 a, int b)
{
    return make_int4(a.x - b, a.y - b, a.z - b,  a.w - b);
}
inline __host__ __device__ int4 operator-(int b, int4 a)
{
    return make_int4(b - a.x, b - a.y, b - a.z, b - a.w);
}
inline __host__ __device__ void operator-=(int4 &a, int b)
{
    a.x -= b; a.y -= b; a.z -= b; a.w -= b;
}
inline __host__ __device__ uint4 operator-(uint4 a, uint4 b)
{
    return make_uint4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline __host__ __device__ void operator-=(uint4 &a, uint4 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}
inline __host__ __device__ uint4 operator-(uint4 a, uint b)
{
    return make_uint4(a.x - b, a.y - b, a.z - b,  a.w - b);
}
inline __host__ __device__ uint4 operator-(uint b, uint4 a)
{
    return make_uint4(b - a.x, b - a.y, b - a.z, b - a.w);
}
inline __host__ __device__ void operator-=(uint4 &a, uint b)
{
    a.x -= b; a.y -= b; a.z -= b; a.w -= b;
}
////////////////////////////////////////////////////////////////////////////////
// multiply
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double2 operator*(double2 a, double2 b)
{
    return make_double2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ void operator*=(double2 &a, double2 b)
{
    a.x *= b.x; a.y *= b.y;
}
inline __host__ __device__ double2 operator*(double2 a, double b)
{
    return make_double2(a.x * b, a.y * b);
}
inline __host__ __device__ double2 operator*(double b, double2 a)
{
    return make_double2(b * a.x, b * a.y);
}
inline __host__ __device__ void operator*=(double2 &a, double b)
{
    a.x *= b; a.y *= b;
}
inline __host__ __device__ int2 operator*(int2 a, int2 b)
{
    return make_int2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ void operator*=(int2 &a, int2 b)
{
    a.x *= b.x; a.y *= b.y;
}
inline __host__ __device__ int2 operator*(int2 a, int b)
{
    return make_int2(a.x * b, a.y * b);
}
inline __host__ __device__ int2 operator*(int b, int2 a)
{
    return make_int2(b * a.x, b * a.y);
}
inline __host__ __device__ void operator*=(int2 &a, int b)
{
    a.x *= b; a.y *= b;
}
inline __host__ __device__ uint2 operator*(uint2 a, uint2 b)
{
    return make_uint2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ void operator*=(uint2 &a, uint2 b)
{
    a.x *= b.x; a.y *= b.y;
}
inline __host__ __device__ uint2 operator*(uint2 a, uint b)
{
    return make_uint2(a.x * b, a.y * b);
}
inline __host__ __device__ uint2 operator*(uint b, uint2 a)
{
    return make_uint2(b * a.x, b * a.y);
}
inline __host__ __device__ void operator*=(uint2 &a, uint b)
{
    a.x *= b; a.y *= b;
}
inline __host__ __device__ double3 operator*(double3 a, double3 b)
{
    return make_double3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ void operator*=(double3 &a, double3 b)
{
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
}
inline __host__ __device__ double3 operator*(double3 a, double b)
{
    return make_double3(a.x * b, a.y * b, a.z * b);
}
inline __host__ __device__ double3 operator*(double b, double3 a)
{
    return make_double3(b * a.x, b * a.y, b * a.z);
}
inline __host__ __device__ void operator*=(double3 &a, double b)
{
    a.x *= b; a.y *= b; a.z *= b;
}
inline __host__ __device__ int3 operator*(int3 a, int3 b)
{
    return make_int3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ void operator*=(int3 &a, int3 b)
{
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
}
inline __host__ __device__ int3 operator*(int3 a, int b)
{
    return make_int3(a.x * b, a.y * b, a.z * b);
}
inline __host__ __device__ int3 operator*(int b, int3 a)
{
    return make_int3(b * a.x, b * a.y, b * a.z);
}
inline __host__ __device__ void operator*=(int3 &a, int b)
{
    a.x *= b; a.y *= b; a.z *= b;
}
inline __host__ __device__ uint3 operator*(uint3 a, uint3 b)
{
    return make_uint3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ void operator*=(uint3 &a, uint3 b)
{
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
}
inline __host__ __device__ uint3 operator*(uint3 a, uint b)
{
    return make_uint3(a.x * b, a.y * b, a.z * b);
}
inline __host__ __device__ uint3 operator*(uint b, uint3 a)
{
    return make_uint3(b * a.x, b * a.y, b * a.z);
}
inline __host__ __device__ void operator*=(uint3 &a, uint b)
{
    a.x *= b; a.y *= b; a.z *= b;
}
inline __host__ __device__ double4 operator*(double4 a, double4 b)
{
    return make_double4(a.x * b.x, a.y * b.y, a.z * b.z,  a.w * b.w);
}
inline __host__ __device__ void operator*=(double4 &a, double4 b)
{
    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;
}
inline __host__ __device__ double4 operator*(double4 a, double b)
{
    return make_double4(a.x * b, a.y * b, a.z * b,  a.w * b);
}
inline __host__ __device__ double4 operator*(double b, double4 a)
{
    return make_double4(b * a.x, b * a.y, b * a.z, b * a.w);
}
inline __host__ __device__ void operator*=(double4 &a, double b)
{
    a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}
inline __host__ __device__ int4 operator*(int4 a, int4 b)
{
    return make_int4(a.x * b.x, a.y * b.y, a.z * b.z,  a.w * b.w);
}
inline __host__ __device__ void operator*=(int4 &a, int4 b)
{
    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;
}
inline __host__ __device__ int4 operator*(int4 a, int b)
{
    return make_int4(a.x * b, a.y * b, a.z * b,  a.w * b);
}
inline __host__ __device__ int4 operator*(int b, int4 a)
{
    return make_int4(b * a.x, b * a.y, b * a.z, b * a.w);
}
inline __host__ __device__ void operator*=(int4 &a, int b)
{
    a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}
inline __host__ __device__ uint4 operator*(uint4 a, uint4 b)
{
    return make_uint4(a.x * b.x, a.y * b.y, a.z * b.z,  a.w * b.w);
}
inline __host__ __device__ void operator*=(uint4 &a, uint4 b)
{
    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;
}
inline __host__ __device__ uint4 operator*(uint4 a, uint b)
{
    return make_uint4(a.x * b, a.y * b, a.z * b,  a.w * b);
}
inline __host__ __device__ uint4 operator*(uint b, uint4 a)
{
    return make_uint4(b * a.x, b * a.y, b * a.z, b * a.w);
}
inline __host__ __device__ void operator*=(uint4 &a, uint b)
{
    a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}
////////////////////////////////////////////////////////////////////////////////
// divide
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double2 operator/(double2 a, double2 b)
{
    return make_double2(a.x / b.x, a.y / b.y);
}
inline __host__ __device__ void operator/=(double2 &a, double2 b)
{
    a.x /= b.x; a.y /= b.y;
}
inline __host__ __device__ double2 operator/(double2 a, double b)
{
    return make_double2(a.x / b, a.y / b);
}
inline __host__ __device__ void operator/=(double2 &a, double b)
{
    a.x /= b; a.y /= b;
}
inline __host__ __device__ double2 operator/(double b, double2 a)
{
    return make_double2(b / a.x, b / a.y);
}
inline __host__ __device__ double3 operator/(double3 a, double3 b)
{
    return make_double3(a.x / b.x, a.y / b.y, a.z / b.z);
}
inline __host__ __device__ void operator/=(double3 &a, double3 b)
{
    a.x /= b.x; a.y /= b.y; a.z /= b.z;
}
 
inline __host__ __device__ double3 operator/(double3 a, double b)
{
    return make_double3(a.x / b, a.y / b, a.z / b);
}
inline __host__ __device__ void operator/=(double3 &a, double b)
{
    a.x /= b; a.y /= b; a.z /= b;
}
inline __host__ __device__ double3 operator/(double b, double3 a)
{
    return make_double3(b / a.x, b / a.y, b / a.z);
}
inline __host__ __device__ double4 operator/(double4 a, double4 b)
{
    return make_double4(a.x / b.x, a.y / b.y, a.z / b.z,  a.w / b.w);
}
inline __host__ __device__ void operator/=(double4 &a, double4 b)
{
    a.x /= b.x; a.y /= b.y; a.z /= b.z; a.w /= b.w;
}
inline __host__ __device__ double4 operator/(double4 a, double b)
{
    return make_double4(a.x / b, a.y / b, a.z / b,  a.w / b);
}
inline __host__ __device__ void operator/=(double4 &a, double b)
{
    a.x /= b; a.y /= b; a.z /= b; a.w /= b;
}
inline __host__ __device__ double4 operator/(double b, double4 a){
    return make_double4(b / a.x, b / a.y, b / a.z, b / a.w);
}
////////////////////////////////////////////////////////////////////////////////
// min
////////////////////////////////////////////////////////////////////////////////
inline  __host__ __device__ double2 fminf(double2 a, double2 b)
{
        return make_double2(fminf(a.x,b.x), fminf(a.y,b.y));
}
inline __host__ __device__ double3 fminf(double3 a, double3 b)
{
        return make_double3(fminf(a.x,b.x), fminf(a.y,b.y), fminf(a.z,b.z));
}
inline  __host__ __device__ double4 fminf(double4 a, double4 b)
{
        return make_double4(fminf(a.x,b.x), fminf(a.y,b.y), fminf(a.z,b.z), fminf(a.w,b.w));
}
inline __host__ __device__ int2 min(int2 a, int2 b)
{
    return make_int2(min(a.x,b.x), min(a.y,b.y));
}
inline __host__ __device__ int3 min(int3 a, int3 b)
{
    return make_int3(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z));
}
inline __host__ __device__ int4 min(int4 a, int4 b)
{
    return make_int4(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z), min(a.w,b.w));
}
inline __host__ __device__ uint2 min(uint2 a, uint2 b)
{
    return make_uint2(min(a.x,b.x), min(a.y,b.y));
}
inline __host__ __device__ uint3 min(uint3 a, uint3 b)
{
    return make_uint3(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z));
}
inline __host__ __device__ uint4 min(uint4 a, uint4 b)
{
    return make_uint4(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z), min(a.w,b.w));
}
////////////////////////////////////////////////////////////////////////////////
// max
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double2 fmaxf(double2 a, double2 b)
{
        return make_double2(fmaxf(a.x,b.x), fmaxf(a.y,b.y));
}
inline __host__ __device__ double3 fmaxf(double3 a, double3 b)
{
        return make_double3(fmaxf(a.x,b.x), fmaxf(a.y,b.y), fmaxf(a.z,b.z));
}
inline __host__ __device__ double4 fmaxf(double4 a, double4 b)
{
        return make_double4(fmaxf(a.x,b.x), fmaxf(a.y,b.y), fmaxf(a.z,b.z), fmaxf(a.w,b.w));
}
inline __host__ __device__ int2 max(int2 a, int2 b)
{
    return make_int2(max(a.x,b.x), max(a.y,b.y));
}
inline __host__ __device__ int3 max(int3 a, int3 b)
{
    return make_int3(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z));
}
inline __host__ __device__ int4 max(int4 a, int4 b)
{
    return make_int4(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z), max(a.w,b.w));
}
inline __host__ __device__ uint2 max(uint2 a, uint2 b)
{
    return make_uint2(max(a.x,b.x), max(a.y,b.y));
}
inline __host__ __device__ uint3 max(uint3 a, uint3 b)
{
    return make_uint3(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z));
}
inline __host__ __device__ uint4 max(uint4 a, uint4 b)
{
    return make_uint4(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z), max(a.w,b.w));
}
////////////////////////////////////////////////////////////////////////////////
// lerp
// - linear interpolation between a and b, based on value t in [0, 1] range
////////////////////////////////////////////////////////////////////////////////
inline __device__ __host__ double lerp(double a, double b, double t)
{
    return a + t*(b-a);
}
inline __device__ __host__ double2 lerp(double2 a, double2 b, double t)
{
    return a + t*(b-a);
}
inline __device__ __host__ double3 lerp(double3 a, double3 b, double t)
{
    return a + t*(b-a);
}
inline __device__ __host__ double4 lerp(double4 a, double4 b, double t)
{
    return a + t*(b-a);
}
////////////////////////////////////////////////////////////////////////////////
// clamp
// - clamp the value v to be in the range [a, b]
////////////////////////////////////////////////////////////////////////////////
inline __device__ __host__ double clamp(double f, double a, double b)
{
    return fmaxf(a, fminf(f, b));
}
inline __device__ __host__ int clamp(int f, int a, int b)
{
    return max(a, min(f, b));
}
inline __device__ __host__ uint clamp(uint f, uint a, uint b)
{
    return max(a, min(f, b));
}
inline __device__ __host__ double2 clamp(double2 v, double a, double b)
{
    return make_double2(clamp(v.x, a, b), clamp(v.y, a, b));
}
inline __device__ __host__ double2 clamp(double2 v, double2 a, double2 b)
{
    return make_double2(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y));
}
inline __device__ __host__ double3 clamp(double3 v, double a, double b)
{
    return make_double3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}
inline __device__ __host__ double3 clamp(double3 v, double3 a, double3 b)
{
    return make_double3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}
inline __device__ __host__ double4 clamp(double4 v, double a, double b)
{
    return make_double4(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b), clamp(v.w, a, b));
}
inline __device__ __host__ double4 clamp(double4 v, double4 a, double4 b)
{
    return make_double4(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z), clamp(v.w, a.w, b.w));
}
inline __device__ __host__ int2 clamp(int2 v, int a, int b)
{
    return make_int2(clamp(v.x, a, b), clamp(v.y, a, b));
}
inline __device__ __host__ int2 clamp(int2 v, int2 a, int2 b)
{
    return make_int2(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y));
}
inline __device__ __host__ int3 clamp(int3 v, int a, int b)
{
    return make_int3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}
inline __device__ __host__ int3 clamp(int3 v, int3 a, int3 b)
{
    return make_int3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}
inline __device__ __host__ int4 clamp(int4 v, int a, int b)
{
    return make_int4(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b), clamp(v.w, a, b));
}
inline __device__ __host__ int4 clamp(int4 v, int4 a, int4 b)
{
    return make_int4(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z), clamp(v.w, a.w, b.w));
}
inline __device__ __host__ uint2 clamp(uint2 v, uint a, uint b)
{
    return make_uint2(clamp(v.x, a, b), clamp(v.y, a, b));
}
inline __device__ __host__ uint2 clamp(uint2 v, uint2 a, uint2 b)
{
    return make_uint2(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y));
}
inline __device__ __host__ uint3 clamp(uint3 v, uint a, uint b)
{
    return make_uint3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}
inline __device__ __host__ uint3 clamp(uint3 v, uint3 a, uint3 b)
{
    return make_uint3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}
inline __device__ __host__ uint4 clamp(uint4 v, uint a, uint b)
{
    return make_uint4(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b), clamp(v.w, a, b));
}
inline __device__ __host__ uint4 clamp(uint4 v, uint4 a, uint4 b)
{
    return make_uint4(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z), clamp(v.w, a.w, b.w));
}
////////////////////////////////////////////////////////////////////////////////
// dot product
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double dot(double2 a, double2 b)
{ 
    return a.x * b.x + a.y * b.y;
}
inline __host__ __device__ double dot(double3 a, double3 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline __host__ __device__ double dot(double4 a, double4 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
inline __host__ __device__ int dot(int2 a, int2 b)
{ 
    return a.x * b.x + a.y * b.y;
}
inline __host__ __device__ int dot(int3 a, int3 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline __host__ __device__ int dot(int4 a, int4 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
inline __host__ __device__ uint dot(uint2 a, uint2 b)
{ 
    return a.x * b.x + a.y * b.y;
}
inline __host__ __device__ uint dot(uint3 a, uint3 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline __host__ __device__ uint dot(uint4 a, uint4 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
////////////////////////////////////////////////////////////////////////////////
// length
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double length(double2 v)
{
    return sqrtf(dot(v, v));
}
inline __host__ __device__ double length(double3 v)
{
    return sqrtf(dot(v, v));
}
inline __host__ __device__ double length(double4 v)
{
    return sqrtf(dot(v, v));
}
////////////////////////////////////////////////////////////////////////////////
// normalize
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double2 normalize(double2 v)
{
    double invLen = rsqrtf(dot(v, v));
    return v * invLen;
}
inline __host__ __device__ double3 normalize(double3 v)
{
    double invLen = rsqrtf(dot(v, v));
    return v * invLen;
}
inline __host__ __device__ double4 normalize(double4 v)
{
    double invLen = rsqrtf(dot(v, v));
    return v * invLen;
}
////////////////////////////////////////////////////////////////////////////////
// floor
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double2 floorf(double2 v)
{
    return make_double2(floorf(v.x), floorf(v.y));
}
inline __host__ __device__ double3 floorf(double3 v)
{
    return make_double3(floorf(v.x), floorf(v.y), floorf(v.z));
}
inline __host__ __device__ double4 floorf(double4 v)
{
    return make_double4(floorf(v.x), floorf(v.y), floorf(v.z), floorf(v.w));
}
////////////////////////////////////////////////////////////////////////////////
// frac - returns the fractional portion of a scalar or each vector component
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double fracf(double v)
{
    return v - floorf(v);
}
inline __host__ __device__ double2 fracf(double2 v)
{
    return make_double2(fracf(v.x), fracf(v.y));
}
inline __host__ __device__ double3 fracf(double3 v)
{
    return make_double3(fracf(v.x), fracf(v.y), fracf(v.z));
}
inline __host__ __device__ double4 fracf(double4 v)
{
    return make_double4(fracf(v.x), fracf(v.y), fracf(v.z), fracf(v.w));
}
////////////////////////////////////////////////////////////////////////////////
// fmod
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double2 fmodf(double2 a, double2 b)
{
    return make_double2(fmodf(a.x, b.x), fmodf(a.y, b.y));
}
inline __host__ __device__ double3 fmodf(double3 a, double3 b)
{
    return make_double3(fmodf(a.x, b.x), fmodf(a.y, b.y), fmodf(a.z, b.z));
}
inline __host__ __device__ double4 fmodf(double4 a, double4 b)
{
    return make_double4(fmodf(a.x, b.x), fmodf(a.y, b.y), fmodf(a.z, b.z), fmodf(a.w, b.w));
}
////////////////////////////////////////////////////////////////////////////////
// absolute value
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double2 fabs(double2 v)
{
        return make_double2(fabs(v.x), fabs(v.y));
}
inline __host__ __device__ double3 fabs(double3 v)
{
        return make_double3(fabs(v.x), fabs(v.y), fabs(v.z));
}
inline __host__ __device__ double4 fabs(double4 v)
{
        return make_double4(fabs(v.x), fabs(v.y), fabs(v.z), fabs(v.w));
}
inline __host__ __device__ int2 abs(int2 v)
{
        return make_int2(abs(v.x), abs(v.y));
}
inline __host__ __device__ int3 abs(int3 v)
{
        return make_int3(abs(v.x), abs(v.y), abs(v.z));
}
inline __host__ __device__ int4 abs(int4 v)
{
        return make_int4(abs(v.x), abs(v.y), abs(v.z), abs(v.w));
}
////////////////////////////////////////////////////////////////////////////////
// reflect
// - returns reflection of incident ray I around surface normal N
// - N should be normalized, reflected vector's length is equal to length of I
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double3 reflect(double3 i, double3 n)
{
        return i - 2.0f * n * dot(n,i);
}
////////////////////////////////////////////////////////////////////////////////
// cross product
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ double3 cross(double3 a, double3 b)
{ 
    return make_double3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); 
}
////////////////////////////////////////////////////////////////////////////////
// smoothstep
// - returns 0 if x < a
// - returns 1 if x > b
// - otherwise returns smooth interpolation between 0 and 1 based on x
////////////////////////////////////////////////////////////////////////////////
inline __device__ __host__ double smoothstep(double a, double b, double x)
{
        double y = clamp((x - a) / (b - a), 0.0f, 1.0f);
        return (y*y*(3.0f - (2.0f*y)));
}
inline __device__ __host__ double2 smoothstep(double2 a, double2 b, double2 x)
{
        double2 y = clamp((x - a) / (b - a), 0.0f, 1.0f);
        return (y*y*(make_double2(3.0f) - (make_double2(2.0f)*y)));
}
 
inline __device__ __host__ double3 smoothstep(double3 a, double3 b, double3 x)
{
        double3 y = clamp((x - a) / (b - a), 0.0f, 1.0f);
        return (y*y*(make_double3(3.0f) - (make_double3(2.0f)*y)));
}
inline __device__ __host__ double4 smoothstep(double4 a, double4 b, double4 x)
{
        double4 y = clamp((x - a) / (b - a), 0.0f, 1.0f);
        return (y*y*(make_double4(3.0f) - (make_double4(2.0f)*y)));
}
#endif 
