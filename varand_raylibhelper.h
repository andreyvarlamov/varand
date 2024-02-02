#ifndef VARAND_RAYLIBHELPER_H
#define VARAND_RAYLIBHELPER_H

#include <varand/varand_types.h>

#include <raylib/raylib.h>
// #include <raylib/raymath.h>

inline Vector2 GetVector2(f32 x, f32 y)
{
    Vector2 result;

    result.x = x;
    result.y = y;

    return result;
}

inline Vector2 GetVector2(f32 value)
{
    Vector2 result = GetVector2(value, value);
    return result;
}

inline Rectangle GetRectangle(f32 x, f32 y, f32 width, f32 height)
{
    Rectangle result;

    result.x = x;
    result.y = y;
    result.width = width;
    result.height = height;

    return result;
}

inline Rectangle GetRectangle(f32 x, f32 y, f32 dim)
{
    Rectangle result = GetRectangle(x, y, dim, dim);
    return result;
}

inline Vector2 GetRectangleMin(Rectangle rectangle)
{
    Vector2 result = GetVector2(rectangle.x, rectangle.y);
    return result;
}

inline Vector3 GetVector3(f32 x, f32 y, f32 z)
{
    Vector3 result;

    result.x = x;
    result.y = y;
    result.z = z;

    return result;
}

inline Color GetColor(int r, int g, int b, int a)
{
    Color result;

    result.r = (u8) r;
    result.g = (u8) g;
    result.b = (u8) b;
    result.a = (u8) a;

    return result;
}

inline Color GetColor(int r, int g, int b)
{
    Color result = GetColor(r, g, b, 255);
    return result;
}

inline Vector3 GetVector3(f32 value)
{
    Vector3 result = GetVector3(value, value, value);
    return result;
}

inline Vector3 GetVector3()
{
    Vector3 result = GetVector3(0.0f);
    return result;
}

inline Vector3 GetVector3(f32 x, f32 y)
{
    Vector3 result = GetVector3(x, y, 0);
    return result;
}

inline Vector2 GetVector2()
{
    Vector2 result = GetVector2(0);
    return result;
}

inline BoundingBox GetBoundingBox(Vector3 min, Vector3 max)
{
    BoundingBox result;

    result.min = min;
    result.max = max;

    return result;
}

#endif
