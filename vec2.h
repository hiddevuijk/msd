#ifndef GUARD_VEC2_H
#define GUARD_VEC2_H

class Vec2 {
public:
  Vec2() : x(), y() {}
  Vec2(double x, double y) : x(x), y(y) {}

  Vec2& operator+=(const Vec2& vRHS)
  { x += vRHS.x;
    y += vRHS.y;
    return *this;
  }

  Vec2& operator-=(const Vec2& vRHS)
  { x -= vRHS.x;
    y -= vRHS.y;
    return *this;
  }

  Vec2& operator*=(const double& RHS)
  { x *= RHS;
    y *= RHS;
    return *this;
  }

  Vec2& operator/=(const double& RHS)
  { x /= RHS;
    y /= RHS;
    return *this;
  }

  Vec2& operator+=(const double& RHS)
  { x += RHS;
    y += RHS;
    return *this;
  }

  Vec2& operator-=(const double& RHS)
  { x -= RHS;
    y -= RHS;
    return *this;
  }

  double x,y;
};


Vec2 operator+(const Vec2& v1, const Vec2& v2)
{ return Vec2(v1.x + v2.x, v1.y + v2.y); }

Vec2 operator-(const Vec2& v1, const Vec2& v2)
{ return Vec2(v1.x - v2.x, v1.y - v2.y); }

Vec2 operator*(const Vec2& v1, const Vec2& v2)
{ return Vec2(v1.x * v2.x, v1.y * v2.y); }
Vec2 operator*(const double& a, const Vec2& v)
{ return Vec2(a*v.x, a*v.y); }
Vec2 operator*(const Vec2& v, const double& a)
{ return Vec2(a*v.x, a*v.y); }


Vec2 operator/(const Vec2& v1, const Vec2& v2)
{ return Vec2(v1.x / v2.x, v1.y / v2.y); }
Vec2 operator/(const Vec2& v, const double& a)
{ return Vec2(v.x / a, v.y /a); }



#endif
