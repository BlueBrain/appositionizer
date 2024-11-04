#pragma once

#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

using t_coordinate = float;

constexpr t_coordinate operator"" _c(long double in) {
    return static_cast<t_coordinate>(in);
}

template <typename T = float>
struct Coord3D {
    T x, y, z;

    using value_type = T;

    inline Coord3D() = default;
    inline Coord3D(T _x, T _y, T _z)
        : x(_x)
        , y(_y)
        , z(_z) {}

    inline Coord3D(const Coord3D&) = default;
    inline Coord3D(Coord3D&&) = default;
    inline Coord3D& operator=(const Coord3D&) = default;
    inline Coord3D& operator=(Coord3D&&) = default;

    inline T& operator[](unsigned int i) {
        assert(i < 3);
        if (i == 0)
            return x;
        if (i == 1)
            return y;
        if (i == 2)
            return z;
        throw std::out_of_range("Index out of range for Coord3D object");
    }

    inline T operator[](unsigned int i) const {
        assert(i < 3);
        if (i == 0)
            return x;
        if (i == 1)
            return y;
        if (i == 2)
            return z;
        throw std::out_of_range("Index out of range for Coord3D object");
    }

    inline bool operator==(const Coord3D& other) const {
        return isStrictEqual(other);
    }

    inline bool operator!=(const Coord3D& other) const {
        return !(*this == other);
    }

    inline Coord3D& operator+=(const Coord3D& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return (*this);
    }

    inline Coord3D& operator-=(const Coord3D& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return (*this);
    }

    inline Coord3D& operator*=(const T n) {
        x *= n;
        y *= n;
        z *= n;
        return (*this);
    }

    inline Coord3D operator+(const Coord3D& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    inline Coord3D operator-(const Coord3D& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    inline Coord3D operator*(T n) const {
        return {x * n, y * n, z * n};
    }

    inline Coord3D operator/(T n) const {
        return {x / n, y / n, z / n};
    }

    inline T operator*(const Coord3D& o) const {
        return x * o.x + y * o.y + z * o.z;
    }

    inline T angle(const Coord3D& o) const {
        return acos((*this * o) / (norm() * o.norm()));
    }

    inline T r2() const {
        return x * x + y * y;
    }

    /**
     * cross product of two vector.
     * (source: http://en.wikipedia.org/wiki/Cross_product#Coordinate_notation)
     * \param u left hand vector
     * \param v right hand vector
     * \return
     */
    inline static Coord3D cross(const Coord3D& u, const Coord3D& v) {
        return {u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
    }

    /**
     * norm of this vector.
     * (source: http://en.wikipedia.org/wiki/Norm_(mathematics)#Euclidean_norm)
     * \return
     */
    inline T norm() const {
        return std::sqrt(norm2());
    }
    inline T norm2() const {
        return x * x + y * y + z * z;
    }

    inline bool isAlmostEqual(const Coord3D& other,
                              T epsilon = std::numeric_limits<T>::epsilon() * 10) const {
        return almost_equal(x, other.x, epsilon) && almost_equal(y, other.y, epsilon) &&
               almost_equal(z, other.z, epsilon);
    }

    inline bool isStrictEqual(const Coord3D& other) const {
        return (x == other.x) && (y == other.y) && (z == other.z);
    }

    inline static T dot(const Coord3D& u, const Coord3D& v) {
        return u.x * v.x + u.y * v.y + u.z * v.z;
    }

    inline static T distance(const Coord3D& u, const Coord3D& v) {
        return std::sqrt((u.x - v.x) * (u.x - v.x) + (u.y - v.y) * (u.y - v.y) +
                         (u.z - v.z) * (u.z - v.z));
    }
};

using Point = Coord3D<t_coordinate>;

static_assert(std::is_trivial<Point>::value, "Point should be trivial");
static_assert(std::is_trivially_copyable<Point>::value, "Point should be trivially copyable");

/**
 * Represents a bounding box as back left bottom, and front right top corner
 */
template <typename T /*= t_coordinate*/>
class Box {
  public:
    Box()
        : min()
        , max() {}

    template <typename S>
    Box(const S& p1, const S& p2) {
        setBounds(p1, p2);
    }

    template <typename S>
    Box(const Box<S>& b)
        : min(b.min)
        , max(b.max) {}

    template <typename S>
    void setBounds(const S& p1, const S& p2) {
        for (size_t i = 0; i < 3; ++i) {
            min[i] = static_cast<T>(std::min(p1[i], p2[i]));
            max[i] = static_cast<T>(std::max(p1[i], p2[i]));
        }
    }

    Coord3D<T> center() const {
        Coord3D<T> p;
        for (size_t i = 0; i < 3; ++i) {
            p[i] = 0.5 * (min[i] + max[i]);
        }
        return p;
    }

    Box& expand(const T r) {
        for (size_t i = 0; i < 3; ++i) {
            min[i] -= r;
            max[i] += r;
        }
        return *this;
    }

    Box& merge(const Box<T>& b) {
        for (size_t i = 0; i < 3; ++i) {
            min[i] = std::min(min[i], b.min[i]);
            max[i] = std::max(max[i], b.max[i]);
        }
        return *this;
    }

    bool overlaps(const Box& o) const {
        for (size_t i = 0; i < 3; ++i) {
            if (not overlap(min[i], max[i], o.min[i], o.max[i]))
                return false;
        }
        return true;
    }

    template <typename S>
    bool contains(const S& p) const {
        for (size_t i = 0; i < 3; ++i) {
            if (p[i] < min[i] or p[i] > max[i])
                return false;
        }
        return true;
    }

    Coord3D<T> min, max;

  private:
    static bool overlap(T p1, T p2, T q1, T q2) {
        return p1 < q2 and q1 < p2;
    }
};
