#ifndef DIRECT_EB_FIELD_VECTOR_H
#define DIRECT_EB_FIELD_VECTOR_H

#include <cmath>

namespace fel
{
  typedef double Double;

  template<typename T>
  class FieldVector
  {
  public:
    FieldVector() : values_{T(), T(), T()} {}
    explicit FieldVector(const T& value)
      : values_{value, value, value} {}

    template<typename U>
    FieldVector(const FieldVector<U>& other)
      : values_{static_cast<T>(other[0]),
                static_cast<T>(other[1]),
                static_cast<T>(other[2])} {}

    const T& operator[](unsigned int component) const
    {
      return values_[component];
    }

    T& operator[](unsigned int component)
    {
      return values_[component];
    }

    T norm2() const
    {
      return values_[0] * values_[0] +
             values_[1] * values_[1] +
             values_[2] * values_[2];
    }

    T norm() const
    {
      using std::sqrt;
      return sqrt(norm2());
    }

    template<typename Scalar, typename U>
    void mv(const Scalar& scalar, const FieldVector<U>& vector)
    {
      for (unsigned int i = 0; i < 3; ++i)
        values_[i] = static_cast<T>(scalar * vector[i]);
    }

    template<typename Scalar, typename U>
    void pmv(const Scalar& scalar, const FieldVector<U>& vector)
    {
      for (unsigned int i = 0; i < 3; ++i)
        values_[i] += static_cast<T>(scalar * vector[i]);
    }

    template<typename U>
    FieldVector& operator=(const FieldVector<U>& other)
    {
      for (unsigned int i = 0; i < 3; ++i)
        values_[i] = static_cast<T>(other[i]);
      return *this;
    }

  private:
    T values_[3];
  };

  template<typename T, typename Scalar>
  FieldVector<T>& operator*=(FieldVector<T>& vector, const Scalar& scalar)
  {
    for (unsigned int i = 0; i < 3; ++i) vector[i] *= scalar;
    return vector;
  }

  template<typename T, typename Scalar>
  FieldVector<T>& operator/=(FieldVector<T>& vector, const Scalar& scalar)
  {
    for (unsigned int i = 0; i < 3; ++i) vector[i] /= scalar;
    return vector;
  }

  template<typename T, typename U>
  FieldVector<T>& operator+=(FieldVector<T>& left,
                             const FieldVector<U>& right)
  {
    for (unsigned int i = 0; i < 3; ++i) left[i] += right[i];
    return left;
  }

  template<typename T, typename U>
  FieldVector<T>& operator-=(FieldVector<T>& left,
                             const FieldVector<U>& right)
  {
    for (unsigned int i = 0; i < 3; ++i) left[i] -= right[i];
    return left;
  }

  template<typename T, typename U>
  T operator*(const FieldVector<T>& left, const FieldVector<U>& right)
  {
    return left[0] * right[0] + left[1] * right[1] +
           left[2] * right[2];
  }

  template<typename T, typename U>
  FieldVector<T> cross(const FieldVector<T>& left,
                       const FieldVector<U>& right)
  {
    FieldVector<T> result(static_cast<T>(0));
    result[0] = left[1] * right[2] - left[2] * right[1];
    result[1] = left[2] * right[0] - left[0] * right[2];
    result[2] = left[0] * right[1] - left[1] * right[0];
    return result;
  }
}

#endif
