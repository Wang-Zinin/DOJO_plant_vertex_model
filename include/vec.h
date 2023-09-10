/*****************/
// Mathematical Operations for vectors
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#ifndef VEC_H
#define VEC_H
#include <stdio.h>
#include <math.h>

//using namespace std;

template <class TEMPLATE> //template refers to the type
class _vec{
  public:
	TEMPLATE x;
	TEMPLATE y;
	TEMPLATE z;
	//�����o�֐�
	//�R���X�g���N�^
	_vec(void);
	_vec(const TEMPLATE&, const TEMPLATE&, const TEMPLATE& );
	_vec(const _vec<TEMPLATE>&);
	
	//�I�y���[�^
	void IN(const TEMPLATE& , const TEMPLATE&, const TEMPLATE& );
	void IN(const _vec<TEMPLATE>&);
	_vec<int> icast(void) const;
	_vec<TEMPLATE>& operator = (const _vec<TEMPLATE>&);
	_vec<TEMPLATE>& operator += (const _vec<TEMPLATE>&);
	_vec<TEMPLATE>& operator -= (const _vec<TEMPLATE>&);
	_vec<TEMPLATE>& operator *= (const _vec<TEMPLATE>&);
	_vec<TEMPLATE>& operator /= (const _vec<TEMPLATE>&);
	_vec<TEMPLATE>& operator *= (const TEMPLATE&);
	_vec<TEMPLATE>& operator /= (const TEMPLATE&);
//	_vec<TEMPLATE>& rot( const TEMPLATE&);
	TEMPLATE norm(void) const;
	TEMPLATE sqr(void) const;
};




//�����o�֐�����
//�R���X�g���N�^
template <class TEMPLATE> inline _vec<TEMPLATE>::_vec(void){
	x = 0;
	y = 0;
	z = 0;
}

template <class TEMPLATE> inline _vec<TEMPLATE>::_vec(const TEMPLATE& a, const TEMPLATE& b, const TEMPLATE& c){
	x = a;
	y = b;
	z = c;
}

template <class TEMPLATE> inline _vec<TEMPLATE>::_vec(const _vec<TEMPLATE>& v){
	x = v.x;
	y = v.y;
	z = v.z;
}

//�I�y���[�^
template <class TEMPLATE> inline void _vec<TEMPLATE>::IN(const TEMPLATE& a, const TEMPLATE& b, const TEMPLATE& c){
	x = a; y = b; z = c;
}

template <class TEMPLATE> inline void _vec<TEMPLATE>::IN(const _vec<TEMPLATE>& v){
	x = v.x; y = v.y; z = v.z;
}

template<class TEMPLATE> inline _vec<int> _vec<TEMPLATE>::icast(void) const
{
  return _vec<int>((int)x, (int)y, (int)z);
}

template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::operator =  ( const _vec<TEMPLATE>& v ){
  x = v.x;    y = v.y; z = v.z;
  return *this;
}

template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::operator +=  ( const _vec<TEMPLATE>& v ){
  x += v.x;    y += v.y; z += v.z;
  return *this;
}

template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::operator -=  ( const _vec<TEMPLATE>& v ){
  x -= v.x;    y -= v.y; z -= v.z;
  return *this;
}

template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::operator *= ( const TEMPLATE& a){
	x *= a; y *= a; z *= a;
	return *this;
}

template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::operator /= ( const TEMPLATE& a){
	x /= a; y /= a; z /= a;
	return *this;
}

template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::operator *= ( const _vec<TEMPLATE>& v){
	x *= v.x; y *= v.y; z *= v.z;
	return *this;
}

template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::operator /= ( const _vec<TEMPLATE>& v){
	x /= v.x; y /= v.y; z /= v.z;
	return *this;
}

//�X�J���΃x�N�g��
template <class TEMPLATE> inline _vec<TEMPLATE> operator*( const TEMPLATE& a, const _vec<TEMPLATE> v  ){
	_vec<TEMPLATE> vec(a*v.x, a*v.y, a*v.z);
	return vec;
}

template <class TEMPLATE> inline _vec<TEMPLATE> operator*( const _vec<TEMPLATE> v, const TEMPLATE& a  ){
	_vec<TEMPLATE> vec(v.x*a, v.y*a, v.z*a);
	return vec;
}


template <class TEMPLATE> inline _vec<TEMPLATE> operator/( const TEMPLATE& a, const _vec<TEMPLATE> v  ){
	_vec<TEMPLATE> vec(a/v.x, a/v.y, a/v.z);
	
	return vec;
}

template <class TEMPLATE> inline _vec<TEMPLATE> operator/(const _vec<TEMPLATE> v, const TEMPLATE& a  ){
	_vec<TEMPLATE> vec(v.x/a,v.y/a, v.z/a);
	return vec;
}

template <class TEMPLATE> inline _vec<TEMPLATE> operator+(const _vec<TEMPLATE> v, const TEMPLATE& a  ){
	_vec<TEMPLATE> vec(v.x+a,v.y+a,v.z+a);
	return vec;
}

template <class TEMPLATE> inline _vec<TEMPLATE> operator+(const TEMPLATE& a ,const _vec<TEMPLATE> v ){
	_vec<TEMPLATE> vec(a+v.x,a+v.y,a+v.z);
	return vec;
}

template <class TEMPLATE> inline _vec<TEMPLATE> operator-(const _vec<TEMPLATE> v, const TEMPLATE& a  ){
	_vec<TEMPLATE> vec(v.x-a,v.y-a,v.z-a);
	return vec;
}

template <class TEMPLATE> inline _vec<TEMPLATE> operator-(const TEMPLATE& a ,const _vec<TEMPLATE> v ){
	_vec<TEMPLATE> vec(a-v.x,a-v.y,a-v.z);
	return vec;
}

//�x�N�g���΃x�N�g��
template <class TEMPLATE> inline TEMPLATE operator*(const _vec<TEMPLATE>& v1, const _vec<TEMPLATE>& v2){
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <class TEMPLATE> inline _vec<TEMPLATE> operator%(const _vec<TEMPLATE>& v1, const _vec<TEMPLATE>& v2){
	return _vec<TEMPLATE>((v1.y * v2.z - v1.z * v2.y), (v1.z * v2.x - v1.x * v2.z), (v1.x * v2.y - v1.y * v2.x));
}

template <class TEMPLATE> inline _vec<TEMPLATE> operator+(const _vec<TEMPLATE>& v1, const _vec<TEMPLATE>& v2){
	
	return _vec<TEMPLATE>(v1.x+v2.x,v1.y+v2.y,v1.z + v2.z);
}

template <class TEMPLATE> inline _vec<TEMPLATE> operator-(const _vec<TEMPLATE>& v1, const _vec<TEMPLATE>& v2){
	return _vec<TEMPLATE>(v1.x-v2.x,v1.y-v2.y,v1.z - v2.z);
}

template <class TEMPLATE> inline _vec<TEMPLATE> operator/(const _vec<TEMPLATE>& v1, const _vec<TEMPLATE>& v2){
	return _vec<TEMPLATE>( (v1.x/v2.x), (v1.y/v2.y), (v1.z/v2.z));
}

//�x�N�g����]
/*
template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::rot( const TEMPLATE& angle )
{
  TEMPLATE c=cos(angle), s=sin(angle);
  TEMPLATE t;
  t = c*x - s*y;
  y = s*x + c*y;
  x = t;

  return *this;  
}
*/
//�x�N�g���m����
template <class TEMPLATE> inline TEMPLATE _vec<TEMPLATE>::norm(void) const{
	TEMPLATE s = x*x + y*y + z*z;
	return sqrt(s);
}


template <class TEMPLATE> inline TEMPLATE _vec<TEMPLATE>::sqr(void) const{
	TEMPLATE s = x*x + y*y + z*z;
	return s;
}
#endif
