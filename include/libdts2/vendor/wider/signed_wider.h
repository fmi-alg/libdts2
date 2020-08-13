#pragma once
#ifndef LIB_DTS2_VENDOR_SIGNED_WIDER_H
#define LIB_DTS2_VENDOR_SIGNED_WIDER_H

#include <libdts2/constants.h>

#include <limits>
#include <libdts2/vendor/wider/wider.h>

namespace LIB_DTS2_NAMESPACE {
	
//A minimal wrapper for Wider to support computation of the determinants
template<class T_UNSIGNED_WIDER>
class SignedWider {
public:
	using unsigned_type = T_UNSIGNED_WIDER;
	using Self = SignedWider<unsigned_type>;
public:
	SignedWider(Self const &) = default;
	SignedWider(Self &&) = default;
	template<class T>
	SignedWider(SignedWider<T> const & other) :
	m_s(other.sign()),
	m_v(other.unsigned_value())
	{}
	template<class T>
	SignedWider(SignedWider<T> && other) :
	m_s(other.sign()),
	m_v(std::move(other.unsigned_value()))
	{}
	SignedWider(uint32_t v) : m_v(uint64_t(v)) {}
	SignedWider(int32_t v) : m_s(v < 0), m_v(uint64_t(std::abs(v))) {}
	SignedWider(uint64_t v) : m_v(v) {}
	SignedWider(int64_t v) : m_s(v < 0), m_v(uint64_t(std::abs(v))) {}
	SignedWider(unsigned_type const & _unsigned_value, bool _sign = false) :
	m_s(_sign),
	m_v(_unsigned_value)
	{
		assert(!m_s || m_v > 0);
	}
	SignedWider(unsigned_type && _unsigned_value, bool _sign = false) :
	m_s(_sign),
	m_v(std::move(_unsigned_value))
	{
		assert(!m_s || m_v > 0);
	}
public:
	Self & operator=(Self const &) = default;
	Self & operator=(Self &&) = default;
public:
	Self operator+(Self const & other) const {
		if (m_s == other.m_s) {
			return Self(m_v + other.m_v, m_s);
		}
		else if (m_s && !other.m_s) {
			if (other.m_v >= m_v) {
				return Self(other.m_v - m_v, false);
			}
			else {
				return Self(m_v - other.m_v, true);
			}
		}
		else { //(!m_s && other.m_s) 
			if (m_v >= other.m_v) {
				return Self(m_v - other.m_v, false);
			}
			else {
				return Self(other.m_v - m_v, true);
			}
		}
	}
	Self operator-(Self const & other) const {
		if (m_s xor other.m_s) {
			return Self(m_v + other.m_v, m_s);
		}
		else if (m_s && other.m_s) {
			if (other.m_v >= m_v) {
				return Self(other.m_v - m_v, false);
			}
			else {
				return Self(m_v - other.m_v, true);
			}
		}
		else { //(!m_s && !other.m_s) 
			if (m_v >= other.m_v) {
				return Self(m_v - other.m_v, false);
			}
			else {
				return Self(other.m_v - m_v, true);
			}
		}
	}
	Self operator*(Self const & other) const {
		return Self(m_v * other.m_v, m_s xor other.m_s);
	}
	Self operator-() const {
		return Self(m_v, !m_s);
	}
	Self & operator*=(int v) {
		m_v *= unsigned_type(uint64_t(std::abs(v)));
		m_s = m_s xor (v < 0);
		return *this;
	}
	Self & operator*=(uint64_t v) {
		m_v *= unsigned_type(v);
		return *this;
	}
public:
	bool operator<(Self const & other) const {
		if (m_s xor other.m_s) { //we are smaller if we are the signed one
			return m_s;
		}
		else if (m_s) { //both are signed, hence we are smaller if we are absolutly greater
			return m_v > other.m_v;
		}
		else { //both are unsigned
			return m_v < other.m_v;
		}
	};
	bool operator<(int other) const {
		if (m_s xor (other < 0)) { //we are smaller if we are the signed one
			return m_s;
		}
		else if (m_s) { //both are signed, henace we are smaller if we are absolutly greater
			return m_v > uint64_t(std::abs(other));
		}
		else { //both are unsigned
			return m_v < uint64_t(other);
		}
	};
public:
	bool operator>(Self const & other) const {
		if (m_s xor other.m_s) { //we are greater if we are the usigned one
			return !m_s;
		}
		else if (m_s) { //both are signed, hence we are greater if we are absolutly smaller
			return other.m_v < m_v;
		}
		else { //both are unsigned
			return m_v > other.m_v;
		}
	}
	bool operator>(int other) const {
		if (m_s xor (other < 0)) { //we are greater if we are the usigned one
			return !m_s;
		}
		else if (m_s) { //both are signed, hence we are greater if we are absolutly smaller
			return uint64_t(abs(other)) < m_v;
		}
		else { //both are unsigned
			return m_v > uint64_t(std::abs(other));
		}
	};
	bool operator==(Self const & other) const {
		if (m_s xor other.m_s) {
			return false;
		}
		else {
			return m_v == other.m_v;
		}
	}
public:
	bool sign() const { return m_s; }
	unsigned_type const & unsigned_value() const { return m_v; }
	unsigned_type & unsigned_value() { return m_v; }
protected:
	bool m_s{false};
	unsigned_type m_v;
};
	
}//end namespace dts

#endif
