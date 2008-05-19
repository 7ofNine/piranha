/***************************************************************************
 *   Copyright (C) 2007, 2008 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef PIRANHA_MPQ_H
#define PIRANHA_MPQ_H

#include <boost/operators.hpp>
#include <gmp.h>
#include <gmpxx.h>
#include <sstream>
#include <string>

#include "exceptions.h"
#include "integer_typedefs.h"

namespace piranha
{
	/// Wrapper for GMP's mpq_class.
	class mpq:
				boost::field_operators < mpq,
				boost::field_operators < mpq, max_fast_int,
				boost::field_operators < mpq, double,
				boost::less_than_comparable< mpq, max_fast_int,
				boost::less_than_comparable< mpq, double,
				boost::less_than_comparable< mpq,
				boost::equality_comparable< mpq, max_fast_int,
				boost::equality_comparable< mpq
				> > > > > > > >
	{
		public:
			explicit mpq(): m_value(0) {}
			explicit mpq(const double &x): m_value(x) {}
			explicit mpq(const max_fast_int &n, const max_fast_int &d): m_value(n, d) {
				m_value.canonicalize();
			}
// 			mpq(const std::string &s): m_value(0) {
// 				try {
// 					m_value = mpq_class(s);
// 					m_value.canonicalize();
// 				} catch (...) {
// 					throw bad_input("Error converting string to rational.");
// 				}
// 			}
			mpq copy() const {return mpq(*this);}
			std::string print_to_string() const {
				std::ostringstream stream;
				stream << m_value;
				std::string retval(stream.str());
				return retval;
			}
			bool operator==(const max_fast_int &n2) const {
				return m_value == n2;
			}
			bool operator==(const mpq &q2) const {
				return m_value == q2.m_value;
			}
			bool operator<(const max_fast_int &n2) const {
				return m_value < n2;
			}
			bool operator<(const double &x2) const {
				return m_value < x2;
			}
			bool operator<(const mpq &q2) const {
				return m_value < q2.m_value;
			}
			bool operator>(const max_fast_int &n2) const {
				return m_value > n2;
			}
			bool operator>(const double &x2) const {
				return m_value > x2;
			}
			mpq &operator+=(const max_fast_int &n) {
				m_value += n;
				return *this;
			}
			mpq &operator+=(const double &x) {
				m_value += x;
				return *this;
			}
			mpq &operator+=(const mpq &m) {
				m_value += m.m_value;
				return *this;
			}
			mpq &operator-=(const max_fast_int &n) {
				m_value -= n;
				return *this;
			}
			mpq &operator-=(const double &x) {
				m_value -= x;
				return *this;
			}
			mpq &operator-=(const mpq &m) {
				m_value -= m.m_value;
				return *this;
			}
			mpq &operator*=(const max_fast_int &n) {
				m_value *= n;
				return *this;
			}
			mpq &operator*=(const double &x) {
				m_value *= x;
				return *this;
			}
			mpq &operator*=(const mpq &m) {
				m_value *= m.m_value;
				return *this;
			}
			mpq &operator/=(const max_fast_int &n) {
				if (n == 0) {
					throw division_by_zero();
				}
				m_value /= n;
				return *this;
			}
			mpq &operator/=(const double &x) {
				if (x == 0) {
					throw division_by_zero();
				}
				m_value /= x;
				return *this;
			}
			mpq &operator/=(const mpq &m) {
				if (m == 0) {
					throw division_by_zero();
				}
				m_value /= m.m_value;
				return *this;
			}
			mpq pow(const max_fast_int &n) const
			{
				mpq retval;
				// If value = 0, then we must make sure we are not dividing by zero.
				if (m_value == 0) {
					if (n == 0) {
						retval.m_value = 1;
					} else if (n > 0) {
						retval.m_value = 0;
					} else {
						throw division_by_zero();
					}
				} else {
					if (n < 0) {
						mpz_pow_ui(mpq_denref(retval.m_value.get_mpq_t()), mpq_numref(m_value.get_mpq_t()), (size_t)(-n));
						mpz_pow_ui(mpq_numref(retval.m_value.get_mpq_t()), mpq_denref(m_value.get_mpq_t()), (size_t)(-n));
					} else {
						mpz_pow_ui(mpq_numref(retval.m_value.get_mpq_t()), mpq_numref(m_value.get_mpq_t()), (size_t)n);
						mpz_pow_ui(mpq_denref(retval.m_value.get_mpq_t()), mpq_denref(m_value.get_mpq_t()), (size_t)n);
					}
				}
				return retval;
			}
			mpq pow(const double &x) const
			{
				mpq retval;
				if (m_value == 1) {
					retval.m_value = 1;
				} else if (m_value == 0) {
					if (x == 0) {
						retval.m_value = 1;
					} else if (x > 0) {
						retval.m_value = 0;
					} else {
						throw division_by_zero();
					}
				} else {
					throw unsuitable("Cannot raise rational number to real power.");
				}
				return retval;
			}
		private:
			mpq_class	m_value;
	};
}

#endif
