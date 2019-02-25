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

// 
// This class is currently nowhere used and would not compile !!
//

#ifndef PIRANHA_MPF_CF_H
#define PIRANHA_MPF_CF_H

#include <complex>
#include <gmp.h>
#include <gmpxx.h>

#include "../base_classes/numerical_container.h"
#include "../common_typedefs.h"
#include "../math.h"

namespace piranha
{
	/// Mpf numerical coefficient.
	/**
	 * This class can be used as coefficient in Poisson series. It encapsulate a GMP mpf_class
	 * value and provides access to it.
	 *
	 * A set of operators is provided to enable interoperability with basic numerical data types.
	 */
	class mpf_cf: public NumericalContainer<mpf_class, mpf_cf>
	{
			/// Alias for self.
			typedef mpf_cf self;
			/// Alias for the parent class.
			typedef NumericalContainer<mpf_class, mpf_cf> ancestor;

		public:

			// Start implementation of basic pseries coefficient interface.
			//------------
			// Ctors and dtor.
			/// Empty constructor.
			explicit mpf_cf(): ancestor::NumericalContainer() {}

			/// Constructor from string.
			template <class ArgsTuple>
			explicit mpf_cf(const std::string &s, const ArgsTuple &a): ancestor::NumericalContainer(s, a) {}

			/// Constructor from integer.
			explicit mpf_cf(int val): ancestor::NumericalContainer(val) {}

			/// Constructor from double.
			explicit mpf_cf(const double &val): ancestor::NumericalContainer(val) {}

			/// Constructor from Psym.
			template <class ArgsTuple>
			explicit mpf_cf(const psym_p &p, const int &n, const ArgsTuple &a): ancestor::NumericalContainer(p, n, a) {}

			// Override norm and evaluation.
			template <class ArgsTuple>
			double norm_(const ArgsTuple &) const 
            {
				return std::abs(g_value().get_d());
			}

			template <class ArgsTuple>
			double t_eval(const double &, const ArgsTuple &) const 
            {
				return g_value().get_d();
			}
			// End implementation of basic pseries coefficient interface.
			//------------

			// Start implementation of trigonometric pseries coefficient interface.
			// Used in:
			// - trigonometric toolbox,
			//------------
			template <class ArgsTuple>
			self besselJ(int n, const ArgsTuple &) const 
            {
				self retval;
				retval.s_value() = math::besselJ(n, g_value().get_d());
				return retval;
			}

			// End implementation of trigonometric pseries coefficient interface.
			//------------
			// Start implementation of power-enabled pseries coefficient interface.
			template <class ArgsTuple>
			self pow(const double &y, const ArgsTuple &) const 
            {
				self retval;
				retval.s_value() = std::pow(g_value().get_d(), y);
				return retval;
			}
	};
}


namespace std
{
	template <>
	struct complex<piranha::mpf_cf>:
				public piranha::NumericalContainer<complex<mpf_class>, complex<piranha::mpf_cf> >,
				public piranha::NumericalContainerComplexToolbox<piranha::mpf_cf> {
private:
		typedef piranha::NumericalContainer<complex<mpf_class>, complex<piranha::mpf_cf> > ancestor;
		typedef piranha::NumericalContainerComplexToolbox<piranha::mpf_cf> complex_toolbox;
		typedef complex self;
		friend class piranha::NumericalContainerComplexToolbox<piranha::mpf_cf>;
public:
		typedef piranha::mpf_cf value_type;
//       using ancestor::swap;
//       using ancestor::print_plain;
//       using ancestor::print_latex;
//       using ancestor::checkup;
//       using ancestor::invert_sign;
//       using ancestor::is_ignorable;
//       using ancestor::is_insertable;
//       using ancestor::needs_padding;
//       using ancestor::pad_right;
//       using ancestor::apply_layout;
//       using ancestor::add;
//       using ancestor::subtract;
//       using ancestor::mult_by;
//       using ancestor::mult_by_self;
//       using ancestor::divide_by;
//       using complex_toolbox::mult_by;
//       using complex_toolbox::divide_by;
//       using complex_toolbox::mult_by_self;
//       using complex_toolbox::real;
//       using complex_toolbox::imag;
//       using complex_toolbox::set_real;
//       using complex_toolbox::set_imag;
		// Start implementation of basic pseries coefficient interface.
		//------------
		// Basic ctors and dtor.
		explicit complex(): ancestor::NumericalContainer() {}
		template <class ArgsTuple>
		explicit complex(const std::string &s, const ArgsTuple &a): ancestor::NumericalContainer(s, a) {}
		explicit complex(int n): ancestor::NumericalContainer(n) {}
		explicit complex(const double &x): ancestor::NumericalContainer(x) {}
		complex(const complex &c): ancestor::NumericalContainer(c), complex_toolbox::NumericalContainerComplexToolbox(c) {}
		// Complex specific contructors.
		explicit complex(int r, int i): complex_toolbox::NumericalContainerComplexToolbox(r, i) {}
		explicit complex(const std::complex<int> &c): complex_toolbox::NumericalContainerComplexToolbox(c) {}
		explicit complex(const double &r, const double &i): complex_toolbox::NumericalContainerComplexToolbox(r, i) {}
		explicit complex(const std::complex<double> &c): complex_toolbox::NumericalContainerComplexToolbox(c) {}
		explicit complex(const value_type &r): complex_toolbox::NumericalContainerComplexToolbox(r) {}
		explicit complex(const value_type &r, const value_type &i):
				complex_toolbox::NumericalContainerComplexToolbox(r, i) {}
		// Operators.
		self &operator=(const self &val2) {
			return assign_self(val2);
		}
		self &operator=(const value_type &r2) {
			return assign_self(r2);
		}
		// Probing.
		template <class ArgsTuple>
		double norm_(const ArgsTuple &) const {
			// NOTICE: the success of this probably depends upon std::complex implementation...
			return std::abs(g_value()).get_d();
		}
		// Override evaluation.
		template <class ArgsTuple>
		std::complex<double> t_eval(const double &, const ArgsTuple &) const {
			return std::complex<double>(g_value().real().get_d(), g_value().imag().get_d());
		}
		// End implementation of complex basic pseries coefficient interface.
		//------------
	};
}
#endif
