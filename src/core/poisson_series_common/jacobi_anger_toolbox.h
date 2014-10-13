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

#ifndef PIRANHA_JACOBI_ANGER_TOOLBOX_H
#define PIRANHA_JACOBI_ANGER_TOOLBOX_H

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/type_traits/is_same.hpp>
#include <cstddef>
#include <complex>
#include <string>
#include <vector>

#include "../config.h"
#include "../exceptions.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <int TrigPos, class Derived>
	class jacobi_anger
	{
			PIRANHA_STATIC_CHECK(TrigPos >= 0, "Wrong trigonometric position in Jacobi-Anger toolbox.");
		public:

			template <class Term, class ArgsTuple>
			static void jacang(const std::vector<Term const *> &v,
				const typename std::vector<Term const *>::const_iterator &it_avoid, 
				std::complex<Derived> &retval, const ArgsTuple &argsTuple)
			{
				typedef typename std::vector<Term const *>::const_iterator const_iterator;
				PIRANHA_ASSERT(retval.empty());
				retval.base_add(1, argsTuple);
				if (v.empty()) 
				{
					return;
				}
				const const_iterator it_f = v.end();
				for (const_iterator it = v.begin(); it != it_f; ++it) 
				{
					// Skip the iterator we want to avoid.
					if (it != it_avoid) 
					{
						retval.base_mult_by(jacang_term(it, argsTuple), argsTuple);
					}
				}
			}


			template <class Iterator, class ArgsTuple>
			static std::complex<Derived> jacang_term(Iterator it, const ArgsTuple &argsTuple)
			{
				typedef typename std::complex<Derived>::TermType complex_term_type;
				// Let's determine the limit of the Jacobi-Anger development from the truncator of the series.
				// The Jacobi-Anger development is a development into Bessel functions of the first kind starting
				// from zero and increasing in unity steps, hence the psi function can be used
				// straightforwardly.
				std::size_t n_;
				{
					std::complex<Derived> tmp;
					complex_term_type tmp_term;
					tmp_term.cf.set_real((*it)->cf, argsTuple);
					tmp.insert(tmp_term, argsTuple);
					try {
						n_ = tmp.psi_(0,1,argsTuple);
					} catch (const value_error &ve) 
					{
						PIRANHA_THROW(value_error,std::string("unable to determine the limit of the Jacobi-Anger development. "
							"The reported error was:\n")+ve.what());
					}
				}
				const int n = boost::numeric_cast<int>(n_);
				std::complex<Derived> retval;
				{
					complex_term_type tmp_term;
					tmp_term.cf.set_real((*it)->cf.besselJ(0, argsTuple), argsTuple);
					retval.insert(tmp_term, argsTuple);
				}
				const std::size_t w = argsTuple.template get<TrigPos>().size();
				std::vector<typename std::complex<Derived>::TermType::key_type::value_type> tmp_trig_mults(w);
				std::complex<double> cos_multiplier(0, 2);
				for (int i = 1; i < n; ++i) 
				{
					complex_term_type tmp_term;
					tmp_term.cf.set_real((*it)->cf.besselJ(i, argsTuple), argsTuple);
					std::copy((*it)->key.begin(),(*it)->key.end(),tmp_trig_mults.begin());
					for (std::size_t j = 0; j < w; ++j) 
					{
						tmp_trig_mults[j] *= i;
					}

					tmp_term.key.resize(boost::numeric_cast<typename std::complex<Derived>::TermType::key_type::size_type>(tmp_trig_mults.size()));
					std::copy(tmp_trig_mults.begin(),tmp_trig_mults.end(),tmp_term.key.begin());
					if ((*it)->key.getFlavour())
					{
						tmp_term.cf.mult_by(cos_multiplier, argsTuple);
					} else 
					{
						if (i % 2 == 0) 
						{
							tmp_term.cf.mult_by(2, argsTuple);
						} else 
						{
							tmp_term.cf.mult_by(std::complex<double>(0, 2), argsTuple);
							tmp_term.key.setFlavour(false);
						}
					}

					retval.insert(tmp_term, argsTuple);
					// Update the multiplier for cosine terms.
					cos_multiplier *= std::complex<double>(0, 1);
				}
				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
