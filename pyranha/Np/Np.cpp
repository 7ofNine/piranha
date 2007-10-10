/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
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

#include "../pyranha.h"

BOOST_PYTHON_MODULE(_Np)
{
  class_<np> inst=ps_basic_instantiation<np>("np","Numerical Poisson series class.");
  ps_instantiate_differential_specifics(inst);
  ps_instantiate_real_specifics(inst);
  def("kep_cosE",&astro::kep_cosE<np>,"Solve Kepler's equation for cosE.");
  def("Pnm",&math::Pnm<np>,"Legendre function of the first kind - Pnm(cos(theta)).");
  def("Ynm",&math::Ynm<np>,"Non-normalized spherical harmonic.");
  def("wig_rot",&math::wig_rot<np>,"Wigner rotation theorem for spherical harmonics.");
  instantiate_tass17<lnp>();
}
