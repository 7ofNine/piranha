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

#include "../pyranha.h"

BOOST_PYTHON_MODULE(_Mpfs)
{
  class_<mpfs> inst=ps_basic_instantiation<mpfs>("mpfs","Multiple precision Fourier series class.");
  ps_instantiate_differential_specifics(inst);
  ps_instantiate_real_specifics(inst);
  fourier_specifics(inst);
  def("kep_cosE",&astro::kep_cosE<mpfs>,"Solve Kepler's equation for cosE.");
  def("Pnm",&math::Pnm<mpfs>,"Legendre function of the first kind - Pnm(cos(theta)).");
  def("Ynm",&math::Ynm<mpfs>,"Non-normalized spherical harmonic.");
  def("wig_rot",&math::wig_rot<mpfs>,"Wigner rotation theorem for spherical harmonics.");
  instantiate_tass17<mpfs>();
}
