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

#ifndef BASE_PSERIES_TA_MACROS_H
#define BASE_PSERIES_TA_MACROS_H

/// Template parameters for piranha::base_pseries.
#define __PIRANHA_BASE_PS_TP Cf,Trig,Term,I,Derived,Allocator
/// Template parameters for piranha::base_pseries (declaration form).
#define __PIRANHA_BASE_PS_TP_DECL class Cf, class Trig, template <class, class> class Term, template <class, class, \
  template <class, class> class> class I, class Derived, class Allocator

#endif
