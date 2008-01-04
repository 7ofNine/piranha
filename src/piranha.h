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

#ifndef PIRANHA_PIRANHA_H
#define PIRANHA_PIRANHA_H

/// Piranha top-level namespace.
namespace piranha {}

// Include all piranha manipulators.
#include "manipulators/ffs.h"
#include "manipulators/ffsc.h"
#include "manipulators/cfs.h"
#include "manipulators/cfsc.h"
#include "manipulators/fs.h"
#include "manipulators/fsc.h"
#include "manipulators/mpfs.h"
#include "manipulators/mpfsc.h"

// Include ipoly for now.
// TODO: remove it later.
#include "bits/ipoly.h"

// Include TASS.
#include "tass17/tass17.h"

#endif
