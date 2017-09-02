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

#ifndef PIRANHA_CONFIG_H
#define PIRANHA_CONFIG_H

// Platform switches.
#ifdef _PIRANHA_WIN32
#ifdef _PIRANHA_DLL_EXPORT_API
#define PIRANHA_VISIBLE __declspec(dllexport)
#elif defined ( _PIRANHA_DLL_IMPORT_API )
#define PIRANHA_VISIBLE __declspec(dllimport)
#else
#define PIRANHA_VISIBLE
#endif
#else
#define PIRANHA_VISIBLE __attribute__ ((visibility("default")))
#endif

#define __PIRANHA_MAX_ECHELON_LEVEL (4)


#endif
