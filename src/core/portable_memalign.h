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

/*****************************************************************************
 *
 *  XVID MPEG-4 VIDEO CODEC
 *  - Aligned Memory Allocator header  -
 *
 *  Copyright(C) 2002-2003 Edouard Gomez <ed.gomez@free.fr>
 *
 *  This program is free software ; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation ; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY ; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program ; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
 *
 * $Id: mem_align.h,v 1.8 2004/03/22 22:36:24 edgomez Exp $
 *
 ****************************************************************************/

#ifndef PIRANHA_PORTABLE_MEMALIGN_H
#define PIRANHA_PORTABLE_MEMALIGN_H

#include <cstdlib>
#include <stdint.h>

namespace piranha
{
	template <int Alignment>
	inline void *portable_memalign(const size_t &size)
	{
		uint8_t *mem_ptr;
		if (!Alignment) {
			// We have not to satisfy any alignment
			if ((mem_ptr = (uint8_t *)malloc(size + 1)) != NULL) {
				// Store (mem_ptr - "real allocated memory") in *(mem_ptr-1)
				*mem_ptr = (uint8_t)1;
				// Return the mem_ptr pointer
				return ((void *)(mem_ptr+1));
			}
		} else {
			uint8_t *tmp;
			// Allocate the required size memory + alignment so we
			// can realign the data if necessary
			if ((tmp = (uint8_t *)malloc(size + Alignment)) != NULL) {
				// Align the tmp pointer
				mem_ptr = (uint8_t *)((size_t)(tmp + Alignment - 1) & (~(size_t)(Alignment - 1)));
				/* Special case where malloc have already satisfied the alignment
				* We must add alignment to mem_ptr because we must store
				* (mem_ptr - tmp) in *(mem_ptr-1)
				* If we do not add alignment to mem_ptr then *(mem_ptr-1) points
				* to a forbidden memory space */
				if (mem_ptr == tmp) {
					mem_ptr += Alignment;
				}
				// (mem_ptr - tmp) is stored in *(mem_ptr-1) so we are able to retrieve
				// the real malloc block allocated and free it in portable_aligned_free
				*(mem_ptr - 1) = (uint8_t)(mem_ptr - tmp);
				// Return the aligned pointer
				return ((void *)mem_ptr);
			}
		}
		return(NULL);
	}

	inline void portable_aligned_free(void *mem_ptr)
	{
		uint8_t *ptr;
		if (mem_ptr == NULL) {
			return;
		}
		// Aligned pointer
		ptr = (uint8_t *)mem_ptr;
		// *(ptr - 1) holds the offset to the real allocated block
		// we sub that offset os we free the real pointer
		ptr -= *(ptr - 1);
		// Free the memory
		free(ptr);
	}
}

#endif
