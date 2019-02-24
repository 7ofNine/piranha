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

#ifndef PIRANHA_STATS_H
#define PIRANHA_STATS_H



#include <boost/lexical_cast.hpp>
#include <map>
#include <sstream>
#include <string>

#include "config.h"
#include "exceptions.h"

#pragma warning (push)
#pragma warning (disable: 4251)

namespace piranha
{

class PIRANHA_VISIBLE stats
{
	public:
		static void set(const std::string &key, const std::string &field)
		{
			m_container[key] = field;
		}


		static std::string get(const std::string &key)
		{
			const std::map<std::string,std::string>::const_iterator it = m_container.find(key);
			if (it == m_container.end()) 
			{
				PIRANHA_THROW(value_error,"stats field not found");
			}
			return it->second;
		}


		template <class T, class Functor>
		static void trace_stat(const std::string &key, const T &initial, const Functor &f)
		{
			T cur_value(initial);
			try {
				cur_value = boost::lexical_cast<T>(get(key));
			} catch (const value_error &) 
			{
				set(key,boost::lexical_cast<std::string>(initial));
			}
			cur_value = f(cur_value);
			set(key,boost::lexical_cast<std::string>(cur_value));
		}


		static std::string dump()
		{
			std::ostringstream oss;
			for  (std::map<std::string,std::string>::const_iterator it = m_container.begin(); it != m_container.end(); ++it) {
				oss << it->first << ": " << it->second << '\n';
			}
			return oss.str();
		}

	private:

		static std::map<std::string, std::string> m_container;
};

}

#pragma warning (pop)
#endif
