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

#ifndef PIRANHA_DIFFERENTIAL_TOOLBOX_H
#define PIRANHA_DIFFERENTIAL_TOOLBOX_H

namespace piranha
{
/// Toolbox for differential calculus.
  template <class Derived>
    class differential_toolbox
  {
    public:
/// Partial derivative.
/**
 * Calculate partial derivative with respect to psymbol "name".
 * @param[in] name std::string, psymbol's name.
 */
// TODO: deal with linargs.
// TODO: special handling of "time" symbol.
// TODO: hinted insertion.
      Derived partial(const std::string &name) const
      {
        typedef typename Derived::ancestor::cf_type cf_type;
        typedef typename Derived::ancestor::it_h_index it_h_index;
        typedef typename Derived::ancestor::term_type term_type;
        Derived retval;
        retval.merge_args(*static_cast<Derived const *>(this));
        const std::pair<bool,size_t> cf_s_index=static_cast<Derived const *>(this)->cf_arg_index(name),
          trig_s_index=static_cast<Derived const *>(this)->trig_arg_index(name);
// If "name" is not a symbol of the series, return empty series.
        if (!cf_s_index.first and !trig_s_index.first)
        {
          std::cout << "No psymbol named '" << name << "', returning empty series in partial derivation." << std::endl;
          return retval;
        }
        const it_h_index it_f=static_cast<Derived const *>(this)->g_h_index().end();
        term_type tmp_term;
        for (it_h_index it=static_cast<Derived const *>(this)->g_h_index().begin();it!=it_f;++it)
        {
// First part of the derivation of the product coefficient * trigonometric_part.
          if (cf_s_index.first)
          {
            it->cf().partial(cf_s_index.second,tmp_term.cf());
            tmp_term.trig()=it->trig();
            tmp_term.trig().flavour()=it->trig().flavour();
            retval.insert_with_checks(tmp_term);
          }
// Second part of the derivation.
// NOTICE: this may be placed somewhere inside trig classes, but probably to do this
// we need to make trig_args() aware of flavour (i.e., move flavour outside Term class).
          if (trig_s_index.first)
          {
            tmp_term.cf()=it->cf();
            switch (it->trig().flavour())
            {
              case true:
                tmp_term.cf().mult_by(-it->trig()[trig_s_index.second]);
                tmp_term.trig().flavour()=false;
                break;
              case false:
                tmp_term.cf().mult_by(it->trig()[trig_s_index.second]);
                tmp_term.trig().flavour()=true;
            }
// Perform this check since if we already assigned trig_args above we don't need to
// do it again now.
            if (!cf_s_index.first)
            {
              tmp_term.trig()=it->trig();
            }
            retval.insert_with_checks(tmp_term);
          }
        }
        return retval;
      }
  };
}
#endif
