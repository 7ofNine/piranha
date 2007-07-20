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

#ifndef PIRANHA_PS_TERM_DEF_H
#define PIRANHA_PS_TERM_DEF_H

namespace piranha
{
/// Poisson series term class.
/**
 * Accepts coefficient and trigonometric part as template parameters.
 */
  template <class Cf, class Trig>
    class ps_term
  {
    public:
/// Alias for self.
      typedef ps_term self;
/// Alias for coefficient type.
      typedef Cf cf_type;
/// Alias for trigonometric type.
      typedef Trig trig_type;
/// Functor to update the coefficient.
      struct modifier_update_cf
      {
        modifier_update_cf(const cf_type &new_cf): new_cf_(new_cf)
          {}
        ~modifier_update_cf()
          {}
        void operator()(ps_term &term)
        {
          term.s_cf() = new_cf_;
        }
// NOTICE: evaluate the impact of using const & here, esp. when using gmp
        const cf_type &new_cf_;
      };
      explicit ps_term();
// FIXME: replace bool with enum.
      explicit ps_term(const cf_type &, bool flavour=true);
/// Copy constructor from term with different template parameters.
      template <class Cf2>
        explicit ps_term(const ps_term<Cf2,trig_type> &term):
      private_flavour_(term.g_flavour()),private_cf_(term.g_cf()),private_trig_(term.g_trig())
      {                                           /*BOOST_STATIC_ASSERT(sizeof(U)==0);*/
      }
// Getters
/// Get coefficient reference.
      cf_type &s_cf()
      {
        return private_cf_;
      }
      const cf_type &g_cf() const
      {
        return private_cf_;
      }
/// Get reference to trigonometric part.
      trig_type &s_trig()
      {
        return private_trig_;
      }
      const trig_type &g_trig() const
      {
        return private_trig_;
      }
/// Get flavour.
      bool &s_flavour()
      {
        return private_flavour_;
      }
      const bool &g_flavour() const
      {
        return private_flavour_;
      }
      double freq(const vector_psym_p &) const;
      double phase(const vector_psym_p &) const;
      bool checkup(const size_t &, const size_t &) const;
      size_t footprint() const;
      bool is_ignorable(const vector_psym_p &) const;
      bool operator<(const ps_term &) const;
// Manipulation.
      void increase_size(const size_t &, const size_t &);
      void invert_trig_args();
// Maths.
/// Assignment operator.
      ps_term &operator=(const ps_term &t2)
      {
        if (this==&t2)
        {
          return *this;
        }
        s_cf()=t2.g_cf();
        s_trig()=t2.g_trig();
        s_flavour()=t2.g_flavour();
        return *this;
      }
// I/O.
      void print_plain(std::ostream &, const vector_psym_p &, const vector_psym_p &) const;
      void print_latex(std::ostream &, const vector_psym_p &, const vector_psym_p &) const;
/// Print to stream.
/**
 * Print format (i.e., plain or latex) is set in piranha::stream_manager.
 */
      void print(std::ostream &out_stream, const vector_psym_p &cv, const vector_psym_p &tv) const
      {
        switch (stream_manager::format())
        {
          case stream_manager::plain:
            print_plain(out_stream,cv,tv);
            break;
          case stream_manager::latex:
            print_latex(out_stream,cv,tv);
        }
      }
/// Print to screen.
      void put(const vector_psym_p &cv, const vector_psym_p &tv) const
      {
        print(std::cout,cv,tv);
      }

      typename cf_type::eval_type t_eval(double, const vector_psym_p &, const vector_psym_p &) const;
      double norm(const vector_psym_p &) const;
    private:
// Probing.
      int trig_sign() const;
    private:
// Data members
      bool        private_flavour_;
      cf_type     private_cf_;
      trig_type   private_trig_;
  };

/// Default constructor.
  template <class Cf, class Trig>
    inline ps_term<Cf,Trig>::ps_term():
  private_flavour_(true),private_cf_(),private_trig_()
    {}

/// Constructor from coefficient and flavour.
  template <class Cf, class Trig>
    inline ps_term<Cf,Trig>::
    ps_term(const cf_type &c, bool flavour):
  private_flavour_(flavour),private_cf_(c),private_trig_()
    {}
}
#endif
