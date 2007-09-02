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

#ifndef PIRANHA_BASE_PSERIES_IO_H
#define PIRANHA_BASE_PSERIES_IO_H

#include "math.h"                                 // math::min.
#include "utils.h"                                // str_to_vector_double.

namespace piranha
{
/// Read coefficient argument.
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::read_cf_arg(std::ifstream &inf)
  {
    std::string temp, temp_name;
    vector_double temp_vdouble;
    std::streampos cur_pos=inf.tellg();
    while (utils::get_valid_string(inf,temp)==0)
    {
      if (temp[0]=='[')
      {
        std::cout << "Finished parsing cf_arg." << std::endl;
        inf.seekg(cur_pos);
        append_cf_args(vector_psym_p(1,psymbol_manager::get_pointer(psymbol(temp_name,temp_vdouble))));
        return;
      }
      deque_string split_v;
      boost::split(split_v,temp,boost::is_any_of("="));
      if (split_v.size()!=2)
      {
        std::cout << "Invalid line in cf_arg section." << std::endl;
      }
      else if (split_v[0]=="name")
      {
        std::cout << "name=" << split_v[1] << std::endl;
        temp_name=split_v[1];
      }
      else if (split_v[0]=="poly_eval")
      {
        std::cout << "poly_eval=" << split_v[1] << std::endl;
        temp_vdouble=utils::str_to_vector_double(split_v[1]);
      }
      else
      {
        std::cout << "Unknown field in cf_arg section." << std::endl;
      }
      cur_pos=inf.tellg();
    }
  }

/// Read trigonometric argument.
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::read_trig_arg(std::ifstream &inf)
  {
    std::string temp, temp_name;
    vector_double temp_vdouble;
    std::streampos cur_pos=inf.tellg();
    while (utils::get_valid_string(inf,temp)==0)
    {
      if (temp[0]=='[')
      {
        std::cout << "Finished parsing trig_arg." << std::endl;
        inf.seekg(cur_pos);
        append_trig_args(vector_psym_p(1,psymbol_manager::get_pointer(psymbol(temp_name,temp_vdouble))));
        return;
      }
      deque_string split_v;
      boost::split(split_v,temp,boost::is_any_of("="));
      if (split_v.size()!=2)
      {
        std::cout << "Invalid line in trig_arg section." << std::endl;
      }
      else if (split_v[0]=="name")
      {
        std::cout << "name=" << split_v[1] << std::endl;
        temp_name=split_v[1];
      }
      else if (split_v[0]=="poly_eval")
      {
        std::cout << "poly_eval=" << split_v[1] << std::endl;
        temp_vdouble=utils::str_to_vector_double(split_v[1]);
      }
      else
      {
        std::cout << "Unknown field in trig_arg section." << std::endl;
      }
      cur_pos=inf.tellg();
    }
  }

// Read linear arguments
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::read_lin_args(std::ifstream &inf)
  {
    std::string temp;
    deque_string split_v;
    std::streampos cur_pos=inf.tellg();
    while (utils::get_valid_string(inf,temp)==0)
    {
      if (temp[0]=='[')
      {
        std::cout << "Finished parsing lin_args." << std::endl;
        inf.seekg(cur_pos);
        return;
      }
      boost::split(split_v,temp,boost::is_any_of(stream_manager::data_separator()));
      if (split_v.size()>trig_width())
      {
        std::cout << "WARNING: lin_args is wider than series, arguments in excess" <<
          " will be cropped" << std::endl;
      }
      else if (split_v.size()<trig_width())
      {
        std::cout << "lin_args is smaller than series" << std::endl;
      }
      size_t min_w=math::min(trig_width(),split_v.size());
      for (size_t j=0;j<min_w;++j)
      {
        lin_args_[j]=utils::lexical_converter<int>(split_v[j]);
      }
      cur_pos=inf.tellg();
    }
  }

// Read terms.
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::read_terms(std::ifstream &inf,const std::string &filename)
  {
    const std::string pl_name=filename+".phl";
    phase_list pl=phase_list(pl_name);
    phase_list::const_iterator it_pl=pl.begin();
    deque_string split_v;
    std::string temp;
    it_s_index it_hint=end();
    while (utils::get_valid_string(inf,temp)==0)
    {
      boost::split(split_v,temp,boost::is_any_of(stream_manager::data_separator()));
      if (split_v.size()<2)
      {
        std::cout << "Warning: not enough elements in term string, ignoring term." << std::endl;
        continue;
      }
      term_type new_term;
// Read cf.
      cf_type cf(split_v[0]);
      while (cf_width()<cf.actual_width())
      {
        std::cout << "Warning: cf width is larger than expected, assigning 'null' to extra cf args."
          << std::endl;
        append_cf_args(vector_psym_p(1,psymbol_manager::get_pointer(psymbol())));
      }
      *new_term.s_cf()=cf;
// Ditch out first element of string vector, now that we read it.
      split_v.pop_front();
// Read trigonometric part.
      trig_type trig(split_v);
// TODO: see if it is possible to group the addition of symbols instead of doing it once at a time.
      while (trig_width()<trig.actual_width())
      {
        std::cout << "Warning: trig width is larger than expected, assigning 'null' to extra trig args."
          << std::endl;
        append_trig_args(vector_psym_p(1,psymbol_manager::get_pointer(psymbol())));
      }
      *new_term.s_trig()=trig;
// Deal with phases.
      if (it_pl==pl.end())
      {
        it_hint=insert(new_term,true,&it_hint);
      }
      else
      {
        term_type tmp_term;
        switch (pl.operation())
        {
          case phase_list::add
            :
          add_phase_to_term(*it_pl,new_term,tmp_term,*this);
          break;
          default:
            add_phase_to_term(*it_pl-new_term.g_trig()->phase(trig_s_vec_),new_term,tmp_term,*this);
        }
        ++it_pl;
      }
    }
  }

// Identify section
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::identify_sections(std::ifstream &inf,const std::string &filename)
  {
    std::string temp;
    while (utils::get_valid_string(inf,temp)==0)
    {
      if (temp[0]=='[')
      {
        std::cout << "New section found: " << temp << std::endl;
        if (temp=="[cf_arg]")
        {
          read_cf_arg(inf);
        }
        else if (temp=="[trig_arg]")
        {
          read_trig_arg(inf);
        }
        else if (temp=="[data]")
        {
          read_terms(inf,filename);
        }
        else if (temp=="[lin_args]")
        {
          read_lin_args(inf);
        }
        else
        {
          std::cout << "Found unknown section \"" << temp << "\", ignoring." << std::endl;
        }
      }
      else
      {
        std::cout << "Found string not belonging to any section: " << temp << std::endl;
      }
    }
  }

// Read data
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::read_data_from_file(std::ifstream &inf,
    const std::string &filename)
  {
    identify_sections(inf,filename);
    std::cout << "EOF" << std::endl;
  }

/// Load series from file.
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::load_from(const std::string &fn)
  {
    std::ifstream inf;
    std::string filename=utils::open_file(fn,inf);
// Read from file
    if (inf.is_open())
    {
      read_data_from_file(inf,filename);
// Close file
      inf.close();
    }
  }

/// Print the first n terms in plain format.
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::print_terms_plain(std::ostream &out_stream,
    int n) const
  {
    stream_manager::setup_print(out_stream);
    size_t j=0, lim;
    if (n<0)
    {
      lim=length();
    }
    else
    {
      lim=(size_t)n;
    }
    for (it_s_index it=g_s_index().begin();it!=g_s_index().end();++it)
    {
      if (j==lim)
      {
        break;
      }
      it->print_plain(out_stream,cf_s_vec_,trig_s_vec_);
      out_stream << std::endl;
      ++j;
    }
  }

/// Print the first n terms in LaTeX format.
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::print_terms_latex(std::ostream &out_stream,
    int n) const
  {
    stream_manager::setup_print(out_stream);
    size_t i=0, lim;
    if (n<0)
    {
      lim=length();
    }
    else
    {
      lim=(size_t)n;
    }
    for (it_s_index it=g_s_index().begin();it!=g_s_index().end();++it)
    {
      if (i==lim)
      {
        break;
      }
// Increase i now, so we print correct term number (as opposed to term index).
      ++i;
      out_stream << i << "&";
      it->print_latex(out_stream,cf_s_vec_,trig_s_vec_);
      out_stream << "\\\\" << std::endl;
    }
  }

/// Print the series in plain format up to term number n.
/**
 * Plain format is the same format used in formatted input from file.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::print_plain(std::ostream &out_stream,
    int limit) const
  {
    stream_manager::setup_print(out_stream);
    size_t j;
    for (j=0;j<cf_s_vec_.size();++j)
    {
      out_stream << "[cf_arg]" << std::endl;
      cf_s_vec_[j]->print(out_stream);
    }
    for (j=0;j<trig_s_vec_.size();++j)
    {
      out_stream << "[trig_arg]" << std::endl;
      trig_s_vec_[j]->print(out_stream);
    }
    out_stream << "[lin_args]" << std::endl;
    for (j=0;j<lin_args_.size();++j)
    {
      out_stream << lin_args_[j];
      if (j==lin_args_.size()-1)
      {
        out_stream << std::endl;
      }
      else
      {
        out_stream << stream_manager::data_separator();
      }
    }
    out_stream << "[data]" << std::endl;
    print_terms_plain(out_stream,limit);
  }

/// Print the series in LaTeX format up to term number n.
/**
 * Write to a LaTeX tabular environment. The output can then be included in a
 * LaTeX document.
 */
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::print_latex(std::ostream &out_stream, int limit) const
  {
    stream_manager::setup_print(out_stream);
    out_stream << "\\begin{xtabular}{rlrrrl}" << std::endl;
    print_terms_latex(out_stream,limit);
    out_stream << "\\end{xtabular}" << std::endl;
  }

/// Save series to file.
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::save_to(const std::string &filename) const
  {
    std::ofstream outf(filename.c_str(),std::ios_base::trunc);
    if (outf.fail())
    {
      std::cout << "Error saving to file " << filename << "." << std::endl;
      outf.close();
      return;
    }
    print(outf);
    outf.close();
  }

/// Print to screen frequencies and phases.
  template <class Cf, class Trig, template <class, class> class Term, template <class, class, template <class, class> class> class I, class Derived>
    inline void base_pseries<Cf, Trig, Term, I, Derived>::put_phases_freqs(int limit) const
  {
    stream_manager::setup_print(std::cout);
    size_t i=0, lim;
    if (limit==0)
    {
      lim=length();
    }
    else
    {
      lim=limit;
    }
    for (iterator it=begin();it!=end();++it)
    {
      it->g_cf()->print_plain(std::cout,cf_s_vec_);
      std::cout << stream_manager::data_separator() << it->g_trig()->phase(trig_s_vec_) <<
        stream_manager::data_separator() << it->g_trig()->freq(trig_s_vec_) << std::endl;
      ++i;
      if (i==lim)
      {
        break;
      }
    }
  }
}
#endif
