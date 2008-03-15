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

#ifndef PIRANHA_NAMED_SERIES_IO_H
#define PIRANHA_NAMED_SERIES_IO_H

#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>

#include "../exceptions.h"

namespace piranha
{
  // TMP for series printing.
  template <class ArgsDescr>
    inline void named_series_print_plain(std::ostream &stream,
    const typename ntuple<vector_psym_p,boost::tuples::length<ArgsDescr>::value>::type &args_tuple)
  {
    for (size_t i=0; i < args_tuple.get_head().size(); ++i)
    {
      stream << "[" << ArgsDescr::head_type::name << "_arg]" << std::endl;
      args_tuple.get_head()[i]->print(stream);
    }
    named_series_print_plain<typename ArgsDescr::tail_type>(stream,args_tuple.get_tail());
  }

  template <>
    inline void named_series_print_plain<boost::tuples::null_type>(std::ostream &,
    const ntuple<vector_psym_p,boost::tuples::length<boost::tuples::null_type>::value>::type &)
  {}

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::print_plain(std::ostream &stream, int limit) const
  {
    named_series_print_plain<arguments_description>(stream,m_arguments);
    stream << "[terms]" << std::endl;
    derived_const_cast->print_terms_plain(stream,m_arguments,limit);
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::print_latex(std::ostream &, int) const
  {
    // TODO: implement.
  }

  /// Print series to stream.
  /**
   * Print first "limit" terms. If limit is negative, print all terms. The output format is read
   * from the piranha::stream_manager class.
   */
  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::print(std::ostream &out_stream, int limit) const
  {
    switch (stream_manager::format())
    {
      case stream_manager::plain:
        print_plain(out_stream,limit);
        break;
      case stream_manager::latex:
        print_latex(out_stream,limit);
    }
  }

  /// Construct from file.
  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::construct_from_file(const std::string &fn)
  {
    std::ifstream inf;
    std::string filename=utils::open_file(fn,inf);
    // Read from file
    if (inf.is_open())
    {
      // Clear the stack of unknown data.
      unknown_data.clear();
      read_sections(inf);
      std::cout << "EOF" << std::endl;
      // Close file
      inf.close();
    }
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::read_sections(std::ifstream &inf)
  {
    std::string temp;
    while (utils::get_valid_string(inf,temp))
    {
      if (temp.size() > 2 and temp[0] == '[' and temp[temp.size()-1] == ']')
      {
        std::string sec_name = temp;
        boost::trim_if(sec_name,boost::is_any_of("[]"));
        std::cout << "New section found: " << sec_name << std::endl;
        std::vector<std::string> split_v;
        boost::split(split_v,sec_name,boost::is_any_of("_"));
        if (split_v.size() == 2 and split_v[1] == "arg")
        {
          read_arg(inf,split_v[0]);
        }
        else if (sec_name == "terms")
        {
          read_terms(inf);
          // If we found the data, then we don't want any more sections.
          return;
        }
        else
        {
          std::cout << "Found unknown section '" << sec_name << "', ignoring." << std::endl;
          unknown_data.push_back(temp);
        }
      }
      else
      {
        std::cout << "Found string not belonging to any (known) section: " << temp << std::endl;
        unknown_data.push_back(temp);
      }
    }
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::read_arg(std::ifstream &inf, const std::string &name)
  {
    // Temporary attributes for the argument.
    std::string temp_name, temp_time_eval;
    // Record stream position, so we can rewind once we finish parsing the argument.
    std::streampos cur_pos=inf.tellg();
    // Temporary storage.
    std::string temp;
    while (utils::get_valid_string(inf,temp))
    {
      // If we found a new section, step back the cursor before exiting.
      if (temp.size() > 2 and temp[0] == '[' and temp[temp.size()-1] == ']')
      {
        std::cout << "Finished parsing " << name << " argument." << std::endl;
        inf.seekg(cur_pos);
        append_arg(name,psymbol_manager::get_pointer(psymbol(temp_name,temp_time_eval)).second);
        return;
      }
      std::vector<std::string> split_v;
      boost::split(split_v,temp,boost::is_any_of("="));
      if (split_v.size() != 2)
      {
        std::cout << "Invalid line in "<< name << " argument section." << std::endl;
      }
      else if (split_v[0] == "name")
      {
        std::cout << "name=" << split_v[1] << std::endl;
        temp_name=split_v[1];
      }
      else if (split_v[0] == "time_eval")
      {
        std::cout << "time_eval=" << split_v[1] << std::endl;
        temp_time_eval=split_v[1];
      }
      else
      {
        std::cout << "Unknown field in " << name << " argument section." << std::endl;
      }
      cur_pos=inf.tellg();
    }
  }

  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::read_terms(std::ifstream &inf)
  {
    typedef typename Derived::term_type term_type;
    typedef typename Derived::const_sorted_iterator const_sorted_iterator;
    std::string temp;
    const_sorted_iterator it_hint = derived_const_cast->template nth_index<0>().end();
    while (!inf.eof())
    {
      getline(inf,temp,derived_const_cast->separator);
      boost::trim(temp);
      // Ignore empty lines.
      if (temp.empty())
      {
        continue;
      }
      try
      {
        term_type term(temp,m_arguments);
        if (!term.is_insertable(m_arguments))
        {
          throw bad_input("Term not insertable in named series.");
        }
        it_hint = derived_cast->insert(term,derived_const_cast->m_arguments,it_hint);
      }
      catch (bad_input &b)
      {
        std::cout << b.what() << std::endl;
      }
    }
  }

  /// Save series to file.
  template <__PIRANHA_NAMED_SERIES_TP_DECL>
    inline void named_series<__PIRANHA_NAMED_SERIES_TP>::save_to(const std::string &filename) const
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
}

#endif
