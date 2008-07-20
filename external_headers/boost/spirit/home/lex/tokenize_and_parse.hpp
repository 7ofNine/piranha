/*=============================================================================
    Copyright (c) 2001-2008 Hartmut Kaiser

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
#if !defined(BOOST_SPIRIT_LEXER_PARSE_NOV_17_2007_0246PM)
#define BOOST_SPIRIT_LEXER_PARSE_NOV_17_2007_0246PM

#include <boost/spirit/home/qi/meta_grammar.hpp>
#include <boost/spirit/home/qi/skip.hpp>
#include <boost/spirit/home/qi/nonterminal/grammar.hpp>
#include <boost/spirit/home/support/unused.hpp>
#include <boost/spirit/home/lex/lexer.hpp>
#include <boost/mpl/assert.hpp>

namespace boost { namespace spirit { namespace lex
{
    ///////////////////////////////////////////////////////////////////////////
    //
    //  The tokenize_and_parse() function is one of the main Spirit API 
    //  functions. It simplifies using a lexer as the underlying token source
    //  while parsing a given input sequence.
    //
    //  The function takes a pair of iterators spanning the underlying input 
    //  stream to parse, the lexer object (built from the token definitions) 
    //  and a parser object (built from the parser grammar definition).
    //
    //  The second version of this function additionally takes an attribute to 
    //  be used as the top level data structure instance the parser should use 
    //  to store the recognized input to.
    //
    //  The function returns true if the parsing succeeded (the given input
    //  sequence has been successfully matched by the given grammar).
    //
    //  first, last:    The pair of iterators spanning the underlying input 
    //                  sequence to parse. These iterators must at least 
    //                  conform to the requirements of the std::intput_iterator 
    //                  category.
    //                  On exit the iterator 'first' will be updated to the 
    //                  position right after the last successfully matched 
    //                  token. 
    //  lex:            The lexer object (encoding the token definitions) to be
    //                  used to convert the input sequence into a sequence of 
    //                  tokens. This token sequence is passed to the parsing 
    //                  process. The LexerExpr type must conform to the 
    //                  lexer interface described in the corresponding section
    //                  of the documentation.
    //  xpr:            The grammar object (encoding the parser grammar) to be
    //                  used to match the token sequence generated by the lex 
    //                  object instance. The ParserExpr type must conform to 
    //                  the grammar interface described in the corresponding 
    //                  section of the documentation.
    //  attr:           The top level attribute passed to the parser. It will 
    //                  be populated during the parsing of the input sequence.
    //                  On exit it will hold the 'parser result' corresponding 
    //                  to the matched input sequence.
    //
    ///////////////////////////////////////////////////////////////////////////
    template <typename Iterator, typename LexerExpr, typename ParserExpr>
    inline bool
    tokenize_and_parse(Iterator& first, Iterator last, LexerExpr const& lex,
        ParserExpr const& xpr)
    {
        typedef typename LexerExpr::iterator_type iterator_type;
        typedef spirit::traits::is_component<qi::domain, ParserExpr> 
            is_component;

        // report invalid expression error as early as possible
        BOOST_MPL_ASSERT_MSG(
            is_component::value,
            xpr_is_not_convertible_to_a_parser, (iterator_type, ParserExpr));

        typedef typename 
            result_of::as_component<qi::domain, ParserExpr>::type 
        component;
        typedef typename component::director director;
        component c = spirit::as_component(qi::domain(), xpr);

        iterator_type iter = lex.begin(first, last);
        return director::parse(c, iter, lex.end(), unused, unused, unused);
    }

    ///////////////////////////////////////////////////////////////////////////
    template <typename Iterator, typename LexerExpr, typename ParserExpr,
        typename Attribute>
    inline bool
    tokenize_and_parse(Iterator& first, Iterator last, LexerExpr const& lex,
        ParserExpr const& xpr, Attribute& attr)
    {
        typedef typename LexerExpr::iterator_type iterator_type;
        typedef spirit::traits::is_component<qi::domain, ParserExpr> 
            is_component;

        // report invalid expression error as early as possible
        BOOST_MPL_ASSERT_MSG(
            is_component::value,
            xpr_is_not_convertible_to_a_parser, (iterator_type, ParserExpr));

        typedef typename 
            result_of::as_component<qi::domain, ParserExpr>::type 
        component;
        typedef typename component::director director;
        component c = spirit::as_component(qi::domain(), xpr);

        iterator_type iter = lex.begin(first, last);
        return director::parse(c, iter, lex.end(), unused, unused, attr);
    }

    ///////////////////////////////////////////////////////////////////////////
    //
    //  The tokenize_and_phrase_parse() function is one of the main Spirit API 
    //  functions. It simplifies using a lexer as the underlying token source
    //  while phrase parsing a given input sequence.
    //
    //  The function takes a pair of iterators spanning the underlying input 
    //  stream to parse, the lexer object (built from the token definitions) 
    //  and a parser object (built from the parser grammar definition). The 
    //  additional skipper parameter will be used as the skip parser during
    //  the parsing process.
    //
    //  The second version of this function additionally takes an attribute to 
    //  be used as the top level data structure instance the parser should use 
    //  to store the recognized input to.
    //
    //  The function returns true if the parsing succeeded (the given input
    //  sequence has been successfully matched by the given grammar).
    //
    //  first, last:    The pair of iterators spanning the underlying input 
    //                  sequence to parse. These iterators must at least 
    //                  conform to the requirements of the std::intput_iterator 
    //                  category.
    //                  On exit the iterator 'first' will be updated to the 
    //                  position right after the last successfully matched 
    //                  token. 
    //  lex:            The lexer object (encoding the token definitions) to be
    //                  used to convert the input sequence into a sequence of 
    //                  tokens. This token sequence is passed to the parsing 
    //                  process. The LexerExpr type must conform to the 
    //                  lexer interface described in the corresponding section
    //                  of the documentation.
    //  xpr:            The grammar object (encoding the parser grammar) to be
    //                  used to match the token sequence generated by the lex 
    //                  object instance. The ParserExpr type must conform to 
    //                  the grammar interface described in the corresponding 
    //                  section of the documentation.
    //  attr:           The top level attribute passed to the parser. It will 
    //                  be populated during the parsing of the input sequence.
    //                  On exit it will hold the 'parser result' corresponding 
    //                  to the matched input sequence.
    //  skipper_:       The skip parser to be used while parsing the given 
    //                  input sequence. Note, the skip parser will have to 
    //                  act on the same token sequence as the main parser 
    //                  'xpr'.
    //
    ///////////////////////////////////////////////////////////////////////////
    template <
        typename Iterator, typename LexerExpr, typename ParserExpr, 
        typename Skipper
    >
    inline bool
    tokenize_and_phrase_parse(Iterator& first, Iterator last, 
        LexerExpr const& lex, ParserExpr const& xpr, Skipper const& skipper_)
    {
        typedef typename LexerExpr::iterator_type iterator_type;
        typedef spirit::traits::is_component<qi::domain, ParserExpr> 
            expr_is_component;
        typedef spirit::traits::is_component<qi::domain, Skipper> 
            skipper_is_component;

        // report invalid expressions error as early as possible
        BOOST_MPL_ASSERT_MSG(
            expr_is_component::value,
            xpr_is_not_convertible_to_a_parser, 
            (iterator_type, ParserExpr, Skipper));

        BOOST_MPL_ASSERT_MSG(
            skipper_is_component::value,
            skipper_is_not_convertible_to_a_parser, 
            (iterator_type, ParserExpr, Skipper));

        typedef typename 
            result_of::as_component<qi::domain, ParserExpr>::type 
        component;
        typedef typename component::director director;
        component c = spirit::as_component(qi::domain(), xpr);

        typename result_of::as_component<qi::domain, Skipper>::type
            skipper = spirit::as_component(qi::domain(), skipper_);

        iterator_type iter = lex.begin(first, last);
        if (!director::parse(c, iter, lex.end(), unused, skipper, unused))
            return false;

        // do a final post-skip
        skip(iter, lex.end(), skipper);
        return true;
    }

    ///////////////////////////////////////////////////////////////////////////
    template <
        typename Iterator, typename LexerExpr, typename ParserExpr, 
        typename Attribute, typename Skipper
    >
    inline bool
    tokenize_and_phrase_parse(Iterator& first, Iterator last, 
        LexerExpr const& lex, ParserExpr const& xpr, Attribute& attr,
        Skipper const& skipper_)
    {
        typedef typename LexerExpr::iterator_type iterator_type;
        typedef spirit::traits::is_component<qi::domain, ParserExpr> 
            expr_is_component;
        typedef spirit::traits::is_component<qi::domain, Skipper> 
            skipper_is_component;

        // report invalid expressions error as early as possible
        BOOST_MPL_ASSERT_MSG(
            expr_is_component::value,
            xpr_is_not_convertible_to_a_parser, 
            (iterator_type, ParserExpr, Skipper));

        BOOST_MPL_ASSERT_MSG(
            skipper_is_component::value,
            skipper_is_not_convertible_to_a_parser, 
            (iterator_type, ParserExpr, Skipper));

        typedef typename 
            result_of::as_component<qi::domain, ParserExpr>::type 
        component;
        typedef typename component::director director;
        component c = spirit::as_component(qi::domain(), xpr);

        typename result_of::as_component<qi::domain, Skipper>::type
            skipper = spirit::as_component(qi::domain(), skipper_);

        iterator_type iter = lex.begin(first, last);
        if (!director::parse(c, iter, lex.end(), unused, skipper, attr))
            return false;

        // do a final post-skip
        skip(iter, lex.end(), skipper);
        return true;
    }

    ///////////////////////////////////////////////////////////////////////////
    //
    //  The tokenize() function is one of the main Spirit API functions. It 
    //  simplifies using a lexer to tokenize a given input sequence. It's main
    //  purpose is to use the lexer to tokenize all the input. 
    
    //  The second version below discards all generated tokens afterwards. 
    //  This is useful whenever all the needed functionality has been 
    //  implemented directly inside the lexer semantic actions, which are being 
    //  executed while the tokens are matched. 
    //
    //  The function takes a pair of iterators spanning the underlying input 
    //  stream to scan, the lexer object (built from the token definitions),
    //  and a (optional) functor being call for each of the generated tokens. 
    //
    //  The function returns true if the scanning of the input succeeded (the 
    //  given input sequence has been successfully matched by the given token
    //  definitions).
    //
    //  first, last:    The pair of iterators spanning the underlying input 
    //                  sequence to parse. These iterators must at least 
    //                  conform to the requirements of the std::intput_iterator 
    //                  category.
    //                  On exit the iterator 'first' will be updated to the 
    //                  position right after the last successfully matched 
    //                  token. 
    //  lex:            The lexer object (encoding the token definitions) to be
    //                  used to convert the input sequence into a sequence of 
    //                  tokens. The LexerExpr type must conform to the 
    //                  lexer interface described in the corresponding section
    //                  of the documentation.
    //  f:              A functor (callable object) taking a single argument of
    //                  the token type and returning a bool, indicating whether
    //                  the tokenization should be canceled.
    //
    ///////////////////////////////////////////////////////////////////////////
    template <typename Iterator, typename LexerExpr, typename F>
    inline bool
    tokenize(Iterator& first, Iterator last, LexerExpr const& lex, F f)
    {
        typedef typename LexerExpr::iterator_type iterator_type;

        iterator_type end = lex.end();
        for (iterator_type iter = lex.begin(first, last); iter != end; ++iter) 
        {
            if (!f(*iter))
                return false;
        }
        return true;
    }

    ///////////////////////////////////////////////////////////////////////////
    template <typename Iterator, typename LexerExpr>
    inline bool
    tokenize(Iterator& first, Iterator last, LexerExpr const& lex)
    {
        typedef typename LexerExpr::iterator_type iterator_type;

        iterator_type iter = lex.begin(first, last);
        iterator_type end = lex.end();

        while (iter != end && token_is_valid(*iter))
            ++iter;
            
        return (iter == end) ? true : false;
    }

}}}

#endif
