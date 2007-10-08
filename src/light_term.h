#ifndef PIRANHA_LIGHT_TERM_H
#define PIRANHA_LIGHT_TERM_H

namespace piranha
{
/// Lightweight term to be used in series multiplications.
    template <class Cf, class Trig>
      struct light_term
    {
/// Default ctor.
/**
 * Won't initialize anything.
 */
      light_term()
        {}
      ~light_term()
        {}
// Manipulation.
      void invert_trig_args()
      {
        trig.invert_sign();
        if (!trig.g_flavour())
        {
// FIXME: maybe here a invert_sign function for the coefficient should be used?
          cf*=-1;
        }
      }
      bool operator==(const light_term &t) const
      {
        return (trig == t.trig);
      }
// Getters and setters.
      const Cf *g_cf() const
      {
        return &cf;
      }
      const Trig *g_trig() const
      {
        return &trig;
      }
      Cf *s_cf()
      {
        return &cf;
      }
      Trig *s_trig()
      {
        return &trig;
      }
      size_t hasher() const
      {
        return trig.hasher();
      }
// Data members.
      mutable Cf    cf;
      mutable Trig  trig;
    };
}

#endif
