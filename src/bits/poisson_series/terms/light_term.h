#ifndef PIRANHA_LIGHT_TERM_H
#define PIRANHA_LIGHT_TERM_H

namespace piranha
{
  /// Lightweight term to be used in series multiplications.
  template <class Cf, class Trig>
    struct light_term
  {
    struct hasher
    {
      size_t operator()(const light_term &t) const
      {
        return t.m_trig.hash_value();
      }
    };
    /// Default ctor.
    /**
     * Won't initialize anything.
     */
    light_term() {}
    ~light_term() {}
    // Manipulation.
    void invert_trig_args()
    {
      m_trig.invert_sign();
      if (!m_trig.flavour())
      {
        m_cf.invert_sign();
      }
    }
    bool operator==(const light_term &t) const
    {
      return (m_trig == t.m_trig);
    }
    // Getters and setters.
    const Cf &cf() const
    {
      return m_cf;
    }
    const Trig &trig() const
    {
      return m_trig;
    }
    Cf &cf()
    {
      return m_cf;
    }
    Trig &trig()
    {
      return m_trig;
    }
    // Data members.
    mutable Cf    m_cf;
    mutable Trig  m_trig;
  };
}
#endif
