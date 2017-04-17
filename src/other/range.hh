namespace range {
    /* free functions
     *  front()
     *  advance()
     *  empty()
     *  begin()/end()
     */

    template <typename Ib_t, typename Ie_t>
    struct range_from_begin_end_t {
        Ib_t m_b;
        Ie_t m_e;
        auto begin() const { return m_b; }
        auto end  () const { return m_e; }
        bool empty() const { return m_b == m_e; }
        void advance()     { ++m_b; }
        auto current_it() const { return m_b; }
    };

    template <typename Ib_t, typename Ie_t>
    auto range_from_begin_end(Ib_t b, Ie_t e) {
        return range_from_begin_end_t<Ib_t,Ie_t>{move(b),move(e)};
    }
} // namespace range
