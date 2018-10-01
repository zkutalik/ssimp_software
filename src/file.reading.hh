#ifndef HH_FILE_READING_HH
#define HH_FILE_READING_HH

#include <vector>
#include <string>
#include <memory>
#include <ostream>
#include <unordered_map>

namespace file_reading {
char   decide_delimiter( std:: string      const & header_line ); // also useful while parsing --sample.names

struct chrpos {
    int chr;
    int pos;
    bool operator< (chrpos const & other) const {
        if(chr < other.chr) return true;
        if(chr > other.chr) return false;

        return pos < other.pos;
    }
    bool operator> (chrpos const & other) const {
        if(chr > other.chr) return true;
        if(chr < other.chr) return false;

        return pos > other.pos;
    }
    bool operator==(chrpos const & other) const {
        return(chr == other.chr
            && pos == other.pos);
    }
    bool operator!=(chrpos const & other) const {
        return !(*this == other);
    }
    bool operator<=(chrpos const & other) const {
        return !(*this > other);
    }
    bool operator>=(chrpos const & other) const {
        return !(*this < other);
    }
};
} // end that namespace, in order to define std:: hash

namespace std {
template<>
struct hash < file_reading:: chrpos >
{
    auto   operator() (file_reading:: chrpos const &crps) const
        noexcept(noexcept(
               hash< decltype(crps.chr) >{}(crps.chr)
             + hash< decltype(crps.pos) >{}(crps.pos)
                ))
        -> decltype(
               hash< decltype(crps.chr) >{}(crps.chr)
             + hash< decltype(crps.pos) >{}(crps.pos)
                )
    {
        return hash< decltype(crps.chr) >{}(crps.chr)
             + hash< decltype(crps.pos) >{}(crps.pos);
    }
};
}


namespace file_reading { // reopen this namespace

inline
std:: ostream& operator<<(std:: ostream &o, chrpos const &c) {
    o
        << "chr"
        << c.chr
        << ':'
        << c.pos;
    return o;
}

struct AnyFile_I {
    virtual int         number_of_snps     ()        const = 0;
    virtual std::string get_SNPname        (int)     const = 0;
    virtual chrpos      get_chrpos         (int)     const = 0;
    virtual std::string get_allele_ref     (int)     const = 0;
    virtual std::string get_allele_alt     (int)     const = 0;
    //virtual char        get_delimiter      ()       const = 0;
};

struct Effects_I : public AnyFile_I {
    virtual void        set_chrpos         (int, chrpos)  = 0; // so that we can fill them in from the ref data
    virtual void        sort_my_entries    ()             = 0;
    virtual void        delete_snps_with_no_position()    = 0;
    virtual int         delete_snps_with_identical_alleles()    = 0;
    virtual double      get_z              (int) const    = 0;
    virtual double      get_N              (int) const    = 0;
    virtual std::string get_column_name_allele_ref () const       = 0;
    virtual std::string get_column_name_allele_alt () const       = 0;
};

using GwasFileHandle_NONCONST     = std:: shared_ptr<Effects_I>;

struct SNPiterator
: public std::iterator<std:: random_access_iterator_tag, chrpos>
{
    using GWASorREF = GwasFileHandle_NONCONST; // TODO: get rid of this

    GWASorREF m_gfh;
    int                m_line_number; // 0 means the first SNP that was read, 1 the second, and so on


    // Constructor
                        SNPiterator(GWASorREF gfh, int line_number) : m_gfh(gfh), m_line_number(line_number) {}


    // operators
    SNPiterator &       operator++();
    bool                operator==(SNPiterator const & other) const;
    bool                operator!=(SNPiterator const & other) const;
    bool                operator< (SNPiterator const & other) const;
    bool                operator>=(SNPiterator const & other) const;
    int                 operator- (SNPiterator const & other) const;
    chrpos              operator* ()                          const;
    SNPiterator &       operator+=(long int            ran);


    chrpos get_chrpos() const {
        return m_gfh->get_chrpos(m_line_number);
    }
    std:: string      get_SNPname() const {
        return m_gfh->get_SNPname(m_line_number);
    }
    std:: string      get_allele_ref() const {
        return m_gfh->get_allele_ref(m_line_number);
    }
    std:: string      get_allele_alt() const {
        return m_gfh->get_allele_alt(m_line_number);
    }

    template<typename T = void>
    auto   set_chrpos(chrpos crps)
        -> decltype ( (std::declval<T>(),m_gfh) -> set_chrpos(m_line_number, crps) )
    {
        return                           m_gfh  -> set_chrpos(m_line_number, crps);
    }
    template<typename T = void>
    auto   get_z()  const
        -> decltype ( (std::declval<T>(),m_gfh) -> get_z(m_line_number) )
    {
        return                           m_gfh  -> get_z(m_line_number);
    }
    template<typename T = void>
    auto   get_N()  const
        -> decltype ( (std::declval<T>(),m_gfh) -> get_N(m_line_number) )
    {
        return                           m_gfh  -> get_N(m_line_number);
    }
};

inline
SNPiterator begin_from_file(GwasFileHandle_NONCONST fh) {
    return {fh, 0};
}
inline
SNPiterator   end_from_file(GwasFileHandle_NONCONST fh) {
    return {fh, fh->number_of_snps()};
}

GwasFileHandle_NONCONST read_in_a_gwas_file(std:: string file_name);

} // namespace file_reading
#endif
