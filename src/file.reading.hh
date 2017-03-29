#include <vector>
#include <string>
#include <memory>
#include <ostream>

namespace file_reading {

struct chrpos {
    int chr;
    int pos;
    bool operator< (chrpos const & other) {
        if(chr < other.chr) return true;
        if(chr > other.chr) return false;

        return pos < other.pos;
    }
};

inline
std:: ostream& operator<<(std:: ostream &o, chrpos const &c) {
    o
        << "chr"
        << c.chr
        << ':'
        << c.pos;
    return o;
}

struct Genotypes_I {
    virtual int         number_of_snps     ()        const = 0;
    virtual chrpos      get_chrpos         (int)     const = 0;
    virtual std::string get_SNPname        (int)     const = 0;
    virtual std::string get_allele_ref     (int)     const = 0;
    virtual std::string get_allele_alt     (int)     const = 0;
};

using GenotypeFileHandle = std:: shared_ptr<Genotypes_I const>;

struct SNPiterator
: public std::iterator<std:: random_access_iterator_tag, chrpos>
{
    GenotypeFileHandle m_gfh;
    int                m_line_number; // 0 means the first SNP that was read, 1 the second, and so on

    // Constructor
                        SNPiterator(GenotypeFileHandle gfh, int line_number) : m_gfh(gfh), m_line_number(line_number) {}


    // operators
    SNPiterator &       operator++();
    bool                operator==(SNPiterator const & other);
    bool                operator!=(SNPiterator const & other);
    bool                operator< (SNPiterator const & other);
    bool                operator>=(SNPiterator const & other);
    int                 operator- (SNPiterator const & other);
    chrpos              operator* ()                           const;
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

    // Maybe I shouldn't have these static methods after all, might be confusing.
    static SNPiterator begin_from_file(GenotypeFileHandle gfh) {
        return {gfh, 0};
    }
    static SNPiterator   end_from_file(GenotypeFileHandle gfh) {
        return {gfh, gfh->number_of_snps()};
    }
};

GenotypeFileHandle      read_in_a_raw_ref_file(std:: string file_name);

} // namespace file_reading
