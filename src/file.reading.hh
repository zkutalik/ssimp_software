#include <vector>
#include <string>
#include <memory>
#include <ostream>

namespace file_reading {

struct chrpos {
    int chr;
    int pos;
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
};

using GenotypeFileHandle = std:: shared_ptr<Genotypes_I const> const ;

struct SNPiterator {
    GenotypeFileHandle m_gfh;
    int                m_line_number; // 0 means the first SNP that was read, 1 the second ...

    SNPiterator &       operator++()        ;

    chrpos get_chrpos() const {
        return m_gfh->get_chrpos(m_line_number);
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
