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

GenotypeFileHandle      read_in_a_raw_ref_file(std:: string file_name);

} // namespace file_reading
