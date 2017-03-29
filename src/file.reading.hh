#include <vector>
#include <string>
#include <memory>

namespace file_reading {

struct Genotypes_I {
    virtual int         number_of_snps() const = 0;
};

using GenotypeFileHandle = std:: shared_ptr<Genotypes_I const> const ;

GenotypeFileHandle      read_in_a_raw_ref_file(std:: string file_name);

} // namespace file_reading
