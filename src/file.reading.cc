#include "file.reading.hh"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <cassert>

#include <gzstream.h>

#include "other/DIE.hh"
#include "other/PP.hh"
#include "other/utils.hh"
#include "other/range_view.hh"

using std:: ifstream;
using std:: string;
using std:: vector;
using std:: map;
using utils:: ssize;
using utils:: tokenize;

#define LOOKUP(hd, fieldname, vec) lookup(hd . fieldname, vec, #fieldname)

namespace file_reading {
    using utils:: operator<<; // to print vectors

    template<typename G>
    SNPiterator<G> &       SNPiterator<G>:: operator++()        {
        assert(m_line_number < m_gfh->number_of_snps());
        ++m_line_number;
        return *this;
    }
    template<typename G>
    bool                SNPiterator<G>:: operator==(SNPiterator<G> const & other) const {
        assert(m_gfh == other.m_gfh);
        return m_line_number == other.m_line_number;
    }
    template<typename G>
    bool                SNPiterator<G>:: operator!=(SNPiterator<G> const & other) const {
        return !(*this == other);
    }
    template<typename G>
    bool                SNPiterator<G>:: operator< (SNPiterator<G> const & other) const {
        assert(m_gfh == other.m_gfh);
        return m_line_number < other.m_line_number;
    }
    template<typename G>
    bool                SNPiterator<G>:: operator>=(SNPiterator<G> const & other) const { return !(*this < other); }
    template<typename G>
    int                 SNPiterator<G>:: operator- (SNPiterator<G> const & other) const {
        assert(m_gfh == other.m_gfh);
        return m_line_number - other.m_line_number;
    }
    template<typename G>
    chrpos              SNPiterator<G>:: operator* () const {
        return get_chrpos();
    }
    template<typename G>
    SNPiterator<G> &       SNPiterator<G>:: operator+=(long int            ran)   {
        m_line_number += ran;
        return *this;
    }

    template struct
    SNPiterator<GenotypeFileHandle>;
    template struct
    SNPiterator<GwasFileHandle>;
    template struct
    SNPiterator<GwasFileHandle_NONCONST>;

struct header_details {
    struct offset_and_name {
        int                  m_offset;
        string               m_name;

        offset_and_name() : m_offset(-1)
        {}
        offset_and_name(int offset, string name) : m_offset(offset), m_name(name)
        {}
        offset_and_name(offset_and_name const &) = default;
        offset_and_name(offset_and_name      &&) = default;

        void operator=(offset_and_name const & from) {
            if(m_offset != -1) {
                WARNING("Two fields with a similar name, using [" << m_name << "] and ignoring [" << from.m_name << "].");
                return;
            }
            assert(m_offset == -1);
            m_offset = from.m_offset;
            m_name   = from.m_name;
        }
        void operator=(offset_and_name && from) {
            return *this = from;
        }
    };

    char                     m_delimiter;
    offset_and_name    SNPname;
    offset_and_name    chromosome;
    offset_and_name    position;
    offset_and_name    allele_ref;
    offset_and_name    allele_alt;
    offset_and_name    qual;
    offset_and_name    filter;
    offset_and_name    info;
    offset_and_name    format;
    offset_and_name    effect_z;
    vector<offset_and_name> unaccounted;
};


// Some forward declarations
static
GenotypeFileHandle      read_in_a_raw_ref_file_as_VCF(std:: string file_name);
static
std::pair<vector<uint8_t>,vector<uint8_t>> parse_many_calls (vector<string> const & calls_as_strings, int);
static
std:: pair<int,int> parse_call_pair (string const & call_for_this_person, int column_number, int i);

static
string          lookup( header_details:: offset_and_name const &on
                      , vector<string> const & all_split_up
                      , const char *fieldname
                      );
static
file_reading::
header_details   parse_header( string      const & header_line );
static bool is_in_this_list(string const & s, std:: initializer_list<char const *> candidates) ;
static
GwasFileHandle_NONCONST      read_in_a_gwas_file_simple(std:: string file_name);
static
char   decide_delimiter( string      const & header_line );
static
vector<int>                actual_lookup_one_ref_get_calls(SNPiterator<GenotypeFileHandle> it);

using call_type = std::pair<uint8_t,uint8_t>;
struct OneLineSummary {
    int                      m_simple_line_number;
    string                   m_SNPname;
    int                      m_chromosome;
    int                      m_position;
    string                   m_allele_alt;
    string                   m_allele_ref;

    vector<bool>             m_calls_compressed;
    vector<call_type>        m_calls_compressed_codes;
};
struct GwasLineSummary {
    string                   m_SNPname;
    chrpos                   m_chrpos;
    string                   m_allele_alt;
    string                   m_allele_ref;
    double                   m_z;
};

struct PlainVCFfile : public file_reading:: Genotypes_I
{
    vector<OneLineSummary> m_each_SNP_and_its_offset;
    string                 m_underlying_file_name;
    header_details         m_header_details;

    virtual int         number_of_snps() const {
        return m_each_SNP_and_its_offset.size();
    }
    virtual chrpos      get_chrpos         (int i)     const {
        OneLineSummary const & ols = get_ols(i);
        return chrpos{ols.m_chromosome, ols.m_position};
    }
    virtual std::string get_SNPname        (int i)     const {
        OneLineSummary const & ols = get_ols(i);
        return ols.m_SNPname;
    }
    virtual std::string get_allele_ref        (int i)     const {
        OneLineSummary const & ols = get_ols(i);
        return ols.m_allele_ref;
    }
    virtual std::string get_allele_alt        (int i)     const {
        OneLineSummary const & ols = get_ols(i);
        return ols.m_allele_alt;
    }
    virtual std::pair<
             std::vector<uint8_t>
            ,std::vector<uint8_t>
        > get_calls          (int i)                   const {
        OneLineSummary const & ols = get_ols(i);

        auto z_out = zip_val(  range:: from_vector( vector<uint8_t>{} )
                            ,  range:: from_vector( vector<uint8_t>{} ) );


        int current_run_of_zeros = 0;
        for(bool bit : ols.m_calls_compressed) {
            if(bit) {
                call_type call = ols.m_calls_compressed_codes.at(current_run_of_zeros);
                current_run_of_zeros = 0;
                z_out.push_back( call );
            } else {
                ++ current_run_of_zeros;
            }
        }
        auto all_calls_decoded_pair =
            utils:: apply
            (   [](auto && ...xs) {
                    return std:: make_pair( std::forward<decltype(xs)>(xs)...);
                }
            ,   move(z_out) | view:: unzip_collect_transpose
            );

        return all_calls_decoded_pair; // use the decompressed version
    }

    OneLineSummary  get_ols         (int i)     const {
        assert(i>=0);
        assert(i<number_of_snps());
        return m_each_SNP_and_its_offset.at(i);
    }

};

GenotypeFileHandle      read_in_a_raw_ref_file(std:: string file_name) {
    // I really should detect the file-type.
    // But for now, we'll just assume a plain (non-gzipped) vcf file.
    return read_in_a_raw_ref_file_as_VCF(file_name);
}

void update_positions_by_comparing_to_another_set( GwasFileHandle_NONCONST gwas, std:: unordered_map<std:: string, file_reading:: chrpos> const & m ) {
    auto       b_gwas = begin_from_file(gwas);
    auto const e_gwas =   end_from_file(gwas);
    for(;b_gwas < e_gwas; ++b_gwas) {
        auto nm= b_gwas.get_SNPname();
        auto chrpos_in_gwas = b_gwas.get_chrpos();
        auto try_to_find_in_ref = m.find(nm);

        if(try_to_find_in_ref != m.end()) {
            auto chrpos_in_ref = try_to_find_in_ref->second;

            // If the gwas doesn't have chrpos, copy it from the ref panel,
            // otherwise, check that they agree on the chrpos:
            if(chrpos_in_gwas == chrpos{-1,-1}) {
                    b_gwas.set_chrpos(chrpos_in_ref);
            }
            else {
                (chrpos_in_gwas == chrpos_in_ref) || DIE("position mismatch between GWAS and RefPanel. "
                        << "What should we do?. "
                        << '[' << nm << ']'
                        << ' ' << chrpos_in_ref
                        << ' ' << chrpos_in_gwas
                        );
            }
        }
    }
}

GwasFileHandle_NONCONST read_in_a_gwas_file(std:: string file_name) {
    // I really should detect the file-type.
    // But for now, we'll just assume a plain (non-gzipped) vcf file.
    return read_in_a_gwas_file_simple(file_name);
}


static
GenotypeFileHandle      read_in_a_raw_ref_file_as_VCF(std:: string file_name) {
    PP(file_name);

    // using gzstream to read gzipped input data. Thanks for http://www.cs.unc.edu/Research/compgeom/gzstream/#doc
    // It's LGPL license, we should remember to document this
    gz:: igzstream f(file_name.c_str());
    (f.rdbuf() && f.rdbuf()->is_open()) || DIE("Can't find file [" << file_name << ']');
    string current_line;

    int simple_line_number = 0;


    // Skip past the '##' lines
    while(getline(f, current_line)) {
        ++ simple_line_number;
        if(current_line.at(0) == '#' && current_line.at(1) == '#')
            continue; // skip these ## lines
        break;
    }

    // current_line is now the first line, i.e. the header


    auto p = std:: make_shared<PlainVCFfile>();
    p->m_underlying_file_name = file_name;
    p->m_header_details       = parse_header(current_line);

    header_details &hd = p->m_header_details;

    while(1) {
        OneLineSummary ols;
        getline(f, current_line);
        ++ simple_line_number;
        ols.m_simple_line_number = simple_line_number;
        // TODO: Should store a real 'line number' field directly in the 'ols' object

        if(!f) {
            f.eof() || DIE("Error before reaching eof() in this file");
            break;
        }
        auto all_split_up = tokenize(current_line, hd.m_delimiter);
        try {
            ols.m_SNPname    =                           LOOKUP(hd, SNPname, all_split_up);
            ols.m_chromosome = utils:: lexical_cast<int>( LOOKUP(hd, chromosome, all_split_up) );
            ols.m_position   = utils:: lexical_cast<int>( LOOKUP(hd, position, all_split_up) );
            ols.m_allele_alt =                           LOOKUP(hd, allele_alt, all_split_up);
            ols.m_allele_ref =                           LOOKUP(hd, allele_ref, all_split_up);

            { // Now, store a 'compressed' version of the 'unaccounted' fields

                int const N_ref = hd.unaccounted.size();
                vector<string> calls_as_strings;
                range:: ints(N_ref) |view:: foreach| [&](int person) {
                    auto column_number = hd.unaccounted.at(person).m_offset;
                    assert(column_number-9 == person);
                    string & call_for_this_person = all_split_up.at(column_number);
                    calls_as_strings.push_back(call_for_this_person);
                };

                auto lefts_and_rights = parse_many_calls(calls_as_strings, -1);

                // Count each distinct observed call type
                map<call_type, int> count_the_call_types;

                zip_val( range:: range_from_begin_end(lefts_and_rights.first)
                       , range:: range_from_begin_end(lefts_and_rights.second))
                |view:: foreach|
                [&](auto x) {
                    ++ count_the_call_types[ std:: make_pair(std::get<0>(x)
                                                            ,std::get<1>(x)
                                                            ) ];
                };

                int const number_of_distinct_call_types = count_the_call_types.size();

                struct count_and_call_t {
                    int count;
                    int left;
                    int right;
                    call_type as_pair() const { return {left,right}; }
                };

                vector<count_and_call_t> most_popular_first;
                for(auto & abc : count_the_call_types) {
                    most_popular_first.push_back( count_and_call_t{ abc.second
                    ,abc.first.first
                    ,abc.first.second}
                    );
                }
                range:: sort(range:: range_from_begin_end(most_popular_first), [](auto &l, auto &r) {
                        return l.count > r.count;
                });

                map<call_type, int> map_call_to_contiguous_ids;
                view:: enumerate_vector(most_popular_first)
                |view::unzip_foreach|
                [&](int i, count_and_call_t cac){
                    map_call_to_contiguous_ids[ cac.as_pair() ] = i;
                };
                assert(utils::ssize(map_call_to_contiguous_ids) == number_of_distinct_call_types);
                assert(utils::ssize(most_popular_first)         == number_of_distinct_call_types);
                assert(utils::ssize(count_the_call_types)       == number_of_distinct_call_types);

                vector<bool> bits_for_this_SNP;

                zip_val( range:: range_from_begin_end(lefts_and_rights.first)
                       , range:: range_from_begin_end(lefts_and_rights.second))
                |view:: unzip_foreach|
                [&](int left, int right) {
                    int code = map_call_to_contiguous_ids.at( std:: make_pair(left,right) );
                    assert(code >= 0);
                    assert(code < number_of_distinct_call_types);
                    for(int rep : range::ints(code)) {
                        (void)rep;
                        bits_for_this_SNP.push_back(false);
                    }
                    bits_for_this_SNP.push_back(true);
                };
                ols.m_calls_compressed = bits_for_this_SNP;

                auto just_calls_ordered_by_popularity =
                range:: from_vector( most_popular_first )
                |view:: map|
                [&](count_and_call_t cac) {
                    return call_type{ cac.left
                                    , cac.right };
                }
                | view:: collect;
                ols.m_calls_compressed_codes = just_calls_ordered_by_popularity;
            }

            p->m_each_SNP_and_its_offset.push_back(ols);
        } catch (std:: invalid_argument &e) {
            WARNING( "Ignoring this line, problem with the chromosome and/or position. ["
                    << "SNPname:" << LOOKUP(hd, SNPname   , all_split_up)
                    << " chr:" << LOOKUP(hd, chromosome, all_split_up)
                    << " pos:" << LOOKUP(hd, position  , all_split_up)
                    << "]" );
        }
    };

    sort( p->m_each_SNP_and_its_offset.begin()
        , p->m_each_SNP_and_its_offset.end()
        , [](auto &l, auto &r) -> bool {
            if(l.m_chromosome < r.m_chromosome) return true;
            if(l.m_chromosome > r.m_chromosome) return false;
            if(l.m_position   < r.m_position  ) return true;
            if(l.m_position   > r.m_position  ) return false;
            return l.m_simple_line_number <  r.m_simple_line_number;
        }
        );

    return p;
}

static
std::pair<vector<uint8_t>,vector<uint8_t>> parse_many_calls (vector<string> const & calls_as_strings, int line_number) {

    auto z_out = zip_val(  range:: from_vector( vector<uint8_t>{} )
                    ,  range:: from_vector( vector<uint8_t>{} ) );

    zip_val( range:: ints(calls_as_strings.size())
           , range:: range_from_begin_end(calls_as_strings) )
    |view:: unzip_foreach|
    [&](int person, string const & call_for_this_person) {
        z_out.push_back( parse_call_pair(call_for_this_person, person, line_number) );
    };
    auto l3 = move(z_out) | view:: unzip_collect_transpose;

    return make_pair(std::get<0>(l3)
                    ,std::get<1>(l3) );
};
static
std:: pair<int,int> parse_call_pair (string const & call_for_this_person, int column_number, int line_number) {
    int d1, d2;
    int n;
    int ret;

    ret = sscanf(call_for_this_person.c_str(), "%d|%d %n", &d1,&d2, &n); // note the space to allow trailing whitespace
    if(   ret == 2 // two digits
       && n == ssize(call_for_this_person) // all characters of input were consumed
       && d1 >= 0
       && d2 >= 0)
    {
        return {d1,d2};
    }

    ret = sscanf(call_for_this_person.c_str(), "%d|%d:%n", &d1,&d2, &n); // note the space to allow trailing whitespace
    if(   ret == 2 // two digits
       // In this case, we don't care what n is
       && d1 >= 0
       && d2 >= 0)
    {
        return {d1,d2};
    }

    ret = sscanf(call_for_this_person.c_str(), "%d/%d:%n", &d1,&d2, &n); // note the space to allow trailing whitespace
    if(   ret == 2 // two digits
       // In this case, we don't care what n is
       && d1 >= 0
       && d2 >= 0)
    {
        return {d1,d2};
    }

    DIE("Couldn't parse \"" << call_for_this_person << "\" in the " << 1+column_number << "th column in the " << 1+line_number << "th SNP in the ref panel");
    return {-1,-1}; // won't reach here - this doesn't matter
};

static
string          lookup( header_details:: offset_and_name const &on
                      , vector<string> const & all_split_up
                      , const char *fieldname
                      ) {
    int offset = on.m_offset;
    (offset >= 0 && offset < ssize(all_split_up)) || DIE("Cannot find ["<<fieldname<<"] field");

    return all_split_up.at(offset);
}

static
file_reading::
header_details   parse_header( string      const & header_line ) {
    header_details hd;

    hd.m_delimiter = decide_delimiter(header_line);
    auto header_names = tokenize(header_line, hd.m_delimiter);

    int field_counter = -1;
    for(auto & one_field_name : header_names) {
        ++field_counter;
        // go through each field name in turn and try to account for it
        if(false) {}
        else if(is_in_this_list(one_field_name, {"ID","rnpid"})) {
            hd.SNPname = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"#CHROM","chr"})) {
            hd.chromosome = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"POS"})) {
            hd.position = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"REF","a1"})) {
            hd.allele_ref = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"ALT","a2"})) {
            hd.allele_alt = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"QUAL"})) {
            hd.qual = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"FILTER"})) {
            hd.filter = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"INFO"})) {
            // actually, ignore this big field
            //hd.info = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"FORMAT"})) {
            hd.format = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"z.from.peff"
                                                ,"binary.mean.beta"
                                                ,"beta"
                                                })) {
            hd.effect_z = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else {
            hd.unaccounted.push_back( header_details:: offset_and_name(field_counter, one_field_name) );
        }
    }
    return hd;
}

static bool is_in_this_list(string const & s, std:: initializer_list<char const *> candidates) {
    for(auto cand : candidates) {
        if(s==cand)
            return true;
    }
    return false;
}


static
char   decide_delimiter( string      const & header_line ) {
    // Which is comma, tab, or space, are most common here?

    size_t commas = std::count(header_line.begin(), header_line.end(), ',');
    size_t tabs   = std::count(header_line.begin(), header_line.end(), '\t');
    size_t spaces = std::count(header_line.begin(), header_line.end(), ' ');
    char delimiter = '\0';
    if(commas >  0 && tabs == 0 && spaces == 0) delimiter = ',';
    if(commas == 0 && tabs >  0 && spaces == 0) delimiter = '\t';
    if(commas == 0 && tabs == 0 && spaces >  0) delimiter = ' ';
    if(delimiter == '\0') {
        DIE("Couldn't decide what delimiter to use."
                <<  " #commas:" << commas
                << ", #tabs:"     << tabs
                << ", #spaces:"     << spaces
           );
    }
    return delimiter;
}

struct SimpleGwasFile : public file_reading:: Effects_I
{
    vector<GwasLineSummary> m_each_SNP_and_its_z;
    string                  m_underlying_file_name;
    char                    m_delimiter;

    virtual int         number_of_snps() const {
        return m_each_SNP_and_its_z.size();
    }
    virtual std::string get_SNPname        (int i)     const {
        auto const & ols = get_gls(i);
        return ols.m_SNPname;
    }
    virtual chrpos      get_chrpos        (int i)     const {
        auto const & ols = get_gls(i);
        return ols.m_chrpos;
    }
    virtual std::string get_allele_ref        (int i)     const {
        auto const & ols = get_gls(i);
        return ols.m_allele_ref;
    }
    virtual std::string get_allele_alt        (int i)     const {
        auto const & ols = get_gls(i);
        return ols.m_allele_alt;
    }
    virtual double      get_z                 (int i)     const {
        auto const & ols = get_gls(i);
        return ols.m_z;
    }
    virtual void        set_chrpos        (int i, chrpos crps)  {
        assert(i>=0);
        assert(i<number_of_snps());
        auto & existing = m_each_SNP_and_its_z.at(i);
        assert(existing.m_chrpos.chr == -1);
        assert(existing.m_chrpos.pos == -1);
        existing.m_chrpos = crps;
    }
    virtual void        sort_my_entries   () {
        sort( m_each_SNP_and_its_z.begin()
            , m_each_SNP_and_its_z.end()
            , [](auto &l, auto &r) {
                return l.m_chrpos < r.m_chrpos;
            }
            );
    }

    GwasLineSummary  get_gls         (int i)     const {
        assert(i>=0);
        assert(i<number_of_snps());
        return m_each_SNP_and_its_z.at(i);
    }

};

static
GwasFileHandle_NONCONST      read_in_a_gwas_file_simple(std:: string file_name) {
    PP(file_name);

    gz:: igzstream f(file_name.c_str());
    (f.rdbuf() && f.rdbuf()->is_open()) || DIE("Can't find file [" << file_name << ']');

    string current_line;

    header_details hd;

    getline(f, current_line);

    // current_line is now the first line, i.e. the header
    hd = parse_header(current_line);

    auto p = std:: make_shared<SimpleGwasFile>();
    p->m_underlying_file_name = file_name;
    p->m_delimiter            = hd.m_delimiter;

    while(1) {
        GwasLineSummary gls;
        getline(f, current_line);
        if(!f) {
            f.eof() || DIE("Error before reaching eof() in this file");
            break;
        }
        auto all_split_up = tokenize(current_line, hd.m_delimiter);
        try {
            gls.m_SNPname    =                           LOOKUP(hd, SNPname, all_split_up);
            gls.m_allele_alt =                           LOOKUP(hd, allele_alt, all_split_up);
            gls.m_allele_ref =                           LOOKUP(hd, allele_ref, all_split_up);
            gls.m_z          = utils:: lexical_cast<double> (LOOKUP( hd, effect_z, all_split_up));

            // Try to read in chromosome and position, but they may not be present
            // If missing, we'll fill them in much later from the reference panel
            gls.m_chrpos.chr = -1;
            gls.m_chrpos.pos = -1;
            if(hd.chromosome.m_offset != -1) {
                gls.m_chrpos.chr = utils:: lexical_cast<int>( LOOKUP(hd, chromosome, all_split_up) );
            }
            if(hd.position.m_offset != -1) {
                gls.m_chrpos.pos = utils:: lexical_cast<int>( LOOKUP(hd, position, all_split_up) );
            }

            p->m_each_SNP_and_its_z.push_back(gls);
        } catch (std:: invalid_argument &e) {
            WARNING( "Ignoring this line, problem with the chromosome and/or position ["
                    << "SNPname:" << LOOKUP(hd, SNPname   , all_split_up)
                    //<< " chr:" << LOOKUP(hd, chromosome, all_split_up)
                    //<< " pos:" << LOOKUP(hd, position  , all_split_up)
                    << "], "
                    << e.what()
                    );
        }
    };

    /* Note: Sort them by position here */

    return p;
}
using utils:: operator<<;
vector<int> CacheOfRefPanelData :: lookup_one_chr_pos(chrpos crps) {

    if(m_cache_of_z12.count(crps) == 1) {
        return m_cache_of_z12[crps];
    }

    auto const b_ref  = begin_from_file(m_rfh);
    auto const e_ref  =   end_from_file(m_rfh);
    auto const    it  = std:: lower_bound(b_ref, e_ref, crps);
    assert(it.get_chrpos() == crps);
    auto pv = it.get_calls();
    assert(!pv.first.empty());
    assert(!pv.second.empty());

    auto has_more_than_one_alt_allele = it.get_allele_alt().find(',') != std::string::npos;
    int fst_max = *max_element(pv.first .begin(), pv.first .end());
    int snd_max = *max_element(pv.second.begin(), pv.second.end());
    if(fst_max>1 || snd_max>1) {
        assert(has_more_than_one_alt_allele);
        return {};
    }

    assert(fst_max <= 1);
    assert(snd_max <= 1);

    assert(pv.first.size() == pv.second.size());
    assert(!has_more_than_one_alt_allele);

    int const N = pv.first.size();

    vector<int> z12;
    for(int i=0; i<N; ++i) {
        z12.push_back( pv.first.at(i) + pv.second.at(i) );
    }

    auto max_z12 = *max_element(z12.begin(), z12.end());
    assert(max_z12 == 2
            || max_z12 == 1
            || max_z12 == 0
            );

    assert(m_cache_of_z12.count(crps) == 0);
    m_cache_of_z12[crps] = z12;
    assert(m_cache_of_z12.count(crps) == 1);
    return lookup_one_chr_pos(crps);
}
vector<int> CacheOfRefPanelData:: lookup_one_ref_get_calls(SNPiterator<GenotypeFileHandle> it) {
    auto cache_key = it.m_line_number;
    if(        m_cache_of_z12_line_number.count(cache_key)==1) {
        return m_cache_of_z12_line_number.at(cache_key);
    }


    m_cache_of_z12_line_number.insert(
                make_pair(cache_key, actual_lookup_one_ref_get_calls(it))
            );

    assert(m_cache_of_z12_line_number.count(cache_key)==1);
    return lookup_one_ref_get_calls(it);
}
static
vector<int>                actual_lookup_one_ref_get_calls(SNPiterator<GenotypeFileHandle> it) {

    auto pv = it.get_calls();
    assert(!pv.first.empty());
    assert(!pv.second.empty());

    auto has_more_than_one_alt_allele = it.get_allele_alt().find(',') != std::string::npos;
    int fst_max = *max_element(pv.first .begin(), pv.first .end());
    int snd_max = *max_element(pv.second.begin(), pv.second.end());
    if(fst_max>1 || snd_max>1) {
        assert(has_more_than_one_alt_allele);
        return {};
    }

    assert(fst_max <= 1);
    assert(snd_max <= 1);

    assert(pv.first.size() == pv.second.size());
    assert(!has_more_than_one_alt_allele);

    int const N = pv.first.size();

    vector<int> z12;
    for(int i=0; i<N; ++i) {
        z12.push_back( pv.first.at(i) + pv.second.at(i) );
    }

    auto max_z12 = *max_element(z12.begin(), z12.end());
    assert(max_z12 == 2
            || max_z12 == 1
            || max_z12 == 0
            );

    return z12;
}

} // namespace file_reading
