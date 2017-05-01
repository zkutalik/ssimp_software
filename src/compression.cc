#include "compression.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include <gzstream.h>

#include "other/DIE.hh"
#include "other/range_action.hh"
#include "other/range_view.hh"
#include "other/PP.hh"

namespace action = range:: action;
namespace view   = range:: view  ;

using std:: ofstream;
using std:: ifstream;
using std:: ostringstream;
using std:: ios_base;
using std:: string;
using std:: vector;
using std:: cout;
using std:: endl;
using std:: unordered_map;

using utils:: operator<<;
using utils:: stdget0;
using utils:: stdget1;
using utils:: print_type;
using utils:: ssize;

using range:: from_vector;

namespace range {
    // put this into the 'range' namespace, to help begin() find it via ADL
    struct read_file_as_a_range_t {
        gz:: igzstream      &   m_f;
        string                  current_line;

        read_file_as_a_range_t(gz:: igzstream & f) : m_f(f) {
            advance();
        }

        bool            empty()     const   { return !m_f; }
        void            advance()           { getline(m_f, current_line); }
        const string&   front_ref() const   { return current_line; }
        string          front_val() const   { return current_line; }
    };
}
namespace compression {
    struct dict_entry {
        string  s;
        int     repetition;
        double  pseudo_count;
        int     code; // to be filled in much later

        string  to_string() const {
            std:: ostringstream oss;
            oss << '(';
            oss << utils::consider_quoting(s);
            oss << ' ' << repetition;
            oss << ',' << pseudo_count;
            oss << ',' << code;
            oss << ')';
            return oss.str();
        }
    };
    static
    auto    decode_from_ints(vector<int> & encoded_as_ints, vector<dict_entry> &dict) {
        vector<string> fields;
        for(int code : encoded_as_ints) {
            assert(code < ssize(dict));
            dict_entry const & d = dict.at(code);
            for(int i : range:: ints(d.repetition)) {
                (void)i;
                fields.push_back(d.s);
            }
        }
        return fields;
    }

    static
    auto build_an_efficient_dictionary_from_a_vector_of_strings(vector<string> const & v) {
        int const N = v.size();
        unordered_map<string, int> frequencies;
        for(string const &s : v) {
            ++frequencies[s];
        }

        vector<dict_entry> dict;

        //PP(N);
        if(frequencies.size() == 1) {
            dict.emplace_back( dict_entry{frequencies.begin()->first, N, double(N), 0} );
            //PP(dict);
            return dict;
        }
        assert(frequencies.size() >= 2);

        range:: range_from_begin_end(frequencies)
        |action::unzip_foreach|
        [&](auto k, auto v) {
            auto rate = double(v)/N;
            auto pseudo_rate = rate;
            //PP(k,v, rate, -std::log2(rate));
            // Entries with very high frequency deserve even more efficient coding
            int repetition = 1;
            while( -std::log2(rate) < 0.5 ) {
                rate = rate * rate;
                pseudo_rate = std:: sqrt(pseudo_rate); // just to ensure it's given high priority in the dictionary
                repetition *= 2;
                //PP(rate, - std::log2(rate) , k, repetition, pseudo_rate*N);
                dict.emplace_back( dict_entry{k, repetition, pseudo_rate*N, -1} );
                assert(repetition <= N);
            }
        };
        //PP(frequencies.size());

        range:: range_from_begin_end(frequencies)
        |action::unzip_foreach|
        [&](auto && k, auto && v) {
            //PP(k,v);
            auto rate = double(v)/N;
            dict.emplace_back( dict_entry{k, 1, rate*N, -1} );
        };

        range:: sort( range:: range_from_begin_end(dict), [&](auto && l, auto && r) {
            if    (l.pseudo_count == r. pseudo_count) {
                if    (l.s == r.s)
                    return l.repetition > r.repetition;
                return l.s > r.s;
            }
            return l.pseudo_count >  r. pseudo_count;
        });
        //for(auto &&d : dict) { PP(d); }

        enumerate(range:: range_from_begin_end(dict) | view:: ref_wraps)
        |action:: unzip_foreach|
        [&](auto &&i , auto && d
                ) {
            d.get().code = i;
        };
        //for(auto &&d : dict) { PP(d); }

        return dict;
    }

    enum class TypesOfCodedOutput : uint8_t {
            code_end_of_file        = 0
        ,   code_null_term_string   = 1 // null-terminated string - when everything else fails

        ,   code_int                = 2
        ,   code_rs_int             = 3
        ,   code_period             = 4 // "." occurs a lot
        ,   code_PASS               = 5 // "PASS" occurs a lot

        ,   code_A                  = 6
        ,   code_T                  = 7
        ,   code_G                  = 8
        ,   code_C                  = 9

        ,   code_dict_start         = 15
        ,   code_0pipe0             = 16
        ,   code_0slsh0             = 17
        ,   code_0pipe1             = 18
        ,   code_0slsh1             = 19
        ,   code_1pipe0             = 20
        ,   code_1slsh0             = 21
        ,   code_1pipe1             = 22
        ,   code_1slsh1             = 23

    };
    struct GTcompressed_output_t {
        ofstream m_f;
        ofstream:: pos_type m_remember_the_begining_position;

        template<typename ...Ts>
        GTcompressed_output_t(Ts&& ...ts)   : m_f( std::forward<Ts>(ts) ... )
                                            , m_remember_the_begining_position( m_f.tellp() )
        {
            // Next two lines to ensure we really are at the start.
            m_f.seekp(0, ios_base::beg);
            assert( m_remember_the_begining_position == m_f.tellp() );
        }

        ~GTcompressed_output_t() { m_f.close(); }

        operator bool() const {
            return !!m_f;
        }

        void output_uint32(uint32_t x) {
            constexpr size_t w = sizeof(x);
            static_assert( w == 4 ,"");
            vector<uint8_t> bytes(w);
            for(int i=w-1; i>=0; --i) {
                bytes.at(i) = x % 256;
                x = x / 256;
            }
            assert(x==0);
            for(int i=0; i<(int)w; ++i) {
                m_f << bytes.at(i);
            }
        }
        void output_uint64(uint64_t x) {
            static_assert( 8 == sizeof(uint64_t) ,"");
            vector<uint8_t> bytes(8);
            for(int i=7; i>=0; --i) {
                bytes.at(i) = x % 256;
                x = x / 256;
            }
            assert(x==0);
            for(int i=0; i<8; ++i) {
                m_f << bytes.at(i);
            }
        }
        void output_string(string const & s) {
            for(auto c : s) {
                static_assert( std:: is_same<decltype(c), char>{} ,"");
                assert(c != '\0');
            }
            m_f << s;
            m_f << '\0';
        }
        void output_smart_string(string const & s) {
            for(auto c : s) {
                static_assert( std:: is_same<decltype(c), char>{} ,"");
                assert(c != '\0');
            }
            try {
                static_assert( sizeof(int) == sizeof(int32_t) ,"");
                int as_int = utils:: lexical_cast<int>(s);
                assert(as_int >= 0); // TODO: Could fail sometime
                output_code(TypesOfCodedOutput:: code_int);
                output_uint32(as_int);
                return;
            } catch (std:: invalid_argument & e) { }
            try {
                if(s.substr(0,2) == "rs") {
                    string s_skip_first_two = s.substr(2, string:: npos);
                    static_assert( sizeof(int) == sizeof(int32_t) ,"");
                    int as_int = utils:: lexical_cast<int>(s_skip_first_two);
                    assert(as_int >= 0); // TODO: Could fail sometime
                    output_code(TypesOfCodedOutput:: code_rs_int);
                    output_uint32(as_int);
                    return;
                }
            } catch (std:: invalid_argument & e) { }
            if(s == "."      ) { output_code(TypesOfCodedOutput:: code_period); return; }
            if(s == "PASS"   ) { output_code(TypesOfCodedOutput:: code_PASS  ); return; }
            if(s == "0|0"   ) { output_code(TypesOfCodedOutput:: code_0pipe0  ); return; }
            if(s == "0/0"   ) { output_code(TypesOfCodedOutput:: code_0slsh0  ); return; }
            if(s == "0|1"   ) { output_code(TypesOfCodedOutput:: code_0pipe1  ); return; }
            if(s == "0/1"   ) { output_code(TypesOfCodedOutput:: code_0slsh1  ); return; }
            if(s == "1|0"   ) { output_code(TypesOfCodedOutput:: code_1pipe0  ); return; }
            if(s == "1/0"   ) { output_code(TypesOfCodedOutput:: code_1slsh0  ); return; }
            if(s == "1|1"   ) { output_code(TypesOfCodedOutput:: code_1pipe1  ); return; }
            if(s == "1/1"   ) { output_code(TypesOfCodedOutput:: code_1slsh1  ); return; }
            if(s == "A"      ) { output_code(TypesOfCodedOutput:: code_A     ); return; }
            if(s == "T"      ) { output_code(TypesOfCodedOutput:: code_T     ); return; }
            if(s == "G"      ) { output_code(TypesOfCodedOutput:: code_G     ); return; }
            if(s == "C"      ) { output_code(TypesOfCodedOutput:: code_C     ); return; }
            // Otherwise, we just need to output a plain null-terminated string
            output_code(TypesOfCodedOutput:: code_null_term_string);
            output_string(s);
        }
        void output_dict(vector<dict_entry>     const   & dict) {
            output_code(TypesOfCodedOutput:: code_dict_start);
            for(dict_entry const & d : dict) {
                assert(d.repetition > 0); // repetition==0 will mark the end
                output_uint32(d.repetition);
                output_smart_string(d.s);
            }
            output_uint32(0); // end of dictionary
        }
        void output_encoded_as_ints(vector<int> const   & encoded_as_ints) {
            for(int i : encoded_as_ints) {
                assert(i>=0);
                output_uint32(i);
            }
        }
        void output_encoded_as_bits(vector<bool> encoded_as_bits
                ) {
            while(encoded_as_bits.size() % 8 != 0)
                encoded_as_bits.push_back(true);
            assert(encoded_as_bits.size() % 8 == 0);
            int num_bytes = encoded_as_bits.size() / 8;
            for(int byte_counter = 0; byte_counter < num_bytes; ++byte_counter) {
                int this_byte_as_an_int = 0;
                for(int bit_counter = 0; bit_counter < 8; ++bit_counter) {
                    if( encoded_as_bits.at(8*byte_counter + bit_counter) ) {
                        //PP( bit_counter, 1 << bit_counter );
                        this_byte_as_an_int += 1 << bit_counter;
                    }
                }
                //PP(this_byte_as_an_int);
                output_uint8(this_byte_as_an_int);
            }
        }
        void output_uint8(uint8_t u) {
                m_f << (uint8_t)u;
        }
        void output_code(TypesOfCodedOutput code) {
                m_f << (uint8_t)code;
        }
    };

    struct binary_reader_t {
        ifstream m_f;

        binary_reader_t(string compressed_out_file_name)
        {
            m_f.open(   compressed_out_file_name
                    //,   ios_base::in
                       //|ios_base::binary
                    );
            m_f || DIE("can't reopen");
        }

        string              read_smart_string() {
            auto code = read_code();
            ostringstream oss;
            switch(code) {
                break; case TypesOfCodedOutput:: code_int: {
                    auto u = read_uint32();
                    //cout << std:: hex << u << endl;
                    oss << u;
                    return oss.str();
                }
                break; case TypesOfCodedOutput:: code_rs_int: {
                    auto u = read_uint32();
                    //cout << std:: hex << u << endl;
                    oss << "rs" << u;
                    return oss.str();
                }
                break; case TypesOfCodedOutput:: code_A     : { return "A"; }
                break; case TypesOfCodedOutput:: code_T     : { return "T"; }
                break; case TypesOfCodedOutput:: code_G     : { return "G"; }
                break; case TypesOfCodedOutput:: code_C     : { return "C"; }
                break; case TypesOfCodedOutput:: code_PASS  : { return "PASS"; }
                break; case TypesOfCodedOutput:: code_0pipe0: { return "0|0"; }
                break; case TypesOfCodedOutput:: code_0slsh0: { return "0/0"; }
                break; case TypesOfCodedOutput:: code_1pipe0: { return "1|0"; }
                break; case TypesOfCodedOutput:: code_1slsh0: { return "1/0"; }
                break; case TypesOfCodedOutput:: code_0pipe1: { return "0|1"; }
                break; case TypesOfCodedOutput:: code_0slsh1: { return "0/1"; }
                break; case TypesOfCodedOutput:: code_1pipe1: { return "1|1"; }
                break; case TypesOfCodedOutput:: code_1slsh1: { return "1/1"; }
                break; case TypesOfCodedOutput:: code_period: { return "."; }
                break; case TypesOfCodedOutput:: code_null_term_string: {
                    string s;
                    getline(m_f, s, '\0');
                    return s;
                }
                break; default: {
                    DIE("unhandled TypesOfCodedOutput [" << (int) code << "]");
                }
            }
            return "?";
        }
        uint32_t read_uint8() {
            uint8_t c = (char) m_f.get();
            return c;
        }
        uint32_t read_uint32() {
            uint32_t x = 0;
            assert(m_f);

            for(int i = 0; i<4; ++i) {
                uint8_t c = (char) m_f.get(); // int -> char -> uint8_t . I think this is the only correct way
                assert(m_f);
                x = 256 * x;
                x += c;
            }
            return x;
        }

        vector<string>  read_dict_and_fully_decode() {
            auto starter = read_code();
            assert(starter == TypesOfCodedOutput:: code_dict_start);

            vector<dict_entry> dict;
            do {
                dict_entry d;
                d.pseudo_count = std::nan("");
                d.code = -1;
                d.repetition = read_uint32();
                if(d.repetition == 0)
                    break;
                else {
                    assert(d.repetition>0);
                    d.s = read_smart_string();
                    dict.push_back(d);
                }
            } while(1);

            // now to read a sequence of bits. Needs care
            struct bit_reader_t {
                int  offset = 7;
                uint8_t current_byte;
                binary_reader_t * m_outer_this;

                bit_reader_t(binary_reader_t * outer_this) : m_outer_this(outer_this) {
                    advance();
                }

                void    advance() {
                    ++offset;
                    if(offset == 8) {
                        offset = 0;
                        current_byte = m_outer_this->read_uint8();
                    }
                }
                bool    get_bit() const {
                    return (current_byte & (1 << offset)) > 0;
                }
                int     get_count_up_to_next_true() {
                    int code_from_bits = 0;
                    while(get_bit() == false) {
                        ++code_from_bits;
                        advance();
                    }
                    assert(get_bit() == true);
                    advance();
                    return code_from_bits;
                }
            } bit_reader{this};

            vector<int> ints_from_bits;
            do {
                int code_from_bits = bit_reader.get_count_up_to_next_true();
                assert(code_from_bits <= ssize(dict));
                if(code_from_bits ==  ssize(dict))
                    break;
                ints_from_bits.push_back(code_from_bits);
                assert(code_from_bits <  ssize(dict));
            }while(1);

            vector<string> decoded;
            for(int code : ints_from_bits) {
                assert(code < ssize(dict));
                dict_entry & d = dict.at(code);
                for(int i = 0 ; i<d.repetition; ++i) {
                    decoded.push_back( d.s );
                }
            }
            return decoded;
        }

        TypesOfCodedOutput  read_code() {
            assert(m_f);
            int i = m_f.get();
            assert(m_f);
            uint8_t x= (char) i;
            assert(m_f);
            return static_cast<TypesOfCodedOutput>(x);
        }
    };

    static
    auto        read_in_one_full_compressed_line( int64_t remember_offset_at_start_of_this_line
            , string compressed_out_file_name
            ) {
        binary_reader_t reopen( compressed_out_file_name );
        reopen.m_f.seekg( remember_offset_at_start_of_this_line
                        , ios_base::beg
                        );
        reopen.m_f || DIE("can't seek the reopen");

        //auto c0 = reopen.read_code();
        //PP( (int)c0 );
        struct stuff_t {
            string field0;
            string field1;
            string field2;
            string field3;
            string field4;
            string field5;
            string field6;
            vector<string> fully_decoded_strings;
        };
        return stuff_t  { reopen.read_smart_string()
                        , reopen.read_smart_string()
                        , reopen.read_smart_string()
                        , reopen.read_smart_string()
                        , reopen.read_smart_string()
                        , reopen.read_smart_string()
                        , reopen.read_smart_string()
                        , reopen.read_dict_and_fully_decode()
        };
    }

    void make_compressed_vcf_file   (std:: string compressed_out_file_name, std:: string vcf_filename) {
        GTcompressed_output_t binary_output (compressed_out_file_name
                                            ,   ios_base::out
                                               |ios_base::binary
                                               |ios_base::trunc
                                            );
        (  binary_output) || DIE("Couldn't open [" << compressed_out_file_name.c_str()
                << "] for writing.");

        gz:: igzstream f(vcf_filename.c_str());
        (f.rdbuf() && f.rdbuf()->is_open()) || DIE("Can't find file [" << vcf_filename << ']');

        binary_output.output_string( "vcf.GT.compressed.0.0.1" ); // later, this will be the offset just past the end
        binary_output.output_uint64( 1234 ); // later, this will be the offset just past the end

        range:: read_file_as_a_range_t vcf_input_range{f};
        auto z = zip( range::ints<int64_t>(), vcf_input_range );

        constexpr int INFO_field_to_skip = 7;
        constexpr int FORMAT_field_to_skip = 8;
        constexpr int first_GT_field = 9;
        int             number_of_snps = 0;

        z |action:: unzip_foreach|
        [&](auto i, std:: string const & s) {
            if(s.substr(0,2) == "##")
                return;
            auto fields = utils:: tokenize(s, '\t');
            //PP(fields);
            if(fields.at(0) == "#CHROM") {
                // #CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT
                assert( fields.at(INFO_field_to_skip) == "INFO" );
                assert( fields.at(FORMAT_field_to_skip) == "FORMAT" );
                assert((
                std::vector<string> {   fields.at(0)
                                    ,   fields.at(1)
                                    ,   fields.at(2)
                                    ,   fields.at(3)
                                    ,   fields.at(4)
                                    ,   fields.at(5)
                                    ,   fields.at(6)
                                    ,   fields.at(7)
                                    ,   fields.at(8)
                                    } ==
                std::vector<string> { "#CHROM" , "POS" , "ID" , "REF", "ALT", "QUAL", "FILTER" , "INFO", "FORMAT"
                                    }));
                return;
            }
            ++ number_of_snps;
            (void)i;
            auto remember_offset_at_start_of_this_line = binary_output.m_f.tellp() - binary_output.m_remember_the_begining_position;
            //PP( i, s.substr(0,90) );
            for(int j=0; j<9; ++j) {
                if(j==INFO_field_to_skip)
                    continue;
                if(j==FORMAT_field_to_skip)
                    continue;
                binary_output.output_smart_string( fields.at(j) );
            }

            /* Next, even though the FORMAT field isn't explicitly encoded, we
             * do need to find 'GT' in it
             */
            auto FORMAT = fields.at(FORMAT_field_to_skip);
            auto where_is_GT =
                view:: enumerate_vector(utils:: tokenize(FORMAT,':'))
                |view:: unzip_filter|
                [](auto /*index*/, auto fmt) {
                        return "GT" == fmt;
                }
                | action:: collect;
            where_is_GT.size() == 1 || DIE("The FORMAT field should have 'GT' exactly once [" << FORMAT << "]");

            int const offset_of_GT_within_FORMAT_field = where_is_GT.at(0) |stdget0;
            assert("GT" ==                              (where_is_GT.at(0) |stdget1));

            auto many_call_pairs_as_strings =
                range:: ints( first_GT_field, fields.size() )
                |view:: map|
                [&](auto index) {
                    return utils::tokenize(fields.at(index),':').at(offset_of_GT_within_FORMAT_field);
                }
                |action:: collect
                ;
            auto dict = build_an_efficient_dictionary_from_a_vector_of_strings(many_call_pairs_as_strings);

            // Now to find long runs of identical things within 'many_call_pairs_as_strings'
            auto r = from_vector(many_call_pairs_as_strings);
            vector<int> encoded_as_ints;
            while(!r.empty()) {
                auto s = front_val(r);
                int  length_of_run = 1;
                r.advance();
                while(!r.empty() && s == front_val(r)) {
                    ++ length_of_run;
                    r.advance();
                }
                //PP(s, length_of_run);
                // To find a suitable set of codings for this
some_more:
                for(dict_entry & d : dict) {
                    if(s == d.s && length_of_run >= d.repetition) {
                        //PP(d);
                        // just use this first one
                        encoded_as_ints.push_back(d.code);
                        assert( &d == &dict.at(d.code) );
                        length_of_run -= d.repetition;
                        assert(length_of_run >= 0);
                        goto some_more;
                    }
                }
                assert(length_of_run==0);
            }

            auto decoded_from_ints = decode_from_ints(encoded_as_ints, dict);
            assert(many_call_pairs_as_strings.size() == decoded_from_ints.size());
            assert(many_call_pairs_as_strings        == decoded_from_ints       );

            binary_output.output_dict(dict);
            encoded_as_ints.push_back( dict.size() );

            vector<bool> encoded_as_bits;
            for(int code : encoded_as_ints) {
                for(int i=0; i<code; ++i) {
                    encoded_as_bits.push_back(false);
                }
                encoded_as_bits.push_back(true);
            }

            // This next line is quite important, to ensure that (when reading)
            // that we don't 'advance' into one byte too far
            encoded_as_bits.push_back(true);

            binary_output.output_encoded_as_bits(std::move(encoded_as_bits));

            binary_output.m_f.flush();

            { // Let's reread everything to check it's OK
                auto stuff = read_in_one_full_compressed_line(   remember_offset_at_start_of_this_line
                                                             ,   compressed_out_file_name);
                assert(stuff.field0 == fields.at(0));
                assert(stuff.field1 == fields.at(1));
                assert(stuff.field2 == fields.at(2));
                assert(stuff.field3 == fields.at(3));
                assert(stuff.field4 == fields.at(4));
                assert(stuff.field5 == fields.at(5));
                assert(stuff.field6 == fields.at(6));
                assert(stuff.fully_decoded_strings == many_call_pairs_as_strings);
            }
        };
        PP( number_of_snps,  z.front_val() | stdget0);

    }
}
