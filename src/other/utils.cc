#include <vector>
#include <string>

using std:: vector;
using std:: string;

namespace utils {

    std:: vector<std:: string>   tokenize       ( std:: string      const & line
                                            , char                delimiter
        ) {
    int pos = 0;
    vector<string> fields;
    while(1) {
        auto next_delim = line.find(delimiter, pos);
        fields.push_back( line.substr(pos, next_delim-pos) );
        if(next_delim== string:: npos)
            break;
        pos=next_delim+1;
    }
    return fields;
}

} // namespace utils
