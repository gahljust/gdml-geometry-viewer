#pragma once
// JSON strings need control-character escaping; std::quoted is only C++ quoting.
#include <string>
namespace viewer_json {
inline std::string quote(const std::string& text) {
    const char* hex="0123456789abcdef";
    std::string out="\"";
    for(unsigned char c:text) {
        if(c=='"'||c=='\\'){out+='\\';out+=c;}
        else if(c<0x20){out+="\\u00";out+=hex[c>>4];out+=hex[c&15];}
        else out+=c;
    }
    return out+'"';
}
}
