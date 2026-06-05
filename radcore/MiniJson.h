// ============================================================================
// MiniJson.h — tiny dependency-free JSON reader for RADiCAL build configs
// ----------------------------------------------------------------------------
// A compact recursive-descent parser sufficient for the config files we author
// ourselves (objects, arrays, strings, numbers, bools, null, nesting). No
// external dependencies, ACLiC-compatible (root -l 'macro.C+'), C++11.
//
// Usage:
//   mj::Value cfg = mj::parseFile("datasets/2023/configs/DSB1.json");
//   std::string build = cfg["build"].str();
//   int drs = cfg["channel_map"]["ends"][0]["hg"][0].asInt();
//   for (auto& e : cfg["channel_map"]["ends"].arr) { ... }
//
// Missing keys / wrong-type access return a shared Null value (never throws on
// access), so callers can probe with .has(key) and .isNull(). Parse errors set
// an error string (mj::parse(text, &err)).
// ============================================================================
#ifndef RADCORE_MINIJSON_H
#define RADCORE_MINIJSON_H

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace mj {

struct Value {
    enum Type { Null, Bool, Num, Str, Arr, Obj };
    Type type = Null;
    bool   b   = false;
    double num = 0.0;
    std::string str;
    std::vector<Value> arr;
    std::map<std::string, Value> obj;

    bool isNull() const { return type == Null; }
    bool isObj()  const { return type == Obj;  }
    bool isArr()  const { return type == Arr;  }

    // object access by key — returns a static Null if absent/not-an-object
    const Value& operator[](const std::string& key) const {
        static const Value kNull;
        if (type != Obj) return kNull;
        auto it = obj.find(key);
        return it == obj.end() ? kNull : it->second;
    }
    // array access by index — returns a static Null if out of range
    const Value& operator[](size_t i) const {
        static const Value kNull;
        if (type != Arr || i >= arr.size()) return kNull;
        return arr[i];
    }
    bool   has(const std::string& key) const { return type == Obj && obj.count(key); }
    size_t size() const { return type == Arr ? arr.size() : (type == Obj ? obj.size() : 0); }

    double      asNum() const { return type == Num ? num : (type == Bool ? (b ? 1.0 : 0.0) : 0.0); }
    int         asInt() const { return (int)std::llround(asNum()); }
    double      asDouble(double def = 0.0) const { return type == Num ? num : def; }
    bool        asBool() const { return type == Bool ? b : (type == Num ? num != 0.0 : false); }
    std::string asStr(const std::string& def = "") const { return type == Str ? str : def; }
};

// ---------------------------------------------------------------------------
class Parser {
public:
    Parser(const std::string& s) : s_(s), i_(0) {}

    Value parse() {
        skipWs();
        Value v = parseValue();
        skipWs();
        if (i_ != s_.size()) fail("trailing characters after JSON value");
        return v;
    }
    const std::string& error() const { return err_; }
    bool ok() const { return err_.empty(); }

private:
    const std::string& s_;
    size_t i_;
    std::string err_;

    void fail(const std::string& m) {
        if (err_.empty()) err_ = m + " at offset " + std::to_string(i_);
        throw std::runtime_error(err_);
    }
    char peek() { return i_ < s_.size() ? s_[i_] : '\0'; }
    char get()  { return i_ < s_.size() ? s_[i_++] : '\0'; }
    void skipWs() {
        while (i_ < s_.size()) {
            char c = s_[i_];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r') { ++i_; continue; }
            if (c == '/' && i_ + 1 < s_.size() && s_[i_ + 1] == '/') {   // // line comment
                i_ += 2; while (i_ < s_.size() && s_[i_] != '\n') ++i_; continue;
            }
            break;
        }
    }

    Value parseValue() {
        skipWs();
        char c = peek();
        switch (c) {
            case '{': return parseObject();
            case '[': return parseArray();
            case '"': { Value v; v.type = Value::Str; v.str = parseString(); return v; }
            case 't': case 'f': return parseBool();
            case 'n': return parseNull();
            default:
                if (c == '-' || (c >= '0' && c <= '9')) return parseNumber();
                fail("unexpected character"); return Value();
        }
    }

    Value parseObject() {
        Value v; v.type = Value::Obj;
        get(); // {
        skipWs();
        if (peek() == '}') { get(); return v; }
        while (true) {
            skipWs();
            if (peek() != '"') fail("expected string key");
            std::string key = parseString();
            skipWs();
            if (get() != ':') fail("expected ':' after key");
            v.obj[key] = parseValue();
            skipWs();
            char c = get();
            if (c == ',') continue;
            if (c == '}') break;
            fail("expected ',' or '}' in object");
        }
        return v;
    }

    Value parseArray() {
        Value v; v.type = Value::Arr;
        get(); // [
        skipWs();
        if (peek() == ']') { get(); return v; }
        while (true) {
            v.arr.push_back(parseValue());
            skipWs();
            char c = get();
            if (c == ',') continue;
            if (c == ']') break;
            fail("expected ',' or ']' in array");
        }
        return v;
    }

    std::string parseString() {
        std::string out;
        get(); // opening "
        while (true) {
            char c = get();
            if (c == '\0') fail("unterminated string");
            if (c == '"') break;
            if (c == '\\') {
                char e = get();
                switch (e) {
                    case '"':  out += '"';  break;
                    case '\\': out += '\\'; break;
                    case '/':  out += '/';  break;
                    case 'n':  out += '\n'; break;
                    case 't':  out += '\t'; break;
                    case 'r':  out += '\r'; break;
                    case 'b':  out += '\b'; break;
                    case 'f':  out += '\f'; break;
                    case 'u': { // \uXXXX -> keep ASCII, approximate others as '?'
                        std::string hex = s_.substr(i_, 4); i_ += 4;
                        long cp = std::strtol(hex.c_str(), nullptr, 16);
                        if (cp < 128) out += (char)cp; else out += '?';
                        break;
                    }
                    default: out += e; break;
                }
            } else {
                out += c;
            }
        }
        return out;
    }

    Value parseNumber() {
        size_t start = i_;
        if (peek() == '-') get();
        while (isdigit((unsigned char)peek())) get();
        if (peek() == '.') { get(); while (isdigit((unsigned char)peek())) get(); }
        if (peek() == 'e' || peek() == 'E') {
            get();
            if (peek() == '+' || peek() == '-') get();
            while (isdigit((unsigned char)peek())) get();
        }
        Value v; v.type = Value::Num;
        v.num = std::strtod(s_.substr(start, i_ - start).c_str(), nullptr);
        return v;
    }

    Value parseBool() {
        Value v; v.type = Value::Bool;
        if (s_.compare(i_, 4, "true") == 0)  { v.b = true;  i_ += 4; }
        else if (s_.compare(i_, 5, "false") == 0) { v.b = false; i_ += 5; }
        else fail("invalid literal");
        return v;
    }

    Value parseNull() {
        if (s_.compare(i_, 4, "null") == 0) i_ += 4; else fail("invalid literal");
        return Value();
    }
};

// ---------------------------------------------------------------------------
inline Value parse(const std::string& text, std::string* errOut = nullptr) {
    Parser p(text);
    try {
        return p.parse();
    } catch (const std::exception& e) {
        if (errOut) *errOut = e.what();
        return Value();
    }
}

inline Value parseFile(const std::string& path, std::string* errOut = nullptr) {
    std::ifstream f(path.c_str());
    if (!f) { if (errOut) *errOut = "cannot open " + path; return Value(); }
    std::stringstream ss; ss << f.rdbuf();
    return parse(ss.str(), errOut);
}

} // namespace mj

#endif // RADCORE_MINIJSON_H
