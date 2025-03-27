//SELF MAKE BIGINT LIBRARY WITH VECTOR<UINT64_T>
//NAIVE MULTIPLICATION
//NO SIGNED OPERATOR

//LE PHUOC THANH
//GITHUB REPOSITORIES (FUNCTION VERSION): https://github.com/TLAquarius/vector-uint64_t-BigPrime-Checking

#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;

class bigInt
{
private:
    uint64_t hexToBinary(const string& hex, int start, int end);
    string decToHex(const string& dec);

public:
    vector<uint64_t> chunk;

    void normalize();
    bigInt();
    bigInt(const uint64_t& num);
    bigInt(const bigInt& other);
    bigInt(const string& num);

    void operator=(const bigInt& other);
    void operator=(const uint64_t& other);
    void operator=(const string& other);
    bool operator== (const bigInt& other) const;
    bool operator!= (const bigInt& other) const;
    bool operator< (const bigInt& other) const;
    bool operator> (const bigInt& other) const;
    bool operator<= (const bigInt& other) const;
    bool operator>= (const bigInt& other) const;

    uint64_t& operator[](size_t index);
    const uint64_t& operator[](size_t index) const;

    bigInt operator<< (const int& bit) const;
    bigInt operator>> (const int& bit) const;

    bigInt operator+(const bigInt& other) const;
    bigInt operator-(const bigInt& other) const;
    bigInt operator*(const bigInt& other) const;
    bigInt operator%(const bigInt& other) const;

    bigInt operator&(const bigInt& other) const;
    bigInt operator|(const bigInt& other) const;
    bigInt operator^(const bigInt& other) const;
    bigInt operator~() const;

    const int bitLength() const;
    const bool getBit(const int& pos) const;
    void setBit(const int& pos);
    int chunk_size() const;
    int leastBit() const;
    string bigIntToHex();
};