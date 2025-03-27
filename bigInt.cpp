//SELF MAKE BIGINT LIBRARY WITH VECTOR<UINT64_T>
//NAIVE MULTIPLICATION
//NO SIGNED OPERATOR

//LE PHUOC THANH
//GITHUB REPOSITORIES (FUNCTION VERSION): https://github.com/TLAquarius/vector-uint64_t-BigPrime-Checking

#include "bigInt.h"

bigInt::bigInt() :chunk(0, 0) {}

bigInt::bigInt(const uint64_t& num) {
    chunk.clear();
    chunk.push_back(num);
}

bigInt::bigInt(const string& num) {
    string temp = num;
    chunk.clear();
    bool isDec = true;
    if (!num.empty())
    {
        for (int i = 0; i < num.size(); i++)
        {
            if (num[i] < 47 || num[i] > 58 || num[i] == 'x')
            {
                isDec = false;
                break;
            }
        }
        if (isDec)
        {
            temp = decToHex(num);
        }
        int start = temp.length();
        int end = start;
        while (end > 0)
        {
            if (end > 15)
                start = end - 16;
            else
                start = 0;

            chunk.push_back(hexToBinary(temp, start, end));
            end = start;
        }
    }
    this->normalize();
}

bigInt::bigInt(const bigInt& other) {
    chunk = other.chunk;
    this->normalize();
}

void bigInt::operator= (const bigInt& other) {

    if (this != &other) {
        chunk = other.chunk;
    }
    this->normalize();
}

void bigInt::operator=(const uint64_t& other) {
    chunk.clear();
    chunk.push_back(other);
    this->normalize();
}

void bigInt::operator=(const string& num) {
    chunk.clear();
    string temp = num;
    bool isDec = true;
    if (!num.empty())
    {
        for (int i = 0; i < num.size(); i++)
        {
            if (num[i] < 47 || num[i] > 58 || num[i] == 'x')
            {
                isDec = false;
                break;
            }
        }
        if (isDec)
        {
            temp = decToHex(num);
        }
        int start = temp.length();
        int end = start;
        while (end > 0)
        {
            if (end > 15)
                start = end - 16;
            else
                start = 0;

            chunk.push_back(hexToBinary(temp, start, end));
            end = start;
        }
    }
    this->normalize();
}

int bigInt::chunk_size() const {
    return int(chunk.size());
}

bool bigInt::operator== (const bigInt& other) const {
    return chunk == other.chunk;
}

bool bigInt::operator!= (const bigInt& other) const {
    return chunk != other.chunk;
}

bool bigInt::operator< (const bigInt& other) const {
    int aSize = chunk.size();
    int otherSize = other.chunk.size();

    if (aSize != otherSize) {
        return aSize < otherSize;
    }

    for (int i = aSize - 1; i >= 0; --i) {
        if (chunk[i] != other.chunk[i])
            return chunk[i] < other.chunk[i];
    }
    return false;
}

bool bigInt::operator<= (const bigInt& other) const {
    return (*this < other) || (*this == other);
}

bool bigInt::operator> (const bigInt& other) const {
    return !(*this <= other);
}

bool bigInt::operator>= (const bigInt& other) const {
    return !(*this < other);
}

uint64_t& bigInt::operator[](size_t index) {
    if (index >= chunk.size()) {
        throw std::out_of_range("Index out of range in BigInt");
    }
    return chunk[index];
}

const uint64_t& bigInt::operator[](size_t index) const {
    if (index >= chunk.size()) {
        throw std::out_of_range("Index out of range in BigInt");
    }
    return chunk[index];
}

bigInt bigInt::operator+(const bigInt& other) const
{
    int aSize = chunk.size();
    int bSize = other.chunk.size();
    int maxSize = max(aSize, bSize);

    bigInt result;
    uint64_t carry = 0;

    for (int i = 0; i < maxSize; ++i) {
        uint64_t chunkA = (i < aSize) ? chunk[i] : 0;
        uint64_t chunkB = (i < bSize) ? other.chunk[i] : 0;

        uint64_t sum = chunkA + chunkB + carry;
        result.chunk.push_back(sum);

        carry = ((chunkA + carry < chunkA) || (chunkB + carry < chunkB) || (chunkA + chunkB < chunkA)) ? 1 : 0;
    }

    if (carry)
        result.chunk.push_back(1);

    result.normalize();

    return result;
}

bigInt bigInt::operator-(const bigInt& other) const
{
    int aSize = chunk.size();
    int bSize = other.chunk_size();

    if (aSize < bSize)
        return 0;

    if (chunk == other.chunk) {
        return 0;
    }

    bool borrow = 0;
    bigInt result;
    for (int i = 0; i < aSize; i++) {
        uint64_t chunkB = (i < bSize) ? other.chunk[i] : 0;
        result.chunk.push_back((chunk[i] - chunkB) - borrow);
        borrow = ((chunkB + borrow < chunkB) || (chunk[i] - borrow > chunk[i]) || (chunk[i] - chunkB > chunk[i])) ? 1 : 0;
    }

    result.normalize();
    return result;
}

bigInt bigInt::operator*(const bigInt& other) const
{
    int aSize = chunk.size();
    int bSize = other.chunk.size();


    if (aSize == 1)
    {
        if (chunk.back() == 0) {
            return bigInt(0);
        }
        if (chunk.back() == 1) {
            return other;
        }
    }
    else if (bSize == 1)
    {
        if (other.chunk.back() == 0) {
            return bigInt(0);
        }
        if (other.chunk.back() == 1) {
            return *this;
        }
    }

    bigInt result;
    result.chunk.resize(aSize + bSize, 0);

    for (int i = 0; i < aSize; ++i) {
        for (int j = 0; j < bSize; ++j) {
            uint64_t a_low = chunk[i] & 0xFFFFFFFFULL;
            uint64_t a_high = chunk[i] >> 32;
            uint64_t b_low = other.chunk[j] & 0xFFFFFFFFULL;
            uint64_t b_high = other.chunk[j] >> 32;

            // Perform 32-bit multiplications
            uint64_t low_low = a_low * b_low;
            uint64_t low_high = a_low * b_high;
            uint64_t high_low = a_high * b_low;
            uint64_t high_high = a_high * b_high;

            // Combine results
            uint64_t carry = (low_low >> 32) + (low_high & 0xFFFFFFFF) + (high_low & 0xFFFFFFFF);
            uint64_t low = (low_low & 0xFFFFFFFF) | (carry << 32);
            uint64_t high = high_high + (low_high >> 32) + (high_low >> 32) + (carry >> 32);

            result[i + j] += low;
            carry = (result[i + j] < low) ? 1 : 0;
            result[i + j + 1] += high + carry;
            carry = result[i + j + 1] < (high + carry) ? 1 : 0;
            int z = i + j + 2;
            while (carry)
            {
                result[z] += 1;
                carry = (result[z] == 0) ? 1 : 0;
                ++z;
            }
        }
    }
    result.normalize();

    return result;
}

bigInt bigInt::operator%(const bigInt& other) const {

    if (*this < other)
        return *this;
    if (*this == other)
        return 0;
    bigInt remainder;
    remainder = *this;

    int otherLength = other.bitLength();
    int remainLength = remainder.bitLength();

    while (remainLength > otherLength) {
        int bitShift = remainLength - otherLength - 1;
        remainder = remainder - (other << bitShift);
        remainLength = remainder.bitLength();
    }

    if (remainder >= other)
        remainder = remainder - other;

    return remainder;
}

bigInt bigInt::operator<<(const int& bit) const
{
    bigInt result = *this;
    result.normalize();
    if (bit == 0)
        return result;
    uint64_t numChunksToShift = bit / 64;
    uint64_t bitShift = bit % 64;

    int size = result.chunk.size();


    if (bitShift == 0) {
        result.chunk.insert(result.chunk.begin(), numChunksToShift, 0);
        return result;
    }

    uint64_t carry = 0;
    uint64_t newCarry = 0;
    for (size_t i = 0; i < size; ++i) {
        newCarry = result.chunk[i] >> (64 - bitShift);
        result.chunk[i] = (result.chunk[i] << bitShift) | carry;
        carry = newCarry;
    }

    if (carry)
        result.chunk.push_back(carry);
    result.chunk.insert(result.chunk.begin(), numChunksToShift, 0);
    return result;
}

bigInt bigInt::operator>>(const int& bit) const
{
    bigInt result = *this;
    result.normalize();
    if (bit == 0)
        return result;
    uint64_t numChunksToShift = bit / 64;
    uint64_t bitShift = bit % 64;

    result.chunk.erase(result.chunk.begin(), result.chunk.begin() + numChunksToShift);
    if (bitShift == 0) {
        return result;
    }
    int size = result.chunk.size();

    uint64_t carry = 0;
    uint64_t newCarry = 0;
    for (int i = size - 1; i >= 0; --i) {
        newCarry = result.chunk[i] << (64 - bitShift);
        result.chunk[i] = (result.chunk[i] >> bitShift) | carry;
        carry = newCarry;
    }

    result.normalize();

    return result;
}

bigInt bigInt::operator&(const bigInt& other) const {
    bigInt result;

    int aSize = chunk.size();
    int bSize = other.chunk.size();
    int minSize = min(aSize, bSize);

    for (int i = 0; i < minSize; i++)
    {
        result.chunk.push_back((*this)[i] & other[i]);
    }
    result.normalize();
    return result;
}

bigInt bigInt::operator|(const bigInt& other) const {
    bigInt result;

    int aSize = chunk.size();
    int bSize = other.chunk.size();
    int maxSize = max(aSize, bSize);

    for (int i = 0; i < maxSize; i++)
    {
        result.chunk.push_back((i < aSize) ? (*this)[i] : 0);
        result.chunk[i] |= (i < bSize) ? other[i] : 0;
    }
    result.normalize();
    return result;
}

bigInt bigInt::operator^(const bigInt& other) const {
    bigInt result;

    int aSize = chunk.size();
    int bSize = other.chunk.size();
    int maxSize = max(aSize, bSize);

    for (int i = 0; i < maxSize; i++)
    {
        result.chunk.push_back((i < aSize) ? (*this)[i] : 0);
        result.chunk[i] ^= (i < bSize) ? other[i] : 0;
    }
    result.normalize();
    return result;
}

bigInt bigInt::operator~() const {
    bigInt result = *this;

    for (int i = 0; i < result.chunk.size(); i++)
    {
        result.chunk[i] = ~result.chunk[i];
    }
    result.normalize();
    return result;
}

const int bigInt::bitLength() const
{
    int pos = 64;
    for (int i = chunk.size() - 1; i >= 0; --i)
    {
        if (chunk[i] != 0) {
            while (((chunk[i] >> --pos) & 1) == 0) {}
            return (64 * i + pos + 1);
        }
    }
    return 0;
}

const bool bigInt::getBit(const int& pos) const
{
    int chunkNum = pos / 64;
    int bitPos = pos % 64;

    return (chunkNum < chunk.size() && ((chunk[chunkNum] >> bitPos) & 1));
}

void bigInt::setBit(const int& pos)
{
    int chunkNum = pos / 64;
    int bitPos = pos % 64;

    chunk[chunkNum] |= (uint64_t(1) << bitPos);
}

int bigInt::leastBit() const
{
    int pos = -1;
    int chunkNum = 0;
    int aSize = chunk.size();
    while (chunkNum < aSize)
    {
        if (chunk[chunkNum] != 0)
        {
            while (((chunk[chunkNum] >> ++pos) & 1) == 0) {}
            return (64 * chunkNum + pos);
            break;
        }
        chunkNum++;
    }
    return -1;
}

void bigInt::normalize()
{
    while (this->chunk.size() > 1 && this->chunk.back() == 0)
        this->chunk.pop_back();
}

uint64_t bigInt::hexToBinary(const string& hex, int start, int end)
{
    uint64_t result = 0;
    for (int i = start; i < end; ++i)
    {
        int value = 0;
        if ('0' <= hex[i] && hex[i] <= '9')
            value = hex[i] - '0';
        else if ('A' <= hex[i] && hex[i] <= 'F')
            value = hex[i] - 'A' + 10;
        else if ('a' <= hex[i] && hex[i] <= 'f')
            value = hex[i] - 'a' + 10;

        result <<= 4;
        result |= value;
    }
    return result;
}

string bigInt::decToHex(const string& dec) {
    string hexResult;
    string num = dec;
    const uint64_t BASE = 16;

    while (!num.empty() && num != "0") {
        uint64_t remainder = 0;
        std::string nextNum;

        for (char digit : num) {
            uint64_t currentDigit = digit - '0';
            uint64_t newValue = remainder * 10 + currentDigit;
            if (!nextNum.empty() || newValue >= BASE)
                nextNum += (newValue / BASE) + '0';
            remainder = newValue % BASE;
        }

        // Convert remainder to hex
        hexResult += "0123456789ABCDEF"[remainder];

        // Remove leading zeros
        num = nextNum.empty() ? "0" : nextNum;
    }

    reverse(hexResult.begin(), hexResult.end());
    return hexResult;
}

string bigInt::bigIntToHex()
{
    string result = "";
    string subString = "";
    int size = chunk.size();
    for (int i = 0; i < size; ++i)
    {
        stringstream stream;
        if (i == size - 1)
            stream << hex << uppercase << chunk[i];
        else
            stream << hex << uppercase << setfill('0') << setw(16) << chunk[i];
        result = stream.str() + result;
    }
    return result;
}