#include <random>
#include <vector>
#include <fstream>
#include <string>

using namespace std;
typedef std::vector<uint64_t> bigInt;
std::random_device rd;  // Hardware-based random seed
std::mt19937_64 generator(rd());

uint64_t hexToBinary(const string& hex, int start, int end) {
    uint64_t result = 0;
    for (int i = start; i < end; ++i) {
        char c = hex[i];
        result = (result << 4) | (isdigit(c) ? (c - '0') : (toupper(c) - 'A' + 10));
    }
    return result;
}

void normalize(bigInt& a)
{
    while (a.size() > 1 && a.back() == 0)
        a.pop_back();
}

bool smallerBigInt(const bigInt& a, const bigInt& b)
{
    int aSize = a.size();
    int bSize = b.size();

    if (aSize != bSize)
        return aSize < bSize;

    for (int i = aSize - 1; i >= 0; --i) {
        if (a[i] != b[i])
            return a[i] < b[i];
    }
    return false;
}

bool smallerNormInt(const bigInt& a, const uint64_t& integer)
{
    return (a[0] < integer && a.size() == 1);
}

bool equalNormInt(const bigInt& a, const uint64_t& integer)
{

    return (a[0] == integer && a.size() == 1);

}

bool greaterNormInt(const bigInt& a, const uint64_t& integer)
{
    return (a[0] > integer && a.size() == 1) || a.size() > 1;
}

//hex string -> bigInt
void hexToBigInt(bigInt& a, const string& hex)
{
    if (!hex.empty())
    {
        int start = hex.length();
        int end = start;

        //Each loop will add one chunk to the bigInt
        while (end > 0)
        {
            if (end > 15)
                start = end - 16;
            else
                start = 0;

            a.push_back(hexToBinary(hex, start, end));
            end = start;
        }
    }

    //Remove leading 0 chunks
    normalize(a);
}

int bitLength(const bigInt& a)
{
    int pos = 64;
    for (int i = a.size() - 1; i >= 0; --i)
    {
        if (a[i] != 0) {
            while (((a[i] >> --pos) & 1) == 0) {}
            return (64 * i + pos + 1);
        }
    }
    return 0;
}

int leastBit(const bigInt& a)
{
    int pos = -1;
    int chunk = 0;
    int aSize = a.size();
    while (chunk < aSize)
    {
        if (a[chunk] != 0)
        {
            while (((a[chunk] >> ++pos) & 1) == 0) {}
            return (64 * chunk + pos);
            break;
        }
        chunk++;
    }
    return -1;
}

const bool getBit(const bigInt& a, const int& pos)
{
    int chunkNum = pos / 64;
    int bitPos = pos % 64;

    return (chunkNum < a.size() && ((a[chunkNum] >> bitPos) & 1) != 0);
}

//shift a to left bits
bigInt shiftLeft(const bigInt& a, const int& bit)
{
    bigInt result = a;
    normalize(result);
    if (bit == 0)
        return result;
    int numChunksToShift = bit / 64;
    int bitShift = bit % 64;

    int size = result.size();

    if (bitShift == 0) {
        result.insert(result.begin(), numChunksToShift, 0);
        return result;
    }

    uint64_t carry = 0;
    uint64_t newCarry = 0;
    for (size_t i = 0; i < size; ++i) {
        newCarry = result[i] >> (64 - bitShift);
        result[i] = (result[i] << bitShift) | carry;
        carry = newCarry;
    }

    if (carry)
        result.push_back(carry);
    result.insert(result.begin(), numChunksToShift, 0);
    return result;
}

//shift a to right bits
void shiftRightSelf(bigInt& a, const int& bit)
{
    normalize(a);
    if (bit == 0)
        return;
    int numChunksToShift = bit / 64;
    int bitShift = bit % 64;

    int size = a.size();

    if (bitShift == 0) {
        a.erase(a.begin(), a.begin() + numChunksToShift);
        return;
    }

    uint64_t carry = 0;
    for (int i = size - 1; i >= 0; --i) {
        uint64_t newCarry = a[i] << (64 - bitShift);
        a[i] = (a[i] >> bitShift) | carry;
        carry = newCarry;
    }

    normalize(a);
}

//result = a + b
void sum(bigInt& result, const bigInt& a, const bigInt& b)
{
    int aSize = a.size();
    int bSize = b.size();
    int maxSize = max(aSize, bSize);
    bool carry = 0;

    if (result != a && result != b)
    {
        result.resize(maxSize + 1);
        for (int i = 0; i < maxSize; ++i) {
            uint64_t chunkA = (i < aSize) ? a[i] : 0;
            uint64_t chunkB = (i < bSize) ? b[i] : 0;
            result[i] = chunkA + chunkB + carry;

            carry = ((chunkA + carry < chunkA) || (chunkB + carry < chunkB) || (chunkA + chunkB < chunkA)) ? 1 : 0;
        }
        if (carry)
            result[maxSize] = 1;
    }
    else
    {
        bigInt temp(maxSize + 1);
        for (int i = 0; i < maxSize; ++i) {
            uint64_t chunkA = (i < aSize) ? a[i] : 0;
            uint64_t chunkB = (i < bSize) ? b[i] : 0;
            temp[i] = chunkA + chunkB + carry;

            carry = ((chunkA + carry < chunkA) || (chunkB + carry < chunkB) || (chunkA + chunkB < chunkA)) ? 1 : 0;
        }
        if (carry)
            temp[maxSize] = 1;
        result = temp;
    }
    normalize(result);
}

//result = result - b
void sub(bigInt& result, const bigInt& a, const bigInt& b)
{
    if (a == b) {
        result = { 0 };
        return;
    }
    if (result != a && result != b)
    {
        int aSize = a.size();
        int bSize = b.size();
        bool borrow = 0;
        result.clear();
        for (int i = 0; i < aSize; i++) {
            uint64_t chunkB = (i < bSize) ? b[i] : 0;
            result.push_back((a[i] - chunkB) - borrow);
            borrow = ((chunkB + borrow < chunkB) || (a[i] - borrow > a[i]) || (a[i] - chunkB > a[i])) ? 1 : 0;
        }
    }
    else
    {
        int aSize = a.size();
        int bSize = b.size();
        bool borrow = 0;
        bigInt temp;
        for (int i = 0; i < aSize; i++) {
            uint64_t chunkB = (i < bSize) ? b[i] : 0;
            temp.push_back((a[i] - chunkB) - borrow);
            borrow = ((chunkB + borrow < chunkB) || (a[i] - borrow > a[i]) || (a[i] - chunkB > a[i])) ? 1 : 0;
        }
        result = temp;
    }

    normalize(result);
}

//result = a * b
void mul(bigInt& result, const bigInt& a, const bigInt& b)
{
    int aSize = a.size();
    int bSize = b.size();
    if (aSize == 1)
    {
        if (a.back() == 0) {
            result = { 0 };
            return;
        }
        if (a.back() == 1) {
            result = b;
            return;
        }
    }
    else if (bSize == 1)
    {
        if (b.back() == 0) {
            result = { 0 };
            return;
        }
        if (b.back() == 1) {
            result = a;
            return;
        }
    }

    if (result != a && result != b)
    {
        result.resize(aSize + bSize);
        for (int i = 0; i < aSize; ++i) {
            for (int j = 0; j < bSize; ++j) {
                uint64_t a_low = a[i] & 0xFFFFFFFFULL;
                uint64_t a_high = a[i] >> 32;
                uint64_t b_low = b[j] & 0xFFFFFFFFULL;
                uint64_t b_high = b[j] >> 32;

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
    }
    else
    {
        bigInt temp(aSize + bSize);
        for (int i = 0; i < aSize; ++i) {
            for (int j = 0; j < bSize; ++j) {
                uint64_t a_low = a[i] & 0xFFFFFFFFULL;
                uint64_t a_high = a[i] >> 32;
                uint64_t b_low = b[j] & 0xFFFFFFFFULL;
                uint64_t b_high = b[j] >> 32;

                // Perform 32-bit multiplications
                uint64_t low_low = a_low * b_low;
                uint64_t low_high = a_low * b_high;
                uint64_t high_low = a_high * b_low;
                uint64_t high_high = a_high * b_high;

                // Combine results
                uint64_t carry = (low_low >> 32) + (low_high & 0xFFFFFFFF) + (high_low & 0xFFFFFFFF);
                uint64_t low = (low_low & 0xFFFFFFFF) | (carry << 32);
                uint64_t high = high_high + (low_high >> 32) + (high_low >> 32) + (carry >> 32);

                temp[i + j] += low;
                carry = (temp[i + j] < low) ? 1 : 0;
                temp[i + j + 1] += high + carry;
                carry = temp[i + j + 1] < (high + carry) ? 1 : 0;
                int z = i + j + 2;
                while (carry)
                {
                    temp[z] += 1;
                    carry = (temp[z] == 0) ? 1 : 0;
                    ++z;
                }
            }
        }
        result = temp;
    }
    normalize(result);
}

//r = a % b
void mod(bigInt& r, const bigInt& a, const bigInt& b)
{
    if (smallerBigInt(a, b)) {
        r = a;
        return;
    }
    if (a == b)
    {
        r = { 0 };
        return;
    }
    if (r != a && r != b)
    {
        r = a;
        int otherLength = bitLength(b);
        int remainLength = bitLength(r);

        while (remainLength > otherLength) {
            int bitShift = remainLength - otherLength - 1;

            sub(r, r, shiftLeft(b, bitShift));
            remainLength = bitLength(r);
        }

        if (!smallerBigInt(r, b)) {
            sub(r, r, b);
        }
    }
    else
    {
        bigInt tempR = a;
        int otherLength = bitLength(b);
        int remainLength = bitLength(tempR);

        while (remainLength > otherLength) {
            int bitShift = remainLength - otherLength - 1;

            sub(tempR, tempR, shiftLeft(b, bitShift));
            remainLength = bitLength(tempR);
        }

        if (!smallerBigInt(tempR, b)) {
            sub(tempR, tempR, b);
        }
        r = tempR;
    }
    normalize(r);
}

// Function to perform modular exponentiation (a^b % mod)
void modExp(bigInt& result, bigInt base, bigInt exp, const bigInt& n) {
    int expSize = bitLength(exp);
    mod(base, base, n); // Handle base larger than mod
    while (expSize > 0) {
        mul(result, result, result);
        mod(result, result, n);
        if (getBit(exp, --expSize) == 1) {
            mul(result, result, base);
            mod(result, result, n);
        }
    }
}

// Miller-Rabin Primality Test
bool miller_rabin(const bigInt& n, int k) {
    // Edge case: Handle small numbers
    if (!greaterNormInt(n, 1)) return false;
    if (equalNormInt(n, 2) || equalNormInt(n, 3)) return true;
    if ((n[0] & 1) == 0)
        return false;

    // Find d and r such that n - 1 = d * 2^r
    bigInt d = { 1 };
    sub(d, n, d);
    bigInt tempN1 = d;
    int r = leastBit(d);
   
    shiftRightSelf(d, r);

    // Pre-determined bases for deterministic Miller-Rabin up to 64-bit numbers
    const uint64_t bases[] = { generator(),generator(),generator(),generator(),generator(),2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
    for (int i = 0; i < k; ++i) {
        bigInt base(1, (i < 17) ? bases[i] : generator());
        if (!smallerBigInt(base, n)) 
            continue;

        bigInt x = { 1 };
        modExp(x, base, d, n);
        if (equalNormInt(x, 1) || x == tempN1)
            continue;
        bool continueTest = false;
        for (int j = 0; j < r; ++j) {
            mul(x, x, x);
            mod(x, x, n);
            
            if (x == tempN1) {
                continueTest = true;
                break;
            }
        }
        if (!continueTest)
            return false;
    }
    return true;
}

bigInt readInput(const string& filename)
{
    ifstream file(filename);
    if (!file) {
        throw std::ios_base::failure("Failed to open file");
    }

    string fileContent = "";
    getline(file, fileContent);

    bigInt bigInt;
    hexToBigInt(bigInt, fileContent);
    file.close();
    return bigInt;
}

void writeOutput(const string& filename, const bigInt& n)
{
    ofstream file(filename);
    if (!file) {
        throw std::ios_base::failure("Failed to open file");
    }
    int nSize = n.size();
    int k = 0;
    if (nSize == 1)
        k = 17;
    else if (nSize <= 3)
        k = 32;
    else if (nSize <= 7)
        k = 32;
    else if (nSize <= 11)
        k = 30;
    else
        k = 5;

    if (miller_rabin(n, k)) {
        file << 1;
    }
    else {
        file << 0;
    }
    file.close();
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    if (argc != 3) {
        return 1;
    }
    bigInt n = readInput(argv[1]);
    writeOutput(argv[2], n);
    return 0;
}
