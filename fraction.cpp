#include<iostream>
using namespace std;

class Fraction {
    int numerator;
    int dominator;
    int sign;

    static int getGCD(int a, int b) {
        if (a < 0) a = -a;
        if (b < 0) b = -b;
        if (a < b) swap(a, b);
        if (b == 0) return a;

        // a > b
        // 24 15
        // 24 / 15 = 1 ... 9
        // 15 / 9  = 1 ... 6
        // 9  / 6  = 1 ... 3
        // 6  / 3  = 2 ... 0

        while (b) {
            int temp = b;
            b = a % b;
            a = temp;
        }

        return a;
    }

public:
    Fraction() {
        this->sign = 1;
        this->numerator = 0;
        this->dominator = 1;
    }
    Fraction(int num) {
        this->sign = 1;
        if (num < 0) {
            this->sign = -1;
            num = -num;
        }
        this->numerator = num;
        this->dominator = 1;
    }
    Fraction(int numerator, int dominator) {
        this->sign = 1;
        if ((numerator < 0 && dominator > 0) || (numerator > 0 && dominator < 0)) {
            this->sign = -1;
        }
        if (numerator < 0) numerator = -numerator;
        if (dominator < 0) dominator = -dominator;
        this->numerator = numerator;
        this->dominator = dominator;
    }
    Fraction(const string& serialized, int begin, int end) {
        if (end - begin <= 0) {
            this->numerator = 0;
            this->dominator = 1;
            return;
        }
        if (serialized[begin] == '-') {
            this->sign = -1;
            ++begin;
        }
        else if (serialized[begin] == '+') {
            this->sign = 1;
            ++begin;
        }

        int num = 0;
        while (begin < end && serialized[begin] != '/') {
            num *= 10;
            num += serialized[begin] - '0';
            ++begin;
        }

        this->numerator = num;
        num = 0;

        ++begin;
        if (begin >= end) this->dominator = 1;

        while (begin < end) {
            num *= 10;
            num += serialized[begin] - '0';
            ++begin;
        }
        this->dominator = num;
    }
    Fraction(const Fraction& obj) {
        this->numerator = obj.numerator;
        this->dominator = obj.dominator;
        this->sign = obj.sign;
    }

    Fraction& standardize() {
        int gcd = getGCD(this->numerator, this->dominator);
        this->numerator /= gcd;
        this->dominator /= gcd;
        return *this;
    }

    string serialize() {
        return to_string(sign*numerator) + '/' + to_string(dominator);
    }

    Fraction& add(const Fraction& obj) {

        int gcd = getGCD(this->dominator, obj.dominator);        
        int factor_obj  = this->dominator / gcd;
        int factor_this = obj.dominator / gcd;

        this->numerator = this->sign * this->numerator * factor_this + obj.sign * obj.numerator * factor_obj;
        this->sign = 1;
        if (this->numerator < 0) {
            this->numerator = -this->numerator;
            this->sign = -1;
        }

        this->dominator *= factor_this;
        return *this;
    }

    Fraction& multiply(const Fraction& obj) {
        
        this->sign = this->sign * obj.sign;

        this->numerator *= obj.numerator;
        this->dominator *= obj.dominator;

        return *this;
    }
    Fraction& reverse() {
        swap(this->numerator, this->dominator);
        return *this;
    }

    Fraction& setNumerator(int num) {
        this->numerator = abs(num);
        return *this;
    }
    Fraction& setDominator(int num) {
        this->dominator = abs(num);
        return *this;
    }
    int getNumerator() {
        return this->numerator;
    }
    int getDominator() {
        return this->dominator;
    }

    // >0:+  <0:-  =0:reverse
    Fraction& setSign(int sign = 0) {
        if (sign == 0) this->sign = -this->sign;
        else if (sign > 0) this->sign = 1;
        else this->sign = -1;

        return *this;
    }
    int getSign() {
        return this->sign;
    }

    Fraction& operator+=(const Fraction& obj) {
        this->add(obj);
        return *this;
    }
    Fraction& operator-=(const Fraction& obj) {
        this->add(Fraction(obj).setSign());
        return *this;
    }
    Fraction& operator=(const Fraction& obj) {
        this->numerator = obj.numerator;
        this->dominator = obj.dominator;
        this->sign = obj.sign;
        return *this;
    }
    Fraction operator+(const Fraction& obj) {
        return Fraction(*this) += obj;
    }
    Fraction operator-(const Fraction& obj) {
        return Fraction(*this) -= obj;
    }
    Fraction& operator*=(const Fraction& obj) {
        this->multiply(obj);
        return *this;
    }
    Fraction& operator/=(const Fraction& obj) {
        this->multiply(Fraction(obj).reverse());
        return *this;
    }
    Fraction operator/(const Fraction& obj) {
        return Fraction(*this) /= obj;
    }
    Fraction operator*(const Fraction& obj) {
        return Fraction(*this) *= obj;
    }

    Fraction& operator+=(const int num) {
        this->add(Fraction(num));
        return *this;
    }
    Fraction& operator-=(const int num) {
        this->add(Fraction(num).setSign());
        return *this;
    }
    Fraction& operator=(const int num) {
        *this = Fraction(num);
        return *this;
    }
    Fraction operator+(const int num) {
        return Fraction(*this) += Fraction(num);
    }
    Fraction operator-(const int num) {
        return Fraction(*this) -= Fraction(num);
    }
    Fraction& operator*=(const int num) {
        this->multiply(Fraction(num));
        return *this;
    }
    Fraction& operator/=(const int num) {
        this->multiply(Fraction(num).reverse());
        return *this;
    }
    Fraction operator/(const int num) {
        return Fraction(*this) /= Fraction(num);
    }
    Fraction operator*(const int num) {
        return Fraction(*this) *= Fraction(num);
    }

    void print() {
        if (sign == -1) cout << '-';
        cout << this->numerator;
        cout << " / ";
        cout << this->dominator;
        cout << endl;
    }
};