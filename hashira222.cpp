#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace std;
using namespace boost::multiprecision;

// Big integer type
using BigInt = cpp_int;
// High precision float (100 decimal digits)
using BigFloat = cpp_dec_float_100;

// Convert string in given base (2–16) to big integer
BigInt toDecimal(const string &val, int base) {
    BigInt result = 0;
    for (char ch : val) {
        int digit;
        if (ch >= '0' && ch <= '9') digit = ch - '0';
        else digit = 10 + (tolower(ch) - 'a'); // for bases > 10
        result = result * base + digit;
    }
    return result;
}

// Gaussian elimination with high precision
vector<BigFloat> gaussSolve(vector<vector<BigFloat>> A, vector<BigFloat> b) {
    int n = A.size();
    for (int i = 0; i < n; i++) A[i].push_back(b[i]);

    for (int i = 0; i < n; i++) {
        // Pivot row
        int pivot = i;
        for (int j = i+1; j < n; j++) {
            if (abs(A[j][i]) > abs(A[pivot][i])) pivot = j;
        }
        swap(A[i], A[pivot]);

        // Normalize pivot row
        BigFloat div = A[i][i];
        for (int k = i; k <= n; k++) A[i][k] /= div;

        // Eliminate below
        for (int j = i+1; j < n; j++) {
            BigFloat factor = A[j][i];
            for (int k = i; k <= n; k++) {
                A[j][k] -= factor * A[i][k];
            }
        }
    }

    // Back substitution
    vector<BigFloat> x(n);
    for (int i = n-1; i >= 0; i--) {
        x[i] = A[i][n];
        for (int j = i+1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
    }
    return x;
}

int main() {
    // JSON-like data
    struct Node { int x, base; string value; };
    vector<Node> data = {
        {1, 6,  "13444211440455345511"},
        {2, 15, "aed7015a346d635"},
        {3, 15, "6aeeb69631c227c"},
        {4, 16, "e1b5e05623d881f"},
        {5, 8,  "316034514573652620673"},
        {6, 3,  "2122212201122002221120200210011020220200"},
        {7, 3,  "20120221122211000100210021102001201112121"},
        {8, 6,  "20220554335330240002224253"},
        {9, 12, "45153788322a1255483"},
        {10, 7, "1101613130313526312514143"}
    };

    int n = 10, k = 7, m = k-1;

    // Convert to points (x,y)
    vector<pair<int, BigInt>> points;
    for (auto &d : data) {
        BigInt y = toDecimal(d.value, d.base);
        points.push_back({d.x, y});
    }

    // Build system (first k points)
    vector<vector<BigFloat>> A(k, vector<BigFloat>(k));
    vector<BigFloat> b(k);
    for (int i = 0; i < k; i++) {
        int x = points[i].first;
        BigFloat power = 1;
        for (int j = 0; j < k; j++) {
            A[i][j] = power;
            power *= x;
        }
        b[i] = BigFloat(points[i].second); // convert BigInt → BigFloat safely
    }

    // Solve coefficients
    vector<BigFloat> coeff = gaussSolve(A, b);

    // Print polynomial
    cout << "Polynomial (degree " << m << "): f(x) = ";
    for (int i = m; i >= 0; i--) {
        cout << coeff[i];
        if (i > 0) cout << " * x^" << i << " + ";
    }
    cout << "\n\n";

    // Verification with tolerance
    const BigFloat TOL = BigFloat("1e-30"); // very strict tolerance
    cout << "Verification (tolerance " << TOL << "):\n";
    for (auto &p : points) {
        BigFloat fx = 0, power = 1;
        for (int j = 0; j < coeff.size(); j++) {
            fx += coeff[j] * power;
            power *= p.first;
        }
        BigFloat expected = BigFloat(p.second);
        BigFloat error = abs(fx - expected);

        cout << "f(" << p.first << ") = " << fx 
             << ", expected = " << expected 
             << " -> " << (error < TOL ? "OK ✅" : "Mismatch ❌")
             << " (error = " << error << ")\n";
    }
}

