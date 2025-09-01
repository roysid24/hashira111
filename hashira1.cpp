#include <bits/stdc++.h>
using namespace std;

// Convert a string number from a given base to decimal
int toDecimal(const string &val, int base) {
    int result = 0;
    for (char ch : val) {
        int digit = ch - '0';
        result = result * base + digit;
    }
    return result;
}

// Gaussian elimination solver for Ax = b
vector<double> gaussSolve(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        A[i].push_back(b[i]); // augment matrix
    }

    // Forward elimination
    for (int i = 0; i < n; i++) {
        // Pivot
        int pivot = i;
        for (int j = i+1; j < n; j++) {
            if (fabs(A[j][i]) > fabs(A[pivot][i])) pivot = j;
        }
        swap(A[i], A[pivot]);

        // Normalize pivot row
        double div = A[i][i];
        for (int k = i; k <= n; k++) A[i][k] /= div;

        // Eliminate below
        for (int j = i+1; j < n; j++) {
            double factor = A[j][i];
            for (int k = i; k <= n; k++) {
                A[j][k] -= factor * A[i][k];
            }
        }
    }

    // Back substitution
    vector<double> x(n);
    for (int i = n-1; i >= 0; i--) {
        x[i] = A[i][n];
        for (int j = i+1; j < n; j++) {
            x[i] -= A[i][j]*x[j];
        }
    }
    return x;
}

int main() {
    // Given JSON-like data
    struct Node { int x, base; string value; };
    vector<Node> data = {
        {1, 10, "4"},
        {2, 2, "111"},
        {3, 10, "12"},
        {6, 4, "213"}
    };

    // Convert to (x,y) points
    vector<pair<int,int>> points;
    for (auto &d : data) {
        int y = toDecimal(d.value, d.base);
        points.push_back({d.x, y});
    }

    int k = 3; // minimum required = degree+1
    int m = k-1;

    // Build linear system
    vector<vector<double>> A(k, vector<double>(k));
    vector<double> b(k);
    for (int i = 0; i < k; i++) {
        int x = points[i].first;
        int y = points[i].second;
        double power = 1;
        for (int j = 0; j < k; j++) {
            A[i][j] = power;
            power *= x;
        }
        b[i] = y;
    }

    // Solve for coefficients
    vector<double> coeff = gaussSolve(A, b);

    // Print polynomial
    cout << "Polynomial: f(x) = ";
    for (int i = m; i >= 0; i--) {
        cout << coeff[i];
        if (i > 0) cout << "x^" << i << " + ";
    }
    cout << "\n\n";

    // Verify all points
    cout << "Verification:\n";
    for (auto &p : points) {
        double fx = 0, power = 1;
        for (int j = 0; j < coeff.size(); j++) {
            fx += coeff[j]*power;
            power *= p.first;
        }
        cout << "f(" << p.first << ") = " << fx 
             << " (expected " << p.second << ")\n";
    }

    return 0;
}
