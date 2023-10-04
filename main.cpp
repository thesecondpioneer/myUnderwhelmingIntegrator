#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

vector<vector<double>> mmultiply(vector<vector<double>> A, vector<vector<double>> B) { //matrix multiplication
    if (A[0].size() == B.size()) {
        vector<vector<double>> C(A.size());
        for (int i = 0; i < C.size(); i++) {
            C[i].resize(B[0].size());
        }
        for (int i = 0; i < A.size(); i++) {
            for (int j = 0; j < B[0].size(); j++) {
                for (int k = 0; k < A[0].size(); k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return C;
    }
}

pair<vector<vector<double>>, vector<vector<double>>>
LUP(vector<vector<double>> A) { //returns a matrix that stores L,U and a matrix P (transposition matrix)
    vector<vector<double>> P(A.size());
    for (int i = 0; i < P.size(); i++) {
        P[i].resize(A.size(), 0);
        P[i][i] = 1;
    }
    for (int i = 0; i < A.size() - 1; i++) {
        double lead = INT64_MIN;
        double nlead = -1;
        for (int j = i; j < A.size(); j++) {
            if (abs(A[j][i]) > lead) {
                lead = abs(A[j][i]);
                nlead = j;
            }
        }
        swap(A[i], A[nlead]);
        swap(P[i], P[nlead]);
        for (int j = i + 1; j < A.size(); j++) {
            A[j][i] = A[j][i] / A[i][i];
            for (int k = i + 1; k < A.size(); k++) {
                A[j][k] = A[j][k] - A[j][i] * A[i][k];
            }
        }
    }
    return make_pair(A, P);
}

vector<double> LUPsolve(vector<vector<double>> A, vector<vector<double>> b) {
    //solves the equation system by using the results of LUP function
    pair<vector<vector<double>>, vector<vector<double>>> LpUaP = LUP(A);
    vector<vector<double>> LU = LpUaP.first;
    b = mmultiply(LpUaP.second, b);
    vector<double> y(b.size());
    for (int i = 0; i < b.size(); i++) {
        y[i] = b[i][0];
    }
    for (int i = 0; i < A.size(); i++) {
        for (int k = 0; k < i; k++) {
            y[i] -= LU[i][k] * y[k];
        }
    }
    vector<double> x(b.size());
    for (int i = b.size() - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int k = i + 1; k < b.size(); k++) {
            x[i] -= LU[i][k] * x[k];
        }
        x[i] = x[i] / LU[i][i];
    }
    return x;
}

vector<double> getcolumn(vector<vector<double>> A, int k) { //get a column of number k as a vector from a matrix
    vector<double> result;
    for (int i = 0; i < A.size(); i++) {
        result.push_back(A[i][k]);
    }
    return result;
}

vector<vector<double>> transpose(vector<vector<double>> A) { //returns the transposed matrix
    vector<vector<double>> result(A[0].size());
    for (int i = 0; i < result.size(); i++) {
        result[i] = getcolumn(A, i);
    }
    return result;
}

double f(double x) {
    return (2 * cos(2.5 * x + 3.75) * exp(x / 3 + 0.5) + 4 * sin(3.5 * x + 5.25) * exp(-3 * x - 4.5) + x + 1.5);
}

double F(double x) {
    return (f(x) / pow(x, 1.0 / 3));
}

vector<double> eq_dist(int n, double a, double b) {
    vector<double> result(0);
    double step = (b - a) / (n + 1);
    for (int i = 0; i < n; i++) {
        a += step;
        result.push_back(a);
    }
    return result;
}

double darbSum(vector<double> (*distFunc)(int, double, double), double (*function)(double), double a, double b, int n) {
    double result = 0;
    vector<double> table = distFunc(n, a, b);
    for (int i = 0; i < table.size() - 1; i++) {
        result += function((table[i + 1] + table[i]) / 2) * (table[i + 1] - table[i]);
    }
    return result;
}

vector<vector<double>> moments(int n, const vector<double> &z) {
    vector<vector<double>> result(z.size() - 1, vector<double>(n));
    for (int i = 0; i < z.size() - 1; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = pow(z[i + 1], j + 2.0 / 3) / (j + 2.0 / 3) - pow(z[i], j + 2.0 / 3) / (j + 2.0 / 3);
        }
    }
    return result;
}

vector<double> qfsum(vector<double> (*distFunc)(int, double, double), int n, vector<double> &z) {
    vector<double> knots, table(0), qfcoeffs;
    double result = 0;
    for (int i = 0; i < z.size() - 1; i++) {
        knots = distFunc(n, z[i], z[i + 1]);
        table.insert(table.end(), knots.begin(), knots.end());
    }
    vector<vector<double>> T(n, vector<double>(n, 1));
    if (n > 1) {
        for (int i = 0; i < n; i++) {
            T[1][i] = table[i];
        }
    }
    for (int i = 2; i < n; i++) {
        for (int j = 0; j < n; j++) {
            T[i][j] = pow(table[j], i);
        }
    }
    qfcoeffs = LUPsolve(T, transpose(moments(n, z)));
    double summod;
    for (int i = 0; i < n; i++) {
        result += qfcoeffs[i] * f(table[i]);
        summod+=abs(qfcoeffs[i]);
        //cout << qfcoeffs[i] << endl;
    }
    return {result, summod};
}

int main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    cout.precision(16);
    double prev = 0, cur, precise = 7.077031437995793610263911711602477164432, summod;
    vector<double> z(2), is;
    vector<double> plot1, plot2, plot3, plot4;
    z[0] = 0;
    z[1] = 1.8;
    for (int n = 1; n < 56; n++) {
        is = qfsum(eq_dist, n, z);
        cur = is[0];
        summod = is[1];
        plot1.push_back(cur);
        plot2.push_back(abs(cur-prev));
        plot3.push_back(abs(precise-cur));
        plot4.push_back(summod);
        cout << n << ' ' << fixed << cur << ' ' << fixed << abs(cur - prev) << ' ' << abs(precise - cur) << ' ' << summod << endl;
        prev = cur;
    }
    /*cout << "[";
    for (int n = 1; n < 56; n++){
        cout << n << ' ';
    }
    cout << "]" << endl;
    cout << "[";
    for (int n = 0; n < 55 ; n++){
        cout << fixed << plot1[n] << ' ';
    }
    cout << "]";
    cout << endl;
    cout << "[";
    for (int n = 1; n < 56; n++){
        cout << n << ' ';
    }
    cout << "]" << endl;
    cout << "[";
    for (int n = 0; n < 55 ; n++){
        cout << fixed << plot2[n] << ' ';
    }
    cout << "]";
    cout << endl;
     /*cout << "[";
    for (int n = 1; n < 56; n++){
        cout << n << ' ';
    }
    cout << "]" << endl;
    cout << "[";
    for (int n = 0; n < 55 ; n++){
        cout << fixed << plot3[n] << ' ';
    }
    cout << "]";
    cout << endl;*/
    /*cout << "[";
    for (int n = 1; n < 56; n++){
        cout << n << ' ';
    }
    cout << "]" << endl;
    cout << "[";
    for (int n = 0; n < 55 ; n++){
        cout << fixed << plot4[n] << ' ';
    }
    cout << "]";
    cout << endl;*/
}
