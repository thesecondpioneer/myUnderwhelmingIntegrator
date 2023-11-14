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

vector<double> durand_kerner_real(vector<double> coefficients, int max_iterations, double epsilon) {
    int degree = coefficients.size() - 1;
    vector<double> roots(degree);
    double leadcoef = coefficients[coefficients.size() - 1];
    for (int i = 0; i < coefficients.size(); ++i) {
        coefficients[i] /= leadcoef;
    }
    for (int i = 0; i < degree; ++i) {
        double angle = M_PI * i / degree;
        roots[i] = cos(angle);
    }
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        vector<double> new_roots(degree);
        for (int i = 0; i < degree; ++i) {
            double denominator = 1.0;
            for (int j = 0; j < degree; ++j) {
                if (i != j) {
                    denominator *= (roots[i] - roots[j]);
                }
            }
            double f_x = 0.0;
            for (int k = 0; k <= degree; ++k) {
                f_x += coefficients[k] * pow(roots[i], k);
            }
            new_roots[i] = roots[i] - f_x / denominator;
        }
        bool converged = true;
        for (int i = 0; i < degree; ++i) {
            if (abs(new_roots[i] - roots[i]) > epsilon) {
                converged = false;
                break;
            }
        }
        if (converged) {
            return new_roots;
        }
        roots = new_roots;
    }
    return roots;
}

double f(double x) {
    return (2 * cos(2.5 * x + 3.75) * exp(x / 3 + 0.5) + 4 * sin(3.5 * x + 5.25) * exp(-3 * x - 4.5) + x + 1.5);
}

double F(double x) {
    return (f(x) / pow(x, 1.0 / 3));
}

vector<double> eq_dist(double m, double a, double b) {
    int n = ceil(m);
    vector<double> result(1, a);
    double step = (b - a) / (n-1);
    for (int i = 1; i < n; i++) {
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

vector<double> ncqfsum(vector<double> (*distFunc)(double, double, double), int n, vector<double> &z) {
    vector<double> knots, qfcoeffs;
    double result = 0, summod = 0;
    for (int i = 0; i < z.size() - 1; i++) {
        vector<double> table(0);
        knots = distFunc(n, z[i], z[i + 1]);
        table.insert(table.end(), knots.begin(), knots.end());
        vector<vector<double>> T(n, vector<double>(n, 1));
        if (n > 1) {
            for (int j = 0; j < n; j++) {
                T[1][j] = table[j];
            }
        }
        for (int j = 2; j < n; j++) {
            for (int k = 0; k < n; k++) {
                T[j][k] = pow(table[k], j);
            }
        }
        qfcoeffs = LUPsolve(T, transpose(moments(n, {z[i], z[i + 1]})));
        for (int j = 0; j < n; j++) {
            result += qfcoeffs[j] * f(table[j]);
            summod += abs(qfcoeffs[j]);
        }
    }
    return {result, summod};
}

vector<double> hqfsum(int n, vector<double> &z) {
    double result = 0, summod = 0;
    for (int i = 0; i < z.size() - 1; i++) {
        vector<vector<double>> mu = moments(2 * n, {z[i], z[i + 1]}), mu2(1), b(n, vector<double>(1));
        vector<vector<double>> M(n, vector<double>(n, 1));
        vector<double> table(0), qfcoeffs;
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                M[k][j] = mu[0][j + k];
            }
            b[k][0] = -mu[0][n + k];
        }
        vector<double> s = LUPsolve(M, b);
        s.push_back(1.0);
        vector<double> roots = durand_kerner_real(s, 100000, 0.000000000000001);
        table.insert(table.end(), roots.begin(), roots.end());
        vector<vector<double>> T(n, vector<double>(n, 1));
        if (n > 1) {
            for (int j = 0; j < n; j++) {
                T[1][j] = table[j];
            }
        }
        for (int j = 2; j < n; j++) {
            for (int k = 0; k < n; k++) {
                T[j][k] = pow(table[k], j);
            }
        }
        for (int j = 0; j < n; j++) {
            mu2[0].push_back(mu[0][j]);
        }
        qfcoeffs = LUPsolve(T, transpose(mu2));
        for (int j = 0; j < n; j++) {
            result += qfcoeffs[j] * f(table[j]);
            summod += abs(qfcoeffs[j]);
        }
    }
    return {result, summod};
}


int main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    cout.precision(16);
    double prev = 0, cur, precise = 7.077031437995793610263911711602477164432, summod, eps, h, p;
    vector<double> z, is;
    vector<double> plot1, plot2, plot3, plot4;
    double l, m;
    int nn; //L multiplier, the amount of intervals and the amount of nodes;
    string mode;
    cin >> mode;
    if (mode == "m") {
        cin >> m >> nn;
        z = eq_dist(m + 1, 0, 1.8);
        for (int n = 1; n <= nn; n++) {
            is = hqfsum(n, z);
            cur = is[0];
            summod = is[1];
            plot1.push_back(cur);
            plot2.push_back(abs(cur - prev));
            plot3.push_back(abs(precise - cur));
            plot4.push_back(summod);
            cout << n << ' ' << fixed << cur << ' ' << fixed << abs(cur - prev) << ' ' << abs(precise - cur) << ' '
                 << summod << endl;
            prev = cur;
        }
    } else if (mode == "ah") { //Hauss + eitken + runge
        cin >> eps >> nn >> l >> h; //precision, amount of nodes, L multiplier, initial step
        vector<double> s(3);
        m = ceil(1.8 / h);
        for (int i = 0; i < 3; i++) {
            z = eq_dist(m + 1, 0, 1.8);
            s[i] = hqfsum(nn, z)[0];
            m *= l;
        }
        double cm1, cm2;
        while (true) {
            m = ceil(1.8 / h) * l;
            p = -log(abs((s[2] - s[1]) / (s[1] - s[0]))) / log(l);
            cm1 = abs((s[1] - s[0]) / (pow(h, p) * (1 - pow(l, -p))));
            cm2 = abs((s[2] - s[1]) / (pow(h/2, p) * (1 - pow(l, -p))));
            if((static_cast<int>(std::floor(std::log10(std::abs(cm1 - cm2)))) < 0.01) and (p > 1)){
                break;
            }
            s[0] = s[1];
            s[1] = s[2];
            z = eq_dist(m*l*l + 1, 0, 1.8);
            s[2] = hqfsum(nn, z)[0];
            h /= l;
        }
        h *= 0.95 * pow((eps * (1 - pow(l, -p))) / abs(s[1] - s[0]), 1 / p);
        if (ceil(1.8 / h) > m*l) {
            double m = ceil(1.8 / h);
            z = eq_dist(m + 1, 0, 1.8);
            s[2] = hqfsum(nn, z)[0];
        }
        cout << fixed << s[2] << ' ' << abs(precise - s[2]) << endl;
    } else if (mode == "anc") { //Newton-Cotes + eitken + runge
        cin >> eps >> nn >> l >> h; //precision, amount of nodes, L multiplier, initial step
        vector<double> s(3);
        m = ceil(1.8 / h);
        for (int i = 0; i < 3; i++) {
            z = eq_dist(m + 1, 0, 1.8);
            s[i] = ncqfsum(eq_dist,nn, z)[0];
            m *= l;
        }
        double cm1, cm2;
        while (true) {
            m = ceil(1.8 / h) * l;
            p = -log(abs((s[2] - s[1]) / (s[1] - s[0]))) / log(l);
            cm1 = abs((s[1] - s[0]) / (pow(h, p) * (1 - pow(l, -p))));
            cm2 = abs((s[2] - s[1]) / (pow(h/2, p) * (1 - pow(l, -p))));
            if((static_cast<int>(std::floor(std::log10(std::abs(cm1 - cm2)))) < 0.01) and (p > 1)){
                break;
            }
            s[0] = s[1];
            s[1] = s[2];
            z = eq_dist(m*l*l + 1, 0, 1.8);
            s[2] = ncqfsum(eq_dist,nn, z)[0];
            h /= l;
        }
        h *= 0.95 * pow((eps * (1 - pow(l, -p))) / abs(s[1] - s[0]), 1 / p);
        if (ceil(1.8 / h) > m*l) {
            double m = ceil(1.8 / h);
            z = eq_dist(m + 1, 0, 1.8);
            s[2] = hqfsum(nn, z)[0];
        }
        cout << fixed << s[2] << ' ' << abs(precise - s[2]) << endl;
    } else if (mode == "p") {
        for(eps = 0.001; eps >= 0.00000000000001; eps*=0.1) {
            vector<double> s(3);
            nn = 3; //node amount
            l = 2; //L multiplier
            h = 0.9; //initial step
            m = ceil(1.8 / h);
            for (int i = 0; i < 3; i++) {
                z = eq_dist(m + 1, 0, 1.8);
                s[i] = hqfsum(nn, z)[0];
                m *= l;
            }
            double cm1, cm2;
            while (true) {
                m = ceil(1.8 / h) * l;
                p = -log(abs((s[2] - s[1]) / (s[1] - s[0]))) / log(l);
                cm1 = abs((s[1] - s[0]) / (pow(h, p) * (1 - pow(l, -p))));
                cm2 = abs((s[2] - s[1]) / (pow(h/2, p) * (1 - pow(l, -p))));
                if ((static_cast<int>(std::floor(std::log10(std::abs(cm1 - cm2)))) < 0.01) and (p > 1)) {
                    break;
                }
                s[0] = s[1];
                s[1] = s[2];
                z = eq_dist(m * l * l + 1, 0, 1.8);
                s[2] = hqfsum(nn, z)[0];
                h /= l;
            }
            h *= 0.95 * pow((eps * (1 - pow(l, -p))) / abs(s[1] - s[0]), 1 / p);
            if (ceil(1.8 / h) > m * l) {
                double m = ceil(1.8 / h);
                z = eq_dist(m + 1, 0, 1.8);
                s[2] = hqfsum(nn, z)[0];
            }
            plot1.push_back(eps);
            plot2.push_back(abs(precise-s[2]));
            plot3.push_back(s[2]);
            }
        //for Matlab plotting
        /*cout << "[";
        for (int n = 1; n < 56; n++){
            cout << n << ' ';
        }
        cout << "]" << endl;*/
        cout << "[";
        for (int n = 0; n < plot1.size() ; n++){
            cout << fixed << plot1[n] << ' ';
        }
        cout << "]" << endl;
        /*cout << endl;
        cout << "[";
        for (int n = 1; n < 56; n++){
            cout << n << ' ';
        }
        cout << "]" << endl;*/
        cout << "[";
        for (int n = 0; n < plot2.size() ; n++){
            cout << fixed << plot2[n] << ' ';
        }
        cout << "]" << endl;
        /*cout << "[";
        for (int n = 1; n < 56; n++){
            cout << n << ' ';
        }
        cout << "]" << endl;*/
        cout << "[";
        for (int n = 0; n < plot3.size() ; n++){
            cout << fixed << plot3[n] << ' ';
        }
        cout << "]" << endl;
        /*
        cout << "[";
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
    return 0;
}
