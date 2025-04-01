#include"Kepler.h"

void KeplerX(Matrix<double, 6, 1> q0, double TF, Matrix<double, 6, 1>* q1) {
    Matrix <double, 3, 1> X0{ q0(0, 0), q0(1, 0), q0(2, 0) };
    Matrix <double, 3, 1> V0{ q0(3, 0), q0(4, 0), q0(5, 0) };
    double DMU = 3.986004418E5, eps = 1e-12;
    double DT = TF;
    Matrix <double, 6, 1> X = q0;
    if (DT == 0.0) {
        for (int i = 0; i < 6; i++) {
            (*q1)(i) = q0(i);
        }
    }
    double V2 = pow(V0(0), 2) + pow(V0(1), 2) + pow(V0(2), 2);
    double R0 = sqrt(pow(X0(0), 2) + pow(X0(1), 2) + pow(X0(2), 2));
    double RV0 = X0(0) * V0(0) + X0(1) * V0(1) + X0(2) * V0(2);
    double H = V2 / 2.0 - DMU / R0;
    double A = -DMU / (2.0 * H);
    double DN1 = sqrt(DMU / pow(A, 3));
    double DN2 = sqrt(DMU * A);

    double DE2, CS, SN;
    double DM = DN1 * DT;
    double DIFFER = 1.0;
    double DE = DM;
    int NCOUNT = 0;
    while (DIFFER > eps && NCOUNT != 1000) {
        CS = cos(DE);
        SN = sin(DE);
        DE2 = DE - (DE + RV0 * (1.0 - CS) / DN2 - (1.0 - R0 / A) * SN - DM) / (1.0 + RV0 * SN / DN2 - (1.0 - R0 / A) * CS);
        DIFFER = abs((DE - DE2) / (DE + DE2));
        DE = DE2;
        NCOUNT += 1;
    }
    if (NCOUNT == 1000) cout << "Divergence!" << DIFFER << endl;

    double OCS = 1.0 - CS;
    double D1 = 1.0 - A * OCS / R0;
    double D2 = DT - (DE - SN) / DN1;
    for (int i = 0; i < 3; i++) {
        X(i) = D1 * X0(i) + D2 * V0(i);
    }
    double R = sqrt(pow(X(0), 2) + pow(X(1), 2) + pow(X(2), 2));
    double D3 = -DN2 * SN / (R * R0);
    double D4 = 1.0 - A * OCS / R;
    for (int i = 3; i < 6; i++) {
        X(i) = D3 * X0(i - 3) + D4 * V0(i - 3);
    }
    for (int i = 0; i < 6; i++) {
        (*q1)(i) = X(i);
    }
};

void KeplerXDX(Matrix<double, 6, 1> q0, double TF, Matrix<double, 6, 1>* q1, Matrix<double, 6, 6>* dx) {
    Matrix <double, 3, 1> X0{ q0(0, 0), q0(1, 0), q0(2, 0) };
    Matrix <double, 3, 1> V0{ q0(3, 0), q0(4, 0), q0(5, 0) };
    double DMU = 3.986004418E5, eps = 1e-12;
    double DT = TF;
    Matrix <double, 6, 1> X = q0;
    if (DT == 0.0) {
        Matrix<double, 6, 6> dx_res;
        dx_res.setZero();
        for (int i = 0; i < 6; i++) {
            dx_res(i, i) = 1.0;
        }
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                (*dx)(i, j) = dx_res(i, j);
            }
        }
        for (int i = 0; i < 6; i++) {
            (*q1)(i) = q0(i);
        }
    }
    double V2 = pow(V0(0), 2) + pow(V0(1), 2) + pow(V0(2), 2);
    double R0 = sqrt(pow(X0(0), 2) + pow(X0(1), 2) + pow(X0(2), 2));
    double RV0 = X0(0) * V0(0) + X0(1) * V0(1) + X0(2) * V0(2);
    double H = V2 / 2.0 - DMU / R0;
    double A = -DMU / (2.0 * H);
    double DN1 = sqrt(DMU / pow(A, 3));
    double DN2 = sqrt(DMU * A);

    double DE2, CS, SN;
    double DM = DN1 * DT;
    double DIFFER = 1.0;
    double DE = DM;
    int NCOUNT = 0;
    while (DIFFER > eps && NCOUNT != 1000) {
        CS = cos(DE);
        SN = sin(DE);
        DE2 = DE - (DE + RV0 * (1.0 - CS) / DN2 - (1.0 - R0 / A) * SN - DM) / (1.0 + RV0 * SN / DN2 - (1.0 - R0 / A) * CS);
        DIFFER = abs((DE - DE2) / (DE + DE2));
        DE = DE2;
        NCOUNT += 1;
    }
    if (NCOUNT == 1000) cout << "Divergence!" << DIFFER << endl;

    double OCS = 1.0 - CS;
    double D1 = 1.0 - A * OCS / R0;
    double D2 = DT - (DE - SN) / DN1;
    for (int i = 0; i < 3; i++) {
        X(i) = D1 * X0(i) + D2 * V0(i);
    }
    double R = sqrt(pow(X(0), 2) + pow(X(1), 2) + pow(X(2), 2));
    double D3 = -DN2 * SN / (R * R0);
    double D4 = 1.0 - A * OCS / R;
    for (int i = 3; i < 6; i++) {
        X(i) = D3 * X0(i - 3) + D4 * V0(i - 3);
    }

    // Derivatives

    Matrix <double, 3, 1> DH{ DMU * X0 / (H * pow(R0, 3)) };
    Matrix <double, 3, 1> HLD{ X0 / pow(R0, 2) + DH };
    Matrix <double, 3, 1> CC1{ (V0 + 0.5 * RV0 * DH) / DN2 };
    Matrix <double, 3, 1> CC2{ R0 * HLD / A };
    Matrix <double, 3, 1> DDE{ A * (1.5 * DH * DM - CC1 * OCS - CC2 * SN) / R };
    Matrix <double, 3, 1> DD1{ A * (OCS * HLD - SN * DDE) / R0 };
    Matrix <double, 3, 1> DD2{ (1.5 * DH * (DE - SN) - OCS * DDE) / DN1 };
    Matrix <double, 3, 1> XDQ;
    int ip3, jp3;
    Matrix <double, 6, 6> DX;
    DX.setZero();

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            (DX)(i, j) = X0(i) * DD1(j) + V0(i) * DD2(j);
            if (i == j) (DX)(i, i) = (DX)(i, i) + D1;
        }
    }

    XDQ.setZero();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            XDQ(i) = XDQ(i) + (X)(j) * (DX)(j, i);
        }
    }

    Matrix <double, 3, 1> DD3{ DN2 * (-CS * DDE + SN * (DH / 2.0 + XDQ / pow(R, 2) + X0 / pow(R0, 2))) / (R * R0) };
    Matrix <double, 3, 1> DD4{ A * (-SN * DDE + OCS * (DH + XDQ / pow(R, 2))) / R };

    for (int i = 0; i < 3; i++) {
        ip3 = i + 3;
        for (int j = 0; j < 3; j++) {
            (DX)(ip3, j) = X0(i) * DD3(j) + V0(i) * DD4(j);
            if (i == j) (DX)(ip3, i) = (DX)(ip3, i) + D3;
        }
    }

    DH = V0 / H;
    CC1 = (X0 + 0.5 * RV0 * DH) / DN2;
    CC2 = R0 * DH / A;
    DDE = A * (1.5 * DH * DM - CC1 * OCS - CC2 * SN) / R;
    DD1 = A * (OCS * DH - SN * DDE) / R0;
    DD2 = (1.5 * DH * (DE - SN) - OCS * DDE) / DN1;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            jp3 = j + 3;
            (DX)(i, jp3) = X0(i) * DD1(j) + V0(i) * DD2(j);
            if (i == j) (DX)(i, jp3) = (DX)(i, jp3) + D2;
        }
    }

    XDQ.setZero();
    for (int i = 0; i < 3; i++) {
        ip3 = i + 3;
        for (int j = 0; j < 3; j++) {
            XDQ(i) = XDQ(i) + (X)(j) * (DX)(j, ip3);
        }
    }

    DD3 = DN2 * (-CS * DDE + SN * (DH / 2.0 + XDQ / pow(R, 2))) / (R * R0);
    DD4 = A * (-SN * DDE + OCS * (DH + XDQ / pow(R, 2))) / R;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            jp3 = j + 3;
            ip3 = i + 3;
            (DX)(ip3, jp3) = X0(i) * DD3(j) + V0(i) * DD4(j);
            if (i == j) (DX)(ip3, jp3) = (DX)(ip3, jp3) + D4;
        }
    }
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            (*dx)(i, j) = DX(i, j);
        }
    }
    for (int i = 0; i < 6; i++) {
        (*q1)(i) = X(i);
    }
};
