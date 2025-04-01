#include"Sourse.h"

static Matrix<double, 1, 9> read_settings() {
	std::ifstream file("SETTINGS/CERTAIN_AN/settings.txt");
	std::string line;
	Matrix<double, 1, 9> set;
	int counter = 0;
	while (std::getline(file, line)) {
		set(counter) = stod(line);
		counter++;
	}
	file.close();
	return set;
}

static void read_receivers(MatrixXd* rec_nums) {
	std::ifstream file("SETTINGS/CERTAIN_AN/receivers.txt");
	std::string line;
	int counter = 0;
	while (std::getline(file, line)) {
		(*rec_nums)(counter) = stoi(line);
		counter++;
	}
	file.close();
}

static void read_transmitters(MatrixXd* tr_nums) {
	std::ifstream file("SETTINGS/CERTAIN_AN/transmitters.txt");
	std::string line;
	int counter = 0;
	while (std::getline(file, line)) {
		(*tr_nums)(counter) = stoi(line);
		counter++;
	}
	file.close();
}

static bool visibility_condition(Matrix<double, 6, 1> q1, Matrix<double, 6, 1> q2, double R_e)
{
	double r1, r2, r12, a, b, R, M_PI = 3.14159265358979323846;
	r1 = sqrt(pow(q1(0), 2) + pow(q1(1), 2) + pow(q1(2), 2));
	r2 = sqrt(pow(q2(0), 2) + pow(q2(1), 2) + pow(q2(2), 2));
	r12 = sqrt(pow(q1(0) - q2(0), 2) + pow(q1(1) - q2(1), 2) + pow(q1(2) - q2(2), 2));
	a = acos((pow(r1, 2) + pow(r2, 2) - pow(r12, 2)) / (2 * r1 * r2));
	b = (pow(r1, 2) + pow(r12, 2) - pow(r2, 2)) / (2 * r1 * r12);
	R = r1 * sqrt(1 - pow(b, 2));

	if (R > R_e || a < M_PI / 2) return true;
	return false;
}

static bool goto_next_epoch(MatrixXd receivers_correction, double eps, int size)
{
	for (int i = 0; i < size; i++) {
		if (receivers_correction(i) > eps) return false;
	}
	return true;
}

static void make_observations(MatrixXd transmitters, MatrixXd receivers, map<int, Matrix<double, 6, 6>>* Q, map<int, Matrix<double, 6, 1>>* G, Matrix<double, 1, 9> settings, int id, MatrixXd rec_nums, MatrixXd tr_nums)
{
	bool visibility_condition(Matrix<double, 6, 1>, Matrix<double, 6, 1>, double);
	std::default_random_engine generator(static_cast<unsigned int>(std::time(0)));
	std::normal_distribution<double> distribution(0.0, settings(3));
	Matrix<double, 6, 1> qrt, qrr, qcr;
	Matrix<double, 6, 1> qrt0, qrr0, qcr0;
	Matrix<double, 1, 6> A;
	Matrix<double, 1, 3> d_p;
	Matrix<double, 6, 6> DX;
	Matrix<double, 3, 6> dx_need;
	int start, end;
	int size_t = settings(6);
	int size_r = settings(7);
	int num_thread = settings(8);
	if (size_r % num_thread == 0) {
		start = id * (size_r / num_thread) + 1;
		end = (id + 1) * (size_r / num_thread);
	}
	else {
		int diff = size_r % num_thread;
		if (id < diff) {
			start = id * floor(size_r / num_thread + 1) + 1;
			end = (id + 1) * floor(size_r / num_thread + 1);
		}
		else {
			start = id * floor(size_r / num_thread) + 1 + diff;
			end = (id + 1) * floor(size_r / num_thread) + diff;
		}
	}
	double t, d_o, d_c;
	int tr, re;
	for (int n = 1; n <= settings(1); n++) {
		t = n * settings(0);
		for (int rec = start - 1; rec < end; rec++) {
			re = rec_nums(rec);
			for (int sat = 0; sat < size_t; sat++) {
				tr = tr_nums(sat);
				if (re != tr) {
					for (int i = 0; i < 6; i++) {
						qrt0(i) = transmitters(tr-1, i);
						qrr0(i) = transmitters(re-1, i);
					}
					KeplerX(qrt0, t, &qrt); // Reference transmitter
					KeplerX(qrr0, t, &qrr); // Reference receiver
					
					if (visibility_condition(qrt, qrr, settings(4))) {
						for (int i = 0; i < 6; i++) {
							qcr0(i) = receivers(rec, i);
						}
						KeplerXDX(qcr0, t, &qcr, &DX); // Calculated receiver
						for (int i = 0; i < 3; i++) {
							for (int j = 0; j < 6; j++) {
								dx_need(i, j) = DX(i, j);
							}
						}
						d_o = sqrt(pow(qrt(0) - qrr(0), 2) + pow(qrt(1) - qrr(1), 2) + pow(qrt(2) - qrr(2), 2)) + distribution(generator) / 1000.0;
						d_c = sqrt(pow(qcr(0) - qrt(0), 2) + pow(qcr(1) - qrt(1), 2) + pow(qcr(2) - qrt(2), 2));
						d_p = { qcr(0) - qrt(0), qcr(1) - qrt(1), qcr(2) - qrt(2) };
						A = (d_p / d_c) * dx_need;
						(*Q)[rec] += A.transpose() * A;
						(*G)[rec] += A.transpose() * (d_o - d_c);
					}
				}
			}
		}
	}
}

static void move_the_epoch(MatrixXd*receivers, MatrixXd*transmitters, map<int, Matrix<double, 6, 6>>* Q, map<int, Matrix<double, 6, 1>>* G, Matrix<double, 1, 9> settings)
{	
	int size_r = settings(7);
	int size_t = settings(6);
	Matrix<double, 6, 1> q0, q1;
	Matrix<double, 6, 6> dx, dx_inv;
	for (int i = 0; i < size_t; i++) {
		for (int j = 0; j < 6; j++) {
			q0(j) = (*transmitters)(i, j);
		}
		KeplerX(q0, settings(2), &q1);
		for (int j = 0; j < 6; j++) {
			(*transmitters)(i, j) = q1(j);
		}
	}
	for (int i = 0; i < size_r; i++){
		for (int j = 0; j < 6; j++) {
			q0(j) = (*receivers)(i, j);
		}
		KeplerXDX(q0, settings(2), &q1, &dx);
		dx_inv = dx.inverse();
		(*Q)[i] = dx_inv.transpose() * (*Q)[i] * dx_inv;
		(*G)[i] = dx_inv.transpose() * (*G)[i];
		for (int j = 0; j < 6; j++) {
			(*receivers)(i, j) = q1(j);
		}
	}
}

static void make_corrections(MatrixXd*corrections, map<int, Matrix<double, 6, 6>>Q, map<int, Matrix<double, 6, 1>>G, MatrixXd*receivers, Matrix<double, 1, 9> settings)
{	
	Matrix<double, 1, 6> scales = { 10000, 10000, 10000, 1, 1, 1 };
	map<int, Matrix<double, 6, 6>> Q_norm = Q;
	map<int, Matrix<double, 6, 1>> G_norm = G;
	int size_r = settings(7);
	Matrix<double, 6, 1> d_q;
	for (int n = 0; n < size_r; n++) {
		for (int i = 0; i < 6; i++) {
			G_norm[n](i) *= scales(i);
			for (int j = 0; j < 6; j++) {
				Q_norm[n](i, j) *= scales(i) * scales(j);
			}
		}
		d_q = (Q_norm[n].inverse()) * G_norm[n];
		for (int i = 0; i < 6; i++) {
			d_q(i) *= scales(i);
		}
		for (int j = 0; j < 6; j++) {
			(*receivers)(n, j) += d_q(j);
		}
		(*corrections)(n) = d_q.norm();
	}
}

static int run_epoch(MatrixXd* transmitters, MatrixXd* receivers, map<int, Matrix<double, 6, 6>>*Q, map<int, Matrix<double, 6, 1>>*G, Matrix<double, 1, 9> settings, MatrixXd*corrections)
{
	void make_observations(MatrixXd, MatrixXd, map<int, Matrix<double, 6, 6>>*, map<int, Matrix<double, 6, 1>>*, Matrix<double, 1, 9>, int, MatrixXd, MatrixXd);
	bool goto_next_epoch(MatrixXd, double, int);
	void make_corrections(MatrixXd*, map<int, Matrix<double, 6, 6>>, map<int, Matrix<double, 6, 1>>, MatrixXd*, Matrix<double, 1, 9>);
	void move_the_epoch(MatrixXd*, MatrixXd*, map<int, Matrix<double, 6, 6>>*, map<int, Matrix<double, 6, 1>>*, Matrix<double, 1, 9>);

	double eps = settings(5);
	int size_t = settings(6);
	int size_r = settings(7);
	int n_threads = settings(8);
	map<int, Matrix<double, 6, 6>> Q_old;
	map<int, Matrix<double, 6, 1>> G_old;
	int n_iter = 0;
	MatrixXd rec_nums(0, 0), tr_nums(0, 0);
	rec_nums.resize(1, size_r);
	tr_nums.resize(1, size_t);
	read_receivers(&rec_nums);
	read_transmitters(&tr_nums);
	while (not goto_next_epoch(*corrections, eps, size_r)) {
		std::vector<std::thread> threads;
		Q_old = *Q;
		G_old = *G;
		for (int id = 0; id < n_threads; id++) {
			threads.emplace_back(make_observations, *transmitters, *receivers, &Q_old, &G_old, settings, id, rec_nums, tr_nums);
		}
		for (auto& thread : threads) {
			thread.join();
		}
		n_iter++;
		make_corrections(corrections, Q_old, G_old, receivers, settings);
		if (goto_next_epoch(*corrections, eps, size_r)) {
			(*Q) = Q_old;
			(*G) = G_old;
			move_the_epoch(receivers, transmitters, Q, G, settings);
		}
	}
	return n_iter;
}

extern "C" __declspec(dllexport) int reform_start_data_and_run_epoch(double*p_transmitters, double*p_receivers, double*normal_matrix, double*gradient, double*corrections) 
{	
	Matrix<double, 1, 9> read_settings();
	void read_receivers(MatrixXd*);
	int run_epoch(MatrixXd*, MatrixXd*, map<int, Matrix<double, 6, 6>>*, map<int, Matrix<double, 6, 1>>*, Matrix<double, 1, 9>, MatrixXd*);

	Matrix<double, 1, 9> settings = read_settings();
	int size_t = settings(6);
	int size_r = settings(7);

	MatrixXd receivers(0, 0), transmitters(0, 0);
	receivers.resize(size_r, 6);
	transmitters.resize(size_t, 6);

	for (int i = 0; i < size_t; i++) {
		for (int j = 0; j < 6; j++) {
			transmitters(i, j) = p_transmitters[6 * i + j];
		}
	} // Unpacking transmitters
	for (int i = 0; i < size_r; i++) {
		for (int j = 0; j < 6; j++) {
			receivers(i, j) = p_receivers[6 * i + j];
		}
	} // Unpacking receivers
	map<int, Matrix<double, 6, 6>> Q;
	map<int, Matrix<double, 6, 1>> G;
	for (int i = 0; i < size_r; i++) {
		for (int j = 0; j < 6; j++) {
			for (int k = 0; k < 6; k++) {
				Q[i](j, k) = normal_matrix[36 * i + j * 6 + k];
			}
			G[i](j) = gradient[i * 6 + j];
		}
	} // Unpacking normal matrixes and gradients
	MatrixXd corrections_norm;
	corrections_norm.resize(1, size_r);
	for (int i = 0; i < size_r; i++) {
		corrections_norm(i) = corrections[i];
	} // Unpacking corrections

	int n_iter = run_epoch(&transmitters, &receivers, &Q, &G, settings, &corrections_norm);

	for (int i = 0; i < size_t; i++) {
		for (int j = 0; j < 6; j++) {
			p_transmitters[6 * i + j] = transmitters(i, j);
		}
	} // Collecting transmitters parameters
	for (int i = 0; i < size_r; i++) {
		for (int j = 0; j < 6; j++) {
			p_receivers[6 * i + j] = receivers(i, j);
		}
	} // Collecting receivers parameters
	for (int i = 0; i < size_r; i++) {
		for (int j = 0; j < 6; j++) {
			for (int k = 0; k < 6; k++) {
				normal_matrix[36 * i + j * 6 + k] = Q[i](j, k);
			}
			gradient[i * 6 + j] = G[i](j);
		}
	} // Collecting normal matrixes and gradients
	for (int i = 0; i < size_r; i++) {
		corrections[i] = corrections_norm(i);
	} // Collectiong corrections
	return n_iter;
}
