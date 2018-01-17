#include <array>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>

double rand_Num();
double r;

class Particle {

  public:
    std::array<double, 3> get_pos() const { return m_pos; }
    std::array<double, 3> get_norm() const { return m_norm; }

  protected:
    std::array<double, 3> m_pos;
    std::array<double, 3> m_norm;
    /*double xpos
    double ypos
    double zpos*/
    double theta;
    double phi;
    double nx;
    double ny;
    double nz;
    Particle() /*: m_pos({{0.0, 0.0, 0.0}}), theta(acos(2.0*rand_Num()-1.0)),
                phi(2.0*M_PI*rand_Num()), m_norm({{sin(theta)*cos(phi),
                sin(theta)*sin(phi), cos(theta)}})*/
    {
      // std::cout << "The initial random theta is:	" << theta << "\n";
      // std::cout << "The initial random phi is:	" << phi << "\n";
      // std::cout << "The initial normals are:	" << m_norm[0] << ", " <<
      // m_norm[1] << ", " << m_norm[2] << "\n";  std::cout << "The initial 1
      // position is:	" << m_pos[0] << ", " << m_pos[1] << ", " << m_pos[2] <<
      // "\n";
    }
};

class Photon : public Particle {
    double wavelength;
    double theta;
    double phi;
    std::array<double, 3> phot_norm;
    std::array<double, 3> phot_pos;
    double taumax;
    double L;
    double tau;
    double rad;
    double total;
    int nphoton;
    double lcount;
    int count;

  public:
    std::array<double, 3> get_phot_pos() const { return phot_pos; }
    std::array<double, 3> get_phot_norm() const { return phot_norm; }

    Photon(double t_wavelength, double t_rad)
        : Particle(), phot_norm({{0.0, 0.0, 0.0}}), phot_pos({{0.0, 0.0, 0.0}}),
          wavelength(t_wavelength), rad(t_rad) {}

    void scatter() {
      theta = acos(2.0 * rand_Num() - 1.0);
      phi = 2.0 * M_PI * rand_Num();
      phot_norm = {{sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)}};
    }

    void move(double taumax) {
      tau = -log(rand_Num());
      L = tau * (rad / taumax);
      // lcount += L;
      phot_pos[0] += L * phot_norm[0];
      phot_pos[1] += L * phot_norm[1];
      phot_pos[2] += L * phot_norm[2];

      // std::cout << "In move position:	" << phot_pos[0] << " " << phot_pos[1]
      // << " " << phot_pos[2] << "\n";
    }
};

double rand_Num() {

  double r = (static_cast<double>(rand()) / RAND_MAX);
  // std::cout << r << " " << "\n";
  return r;
}

int main() {
  srand(time(nullptr));
  // std::array<double, 3> m_pos, m_norm;
  int nphoton = 1e6;
  // int total = 0;
  double rad = 1.0;
  double taumax;
  std::array<double, 5> taus({{0.1, 0.5, 1.0, 5.0, 10.0}});

  // double nx, ny, nz;

  for (int j = 0; j < 5; j++) {
    taumax = taus[j];
    int total(0);
    /*std::cout << "Tau max:	" << taumax << "\n";*/

    for (int i = 0; i < nphoton; ++i) {
      /*std::cout << " " << "\n";
      std::cout << " " << "\n";
      std::cout << "Photon number:	" << i << "\n";
      std::cout << " " << "\n";
      std::cout << " " << "\n";*/

      Photon phot(0.0, 1.0);
      phot.scatter();
      std::array<double, 3> pos = phot.get_phot_pos();
      std::array<double, 3> norm = phot.get_phot_norm();
      // std::cout << "Photon initial position:	" << pos[0] << " " <<
      // pos[1] << " " << pos[2] << "\n";  std::cout << "Photon initial
      // direction: " << norm[0] << " " << norm[1] << " " << norm[2] << "\n";

      while (((pos[0] * pos[0]) + (pos[1] * pos[1]) + (pos[2] * pos[2])) <=
          (rad * rad)) {
        phot.move(taumax);
        pos = phot.get_phot_pos();

        // std::cout << "Pos:	" << pos[0] << " " << pos[1] << " " << pos[2]
        // << "\n";

        if (((pos[0] * pos[0]) + (pos[1] * pos[1]) + (pos[2] * pos[2])) >
            (rad * rad)) {
          // std::cout << "Final squared pos:	" << (pos[0]*pos[0]) +
          // (pos[1]*pos[1]) + (pos[2]*pos[2]) << "\n";  std::cout << "EXIT" <<
          // "\n";  std::cout << "Total" << total << "\n";
          break;
        }

        // std::cout << "Total:	" << i << " " << total << "\n";

        total++;
        phot.scatter();
        norm = phot.get_phot_norm();
        // std::cout << "Scatter total:	" << total << "\n";
      }
    }

    /*std::cout << "Average number of scatterings before escaping:	"
              << static_cast<double>(total) / static_cast<double>(nphoton)
              << "\n";
    std::cout
        << "Predicted average number of scatterings before escaping:	"
        << taumax + (taumax * taumax) / 2.0 << "\n";*/
  }
  return 0;
}
