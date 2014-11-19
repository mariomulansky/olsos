/* Ensemble of coupled Roessler attractors
 * straight forward implementation, suffering from memory bandwidth bottleneck
 */


#include <iostream>
#include <vector>
#include <utility>
#include <tuple>
#include <memory>
#include <random>

#include <boost/timer.hpp>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

namespace odeint = boost::numeric::odeint;

typedef boost::timer timer_type;

typedef std::vector<double> state_type;

//---------------------------------------------------------------------------
struct roessler_system {
    const double m_a, m_b, m_c;
    const double m_d;
    static const int dim = 3;

    roessler_system(const double a, const double b, const double c,
                           const double d)
        : m_a(a), m_b(b), m_c(c), m_d(d)
    {}

    void operator()(const state_type &x, state_type &dxdt, const double t) const
    {
        const int N = x.size();
        // left boundary (periodic)
        dxdt[0] = -x[1] - x[2] + m_d * (x[dim] + x[N - dim] - 2 * x[0]);
        dxdt[1] = x[0] + m_a * x[1];
        dxdt[2] = m_b + x[2] * (x[0] - m_c);
        // center loop
        for( int j=1; j<N/dim-1; ++j )
        {
            const int i = j*dim;
            dxdt[i] = -x[i + 1] - x[i + 2] +
                      m_d * (x[i - dim] + x[i + dim] - 2 * x[i]);
            dxdt[i + 1] = x[i] + m_a * x[i + 1];
            dxdt[i + 2] = m_b + x[i + 2] * (x[i] - m_c);
        }
        // right boundary (periodic)
        dxdt[N - dim] = -x[N - dim + 1] - x[N - dim + 2] +
                        m_d * (x[N - 2 * dim] + x[0] - 2 * x[N - dim]);
        dxdt[N - dim + 1] = x[N - dim] + m_a * x[N - dim + 1];
        dxdt[N - dim + 2] = m_b + x[N - dim + 2] * (x[N - dim] - m_c);
    }
};

//---------------------------------------------------------------------------
int main(int argc, char *argv[]) {
if(argc<3)
{
    std::cerr << "Expected size and steps as parameter" << std::endl;
    exit(1);
}
const size_t n = atoi(argv[1]);
const size_t steps = atoi(argv[2]);

const double dt = 0.01;

const double a = 0.2;
const double b = 1.0;
const double c = 9.0;
const double d = 1.0;

// random initial conditions on the device
std::vector<double> x(n), y(n), z(n);
std::default_random_engine generator;
std::uniform_real_distribution<double> distribution_xy(-8.0, 8.0);
std::uniform_real_distribution<double> distribution_z(0.0, 20.0);
auto rand_xy = std::bind(distribution_xy, std::ref(generator));
auto rand_z = std::bind(distribution_z, std::ref(generator));
std::generate(x.begin(), x.end(), rand_xy);
std::generate(y.begin(), y.end(), rand_xy);
std::generate(z.begin(), z.end(), rand_z);

state_type state(3*n);
for(size_t i=0; i<n; ++i)
{
    state[3*i] = x[i];
    state[3*i+1] = y[i];
    state[3*i+2] = z[i];
}

std::cout << "# n: " << n << std::endl;

// Stepper type - use never_resizer for slight performance improvement
odeint::runge_kutta4_classic<state_type, double, state_type, double,
                             odeint::range_algebra,
                             odeint::default_operations,
                             odeint::never_resizer> stepper;
// adjust size manually in the beginning
stepper.adjust_size(state);

roessler_system sys(a, b, c, d);

timer_type timer;

double t = 0.0;

std::cout.precision(16);

//std::cout << state[n/2+4] << std::endl;

// transient
for(int step = 0; step < steps; step++)
{
    stepper.do_step(sys, state, t, dt);
    t += dt;
    // std::cout << state[512+4] << std::endl;
}

std::cout << "Integration finished, runtime for " << steps << " steps: ";
std::cout << timer.elapsed() << " s" << std::endl;

// compute some accumulation to make sure all results have been computed
double s = 0.0;
for(size_t i = 0; i < n; ++i)
{
    s += state[3*i];
    // std::cout << state[i] << " , ";
}
std::cout << std::endl;

std::cout << state[0] << std::endl;
std::cout << s/n << std::endl;

}
