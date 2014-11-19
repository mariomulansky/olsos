/* Ensemble of coupled Roessler attractors
 * divide into clusters with double-calculation at the edges
 * and employ SIMD vectorization
 * performance gain of roughly x3
 */


#include <iostream>
#include <vector>
#include <random>

#include <boost/timer.hpp>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/simd/sdk/simd/pack.hpp>
#include <boost/simd/sdk/simd/io.hpp>
#include <boost/simd/memory/allocator.hpp>
#include <boost/simd/include/functions/splat.hpp>
#include <boost/simd/include/functions/plus.hpp>
#include <boost/simd/include/functions/multiplies.hpp>

namespace odeint = boost::numeric::odeint;
namespace simd = boost::simd;

typedef boost::timer timer_type;

typedef simd::pack<double> simd_pack;
typedef std::vector< simd_pack, simd::allocator<simd_pack> > state_type;

static const size_t dim = 3;  // roessler is 3D
static const size_t pack_size = simd_pack::static_size;

static const int overlap = 4;
typedef boost::array<simd_pack, 2*overlap*3> array_type1;
typedef boost::array<simd_pack, overlap*3> array_type2;

//---------------------------------------------------------------------------
struct roessler_system {
    const double m_a, m_b, m_c;
    const double m_d;

    roessler_system(const double a, const double b, const double c,
                           const double d)
        : m_a(a), m_b(b), m_c(c), m_d(d)
    {}

    template< class State , class Deriv >
    void operator()( const State &x_ , Deriv &dxdt_ , double t ) const
    {
        typename boost::range_iterator< const State >::type x = boost::begin( x_ );
        typename boost::range_iterator< Deriv >::type dxdt = boost::begin( dxdt_ );

        const int N = boost::size(x_);
        // left boundary -- can be neglected as it will be wrong anyways
        dxdt[0] = boost::simd::splat<simd_pack>(0.0);
        dxdt[1] = boost::simd::splat<simd_pack>(0.0);
        dxdt[2] = boost::simd::splat<simd_pack>(0.0);
        // center loop
        for( int j=1; j<N/dim-1; ++j )
        {
            const int i = j*dim;
            dxdt[i] = -1.0*x[i + 1] - x[i + 2] +
                      m_d * (x[i - dim] + x[i + dim] - 2.0 * x[i]);
            dxdt[i + 1] = x[i] + m_a * x[i + 1];
            dxdt[i + 2] = m_b + x[i + 2] * (x[i] - m_c);
        }
        // right boundary -- can be neglected as it will be wrong anyways
        dxdt[N - dim] = boost::simd::splat<simd_pack>(0.0);
        dxdt[N - dim + 1] = boost::simd::splat<simd_pack>(0.0);
        dxdt[N - dim + 2] = boost::simd::splat<simd_pack>(0.0);
    }
};

template< typename Array>
void save_state(const state_type &state, const size_t index, Array &save)
{
    const size_t len = save.size();
    std::copy(state.begin() + index, state.begin() + index + len, save.begin());
}

template< typename Array>
void load_state(state_type &state, const size_t index, const Array &save)
{
    std::copy(save.begin(), save.end(), state.begin() + index);
}

void rotate_simd_boundaries(state_type &state)
{
    const size_t N = state.size() - 2*overlap*dim;
    for(int p=0; p<pack_size; ++p)
    {
        for(int i=0; i<overlap*dim; ++i)
        {
            state[i][p] = state[N+i][(p+pack_size-1) % pack_size];
            state[N+overlap*dim+i][p] = state[overlap*dim+i][(p+1) % pack_size];
        }
    }
}

//---------------------------------------------------------------------------
int main(int argc, char *argv[]) {
if(argc<4)
{
    std::cerr << "Expected size, steps and cluster size as parameter" << std::endl;
    exit(1);
}
const size_t n = atoi(argv[1]);
const int steps = atoi(argv[2]);
const size_t clusters = atoi(argv[3]);

const size_t state_size = n*dim/pack_size;
// we need some extra space at the start and end for the overlap
const size_t work_size = state_size + 2*overlap*dim;
const size_t cluster_size = clusters*dim;
const size_t cluster_count = state_size/cluster_size;

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

state_type state(work_size);
for(size_t i=0; i<n/pack_size; ++i)
{
    for(size_t p=0; p<pack_size; ++p)
    {
        state[(overlap+i)*dim][p] = x[p*n/pack_size+i];
        state[(overlap+i)*dim+1][p] = y[p*n/pack_size+i];
        state[(overlap+i)*dim+2][p] = z[p*n/pack_size+i];
    }
}

std::cout << "Systems: " << n << std::endl;
std::cout << "SIMD pack size: " << pack_size << std::endl;
std::cout << "State size: " << state_size << std::endl;
std::cout << "Work size: " << work_size << std::endl;
std::cout << "Cluster size: " << cluster_size << std::endl;
std::cout << "Cluster count: " << cluster_count << std::endl;

std::cout.precision(16);

std::cout << state[overlap*dim] << std::endl;
std::cout << state[state_size] << std::endl << std::endl;

// Stepper type
odeint::runge_kutta4_classic<state_type, double, state_type, double,
                             odeint::range_algebra, odeint::default_operations,
                             odeint::never_resizer> stepper;
// adjust size manually in the beginning
stepper.adjust_size(state_type(cluster_size + 2*overlap*dim));

roessler_system sys(a, b, c, d);

timer_type timer;

double t = 0.0;

array_type1 saved1;
array_type2 saved2;

for(int step = 0; step < steps; step++)
{
    rotate_simd_boundaries(state);
    save_state(state, cluster_size, saved1);
    stepper.do_step(
        sys, std::make_pair(state.begin(),
                            state.begin() + cluster_size + 2*overlap*dim),
        t, dt);
    save_state(state, cluster_size, saved2);
    for(int c=1; c<cluster_count-1; ++c)
    {
        load_state(state, c*cluster_size, saved1);
        save_state(state, (c+1)*cluster_size, saved1);
        stepper.do_step(
            sys, std::make_pair(state.begin() + c * cluster_size,
                                state.begin() + (c + 1) * cluster_size + 2*overlap*dim),
            t, dt);
        load_state(state, c*cluster_size, saved2);
        save_state(state, (c+1)*cluster_size, saved2);
    }
    load_state(state, (cluster_count-1)*cluster_size, saved1);
    stepper.do_step(
        sys, std::make_pair(state.begin() + (cluster_count-1)*cluster_size,
                            state.end()),
        t, dt);
    load_state(state, (cluster_count-1) * cluster_size, saved2);
    t += dt;
}

std::cout << "Integration finished, runtime for " << steps << " steps: ";
std::cout << timer.elapsed() << " s" << std::endl;

// compute some accumulation to make sure all results have been computed
simd_pack s_pack = 0.0;
for(size_t i = 0; i < n/pack_size; ++i)
{
    s_pack += state[overlap*dim+i*dim];
}
double s = 0.0;
for(size_t p=0; p<pack_size; ++p)
{
    s += s_pack[p];
}

std::cout << state[overlap*dim] << std::endl;
std::cout << s/n << std::endl;

}
