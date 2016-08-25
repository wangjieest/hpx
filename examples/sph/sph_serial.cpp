//  Copyright (c) 2016 Thomas Heller
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <hpx/hpx_init.hpp>

#include <hpx/util/high_resolution_timer.hpp>
#include <hpx/util/zip_iterator.hpp>

#include <boost/program_options.hpp>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <tuple>
#include <vector>

#ifndef SPH_PI
#define SPH_PI 3.14159265358979323846f
#endif

std::size_t box_indicator(float x, float y)
{
    return (x < 0.25f) && (y < 1.f);
}

std::size_t count_particles(float hh)
{
    std::size_t res = 0;

    for (float x = 0.0f; x < 1.0f; x += hh)
    {
        for (float y = 0.0f; y < 1.0f; y += hh)
        {
            res += box_indicator(x, y);
        }
    }
    return res;
}

void place_particles(
    std::vector<float>& pos_x,
    std::vector<float>& pos_y,
    float hh)
{
    std::size_t index = 0;
    for (float x = 0.f; x < 1.f; x += hh) {
        for (float y = 0.f; y < 1.f; y += hh) {
            if (box_indicator(x,y)) {
                pos_x[index] = x;
                pos_y[index] = y;
                ++index;
            }
        }
    }
    HPX_ASSERT(pos_x.size() == index);
    HPX_ASSERT(pos_y.size() == index);
}

float normalize_mass(float mass, std::vector<float> const& rho, float rho0)
{
    float rho_squared_sum = 0.f;
    float rho_sum = 0.f;

    for (float r: rho)
    {
        rho_squared_sum += r * r;
        rho_sum += r;
    }

    return mass * rho0 * rho_sum / rho_squared_sum;
}

void dump_timestep(std::size_t t,
    std::vector<float> const& pos_x, std::vector<float> const& pos_y)
{
    std::stringstream ss;
    ss << "output" << std::setw(10) << std::setfill('0') << t << ".csv";
    std::string filename(ss.str());

    std::ofstream of(filename.c_str());

    float rho = 1.0;

    for (std::size_t i = 0; i < pos_x.size(); ++i)
    {
        of << pos_x[i] << ", " << pos_y[i] << ", 0.0, " << rho << "\n";
    }

    of.flush();
}

std::vector<float> compute_density(
    std::vector<float> const& pos_x,
    std::vector<float> const& pos_y,
    float h,
    float mass)
{

    float const h_squared = h * h;
    float const h_pow_8 = h_squared * h_squared * h_squared * h_squared;
    float const C = 4.f * mass / SPH_PI / h_pow_8;

    std::vector<float> rho(pos_x.size(), 4.f * mass / SPH_PI / h_squared);

    std::size_t count = rho.size();

    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t j = i + 1; j < count; ++j)
        {
            float delta_x = pos_x[i] - pos_x[j];
            float delta_y = pos_y[i] - pos_y[j];
            float dist_squared = delta_x * delta_x + delta_y * delta_y;
            float overlap = h_squared - dist_squared;

            if (overlap > 0.f)
            {
                float rho_ij = C * overlap * overlap * overlap;
                rho[i] += rho_ij;
                rho[j] += rho_ij;
            }
        }
    }

    return rho;
}

std::tuple<std::vector<float>, std::vector<float>> compute_accel(
    std::vector<float> const& pos_x,
    std::vector<float> const& pos_y,
    std::vector<float> const& v_x,
    std::vector<float> const& v_y,
    std::vector<float> const& rho,
    float mass, float g, float h, float k, float rho0, float mu)
{
    std::vector<float> a_x(pos_x.size(), 0.f);
    std::vector<float> a_y(pos_x.size(), -g);

    float const h_squared = h * h;
    float const C_0 = mass / SPH_PI / (h_squared * h_squared);
    float const C_p = 15.f * k;
    float const C_v = -40.f * mu;

    // Now compute interaction forces
    std::size_t count = pos_x.size();

    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t j = i + 1; j < count; ++j)
        {
            float delta_x = pos_x[i] - pos_x[j];
            float delta_y = pos_y[i] - pos_y[j];
            float dist_squared = delta_x * delta_x + delta_y * delta_y;

            if (dist_squared < h_squared)
            {
                float q = std::sqrt(dist_squared) / h;
                float u = 1 - q;
                float w_0 = C_0 * u / rho[i] / rho[j];
                float w_p = w_0 * C_p * (rho[i] + rho[j] - 2 * rho0) * u / q;
                float w_v = w_0 * C_v;
                float delta_v_x = v_x[i] - v_y[j];
                float delta_v_y = v_y[i] - v_y[j];
                a_x[i] += (w_p * delta_x + w_v * delta_v_x);
                a_y[i] += (w_p * delta_y + w_v * delta_v_y);
                a_x[j] -= (w_p * delta_x + w_v * delta_v_x);
                a_y[j] -= (w_p * delta_y + w_v * delta_v_y);
            }
        }
    }

    return std::make_tuple(a_x, a_y);
}

void leapfrog(
    std::vector<float>& pos_x,
    std::vector<float>& pos_y,
    std::vector<float>& v_x,
    std::vector<float>& v_y,
    std::vector<float> const& a_x,
    std::vector<float> const& a_y,
    float dt)
{
    for (std::size_t i = 0; i < pos_x.size(); ++i)
    {
        v_x[i] += a_x[i] * dt;
        v_y[i] += a_y[i] * dt;

        pos_x[i] += v_x[i] * dt;
        pos_y[i] += v_y[i] * dt;
    }
}

void damp_reflect(
    int which,
    float barrier,
    float& pos_x,
    float& pos_y,
    float& v_x,
    float& v_y)
{
    float& v_which   = (which == 0) ? v_x   : v_y;
    float& pos_which = (which == 0) ? pos_x : pos_y;

    // Coefficient of resitiution
    const float DAMP = 0.75;
    // Ignore degenerate cases
    if (std::abs(v_which) <= 1e-3f)
        return;

    // Scale back the distance traveled based on time from collision
    float tbounce = (pos_which - barrier) / v_which;
    pos_x -= v_x*(1-DAMP)*tbounce;
    pos_y -= v_y*(1-DAMP)*tbounce;

    // Reflect the position and velocity
    pos_which = 2 * barrier - pos_which;
    v_which   = -v_which;

    // Damp the velocities
    v_x *= DAMP;
    v_y *= DAMP;
}

void reflect_bc(
    std::vector<float>& pos_x,
    std::vector<float>& pos_y,
    std::vector<float>& v_x,
    std::vector<float>& v_y)
{
    // Boundaries of the computational domain
    float const XMIN = 0.0;
    float const XMAX = 1.0;
    float const YMIN = 0.0;
    float const YMAX = 1.0;

    auto it = hpx::util::make_zip_iterator(
        pos_x.begin(),
        pos_y.begin(),
        v_x.begin(),
        v_y.begin());
    auto end = hpx::util::make_zip_iterator(
        pos_x.end(),
        pos_y.end(),
        v_x.end(),
        v_y.end());

    //for (auto& v: std::make_pair(begin, end))
    for(; it != end; ++it)
    {
        auto& px = hpx::util::get<0>(*it);
        auto& py = hpx::util::get<1>(*it);
        auto& vx = hpx::util::get<2>(*it);
        auto& vy = hpx::util::get<3>(*it);

        if (px < XMIN) {
            damp_reflect(0, XMIN, px, py, vx, vy);
        }
        if (px > XMAX) {
            damp_reflect(0, XMAX, px, py, vx, vy);
        }
        if (py < YMIN) {
            damp_reflect(1, YMIN, px, py, vx, vy);
        }
        if (py > YMAX) {
            damp_reflect(1, YMAX, px, py, vx, vy);
        }
    }
}

int hpx_main(boost::program_options::variables_map& vm)
{
    float dt = vm["dt"].as<float>();
    float h = vm["pitch"].as<float>();
    float rho0 = vm["rho0"].as<float>();
    float k = vm["bulk-modulus"].as<float>();
    float mu = vm["viscosity"].as<float>();
    float g = vm["gravitation"].as<float>();
    float mass = vm["mass"].as<float>();
    std::size_t timesteps = vm["timesteps"].as<std::size_t>();
    std::size_t io_period = vm["io-period"].as<std::size_t>();


    float hh = h / 1.3f;

    std::size_t count = count_particles(hh);

    std::cout << "Running with " << count << " particles\n";

    hpx::util::high_resolution_timer t;
    std::vector<float> pos_x(count);
    std::vector<float> pos_y(count);
    std::vector<float> v_x(count, 0.f);
    std::vector<float> v_y(count, 0.f);

    place_particles(pos_x, pos_y, hh);

    std::vector<float> rho = compute_density(pos_x, pos_y, h, mass);
    mass = normalize_mass(mass, rho, rho0);

    for (std::size_t t = 0; t != timesteps; ++t)
    {
        if ((t % io_period) == 0)
        {
            dump_timestep(t, pos_x, pos_y);
        }

        rho = compute_density(pos_x, pos_y, h, mass);

        std::vector<float> a_x;
        std::vector<float> a_y;

        std::tie(a_x, a_y) = compute_accel(pos_x, pos_y, v_x, v_y, rho,
            mass, g, h, k, rho0, mu);

        leapfrog(pos_x, pos_y, v_x, v_y, a_x, a_y, dt);

        reflect_bc(pos_x, pos_y, v_x, v_y);
    }
    double elapsed = t.elapsed();

    dump_timestep(timesteps, pos_x, pos_y);

    std::cout << "Elapsed seconds: " << elapsed << "\n";

    return hpx::finalize();
}

int main(int argc, char** argv)
{
    using namespace boost::program_options;

    options_description desc_commandline;
    desc_commandline.add_options()
        ("dt", value<float>()->default_value(1e-4),
         "Timestep length")
        ("pitch", value<float>()->default_value(2e-2),
         "Pitch")
        ("rho0", value<float>()->default_value(1000),
         "Target density")
        ("bulk-modulus", value<float>()->default_value(1e3),
         "Bulk modulus")
        ("viscosity", value<float>()->default_value(0.1),
         "Viscosity")
        ("gravitation", value<float>()->default_value(9.8),
         "Gravitational acceleration")
        ("mass", value<float>()->default_value(1.f),
         "Total mass")
        ("timesteps", value<std::size_t>()->default_value(20000),
         "Number of timesteps")
        ("io-period", value<std::size_t>()->default_value(20000),
         "Frequency of IO")
    ;

    // Initialize and run HPX, this example requires to run hpx_main on all
    // localities
    std::vector<std::string> const cfg = {
        "hpx.run_hpx_main!=1"
    };

    return hpx::init(desc_commandline, argc, argv, cfg);
}
