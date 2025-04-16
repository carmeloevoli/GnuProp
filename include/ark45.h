#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <type_traits>
#include <vector>

// ---------------------------------------------------------------------
// A generic adaptive Runge-Kutta 5(4) Dormand–Prince ODE solver
//
// This version is designed for C++14 compatibility, removing usage
// of 'if constexpr'. We use template SFINAE helpers to handle
// scalars vs. vectors.
//
// Usage examples are given under `#ifdef RUN_EXAMPLE_MAIN`.
// ---------------------------------------------------------------------

template <typename State>
class AdaptiveRK45 {
 public:
  using OdeFunction = std::function<State(double, const State&)>;

  AdaptiveRK45(OdeFunction f, double abs_tol = 1e-6, double rel_tol = 1e-6, double min_dt = 1e-12,
               double max_dt = 0.1)
      : f_(std::move(f)), abs_tol_(abs_tol), rel_tol_(rel_tol), min_dt_(min_dt), max_dt_(max_dt) {}

  // Integrate from t0 to t1 with initial value y0.
  State solve(double t0, const State& y0, double t1) {
    State y = y0;

    // Choose initial dt sign based on whether t1 > t0
    double h = (t1 > t0) ? 1e-3 : -1e-3;

    double t = t0;
    while ((t1 - t) * h > 0.0)  // while we have not passed t1
    {
      // If the remaining time < current step, reduce step to hit t1 exactly
      if (std::fabs(h) > std::fabs(t1 - t)) {
        h = t1 - t;
      }

      double new_t;
      State new_y;
      bool step_ok = try_step(t, y, h, new_t, new_y);
      if (step_ok) {
        // Accept the step
        t = new_t;
        y = new_y;
      } else {
        // Step was rejected; try a smaller step-size.
        continue;
      }

      // If we are  orced below min_dt, break
      if (std::fabs(h) < min_dt_) {
        std::cerr << "Warning: minimum step size reached.\n";
        break;
      }

      std::cout << std::scientific << std::fabs(h) << " " << min_dt_ << " " << max_dt_ << "\n";
    }

    return y;
  }

 private:
  // Perform one Dormand–Prince 5(4) step from time t, state y, with step h.
  // Returns whether the step is accepted; if accepted, updates new_t and new_y.
  bool try_step(double t, const State& y, double& h, double& new_t, State& new_y) {
    // Dormand–Prince coefficients
    static const double c2 = 1.0 / 5.0;
    static const double c3 = 3.0 / 10.0;
    static const double c4 = 4.0 / 5.0;
    static const double c5 = 8.0 / 9.0;

    static const double a21 = 1.0 / 5.0;

    static const double a31 = 3.0 / 40.0;
    static const double a32 = 9.0 / 40.0;

    static const double a41 = 44.0 / 45.0;
    static const double a42 = -56.0 / 15.0;
    static const double a43 = 32.0 / 9.0;

    static const double a51 = 19372.0 / 6561.0;
    static const double a52 = -25360.0 / 2187.0;
    static const double a53 = 64448.0 / 6561.0;
    static const double a54 = -212.0 / 729.0;

    static const double a61 = 9017.0 / 3168.0;
    static const double a62 = -355.0 / 33.0;
    static const double a63 = 46732.0 / 5247.0;
    static const double a64 = 49.0 / 176.0;
    static const double a65 = -5103.0 / 18656.0;

    // b-coeffs for the 5th-order solution
    static const double b1 = 35.0 / 384.0;
    static const double b2 = 0.0;
    static const double b3 = 500.0 / 1113.0;
    static const double b4 = 125.0 / 192.0;
    static const double b5 = -2187.0 / 6784.0;
    static const double b6 = 11.0 / 84.0;

    // b*-coeffs for the 4th-order solution
    static const double b1s = 5179.0 / 57600.0;
    static const double b2s = 0.0;
    static const double b3s = 7571.0 / 16695.0;
    static const double b4s = 393.0 / 640.0;
    static const double b5s = -92097.0 / 339200.0;
    static const double b6s = 187.0 / 2100.0;
    // Note: Full Dormand–Prince sometimes has a b7s * k7. For brevity, we skip it here.

    State k1 = f_(t, y);
    State k2 = f_(t + c2 * h, add(y, mul(k1, a21 * h)));
    State k3 = f_(t + c3 * h, add(y, mul(k1, a31 * h), mul(k2, a32 * h)));
    State k4 = f_(t + c4 * h, add(y, mul(k1, a41 * h), mul(k2, a42 * h), mul(k3, a43 * h)));
    State k5 = f_(t + c5 * h,
                  add(y, mul(k1, a51 * h), mul(k2, a52 * h), mul(k3, a53 * h), mul(k4, a54 * h)));
    State k6 = f_(t + h, add(y, mul(k1, a61 * h), mul(k2, a62 * h), mul(k3, a63 * h),
                             mul(k4, a64 * h), mul(k5, a65 * h)));

    // 5th-order estimate
    State y5 = add(y, mul(k1, b1 * h), mul(k2, b2 * h), mul(k3, b3 * h), mul(k4, b4 * h),
                   mul(k5, b5 * h), mul(k6, b6 * h));

    // 4th-order estimate (for error)
    State y4 = add(y, mul(k1, b1s * h), mul(k2, b2s * h), mul(k3, b3s * h), mul(k4, b4s * h),
                   mul(k5, b5s * h), mul(k6, b6s * h));

    State error_vec = sub(y5, y4);
    double error_norm = norm_impl(error_vec);

    // Evaluate tolerance
    double y_norm = std::max(norm_impl(y5), 1e-12);
    double tolerance = abs_tol_ + rel_tol_ * y_norm;

    bool accept_step = (error_norm < tolerance);
    if (accept_step) {
      // We accept the step
      new_t = t + h;
      new_y = y5;
    }

    // Calculate next step-size
    static const double SAFETY = 0.9;
    static const double EXPONENT = 0.2;  // ~ 1/(order+1). 5th order => 1/6, but 0.2 is typical
    if (error_norm < 1e-14) {
      // Avoid zero in denominator
      error_norm = 1e-14;
    }
    double scale = SAFETY * std::pow(tolerance / error_norm, EXPONENT);
    if (scale < 0.2) scale = 0.2;
    if (scale > 5.0) scale = 5.0;

    // Update h
    h *= scale;
    if (std::fabs(h) > max_dt_) {
      h = (h > 0.0) ? max_dt_ : -max_dt_;
    } else if (std::fabs(h) < min_dt_) {
      h = (h > 0.0) ? min_dt_ : -min_dt_;
    }

    return accept_step;
  }

  // ---------------------------------------------------------------
  // Utilities to handle both scalar (double) and vector states
  // in C++14 without 'if constexpr'.
  // ---------------------------------------------------------------

  // Overload to detect arithmetic vs. vector
  // We'll do this by making small helpers with SFINAE:

  // ---------- For arithmetic -----------
  template <typename T = State>
  typename std::enable_if<std::is_arithmetic<T>::value, T>::type add(const T& s1) {
    return s1;
  }

  template <typename T = State, typename... Ts>
  typename std::enable_if<std::is_arithmetic<T>::value, T>::type add(const T& s1,
                                                                     const Ts&... tail) {
    return add_impl_arithmetic(s1, tail...);
  }

  template <typename T, typename... Ts>
  T add_impl_arithmetic(const T& s1, const T& s2, const Ts&... tail) {
    return add_impl_arithmetic(s1 + s2, tail...);
  }

  template <typename T>
  T add_impl_arithmetic(const T& s) {
    return s;
  }

  template <typename T = State>
  typename std::enable_if<std::is_arithmetic<T>::value, T>::type sub(const T& s1, const T& s2) {
    return s1 - s2;
  }

  template <typename T = State>
  typename std::enable_if<std::is_arithmetic<T>::value, T>::type mul(const T& s, double a) {
    return s * a;
  }

  template <typename T = State>
  typename std::enable_if<std::is_arithmetic<T>::value, double>::type norm_impl(const T& s) {
    return std::fabs(s);
  }

  // ---------- For vector -----------
  template <typename T = State>
  typename std::enable_if<!std::is_arithmetic<T>::value, T>::type add(const T& s1) {
    return s1;  // Single argument, just return it
  }

  template <typename T = State, typename... Ts>
  typename std::enable_if<!std::is_arithmetic<T>::value, T>::type add(const T& s1,
                                                                      const Ts&... tail) {
    return add_impl_vector(s1, tail...);
  }

  template <typename T>
  T add_impl_vector(const T& lhs, const T& rhs) {
    assert(lhs.size() == rhs.size());
    T out(lhs.size());
    for (size_t i = 0; i < lhs.size(); ++i) out[i] = lhs[i] + rhs[i];
    return out;
  }

  // Overload that can handle multiple states (s1 + s2 + s3 ...)
  template <typename T, typename... Ts>
  T add_impl_vector(const T& s1, const T& s2, const Ts&... tail) {
    T tmp = add_impl_vector(s1, s2);
    return add_impl_vector(tmp, tail...);
  }

  template <typename T = State>
  typename std::enable_if<!std::is_arithmetic<T>::value, T>::type sub(const T& lhs, const T& rhs) {
    assert(lhs.size() == rhs.size());
    T out(lhs.size());
    for (size_t i = 0; i < lhs.size(); ++i) out[i] = lhs[i] - rhs[i];
    return out;
  }

  template <typename T = State>
  typename std::enable_if<!std::is_arithmetic<T>::value, T>::type mul(const T& s, double a) {
    T out(s.size());
    for (size_t i = 0; i < s.size(); ++i) out[i] = s[i] * a;
    return out;
  }

  template <typename T = State>
  typename std::enable_if<!std::is_arithmetic<T>::value, double>::type norm_impl(const T& s) {
    double sum = 0.0;
    for (auto val : s) sum += val * val;
    return std::sqrt(sum);
  }

 private:
  OdeFunction f_;
  double abs_tol_;
  double rel_tol_;
  double min_dt_;
  double max_dt_;
};

// ---------------------------------------------------------------------
// Example usage (optional)
// ---------------------------------------------------------------------
#ifdef RUN_EXAMPLE_MAIN

int main() {
  // (1) Scalar example: dy/dt = -y, y(0) = 1
  //     Exact solution => y(t) = e^(-t)
  auto ode_scalar = [](double /*t*/, double y) { return -y; };
  AdaptiveRK45<double> solver_scalar(ode_scalar, 1e-8, 1e-8, 1e-12, 0.1);

  double t0 = 0.0, t1 = 5.0, y0 = 1.0;
  double y_final = solver_scalar.solve(t0, y0, t1);
  std::cout << "Scalar ODE:\n";
  std::cout << "  final t=" << t1 << ", numerical y=" << y_final << ", exact=" << std::exp(-t1)
            << "\n\n";

  // (2) 2D system: dx/dt = y, dy/dt = -x
  //     With x(0)=1, y(0)=0 => x(t) = cos t, y(t) = sin t
  auto ode_2d = [](double /*t*/, const std::vector<double>& state) {
    // state[0] = x, state[1] = y
    std::vector<double> dydt(2);
    dydt[0] = state[1];   // dx/dt = y
    dydt[1] = -state[0];  // dy/dt = -x
    return dydt;
  };
  AdaptiveRK45<std::vector<double>> solver_2d(ode_2d, 1e-8, 1e-8, 1e-12, 0.1);

  std::vector<double> y0_2d = {1.0, 0.0};
  auto result = solver_2d.solve(0.0, y0_2d, 2 * M_PI);

  std::cout << "2D ODE:\n";
  std::cout << "  final t = 2*pi\n"
            << "  numerical = (" << result[0] << ", " << result[1] << ")\n"
            << "  exact     = (1, 0)\n";

  return 0;
}

#endif
