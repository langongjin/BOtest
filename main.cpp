// please see the explanation in the documentation
// http://www.resibots.eu/limbo

#include <iostream>
#include <sys/time.h>

#include <Eigen/Core>

// you can also include <limbo/limbo.hpp> but it will slow down the compilation
#include "limbo/bayes_opt/boptimizer.hpp"

using namespace limbo;

using vec_t = Eigen::VectorXd;
using mat_t = Eigen::MatrixXd;

#define ASSIGNMENT_OP <<

#define USE_NLOPT; //installed NLOPT

// Here we define the parameters
struct Params {
    struct bayes_opt_boptimizer : public defaults::bayes_opt_boptimizer {
    };
//    struct opt_nloptnograd : public defaults::opt_nloptnograd {
//    };
// depending on which internal optimizer we use, we need to import different parameters
#ifdef USE_NLOPT
    struct opt_nloptnograd : public defaults::opt_nloptnograd {
    };
#elif defined(USE_LIBCMAES)
    struct opt_cmaes : public defaults::opt_cmaes {
    };
#else
    struct opt_gridsearch : public defaults::opt_gridsearch {
    };
#endif

    struct kernel : public defaults::kernel {
        BO_PARAM(double, noise, 0.000001);
        BO_PARAM(bool, optimize_noise, false);
    };

    struct bayes_opt_bobase : public defaults::bayes_opt_bobase {
        BO_PARAM(bool, stats_enabled, true);
        BO_PARAM(bool, bounded, true); //false

    };

  struct kernel_exp : public defaults::kernel_exp {
    /// @ingroup kernel_defaults
    BO_PARAM(double, sigma_sq, 1);
    BO_PARAM(double, l, 0.5); // the width of the kernel
  };

  struct kernel_squared_exp_ard : public defaults::kernel_squared_exp_ard {
    /// @ingroup kernel_defaults
    BO_PARAM(int, k, 0); //equivalent to the standard exp ARD
    /// @ingroup kernel_defaults
    BO_PARAM(double, sigma_sq, 1); //brochu2010tutorial p.9 without sigma_sq
  };

//  struct kernel_maternthreehalves
//      : public defaults::kernel_maternthreehalves
//  {
//    BO_PARAM(double, sigma_sq, 0.001);
//    BO_PARAM(double, l, 0.1);
//  };

  struct kernel_maternfivehalves : public defaults::kernel_maternfivehalves
  {
    BO_PARAM(double, sigma_sq, 1); //brochu2010tutorial p.9 without sigma_sq
    BO_PARAM(double, l, 0.2); //characteristic length scale
  };


    struct init_randomsampling : public defaults::init_randomsampling {
        BO_PARAM(int, samples, 10);
    };

    struct stop_maxiterations : public defaults::stop_maxiterations {
        BO_PARAM(int, iterations, 390); //not including randomsampling

    };

  // we use the default parameters for acqui_gpucb
  struct acqui_gpucb : public limbo::defaults::acqui_gpucb {
    //UCB(x) = \mu(x) + \kappa \sigma(x).
    BO_PARAM(double, delta, 0.1); // default delta = 0.1
  };

  // we use the default parameters for acqui_ucb
  struct acqui_ucb : public limbo::defaults::acqui_ucb
  {
    //UCB(x) = \mu(x) + \alpha \sigma(x). high alpha have high exploration
    //iterations is high, alpha can be low for high accuracy in enoughiterations.
    // In contrast, the low iterations should have high alpha for high
    // searching in limited iterations, which guarantee to optimal.
    BO_PARAM(double, alpha, 5); // default alpha = 0.5
    //        BO_DECLARE_DYN_PARAM(int, Params::acqui_ucb, alpha);
  };


};

// Here we define the evaluation function
struct eval_func {
    // number of input dimension (x.size())
    BO_PARAM(size_t, dim_in, 5);
    // number of dimenions of the result (res.size())
    BO_PARAM(size_t, dim_out, 1);

    // the function to be optimized
    Eigen::VectorXd operator()(const Eigen::VectorXd& x) const
    {
//        double y = 20 - std::pow(x(0),2) - std::pow(x(1),2) -
//                   std::pow(x(2),2) - std::pow(x(3),2) - std::pow(x(4),2) -
//                   std::pow(x(5),2) - std::pow(x(6),2) - std::pow(x(7),2) -
//                   std::pow(x(8),2) - std::pow(x(9),2) - std::pow(x(10),2) -
//                   std::pow(x(11),2) - std::pow(x(12),2) - std::pow(x(13),2) -
//                   std::pow(x(14),2) - std::pow(x(15),2) - std::pow(x(16),2) -
//                   std::pow(x(17),2) - std::pow(x(18),2) - std::pow(x(19),2);

      /********** 1# ACKLEY function N Dimensions **************/
//      size_t dim_in = 5;
//      auto xx = x;
//      // transfer interval from [0, 1] to [-32.768, 32.768]
//      for (int i = 0; i < dim_in; i++)
//      {
//        xx[i] = 65.536 * x[i] - 32.768;
//      }
//      const double a = 20.;
//      const double b = 0.2;
//      const double c = 2 * M_PI;
//      double sum1 = 0.;
//      double sum2 = 0.;
//      for (size_t i = 0; i < dim_in; i++)
//      {
//        sum1 = sum1 + xx[i] * xx[i];
//        sum2 = sum2 + std::cos(c * xx[i]);
//      }
//      double term1 = -a * std::exp(-b * std::sqrt(sum1 / dim_in));
//      double term2 = -std::exp(sum2 / dim_in);
//      double obj = term1 + term2 + a + std::exp(1);
//      return tools::make_vector(-obj); //max = 0, at (0,...,0)

      /********** 2# SCHWEFEL function N Dimensions **************/
      size_t dim_in = 5; //todo not accurate results
      auto xx = x;
      // transfer interval from [0, 1] to [-500, 500]
      for (int i = 0; i < dim_in; i++)
      {
        xx[i] = 1000. * x[i] - 500.;
      }
      double sum = 0.;
      for (size_t i = 0; i < dim_in; i++)
      {
        sum = sum + xx[i] * sin(sqrt(abs(xx[i])));
      }
      double obj = 418.9829 * dim_in - sum;
      return tools::make_vector(-obj); //maximum = 0 with (420.9687, ...,420.9687)

      /********** 3# Ellipsoid function N Dimensions **************/
//      size_t dim_in = 5;
//      double inner = 0., outer = 0.;
//      for (size_t i = 0; i < dim_in; ++i)
//      {
//        for(size_t j = 0; j < i; j++)
//        {
//          inner = inner + std::pow((131.072 * x[j] - 65.536), 2); //(-65.536, 65.536)
//        }
//        outer = outer + inner;
//      }
//      return tools::make_vector(-outer); //maximum = 0 at (0, ..., 0)

      /********** 4# Sphere function N Dimensions **************/
//      size_t dim_in = 5;
//      double inner = 0.;
//      for (size_t i = 0; i < dim_in; ++i)
//      {
//        inner = inner + std::pow((10. * x[i] - 5.), 2);
//      }
//      return tools::make_vector(-inner); //max = 0 with (0, 0)

      /********** 5# Rosenbrock  function N Dimensions **************/
//      size_t dim_in = 5;  //todo not good results
//      auto xx = x;
//      // transfer interval from [0, 1] to [-5, 10]
//      for (int i = 0; i < dim_in; i++)
//        xx[i] = 15. * x[i] - 5.;
//      double sum = 0.;
//      double term = 0.;
//      double xnext = 0.;
//      for(size_t i = 0; i < (dim_in - 1); i++)
//      {
//        xnext = xx[i + 1];
//        term = 100. * std::pow((xnext - xx[i] * xx[i]), 2.0) + std::pow((xx[i] - 1), 2.0);
//        sum = sum + term;
//      }
//      double obj = 0.001 * sum;
////      double obj = (sum - 382700)/375500.; //rescale
//      return tools::make_vector(-obj); //maximum = 0 with (1,...,1)

      /********* 6# Michalewicz function N = 2/5/10 Dimensions **********/
//      size_t dim_in = 5;
//      auto xx = x;
//      // transfer interval from [0, 1] to [0, pi]
//      for (int i = 0; i < dim_in; i++)
//        xx[i] = M_PI * x[i];
//      double sum = 0.;
//      double term = 0.;
//      double m = 10.;
//      for(size_t i = 0; i < dim_in; i++)
//      {
//        term = std::sin(xx[i]) * std::pow(std::sin(i * xx[i] * xx[i]/M_PI), 2 * m);
//        sum = sum + term;
//      }
//      double obj = sum;
//      return tools::make_vector(obj); //max= -1.8013(2D) at (2.20,1.57)/-4.687658(5D)/-9.66015(10D)

      /********** 7# StyblinskiTang function N Dimensions ****************/
//      size_t dim_in = 5;
//      auto xx = x;
//      // transfer interval from [0, 1] to [-5, 5]
//      for (int i = 0; i < dim_in; i++)
//        xx[i] = 10 * x[i] - 5;
//      double sum = 0.;
//      double term;
//      for(size_t i = 0; i < dim_in; i++)
//      {
//        term = std::pow(xx[i], 4.0) - 16 * xx[i] * xx[i] + 5 * xx[i];
//        sum = sum + term;
//      }
//      double obj = sum/2.0;
//      //max= 39.16599 * d, (5D:195.82995) at (-2.903534,...,-2.903534)
//      return tools::make_vector(-obj);

      /********** 8# Powell function N >= 4 Dimensions ****************/
//      size_t dim_in = 5;
//      auto xx = x;
//      // transfer interval from [0, 1] to [-4, 5]
//      for (int i = 0; i < dim_in; i++)
//        xx[i] = 9 * x[i] - 4;
//      double sum = 0.;
//      double term1, term2, term3, term4;
//      for(size_t i = 0; i < dim_in/4; i++)
//      {
//        term1 = std::pow((xx[4 * i - 3] + 10 * xx[4 * i -2]), 2.0);
//        term2 = 5 * std::pow((xx[4 * i - 1] - xx[4 * i]), 2.0);
//        term3 = std::pow(xx[4 * i - 2] - 2 * xx[4 * i - 1], 4.0);
//        term4 = 10 * std::pow(xx[4 * i - 3] - xx[4 * i], 4.0);
//        sum = sum + term1 + term2 + term3 + term4;
//      }
//      double obj = sum;
//      return tools::make_vector(-obj); //max= 0 at (0,...,0)

      ///********** Branin function 2 Dimensions **************///
//        auto x1 = x(0) * 15 - 5;
//        auto x2 = x(1) * 15;
//        auto term1 = sqr(x2 - (5.1 * sqr(x1) / (4. * sqr(M_PI))) + 5. * x1 / M_PI - 6);
//        auto term2 = (10. - 10. / (8. * M_PI)) * std::cos(x1);
//        //double y = (term1 + term2 - 44.81) / 51.95;
//        auto y = - (term1 + term2 + 10);
//        // we return a 1-dimensional vector
//        return tools::make_vector(y); //maximum=0.397887 with (-/+pi,12.275) and (9.42478,2.475)

        /********** Ellipsoid function 2 Dimensions ***********///
//        size_t dim_in = 2;
//        double inner = 0., outer = 0.;
//        for (size_t i = 0; i < dim_in; ++i)
//        {
//            for(size_t j = 0; j < i; j++)
//            {
//                inner = inner + std::pow((131.072 * x[j] - 65.536), 2); //(-65.536, 65.536)
//            }
//            outer = outer + inner;
//        }
//        return tools::make_vector(-outer); //maximum = 0 with (0, 0)

        /********** Goldsteinprice function 2 Dimensions ******///
//        auto xx = x;
//        // x.size = 2, transfer interval from [0, 1] to [-2, 2]
//        for (int i = 0; i < xx.size(); i++)
//            xx[i] = 4. * x[i] - 2.;
//
//        double fact1a = sqr(xx[0] + xx[1] + 1.);
//        double fact1b = 19. - 14. * xx[0] + 3. * sqr(xx[0]) - 14. * xx[1] + 6. * xx[0] * xx[1] + 3. * sqr(xx[1]);
//        double fact1 = 1. + fact1a * fact1b;
//
//        double fact2a = sqr(2. * xx[0] - 3. * xx[1]);
//        double fact2b = 18. - 32. * xx[0] + 12. * sqr(xx[0]) + 48. * xx[1] - 36. * xx[0] * xx[1] + 27. * sqr(xx[1]);
//        double fact2 = 30. + fact2a * fact2b;
//
//        double obj = 0.01 * fact1 * fact2;
//
//        return tools::make_vector(-obj); //maximum = 0.03 with (0,-1);

        /********** Hartmann3 function 3 Dimensions ***********///
//        mat_t a(4, 3);
//        mat_t p(4, 3);
//        a ASSIGNMENT_OP 3.0, 10., 30.,
//                0.1, 10., 35.,
//                3.0, 10., 30.,
//                0.1, 10., 35.;
//        p ASSIGNMENT_OP 0.3689, 0.1170, 0.2673,
//                0.4699, 0.4387, 0.7470,
//                0.1091, 0.8732, 0.5547,
//                0.0381, 0.5743, 0.8828;
//        vec_t alpha(4);
//        alpha ASSIGNMENT_OP 1.0, 1.2, 3.0, 3.2;
//
//        double res = 0.;
//        for (int i = 0; i < 4; i++) {
//            double s = 0.;
//            for (size_t j = 0; j < 3; j++) {
//                s += a(i, j) * sqr(x[j] - p(i, j));
//            }
//            res += alpha(i) * std::exp(-s);
//        }
//        return tools::make_vector(res); //maximum = -3.86278 with (0.114614,0.555649, 0.852547);

        /********** Hartmann6 function 6 Dimensions ***********///
//        mat_t a(4, 6);
//        mat_t p(4, 6);
//        a ASSIGNMENT_OP 10., 3., 17., 3.5, 1.7, 8.,
//                0.05, 10., 17., 0.1, 8., 14.,
//                3., 3.5, 1.7, 10., 17., 8.,
//                17., 8., 0.05, 10., 0.1, 14.;
//        p ASSIGNMENT_OP 0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886,
//                0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991,
//                0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650,
//                0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381;
//
//        vec_t alpha(4);
//        alpha ASSIGNMENT_OP 1.0, 1.2, 3.0, 3.2;
//
//        double res = 0.;
//        for (int i = 0; i < 4; i++) {
//            double s = 0.;
//            for (size_t j = 0; j < 6; j++) {
//                s += a(i, j) * sqr(x[j] - p(i, j));
//            }
//            res += alpha(i) * std::exp(-s);
//        }
//        return tools::make_vector(res); //maximum = -3.32237 with (0.20169,0.150011,0.476874,0.275332,0.311652, 0.6573);

        /********** Rastrigin function 4 Dimensions ***********///
//        size_t dim_in = 4;
//        auto xx = x;
//        for (int i = 0; i < xx.size(); i++)
//            xx[i] = 2. * x[i] - 1.;
//        double f = 10. * dim_in;
//        for (size_t i = 0; i < dim_in; ++i)
//            f += xx[i] * xx[i] - 10. * std::cos(2 * M_PI * xx[i]);
//        return tools::make_vector(-f); //maximum = 0 with (0, 0, 0, 0);

        /********** Sixhumpcamel function 2 Dimensions ********///
//        double x1 = -3 + 6 * x[0];
//        double x2 = -2 + 4 * x[1];
//        double x1_2 = sqr(x1);
//        double x2_2 = sqr(x2);
//
//        double tmp1 = (4 - 2.1 * x1_2 + sqr(x1_2) / 3.) * x1_2;
//        double tmp2 = x1 * x2;
//        double tmp3 = (-4 + 4 * x2_2) * x2_2;
//        double obj = tmp1 + tmp2 + tmp3;
//        return tools::make_vector(-obj); //maximum = 1.0316 with (0.0898, -0.7126) and (-0.0898, 0.7126);

        /********** Sphere function 2 Dimensions **************///
//        size_t dim_in = 2;
//        double inner = 0.;
//        for (size_t i = 0; i < dim_in; ++i)
//        {
//            inner = inner + std::pow((10. * x[i] - 5.), 2);
//        }
//        return tools::make_vector(-inner); //maximum = 0 with (0, 0)

    }
};

int main()
{
    struct timeval timeStartall, timeEndall;
    double timeDiffall;
    gettimeofday(&timeStartall,NULL);

  /*********** configuration ************/
//  using Kernel_t = kernel::SquaredExpARD<Params>;
  using Kernel_t = kernel::MaternFiveHalves <Params>;
  using Mean_t = mean::Data<Params>;
  using GP_t = model::GP<Params, Kernel_t, Mean_t>;
  using Acqui_t = acqui::GP_UCB<Params, GP_t>;


  int iters = 10; // the number of run
    for(int run = 0; run < iters; run++)
    {
        // we use the default acquisition function / model / stat / etc.
//        bayes_opt::BOptimizer<Params> boptimizer;
      bayes_opt::BOptimizer<Params, modelfun<GP_t>, acquifun<Acqui_t>> boptimizer;
      // run the evaluation
        boptimizer.optimize(eval_func());
        // the best sample found
//        std::cout << "Best sample: " << boptimizer.best_sample().transpose() << " - Best observation: " << boptimizer.best_observation()(0) << std::endl;
    }
    gettimeofday(&timeEndall,NULL);

    timeDiffall = 1000000 * (timeEndall.tv_sec - timeStartall.tv_sec)
                  + timeEndall.tv_usec - timeStartall.tv_usec; //tv_sec: value of second, tv_usec:value of microsecond
    timeDiffall/=1000;
//    std::cout << "timeDiffall: " << timeDiffall << std::endl;

    return 0;
}
